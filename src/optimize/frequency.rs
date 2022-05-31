use psqs::geom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::program::{Job, Template};
use psqs::queue::Queue;

use crate::{build_jobs, setup, takedown};
use nalgebra as na;
use std::rc::Rc;

use super::Optimize;

static DELTA: f64 = 1e-4;

pub struct Frequency {
    pub config: rust_pbqff::config::Config,
    pub intder: rust_pbqff::Intder,
    pub spectro: rust_pbqff::Spectro,
}

pub fn optimize_geometry<Q: Queue<Mopac>>(
    geom: Geom,
    params: &Params,
    queue: &Q,
    dir: &str,
    name: &str,
) -> Geom {
    static OPT_TMPL: Template =
        Template::from("scfcrt=1.D-21 aux(precision=14) PM6");
    let opt = Job::new(
        Mopac::new(
            format!("{}/{}", dir, name),
            Some(Rc::new(params.clone())),
            Rc::new(geom),
            0,
            &OPT_TMPL,
        ),
        0,
    );
    let mut res = vec![Geom::default(); 1];
    queue.optimize(&mut [opt], &mut res);
    res.pop().unwrap()
}

#[derive(Clone)]
struct FreqParts {
    intder: rust_pbqff::Intder,
    taylor: taylor::Taylor,
    taylor_disps: Vec<Vec<isize>>,
    atomic_numbers: Vec<usize>,
}

impl FreqParts {
    fn new(
        intder: rust_pbqff::Intder,
        taylor: taylor::Taylor,
        taylor_disps: Vec<Vec<isize>>,
        atomic_numbers: Vec<usize>,
    ) -> Self {
        Self {
            intder,
            taylor,
            taylor_disps,
            atomic_numbers,
        }
    }
}

impl Frequency {
    fn build_jobs(
        &self,
        geom: Geom,
        params: &Params,
        start_index: usize,
        job_num: usize,
        charge: isize,
    ) -> (FreqParts, Vec<Job<Mopac>>) {
        let mut intder = self.intder.clone();
        // generate pts == moles
        let (moles, taylor, taylor_disps, atomic_numbers) =
            rust_pbqff::coord_type::generate_pts(
                geom,
                &mut intder,
                self.config.step_size,
            );
        // dir created in generate_pts but unused here
        let _ = std::fs::remove_dir_all("pts");
        // call build_jobs like before
        let jobs =
            build_jobs(&moles, params, start_index, 1.0, job_num, charge);
        (
            FreqParts::new(intder, taylor, taylor_disps, atomic_numbers),
            jobs,
        )
    }
}

impl Optimize for Frequency {
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q: Queue<Mopac>>(
        &self,
        _moles: &Vec<Rc<Geom>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DVector<f64> {
        setup();
        let geom = optimize_geometry(
            self.config.geometry.clone(),
            params,
            submitter,
            "inp",
            "opt",
        );
        // start_index and job_num always zero for this since I'm only running
        // one set of jobs at a time
        let (mut freq, mut jobs) = self.build_jobs(geom, params, 0, 0, charge);
        let mut energies = vec![0.0; jobs.len()];
        // drain to get energies
        setup();
        submitter.drain(&mut jobs, &mut energies);
        takedown();
        // convert energies to frequencies and return those
        let _ = std::fs::create_dir("freqs");
        let res = na::DVector::from(
            rust_pbqff::coord_type::freqs(
                "freqs",
                &mut energies,
                &mut freq.intder,
                &freq.taylor,
                &freq.taylor_disps,
                &freq.atomic_numbers,
                &self.spectro,
                &self.config.gspectro_cmd,
                &self.config.spectro_cmd,
                self.config.step_size,
            )
            .corr,
        );
        let _ = std::fs::remove_dir_all("freqs");
        res
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned
    fn num_jac<Q: Queue<Mopac>>(
        &self,
        _moles: &Vec<Rc<Geom>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DMatrix<f64> {
        let start = std::time::SystemTime::now();

        let rows = params.values.len();
        // build all the optimizations
        static OPT_TMPL: Template =
            Template::from("scfcrt=1.D-21 aux(precision=14) PM6");
        let mut opts = Vec::new();
        for row in 0..rows {
            // forward
            {
                let mut pf = params.clone();
                pf.values[row] += DELTA;
                opts.push(Job::new(
                    Mopac::new(
                        format!("inp/opt{row}_fwd"),
                        Some(Rc::new(pf)),
                        Rc::new(self.config.geometry.clone()),
                        charge,
                        &OPT_TMPL,
                    ),
                    2 * row,
                ));
            }

            // backward
            {
                let mut pb = params.clone();
                pb.values[row] -= DELTA;
                opts.push(Job::new(
                    Mopac::new(
                        format!("inp/opt{row}_bwd"),
                        Some(Rc::new(pb)),
                        Rc::new(self.config.geometry.clone()),
                        charge,
                        &OPT_TMPL,
                    ),
                    2 * row + 1,
                ));
            }
        }
        setup();
        let mut geoms = vec![Default::default(); opts.len()];
        assert!(geoms.len() == 2 * rows);
        submitter.optimize(&mut opts, &mut geoms);

        eprintln!(
            "finished optimizations after {:.1} sec",
            start.elapsed().unwrap().as_millis() as f64 / 1000.
        );
        let start = std::time::SystemTime::now();

        // build all the jobs, including optimizations
        let mut idx = 0;
        let mut jobs = Vec::new();
        let mut freqs = Vec::new();
        let mut cols = 0;
        for row in 0..rows {
            let mut pf = params.clone();
            pf.values[row] += DELTA;
            // idx = job_num so use it twice
            let (freq, fwd_jobs) =
                self.build_jobs(geoms[2 * row].clone(), &pf, idx, idx, charge);
            idx += fwd_jobs.len();
            jobs.extend(fwd_jobs);
            freqs.push(freq);

            let mut pb = params.clone();
            pb.values[row] -= DELTA;
            let (freq, bwd_jobs) = self.build_jobs(
                geoms[2 * row + 1].clone(),
                &pb,
                idx,
                idx,
                charge,
            );
            idx += bwd_jobs.len();
            if row == 0 {
                // set this once
                cols = bwd_jobs.len();
            }
            jobs.extend(bwd_jobs);
            freqs.push(freq);
        }
        // run all of the energies
        let mut energies = vec![0.0; idx];
        setup();
        submitter.drain(&mut jobs, &mut energies);
        takedown();

        eprintln!(
            "finished energies after {:.1} sec",
            start.elapsed().unwrap().as_millis() as f64 / 1000.
        );
        let start = std::time::SystemTime::now();

        use rayon::prelude::*;

        // calculate all of the frequencies and assemble the jacobian
        let energies = energies.chunks_exact_mut(cols).map(|c| c.to_vec());
        let pairs: Vec<_> = energies.zip(freqs).collect();
        // TODO more parallel
        let freqs: Vec<_> = pairs
            .par_iter()
            .enumerate()
            .map(|(i, (energy, freq))| {
                let mut energy = energy.clone();
                let mut freq = freq.clone();
                let dir = format!("freqs{}", i);
                let _ = std::fs::create_dir(&dir);
                na::DVector::from(
                    rust_pbqff::coord_type::freqs(
                        &dir,
                        &mut energy,
                        &mut freq.intder,
                        &freq.taylor,
                        &freq.taylor_disps,
                        &freq.atomic_numbers,
                        &self.spectro,
                        &self.config.gspectro_cmd,
                        &self.config.spectro_cmd,
                        self.config.step_size,
                    )
                    .corr,
                )
            })
            .collect();
        // remove the directories created by the iteration above
        for i in 0..pairs.len() {
            let dir = format!("freqs{}", i);
            let _ = std::fs::remove_dir_all(&dir);
        }
        let jac_t: Vec<_> = freqs
            .chunks(2)
            .map(|pair| {
                // use 'let else' if that gets stabilized
                let (fwd, bwd) = if let [fwd, bwd] = pair {
                    (fwd, bwd)
                } else {
                    panic!("unmatched pair");
                };
                (fwd - bwd) / (2. * DELTA)
            })
            .collect();

        eprintln!(
            "finished frequencies after {:.1} sec",
            start.elapsed().unwrap().as_millis() as f64 / 1000.
        );

        na::DMatrix::from_columns(&jac_t)
    }

    fn stat_multiplier(&self) -> f64 {
        1.0
    }
}
