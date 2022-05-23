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
) -> Geom {
    // TODO handle error
    let _ = std::fs::remove_dir_all("opt");
    let _ = std::fs::create_dir("opt");

    static OPT_TMPL: Template =
        Template::from("scfcrt=1.D-21 aux(precision=14) PM6");
    let opt = Job::new(
        Mopac::new(
            "opt/opt".to_string(),
            Some(Rc::new(params.clone())),
            Rc::new(geom),
            0,
            &OPT_TMPL,
        ),
        0,
    );
    queue.optimize(opt)
}

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
    fn build_jobs<Q: Queue<Mopac>>(
        &self,
        params: &Params,
        start_index: usize,
        job_num: usize,
        submitter: &Q,
        charge: isize,
    ) -> (FreqParts, Vec<Job<Mopac>>) {
        let mut intder = self.intder.clone();
        // optimize
        setup();
        // setup because submitter tries to write in inp, not opt
        let geom =
            optimize_geometry(self.config.geometry.clone(), params, submitter);
        // generate pts == moles
        let (moles, taylor, taylor_disps, atomic_numbers) =
            rust_pbqff::generate_pts(geom, &mut intder);
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
        // start_index and job_num always zero for this since I'm only running
        // one set of jobs at a time
        let (mut freq, mut jobs) =
            self.build_jobs(params, 0, 0, submitter, charge);
        let mut energies = vec![0.0; jobs.len()];
        // drain to get energies
        submitter.drain(&mut jobs, &mut energies);
        takedown();
        // convert energies to frequencies and return those
        let res = na::DVector::from(
            rust_pbqff::freqs(
                &mut energies,
                &mut freq.intder,
                &freq.taylor,
                &freq.taylor_disps,
                &freq.atomic_numbers,
                &self.spectro,
                &self.config.gspectro_cmd,
                &self.config.spectro_cmd,
            )
            .corr,
        );
        // this gets made inside `freqs
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
	let rows = params.values.len();
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
                self.build_jobs(&pf, idx, idx, submitter, charge);
            idx += fwd_jobs.len();
            jobs.extend(fwd_jobs);
            freqs.push(freq);

            let mut pb = params.clone();
            pb.values[row] -= DELTA;
            let (freq, bwd_jobs) =
                self.build_jobs(&pb, idx, idx, submitter, charge);
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

        // calculate all of the frequencies and assemble the jacobian
        let mut jac_t = Vec::new();
        let mut energies = energies.chunks_exact_mut(cols);
        let mut freqs = freqs.iter_mut();
	// TODO more parallel
        for _ in 0..rows {
            // need to pull out two slices of the energies and two freqs
            let mut fwd_energies = energies.next().unwrap();
	    let fwd_freq = freqs.next().unwrap();
            let fwd_freqs = na::DVector::from(
                rust_pbqff::freqs(
                    &mut fwd_energies,
                    &mut fwd_freq.intder,
                    &fwd_freq.taylor,
                    &fwd_freq.taylor_disps,
                    &fwd_freq.atomic_numbers,
                    &self.spectro,
                    &self.config.gspectro_cmd,
                    &self.config.spectro_cmd,
                )
                .corr,
            );

            let mut bwd_energies = energies.next().unwrap();
	    let bwd_freq = freqs.next().unwrap();
            let bwd_freqs = na::DVector::from(
                rust_pbqff::freqs(
                    &mut bwd_energies,
                    &mut bwd_freq.intder,
                    &bwd_freq.taylor,
                    &bwd_freq.taylor_disps,
                    &bwd_freq.atomic_numbers,
                    &self.spectro,
                    &self.config.gspectro_cmd,
                    &self.config.spectro_cmd,
                )
                .corr,
            );

            // (fwd_freqs - bwd_freqs) / 2DELTA gives a column of jac_t
            jac_t.extend_from_slice(
                ((fwd_freqs - bwd_freqs) / (2. * DELTA)).as_slice(),
            );
        }
        assert!(energies.into_remainder().is_empty());
        let cols = jac_t.len() / rows;
        na::DMatrix::from_vec(cols, rows, jac_t)
    }

    fn stat_multiplier(&self) -> f64 {
        1.0
    }
}
