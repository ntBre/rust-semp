use na::DVector;
use psqs::geom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::program::{Job, Template};
use psqs::queue::Queue;
use rust_pbqff::coord_type::generate_pts;
use rust_pbqff::Intder;
use symm::{Irrep, Molecule, PointGroup};
use taylor::Taylor;

use crate::{config, setup, sort_irreps, takedown, MOPAC_TMPL};
use nalgebra as na;
use std::fs::{read_to_string, File};
use std::io::Write;
use std::rc::Rc;
use std::sync::Mutex;

use super::Optimize;

static DELTA: f64 = 1e-4;

static DEBUG: bool = true;

fn output_stream() -> Box<dyn Write> {
    if DEBUG {
        Box::new(std::io::stderr())
    } else {
        Box::new(std::io::sink())
    }
}

pub struct Frequency {
    pub config: rust_pbqff::config::Config,
    pub intder: rust_pbqff::Intder,
    pub spectro: rust_pbqff::Spectro,
    pub dummies: Vec<(usize, usize)>,
    pub reorder: bool,
    pub irreps: Vec<Irrep>,
    logger: Mutex<File>,
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

type Dummies = Vec<(usize, usize)>;

impl Frequency {
    pub fn new(
        config: rust_pbqff::config::Config,
        intder: rust_pbqff::Intder,
        spectro: rust_pbqff::Spectro,
        dummies: Dummies,
        reorder: bool,
        irreps: Vec<Irrep>,
    ) -> Self {
        let logger = Mutex::new(
            std::fs::File::create("freqs.log")
                .expect("failed to create 'freqs.log'"),
        );
        Self {
            config,
            intder,
            spectro,
            dummies,
            logger,
            reorder,
            irreps,
        }
    }

    pub fn load_irreps(filename: &str) -> Vec<Irrep> {
        let data = read_to_string(filename).expect("failed to load irrep file");
        data.lines()
            .map(|s| match s.trim().parse() {
                Ok(v) => v,
                Err(_) => panic!("failed to parse irrep {}", s),
            })
            .collect()
    }

    fn build_jobs<W>(
        &self,
        w: &mut W,
        geom: Geom,
        params: &Params,
        start_index: usize,
        job_num: usize,
        charge: isize,
    ) -> (FreqParts, Vec<Job<Mopac>>)
    where
        W: Write,
    {
        let mol = {
            let mut mol = Molecule::new(geom.xyz().unwrap().to_vec());
            mol.normalize();
            if self.reorder {
                mol.reorder();
            }
            mol
        };
        const SYMM_EPS: f64 = 1e-6;
        let pg = mol.point_group_approx(SYMM_EPS);

        // TODO assuming that the intder coordinates max out at C2v, this was
        // the case for ethylene
        let pg = if let PointGroup::D2h { axes, planes } = pg {
            PointGroup::C2v {
                axis: axes[0],
                planes: [planes[1], planes[2]],
            }
        } else {
            pg
        };

        writeln!(w, "Normalized Geometry:\n{:20.12}", mol).unwrap();
        writeln!(w, "Point Group = {}", pg).unwrap();

        let mut intder = self.intder.clone();
        let (moles, taylor, taylor_disps, atomic_numbers) = generate_pts(
            w,
            &mol,
            &pg,
            &mut intder,
            self.config.step_size,
            &self.dummies,
        );

        // dir created in generate_pts but unused here
        let _ = std::fs::remove_dir_all("pts");

        // call build_jobs like before
        let jobs = Mopac::build_jobs(
            &moles,
            Some(&params),
            "inp",
            start_index,
            1.0,
            job_num,
            charge,
            &MOPAC_TMPL,
        );
        (
            FreqParts::new(intder, taylor, taylor_disps, atomic_numbers),
            jobs,
        )
    }
    /// call `rust_pbqff::coord_type::freqs`, but sort the frequencies by irrep and
    /// then frequency and return the result as a DVector
    pub fn freqs<W: std::io::Write>(
        &self,
        w: &mut W,
        dir: &str,
        energies: &mut [f64],
        intder: &mut Intder,
        taylor: &Taylor,
        taylor_disps: &Vec<Vec<isize>>,
        atomic_numbers: &Vec<usize>,
    ) -> DVector<f64> {
        let summary = rust_pbqff::coord_type::freqs(
            w,
            &dir,
            energies,
            intder,
            taylor,
            taylor_disps,
            atomic_numbers,
            &self.spectro,
            &self.config.gspectro_cmd,
            &self.config.spectro_cmd,
            self.config.step_size,
        );
        let freqs = sort_irreps(&summary.corr, &summary.irreps);
        DVector::from(freqs)
    }
}

impl Optimize for Frequency {
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q: Queue<Mopac>>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[config::Molecule],
    ) -> DVector<f64> {
        let mut w = output_stream();

        writeln!(w, "Params:\n{}", params.to_string()).unwrap();

        setup();
        // TODO the Molecules should contain their own geometries
        let geom = optimize_geometry(
            self.config.geometry.clone(),
            params,
            submitter,
            "inp",
            "opt",
        );

        writeln!(w, "Optimized Geometry:\n{:20.12}", geom).unwrap();

        let (mut freq, mut jobs) =
            self.build_jobs(&mut w, geom, params, 0, 0, molecules[0].charge);
        writeln!(
            w,
            "\n{} atoms require {} jobs",
            freq.atomic_numbers.len(),
            jobs.len()
        )
        .unwrap();

        let mut energies = vec![0.0; jobs.len()];
        // drain to get energies
        setup();
        submitter.drain(&mut jobs, &mut energies);
        takedown();
        // convert energies to frequencies and return those
        let _ = std::fs::create_dir("freqs");
        let res = self.freqs(
            &mut w,
            "freqs",
            &mut energies,
            &mut freq.intder,
            &freq.taylor,
            &freq.taylor_disps,
            &freq.atomic_numbers,
        );
        let _ = std::fs::remove_dir_all("freqs");
        res
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned
    fn num_jac<Q: Queue<Mopac>>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[config::Molecule],
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
                        molecules[0].charge,
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
                        molecules[0].charge,
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
            let (freq, fwd_jobs) = self.build_jobs(
                &mut output_stream(),
                geoms[2 * row].clone(),
                &pf,
                idx,
                idx,
                molecules[0].charge,
            );
            idx += fwd_jobs.len();
            jobs.extend(fwd_jobs);
            freqs.push(freq);

            let mut pb = params.clone();
            pb.values[row] -= DELTA;
            let (freq, bwd_jobs) = self.build_jobs(
                &mut output_stream(),
                geoms[2 * row + 1].clone(),
                &pb,
                idx,
                idx,
                molecules[0].charge,
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
                self.freqs(
                    &mut output_stream(),
                    &dir,
                    &mut energy,
                    &mut freq.intder,
                    &freq.taylor,
                    &freq.taylor_disps,
                    &freq.atomic_numbers,
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

    fn log(&self, iter: usize, got: &DVector<f64>, want: &DVector<f64>) {
        let mut logger = self.logger.lock().unwrap();
        write!(logger, "{:5}", iter).unwrap();
        for g in got {
            write!(logger, "{:8.1}", g).unwrap();
        }
        writeln!(logger, "{:8.1}", mae(got, want)).unwrap();
    }
}

fn mae(a: &DVector<f64>, b: &DVector<f64>) -> f64 {
    (a - b).abs().sum() / a.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]

    fn test_mae() {
        let se = na::dvector![
            3130.3, 3107.9, 1590.1, 1227.0, 1117.1, 1002.9, 876.2, 925.6, 773.6
        ];
        let ai = na::dvector![
            3142.6, 3120.8, 1600.9, 1278.8, 1065.1, 970.2, 888.6, 878.6, 772.8
        ];
        approx::assert_abs_diff_eq!(mae(&se, &ai), 25.855555555555547);
    }
}
