use na::DVector;
use psqs::geom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::program::{Job, ProgramResult, Template};
use psqs::queue::Queue;
use rust_pbqff::coord_type::cart::{make_fcs, BigHash};
use rust_pbqff::coord_type::{sic, Cart};
use symm::Molecule;

use crate::utils::sort_irreps;
use crate::{config, utils::setup, utils::takedown};
use nalgebra as na;
use std::fs::File;
use std::io::Write;
use std::rc::Rc;
use std::sync::Mutex;

use super::Optimize;

static DELTA: f64 = 1e-4;

use std::sync::Once;
static mut DEBUG: bool = false;
static INIT: Once = Once::new();

/// check the `SEMP_DUMP_DEBUG` environment variable on first call and from then
/// on report whether or not it was set to `1`
fn is_debug() -> bool {
    unsafe {
        INIT.call_once(|| {
            let v = std::env::var("SEMP_FREQ_DEBUG").unwrap_or_default();
            if v == "1" {
                DEBUG = true;
            }
        });
        DEBUG
    }
}

/// pbqff step size
const STEP_SIZE: f64 = 0.005;

fn output_stream() -> Box<dyn Write> {
    if is_debug() {
        Box::new(std::io::stderr())
    } else {
        Box::new(std::io::sink())
    }
}

pub struct Frequency {
    logger: Mutex<File>,
}

pub fn optimize_geometry<Q: Queue<Mopac> + std::marker::Sync>(
    geom: Geom,
    params: &Params,
    queue: &Q,
    dir: &str,
    name: &str,
    charge: isize,
    template: Template,
) -> Option<ProgramResult> {
    let opt = Job::new(
        Mopac::new_full(
            format!("{}/{}", dir, name),
            Some(params.clone()),
            geom,
            charge,
            template,
        ),
        0,
    );
    let mut res = vec![Default::default(); 1];
    let status = queue.energize(dir, &mut [opt], &mut res);
    if status.is_err() {
        return None;
    }
    Some(res.pop().unwrap())
}

#[derive(Clone, Debug)]
enum FreqParts {
    Sic {
        intder: rust_pbqff::Intder,
        taylor: taylor::Taylor,
        taylor_disps: taylor::Disps,
        atomic_numbers: Vec<usize>,
    },
    Cart {
        fcs: Vec<f64>,
        target_map: BigHash,
        n: usize,
        nfc2: usize,
        nfc3: usize,
        mol: Molecule,
    },
}

impl FreqParts {
    fn sic(
        intder: rust_pbqff::Intder,
        taylor: taylor::Taylor,
        taylor_disps: taylor::Disps,
        atomic_numbers: Vec<usize>,
    ) -> Self {
        Self::Sic {
            intder,
            taylor,
            taylor_disps,
            atomic_numbers,
        }
    }

    fn cart(
        fcs: Vec<f64>,
        target_map: BigHash,
        n: usize,
        nfc2: usize,
        nfc3: usize,
        mol: Molecule,
    ) -> Self {
        Self::Cart {
            fcs,
            target_map,
            n,
            nfc2,
            nfc3,
            mol,
        }
    }
}

impl Frequency {
    pub fn new() -> Self {
        let logger = Mutex::new(
            std::fs::File::create("freqs.log")
                .expect("failed to create 'freqs.log'"),
        );
        Self { logger }
    }

    /// build jobs for a fixed set of Params. instead of passing the params as
    /// part of each separate Mopac, write the parameters once and update the
    /// Template passed to Mopac::new
    fn build_jobs<W>(
        &self,
        w: &mut W,
        geom: ProgramResult,
        params: &Params,
        start_index: usize,
        job_num: usize,
        molecule: &config::Molecule,
    ) -> (FreqParts, Vec<Job<Mopac>>)
    where
        W: Write,
    {
        let mol = {
            let mut mol = Molecule::new(geom.cart_geom.unwrap());
            mol.normalize();
            mol
        };
        const SYMM_EPS: f64 = 1e-6;
        let pg = mol.point_group_approx(SYMM_EPS);

        writeln!(w, "Normalized Geometry:\n{:20.12}", mol).unwrap();
        writeln!(w, "Point Group = {}", pg).unwrap();

        let param_file = Rc::new(format!("tmparam/{}.dat", job_num));
        Mopac::write_params(params, &param_file);
        let mut tmpl = molecule.template.clone();
        use std::fmt::Write;
        write!(tmpl.header, " external={}", param_file).unwrap();

        match &molecule.intder {
            Some(intder) => {
                let mut intder = intder.clone();
                // NOTE: assuming that the intder coordinates max out at C2v,
                // this was the case for ethylene
                let pg = if pg.is_d2h() {
                    pg.subgroup(symm::Pg::C2v).unwrap()
                } else {
                    pg
                };

                let (moles, taylor, taylor_disps, atomic_numbers) =
                    sic::generate_pts(w, &mol, &pg, &mut intder, STEP_SIZE);

                // dir created in generate_pts but unused here
                let _ = std::fs::remove_dir_all("pts");

                // call build_jobs like before
                let jobs = Mopac::build_jobs(
                    &moles,
                    None,
                    "inp",
                    start_index,
                    1.0,
                    job_num,
                    molecule.charge,
                    tmpl,
                );
                (
                    FreqParts::sic(
                        intder,
                        taylor,
                        taylor_disps,
                        atomic_numbers,
                    ),
                    jobs,
                )
            }
            None => {
                let n = 3 * mol.atoms.len();
                let nfc2 = n * n;
                let nfc3 = n * (n + 1) * (n + 2) / 6;
                let nfc4 = n * (n + 1) * (n + 2) * (n + 3) / 24;
                let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];

                let mut target_map = BigHash::new(mol.clone(), pg);

                let geoms = Cart.build_points(
                    Geom::Xyz(mol.atoms.clone()),
                    STEP_SIZE,
                    geom.energy,
                    nfc2,
                    nfc3,
                    &mut fcs,
                    &mut target_map,
                );
                let dir = "inp";
                let mut job_num = start_index;
                let mut jobs = Vec::new();
                for mol in geoms {
                    let filename = format!("{dir}/job.{:08}", job_num);
                    job_num += 1;
                    jobs.push(Job::new(
                        Mopac::new_full(
                            filename,
                            None,
                            mol.geom.clone(),
                            molecule.charge,
                            tmpl.clone(),
                        ),
                        mol.index + start_index,
                    ));
                }
                (FreqParts::cart(fcs, target_map, n, nfc2, nfc3, mol), jobs)
            }
        }
    }

    /// call `rust_pbqff::coord_type::freqs`, but sort the frequencies by irrep
    /// and then frequency and return the result as a DVector
    fn freqs<W: std::io::Write>(
        &self,
        w: &mut W,
        dir: &str,
        energies: &mut [f64],
        freq: &mut FreqParts,
    ) -> DVector<f64> {
        let (_, summary) = match freq {
            FreqParts::Sic {
                intder,
                taylor,
                taylor_disps,
                atomic_numbers,
            } => rust_pbqff::coord_type::sic::freqs(
                w,
                dir,
                energies,
                &mut intder.clone(),
                taylor,
                taylor_disps,
                atomic_numbers,
                STEP_SIZE,
            )
            .unwrap_or_default(),
            FreqParts::Cart {
                fcs,
                target_map,
                n,
                nfc2,
                nfc3,
                mol,
            } => {
                let (fc2, f3, f4) =
                    make_fcs(target_map, energies, fcs, *n, *nfc2, *nfc3, dir);
                rust_pbqff::coord_type::cart::freqs(dir, mol, fc2, f3, f4)
            }
        };
        if is_debug() {
            writeln!(w, "dir={dir}").unwrap();
            writeln!(w, "{}", summary).unwrap();
        }
        let freqs: Vec<f64> = summary
            .corrs
            .into_iter()
            .filter(|x| x.is_finite())
            .collect();
        let freqs = sort_irreps(&freqs, &summary.irreps);
        DVector::from(freqs)
    }

    /// helper method for computing the single-point energies for the jacobian
    #[allow(clippy::type_complexity)]
    fn jac_energies(
        &self,
        rows: usize,
        params: &Params,
        geoms: &[ProgramResult],
        molecules: &[config::Molecule],
    ) -> (Vec<Job<Mopac>>, Vec<Vec<FreqParts>>, Vec<Vec<usize>>) {
        let mut idx = 0;
        let mut jobs = Vec::new();
        let mut freqs = vec![vec![]; molecules.len()];
        // boundaries of each chunk for each molecule
        let mut indices = vec![];
        for (i, molecule) in molecules.iter().enumerate() {
            indices.push(vec![idx]);
            for row in 0..rows {
                let index = 2 * i * rows + 2 * row;
                let mut pf = params.clone();
                pf.values[row] += DELTA;
                // idx = job_num so use it twice
                let (freq, fwd_jobs) = self.build_jobs(
                    &mut output_stream(),
                    geoms[index].clone(),
                    &pf,
                    idx,
                    idx,
                    molecule,
                );
                idx += fwd_jobs.len();
                indices[i].push(idx);
                jobs.extend(fwd_jobs);
                freqs[i].push(freq);

                let mut pb = params.clone();
                pb.values[row] -= DELTA;
                let (freq, bwd_jobs) = self.build_jobs(
                    &mut output_stream(),
                    geoms[index + 1].clone(),
                    &pb,
                    idx,
                    idx,
                    molecule,
                );
                idx += bwd_jobs.len();
                indices[i].push(idx);
                jobs.extend(bwd_jobs);
                freqs[i].push(freq);
            }
        }
        (jobs, freqs, indices)
    }
}

impl Default for Frequency {
    fn default() -> Self {
        Self::new()
    }
}

impl Optimize for Frequency {
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q: Queue<Mopac> + std::marker::Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[config::Molecule],
    ) -> Option<DVector<f64>> {
        let mut w = output_stream();

        writeln!(w, "Params:\n{}", params).unwrap();

        let mut ret = Vec::new();

        for molecule in molecules {
            setup();
            let geom = optimize_geometry(
                molecule.geometry.clone(),
                params,
                submitter,
                "inp",
                "opt",
                molecule.charge,
                molecule.template.clone(),
            )?;

            writeln!(
                w,
                "Optimized Geometry:\n{:20.12}",
                Molecule::new(geom.cart_geom.clone().unwrap())
            )
            .unwrap();

            // have to call this before build_jobs since it writes params now
            setup();

            let (mut freq, mut jobs) =
                self.build_jobs(&mut w, geom, params, 0, 0, molecule);
            writeln!(
                w,
                "\n{} atoms require {} jobs",
                molecule.atom_names.len(),
                jobs.len()
            )
            .unwrap();

            let mut energies = vec![0.0; jobs.len()];
            // drain to get energies
            let status = submitter.drain("inp", &mut jobs, &mut energies);

            if status.is_err() {
                return None;
            }
            takedown();

            // convert energies to frequencies and return those
            let _ = std::fs::create_dir("freqs");
            let res = self.freqs(&mut w, "freqs", &mut energies, &mut freq);
            let _ = std::fs::remove_dir_all("freqs");
            ret.extend(res.iter());
        }
        Some(DVector::from(ret))
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned. `ntrue` is the expected ("true") size of
    /// the frequency vectors in case something goes wrong and a different
    /// number of entries comes out of freqs
    fn num_jac<Q: Queue<Mopac> + std::marker::Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[config::Molecule],
        ntrue: usize,
    ) -> na::DMatrix<f64> {
        let start = std::time::Instant::now();

        let rows = params.len();
        // build all the optimizations
        let mut opts = jac_opt(rows, params, molecules);

        setup();
        let mut geoms = vec![Default::default(); opts.len()];
        submitter
            .energize("inp", &mut opts, &mut geoms)
            .expect("numjac optimizations failed");

        eprintln!(
            "finished optimizations after {:.1} sec",
            start.elapsed().as_millis() as f64 / 1000.
        );
        let start = std::time::Instant::now();

        // jac_energies writes to tmparam so have to setup first
        setup();

        // build all the single-point energies
        let (mut jobs, mut freqs, indices) =
            self.jac_energies(rows, params, &geoms, molecules);

        eprintln!(
            "finished building {} energies after {:.1} sec",
            jobs.len(),
            start.elapsed().as_millis() as f64 / 1000.
        );
        let start = std::time::Instant::now();

        // freqs is nmolecules long, so check that set of freqs is equal to
        // 2*nparam since there is a forward and back for each
        assert!(freqs.iter().all(|x| x.len() == 2 * params.len()));

        // run all of the energies
        let mut energies = vec![0.0; jobs.len()];
        submitter
            .drain("inp", &mut jobs, &mut energies)
            .expect("numjac optimizations failed");
        takedown();

        // reverse the freqs so I can pop instead of cloning out of it
        freqs.reverse();

        eprintln!(
            "finished running energies after {:.1} sec",
            start.elapsed().as_millis() as f64 / 1000.
        );
        let start = std::time::Instant::now();

        use rayon::prelude::*;

        // calculate all of the frequencies and assemble the jacobian
        let mut jacs = Vec::new();

        for m in 0..molecules.len() {
            let mut energy_chunks = Vec::new();
            assert_eq!(indices[m].len(), 2 * params.len() + 1);
            // -1 because I'm going to refer to idx+1 in the loop
            for idx in 0..indices[m].len() - 1 {
                energy_chunks.push(
                    energies[indices[m][idx]..indices[m][idx + 1]].to_vec(),
                );
            }
            // zip it with the FreqParts
            let pairs: Vec<_> = energy_chunks
                .iter()
                .cloned()
                .zip(freqs.pop().unwrap())
                .collect();
            // iterate over the energy/freqparts pairs in parallel
            let freqs: Vec<_> = pairs
                .par_iter()
                .enumerate()
                .map(|(i, (energy, freq))| {
                    let mut energy = energy.clone();
                    let mut freq = freq.clone();
                    let dir = format!("freqs{}_{}", i, m);
                    let _ = std::fs::create_dir(&dir);
                    self.freqs(
                        &mut output_stream(),
                        &dir,
                        &mut energy,
                        &mut freq,
                    )
                })
                .collect();

            // should be a forward and backward step for each parameter
            assert_eq!(freqs.len(), 2 * params.len());

            for i in 0..pairs.len() {
                let dir = format!("freqs{}_{}", i, m);
                let _ = std::fs::remove_dir_all(&dir);
            }
            let jac_t: Vec<_> = freqs
                .chunks(2)
                .map(|pair| {
                    // I hate these clones but I can't figure out how to borrow
                    // them to avoid it. I think chunks is the root of the issue
                    // because it gives me references to freqs when I really
                    // just want to consume freqs
                    let mut fwd = pair[0].clone();
                    let mut bwd = pair[1].clone();
                    let (fl, bl) = (fwd.len(), bwd.len());
                    if fl < bl {
                        bwd.resize_vertically_mut(fl, 0.0);
                    } else if fl > bl {
                        fwd.resize_vertically_mut(bl, 0.0);
                    }
                    let mut ret = (fwd - bwd) / (2. * DELTA);
                    if ret.len() != ntrue {
                        ret.resize_vertically_mut(ntrue, 0.0);
                    }
                    ret
                })
                .collect();
            assert_eq!(jac_t.len(), params.len());
            jacs.push(jac_t);
        }

        // join the columns from different molecules

        let mut cols = Vec::new();
        // the length of each jacobian should be the same for each molecule
        // since it is determined by the number of parameters
        for i in 0..jacs[0].len() {
            let mut tmp = Vec::new();
            // for each molecule
            for jac in &jacs {
                assert_eq!(params.len(), jac.len());
                tmp.extend(jac[i].iter());
            }
            // if the whole column is zeros, make the first element 1 to avoid
            // linear algebra issues
            if tmp.iter().all(|s| *s == 0.0) {
                eprintln!("singular column in jacobian, fixing");
                tmp[0] = 1.0;
            }
            cols.push(DVector::from(tmp));
        }

        eprintln!(
            "finished frequencies after {:.1} sec",
            start.elapsed().as_millis() as f64 / 1000.
        );

        na::DMatrix::from_columns(&cols)
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

/// helper function for generating the optimizations for the jacobian. returns a
/// vec of Jobs and also a vec of indices separating each molecule
fn jac_opt(
    rows: usize,
    params: &Params,
    molecules: &[config::Molecule],
) -> Vec<Job<Mopac>> {
    let mut opts = Vec::new();
    // each row corresponds to one parameter
    for (i, molecule) in molecules.iter().enumerate() {
        for row in 0..rows {
            let index = 2 * i * rows + 2 * row;
            // forward
            {
                let mut pf = params.clone();
                pf.values[row] += DELTA;
                opts.push(Job::new(
                    Mopac::new_full(
                        format!("inp/opt{row}_fwd{i}"),
                        Some(pf),
                        molecule.geometry.clone(),
                        molecule.charge,
                        molecule.template.clone(),
                    ),
                    // this is the index because I have rows entries for each
                    // molecule. for the first molecule i = 0 and this reduces
                    // to the original formula of 2*row
                    index,
                ));
            }

            // backward
            {
                let mut pb = params.clone();
                pb.values[row] -= DELTA;
                opts.push(Job::new(
                    Mopac::new_full(
                        format!("inp/opt{row}_bwd{i}"),
                        Some(pb),
                        molecule.geometry.clone(),
                        molecule.charge,
                        molecule.template.clone(),
                    ),
                    index + 1,
                ));
            }
        }
    }
    opts
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
