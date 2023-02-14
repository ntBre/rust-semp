use na::DVector;
use psqs::geom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::program::{Job, ProgramResult, Template};
use psqs::queue::Queue;
use rust_pbqff::coord_type::cart::FirstPart;
use rust_pbqff::coord_type::findiff::bighash::BigHash;
use rust_pbqff::coord_type::findiff::FiniteDifference;
use rust_pbqff::coord_type::fitting::Fitted;
use rust_pbqff::coord_type::normal::{
    fc3_index, fc4_index, to_qcm, F3qcm, F4qcm, Fc, Normal,
};
use rust_pbqff::coord_type::sic::IntderError;
use rust_pbqff::coord_type::{Cart, Derivative, Sic};
use rust_pbqff::{Output, Spectro};
use symm::{Molecule, PointGroup};
use taylor::{Disps, Taylor};

use crate::config::CoordType;
use crate::utils::sort_irreps;
use crate::BAD_FLOAT;

use crate::{config, utils::setup, utils::takedown};
use nalgebra as na;
use std::cmp::Ordering;
use std::fs::File;
use std::io::Write;
use std::marker::Sync;
use std::rc::Rc;
use std::sync::Mutex;

use super::Optimize;

lazy_static::lazy_static! {
    static ref DEBUG: bool = std::env::var("SEMP_FREQ_DEBUG").is_ok();
}

type Dvec = DVector<f64>;

/// pbqff step size
const STEP_SIZE: f64 = 0.005;

fn output_stream() -> Box<dyn Write> {
    if *DEBUG {
        Box::new(std::io::stderr())
    } else {
        Box::new(std::io::sink())
    }
}

pub struct Frequency {
    train: DVector<f64>,
    cur: Option<DVector<f64>>,
    delta: f64,
    logger: Mutex<File>,
}

pub fn optimize_geometry<Q: Queue<Mopac> + Sync>(
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
            format!("{dir}/{name}"),
            Some(params.clone()),
            geom,
            charge,
            template,
        ),
        0,
    );
    let mut res = vec![Default::default(); 1];
    let status = queue.energize(dir, vec![opt], &mut res);
    if status.is_err() {
        return None;
    }
    Some(res.pop().unwrap())
}

/// the state generated in build_jobs needed to finish running frequencies after
/// the energies are run
#[derive(Clone, Debug)]
#[allow(clippy::large_enum_variant)]
enum FreqParts {
    Sic {
        intder: rust_pbqff::Intder,
        taylor: Taylor,
        taylor_disps: Disps,
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
    Norm {
        normal: Normal,
        taylor: Taylor,
        taylor_disps: Disps,
        spectro: Spectro,
        output: Output,
    },
    NormHarm {
        fcs: Vec<f64>,
        target_map: BigHash,
        n: usize,
        nfc2: usize,
        nfc3: usize,
        mol: Molecule,
        pg: PointGroup,
    },
}

impl FreqParts {
    fn sic(
        intder: rust_pbqff::Intder,
        taylor: Taylor,
        taylor_disps: Disps,
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

    fn norm(
        normal: Normal,
        taylor: Taylor,
        taylor_disps: Disps,
        spectro: Spectro,
        output: Output,
    ) -> Self {
        Self::Norm {
            normal,
            taylor,
            taylor_disps,
            spectro,
            output,
        }
    }
}

#[allow(clippy::large_enum_variant)]
#[derive(Clone)]
enum Builder {
    None,
    Norm {
        norm: Normal,
        s: Spectro,
        o: Output,
        pg: PointGroup,
    },
}

impl Frequency {
    pub fn new(train: DVector<f64>, delta: f64) -> Self {
        let logger = Mutex::new(
            std::fs::File::create("freqs.log")
                .expect("failed to create 'freqs.log'"),
        );
        Self {
            train,
            delta,
            logger,
            cur: None,
        }
    }

    /// build jobs for a fixed set of Params. instead of passing the params as
    /// part of each separate Mopac, write the parameters once and update the
    /// Template passed to Mopac::new. returns FreqParts, which contains
    /// whatever this molecule's coordinate type requires to finish running
    /// frequencies, and the list of jobs to run
    #[allow(clippy::too_many_arguments)]
    fn build_jobs<W>(
        &self,
        w: &mut W,
        geom: ProgramResult,
        params: &Params,
        start_index: usize,
        job_num: usize,
        molecule: &config::Molecule,
        coord_type: &CoordType,
        builder: Builder,
    ) -> Result<(FreqParts, Vec<Job<Mopac>>), IntderError>
    where
        W: Write,
    {
        let mol = {
            let mut mol =
                Molecule::new(geom.cart_geom.expect(
                    "no geom found, try adding GEO-OK to input template",
                ));
            mol.normalize();
            mol
        };
        const SYMM_EPS: f64 = 1e-6;
        let pg = mol.point_group_approx(SYMM_EPS);

        writeln!(w, "Normalized Geometry:\n{mol:20.12}").unwrap();
        writeln!(w, "Point Group = {pg}").unwrap();

        let tmpl = write_params(job_num, params, molecule.template.clone());

        match &coord_type {
            CoordType::Sic(intder) => {
                let intder = intder.clone();
                // NOTE: assuming that the intder coordinates max out at C2v,
                // this was the case for ethylene
                let pg = if pg.is_d2h() {
                    pg.subgroup(symm::Pg::C2v).unwrap()
                } else {
                    pg
                };

                let mut sic = Sic::new(intder);
                let (moles, taylor, taylor_disps, atomic_numbers) =
                    sic.generate_pts(w, &mol, &pg, STEP_SIZE)?;

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
                Ok((
                    FreqParts::sic(
                        sic.intder,
                        taylor,
                        taylor_disps,
                        atomic_numbers,
                    ),
                    jobs,
                ))
            }
            CoordType::Cart => {
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
                    Derivative::Quartic(nfc2, nfc3, nfc4),
                    &mut fcs,
                    &mut target_map,
                    n,
                );
                let dir = "inp";
                let mut job_num = start_index;
                let mut jobs = Vec::new();
                for mol in geoms {
                    let filename = format!("{dir}/job.{job_num:08}");
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
                Ok((FreqParts::cart(fcs, target_map, n, nfc2, nfc3, mol), jobs))
            }
            CoordType::Normal => {
                // there must be a way to tie these types together more smoothly
                let Builder::Norm { mut norm, s, o, pg } = builder else {
		    unreachable!()
		};
                assert!(
                    !norm.findiff,
                    "only fitted normal coordinates are currently supported"
                );
                // adapted from Normal::run_fitted in pbqff
                norm.prep_qff(w, &o);
                let (geoms, taylor, taylor_disps, _atomic_numbers) =
                    norm.generate_pts(w, &o.geom, &pg, STEP_SIZE).unwrap();
                let dir = "inp";
                let jobs = Mopac::build_jobs(
                    &geoms,
                    None,
                    dir,
                    start_index,
                    1.0,
                    job_num,
                    molecule.charge,
                    tmpl,
                );
                Ok((FreqParts::norm(norm, taylor, taylor_disps, s, o), jobs))
            }
            CoordType::NormalHarm => {
                // kinda sloppy but copied from Cart case with just Derivative
                // changed
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
                    Derivative::Harmonic(nfc2),
                    &mut fcs,
                    &mut target_map,
                    n,
                );
                let dir = "inp";
                let mut job_num = start_index;
                let mut jobs = Vec::new();
                for mol in geoms {
                    let filename = format!("{dir}/job.{job_num:08}");
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
                Ok((
                    FreqParts::NormHarm {
                        fcs,
                        target_map,
                        n,
                        nfc2,
                        nfc3,
                        mol,
                        pg,
                    },
                    jobs,
                ))
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
            } => Sic::new(intder.clone())
                .freqs(
                    w,
                    dir,
                    energies,
                    taylor,
                    taylor_disps,
                    atomic_numbers,
                    STEP_SIZE,
                )
                .unwrap_or_else(|_| {
                    eprintln!("SIC freqs failed");
                    Default::default()
                }),
            FreqParts::Cart {
                fcs,
                target_map,
                n,
                nfc2,
                nfc3,
                mol,
            } => {
                let (fc2, f3, f4) = Cart.make_fcs(
                    target_map,
                    energies,
                    fcs,
                    *n,
                    rust_pbqff::coord_type::cart::Derivative::Quartic(
                        *nfc2, *nfc3, 0,
                    ),
                    dir,
                );
                rust_pbqff::coord_type::cart::freqs(dir, mol, fc2, f3, f4)
            }
            FreqParts::Norm {
                normal,
                taylor,
                taylor_disps,
                spectro,
                output,
            } => {
                let (fcs, _) = normal
                    .anpass(None, energies, taylor, taylor_disps, STEP_SIZE, w)
                    .unwrap();
                // needed in case taylor eliminated some of the higher
                // derivatives by symmetry. this should give the maximum, full
                // sizes without resizing
                let n = normal.ncoords;
                let mut f3qcm = vec![0.0; fc3_index(n, n, n) + 1];
                let mut f4qcm = vec![0.0; fc4_index(n, n, n, n) + 1];
                for Fc(i, j, k, l, val) in fcs {
                    // adapted from Intder::add_fc
                    match (k, l) {
                        (0, 0) => {
                            // harmonic, skip it
                        }
                        (_, 0) => {
                            let idx = fc3_index(i, j, k);
                            f3qcm[idx] = val;
                        }
                        (_, _) => {
                            let idx = fc4_index(i, j, k, l);
                            f4qcm[idx] = val;
                        }
                    }
                }
                let (f3, f4) =
                    to_qcm(&output.harms, normal.ncoords, &f3qcm, &f4qcm, 1.0);
                let (o, _) = spectro.finish(
                    DVector::from(output.harms.clone()),
                    F3qcm::new(f3),
                    F4qcm::new(f4),
                    output.irreps.clone(),
                    normal.lxm.as_ref().unwrap().clone(),
                    normal.lx.as_ref().unwrap().clone(),
                );
                (std::mem::take(spectro), o)
            }
            FreqParts::NormHarm { .. } => unimplemented!(),
        };
        if *DEBUG {
            writeln!(w, "dir={dir}").unwrap();
            writeln!(w, "{:?}", summary.corrs).unwrap();
        }
        let s = summary.corrs.len();
        let freqs: Vec<f64> = summary
            .corrs
            .into_iter()
            .filter(|x| x.is_finite())
            .collect();
        let e = freqs.len();
        if s != e {
            eprintln!("filtered out {} non-finite corrs from", s - e);
        }
        let freqs = sort_irreps(&freqs, &summary.irreps);
        DVector::from(freqs)
    }

    /// generate all of the HFF energies for finding normal coordinates
    fn harm_jac_energies<Q: Queue<Mopac> + Sync>(
        &self,
        rows: usize,
        params: &Params,
        geoms: &[ProgramResult],
        molecules: &[config::Molecule],
        submitter: &Q,
    ) -> Vec<Builder> {
        let mut idx = 0;
        let mut jobs = Vec::new();
        let mut freqs = vec![vec![]; molecules.len()];
        // boundaries of each chunk for each molecule
        let mut indices = vec![];
        for (i, molecule) in molecules.iter().enumerate() {
            // push some junk in for non-normals. it'll be skipped later
            indices.push(vec![idx]);
            if !molecule.coord_type.is_normal() {
                continue;
            }
            // at this point we know the coord_type is normal
            for row in 0..rows {
                let index = 2 * i * rows + 2 * row;
                let mut pf = params.clone();
                pf.values[row] += self.delta;
                // idx = job_num so use it twice
                let (freq, fwd_jobs) = self
                    .build_jobs(
                        &mut output_stream(),
                        geoms[index].clone(),
                        &pf,
                        idx,
                        idx,
                        molecule,
                        &CoordType::NormalHarm,
                        Builder::None,
                    )
                    .unwrap();
                // we unwrap here because it's not entirely clear how to recover
                // from one half of one column going off the rails.
                idx += fwd_jobs.len();
                indices[i].push(idx);
                jobs.extend(fwd_jobs);
                freqs[i].push(freq);

                let mut pb = params.clone();
                pb.values[row] -= self.delta;
                let (freq, bwd_jobs) = self
                    .build_jobs(
                        &mut output_stream(),
                        geoms[index + 1].clone(),
                        &pb,
                        idx,
                        idx,
                        molecule,
                        &CoordType::NormalHarm,
                        Builder::None,
                    )
                    .unwrap();
                idx += bwd_jobs.len();
                indices[i].push(idx);
                jobs.extend(bwd_jobs);
                freqs[i].push(freq);
            }
        }

        // at this point, all of the jobs are built for running all of the HFFs,
        // but only the builders for non-normal molecules are built.
        let mut energies = vec![0.0; jobs.len()];
        submitter
            .drain("inp", jobs, &mut energies, 0)
            .expect("numjac optimizations failed");
        // reset for rest of num_jac
        setup();

        use rayon::prelude::*;
        freqs.reverse(); // for popping
        let mut builders = Vec::new();
        for m in 0..molecules.len() {
            if !molecules[m].coord_type.is_normal() {
                builders.extend(vec![Builder::None; 2 * params.len()]);
                continue;
            }
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
            let bs: Vec<_> = pairs
                .par_iter()
                .enumerate()
                .map(|(i, (energy, freq))| {
                    let freq = freq.clone();
                    let dir = format!("freqs{i}_{m}");
                    let _ = std::fs::create_dir(&dir);
                    let norm = Normal::findiff(false);
                    let FreqParts::NormHarm {
                        mut fcs,
                        mut target_map,
                        n,
                        nfc2,
                        nfc3,
                        mol,
                        pg,
                    } = freq else { unimplemented!() } ;
                    let (fc2, _, _) = norm.make_fcs(
                        &mut target_map,
                        energy,
                        &mut fcs,
                        n,
                        Derivative::Quartic(nfc2, nfc3, 0),
                        &dir,
                    );
                    let (s, o) = norm.harm_freqs(&dir, &mol, fc2);
                    Builder::Norm { norm, s, o, pg }
                })
                .collect();
            builders.extend(bs)
        }
        builders
    }

    /// helper method for building the single-point energy jobs for the
    /// jacobian. `builders` should be of length 2 * molecules.len() *
    /// params.len()
    #[allow(clippy::type_complexity)]
    fn jac_energies(
        &self,
        rows: usize,
        params: &Params,
        geoms: &[ProgramResult],
        molecules: &[config::Molecule],
        mut builders: Vec<Builder>,
    ) -> (Vec<Job<Mopac>>, Vec<Vec<FreqParts>>, Vec<Vec<usize>>) {
        let mut idx = 0;
        let mut jobs = Vec::new();
        let mut freqs = vec![vec![]; molecules.len()];
        // boundaries of each chunk for each molecule
        let mut indices = vec![];
        // reverse for popping
        builders.reverse();
        for (i, molecule) in molecules.iter().enumerate() {
            indices.push(vec![idx]);
            for row in 0..rows {
                let index = 2 * i * rows + 2 * row;
                let mut pf = params.clone();
                pf.values[row] += self.delta;
                // idx = job_num so use it twice
                let (freq, fwd_jobs) = self
                    .build_jobs(
                        &mut output_stream(),
                        geoms[index].clone(),
                        &pf,
                        idx,
                        idx,
                        molecule,
                        &molecule.coord_type,
                        builders.pop().unwrap(),
                    )
                    .unwrap();
                // we unwrap here because it's not entirely clear how to recover
                // from one half of one column going off the rails.
                idx += fwd_jobs.len();
                indices[i].push(idx);
                jobs.extend(fwd_jobs);
                freqs[i].push(freq);

                let mut pb = params.clone();
                pb.values[row] -= self.delta;
                let (freq, bwd_jobs) = self
                    .build_jobs(
                        &mut output_stream(),
                        geoms[index + 1].clone(),
                        &pb,
                        idx,
                        idx,
                        molecule,
                        &molecule.coord_type,
                        builders.pop().unwrap(),
                    )
                    .unwrap();
                idx += bwd_jobs.len();
                indices[i].push(idx);
                jobs.extend(bwd_jobs);
                freqs[i].push(freq);
            }
        }
        (jobs, freqs, indices)
    }

    fn rel_diff(&self, mut v: Dvec) -> Dvec {
        for (i, f) in v.iter_mut().enumerate() {
            let t = self.train.get(i).unwrap_or(&1.0);
            *f = 100.0 * (*f - t) / t;
        }
        if v.len() != self.train.len() {
            v = v.resize_vertically(self.train.len(), BAD_FLOAT);
        }
        v
    }
}

/// generate a parameter filename using `job_num`, write the parameters to that
/// file, and return the template with the parameter filename added
fn write_params(
    job_num: usize,
    params: &Params,
    template: Template,
) -> Template {
    let param_file = Rc::new(format!("tmparam/{job_num}.dat"));
    Mopac::write_params(params, &param_file);
    let mut tmpl = template;
    use std::fmt::Write;
    write!(tmpl.header, " external={param_file}").unwrap();
    tmpl
}

impl Default for Frequency {
    fn default() -> Self {
        Self::new(DVector::zeros(0), 1e-4)
    }
}

impl Optimize for Frequency {
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q: Queue<Mopac> + Sync>(
        &mut self,
        params: &Params,
        submitter: &Q,
        molecules: &[config::Molecule],
    ) -> Option<DVector<f64>> {
        let mut w = output_stream();

        writeln!(w, "Params:\n{params}").unwrap();

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

            let builder = if molecule.coord_type.is_normal() {
                // NOTE using fitted normal coordinates
                let norm = Normal::findiff(false);
                // safe to use 0 as job_num because we have to reset directories
                // after the harmonic part anyway
                let tmpl = write_params(0, params, molecule.template.clone());
                let (s, o, _ref_energy, pg, _) = norm.cart_part(
                    &FirstPart {
                        template: tmpl.header,
                        optimize: false,
                        geometry: Geom::Xyz(
                            geom.cart_geom.as_ref().unwrap().clone(),
                        ),
                        charge: molecule.charge,
                        step_size: STEP_SIZE,
                    },
                    submitter,
                    &mut w,
                    "inp",
                );
                setup();
                Builder::Norm { norm, s, o, pg }
            } else {
                Builder::None
            };

            let Ok((mut freq, jobs)) =
		self.build_jobs(&mut w, geom, params, 0, 0, molecule, &molecule.coord_type, builder)
	    else {
		return None
	    };

            writeln!(
                w,
                "\n{} atoms require {} jobs",
                molecule.atom_names.len(),
                jobs.len()
            )
            .unwrap();

            let mut energies = vec![0.0; jobs.len()];
            // drain to get energies
            let status = submitter.drain("inp", jobs, &mut energies, 0);

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
        let ret = DVector::from(ret);
        self.cur = Some(ret.clone());
        Some(self.rel_diff(ret))
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned. `ntrue` is the expected ("true") size of
    /// the frequency vectors in case something goes wrong and a different
    /// number of entries comes out of freqs
    fn num_jac<Q: Queue<Mopac> + Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[config::Molecule],
        ntrue: usize,
    ) -> na::DMatrix<f64> {
        let start = std::time::Instant::now();

        let rows = params.len();
        // build all the optimizations
        let opts = self.jac_opt(rows, params, molecules);

        setup();
        let mut geoms = vec![Default::default(); opts.len()];
        submitter
            .energize("inp", opts, &mut geoms)
            .expect("numjac optimizations failed");

        eprintln!(
            "finished optimizations after {:.1} sec",
            start.elapsed().as_millis() as f64 / 1000.
        );
        let start = std::time::Instant::now();

        // jac_energies writes to tmparam so have to setup first
        setup();

        let builders = if molecules.iter().any(|m| m.coord_type.is_normal()) {
            self.harm_jac_energies(rows, params, &geoms, molecules, submitter)
        } else {
            vec![Builder::None; 2 * molecules.len() * params.len()]
        };

        // build all the single-point energies
        let (jobs, mut freqs, indices) =
            self.jac_energies(rows, params, &geoms, molecules, builders);

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
            .drain("inp", jobs, &mut energies, 0)
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
                    let dir = format!("freqs{i}_{m}");
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
                let dir = format!("freqs{i}_{m}");
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
                    if fl != bl {
                        eprintln!("fl ({fl}) != bl ({bl}), resizing");
                    }
                    match fl.cmp(&bl) {
                        Ordering::Less => bwd.resize_vertically_mut(fl, 0.0),
                        Ordering::Equal => {}
                        Ordering::Greater => fwd.resize_vertically_mut(bl, 0.0),
                    }
                    let fwd = self.rel_diff(fwd);
                    let bwd = self.rel_diff(bwd);
                    (fwd - bwd) / (2. * self.delta)
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
            if tmp.len() != ntrue {
                tmp.resize(ntrue, 0.0);
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

    fn log(&self, iter: usize, _got: &DVector<f64>, _want: &DVector<f64>) {
        let mut logger = self.logger.lock().unwrap();
        write!(logger, "{iter:5}").unwrap();
        let got = self.cur.as_ref().unwrap();
        for g in got {
            write!(logger, "{g:8.1}").unwrap();
        }
        writeln!(logger, "{:8.1}", mae(got, &self.train)).unwrap();
    }
}

impl Frequency {
    /// helper function for generating the optimizations for the jacobian. returns a
    /// vec of Jobs and also a vec of indices separating each molecule
    fn jac_opt(
        &self,
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
                    pf.values[row] += self.delta;
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
                    pb.values[row] -= self.delta;
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
}

/// return the mean absolute error (1/N ∑ᵢᴺ |aᵢ - bᵢ|) where N is the shorter of
/// `a.len()` and `b.len()` to guard against panics
fn mae(a: &DVector<f64>, b: &DVector<f64>) -> f64 {
    let al = a.len();
    let bl = b.len();
    if al != bl {
        eprintln!("len mismatch in mae: {al} vs {bl}");
    }
    let l = al.min(bl);
    let mut sum = 0.0;
    for i in 0..l {
        sum += f64::abs(a[i] - b[i]);
    }
    sum / l as f64
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
