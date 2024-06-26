use super::Optimize;
use crate::{
    config::{self, CoordType},
    driver::{Driver, Params},
    utils::{mae, setup, takedown},
};
use na::DVector;
use nalgebra as na;
use pbqff::{
    coord_type::{
        cart::{self, Derivative, FirstPart},
        findiff::{
            bighash::{BigHash, Target},
            FiniteDifference,
        },
        fitted::Fitted,
        normal::{to_qcm, CartPart, F3qcm, F4qcm, Normal},
        Cart, Sic,
    },
    Output, Spectro,
};
use psqs::{
    geom::Geom,
    program::{Job, ProgramResult, Template},
    queue::{Check, Queue},
};
use std::{
    cmp::Ordering,
    error::Error,
    fmt::Display,
    fs::File,
    io::Write,
    path::Path,
    sync::{LazyLock, Mutex},
};
use symm::{Irrep, Molecule, PointGroup};
use taylor::Taylor;

#[cfg(test)]
mod tests;

static DEBUG: LazyLock<bool> =
    LazyLock::new(|| std::env::var("SEMP_FREQ_DEBUG").is_ok());

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
    delta: f64,
    logger: Mutex<File>,
    sort_fn: fn(&[f64], &[Irrep]) -> Vec<f64>,
}

/// the shared components of a finite difference FreqPart
#[derive(Clone, Debug)]
struct FinDiff {
    fcs: Vec<f64>,
    targets: Vec<Target>,
    n: usize,
    nfc2: usize,
    nfc3: usize,
}

/// the state generated in build_jobs needed to finish running frequencies after
/// the energies are run
#[derive(Clone, Debug)]
enum FreqParts {
    Sic {
        intder: pbqff::Intder,
        taylor: Taylor,
        atomic_numbers: Vec<usize>,
    },
    Cart {
        inner: FinDiff,
        mol: Molecule,
    },
    FitNorm {
        normal: Normal,
        taylor: Taylor,
        spectro: Spectro,
        output: Output,
    },
    FinNorm {
        normal: Normal,
        spectro: Spectro,
        output: Output,
        inner: FinDiff,
    },
    NormHarm {
        inner: FinDiff,
        mol: Molecule,
        pg: PointGroup,
    },
}

#[derive(Debug)]
struct GeomError;

impl Display for GeomError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl Error for GeomError {}

impl FreqParts {
    fn sic(
        intder: pbqff::Intder,
        taylor: Taylor,
        atomic_numbers: Vec<usize>,
    ) -> Self {
        Self::Sic {
            intder,
            taylor,
            atomic_numbers,
        }
    }

    fn norm(
        normal: Normal,
        taylor: Taylor,
        spectro: Spectro,
        output: Output,
    ) -> Self {
        Self::FitNorm {
            normal,
            taylor,
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
        ref_energy: Option<f64>,
    },
}

type Jobs<D> = Vec<Job<D>>;

impl Frequency {
    pub fn new(delta: f64, sort_ascending: bool) -> Self {
        let logger = Mutex::new(
            std::fs::File::create("freqs.log")
                .expect("failed to create 'freqs.log'"),
        );
        let sort_fn = if sort_ascending {
            crate::utils::sort_ascending
        } else {
            crate::utils::sort_irreps
        };
        Self {
            delta,
            logger,
            sort_fn,
        }
    }

    /// build jobs for a fixed set of Params. instead of passing the params as
    /// part of each separate Mopac, write the parameters once and update the
    /// Template passed to Mopac::new. returns FreqParts, which contains
    /// whatever this molecule's coordinate type requires to finish running
    /// frequencies, and the list of jobs to run. At first glance, it appears
    /// redundant to pass both `molecule` and `coord_type` because each
    /// [config::Molecule] contains its own coordinate type, but this is needed
    /// to call `build_jobs` for the harmonic portion of a normal QFF.
    #[allow(clippy::too_many_arguments)]
    fn build_jobs<W, D>(
        &self,
        w: &mut W,
        geom: ProgramResult,
        params: &D::Params,
        start_index: usize,
        job_num: usize,
        molecule: &config::Molecule,
        coord_type: &CoordType,
        builder: Builder,
    ) -> Result<(FreqParts, Jobs<D>), Box<dyn Error>>
    where
        W: Write,
        D: Driver,
    {
        let ProgramResult {
            energy,
            cart_geom: Some(g),
            ..
        } = geom
        else {
            return Err(GeomError.into());
        };

        let mol = {
            let mut mol = Molecule::new(g);
            mol.normalize();
            mol
        };
        const SYMM_EPS: f64 = 1e-6;
        let pg = mol.point_group_approx(SYMM_EPS);

        writeln!(w, "Normalized Geometry:\n{mol:20.12}").unwrap();
        writeln!(w, "Point Group = {pg}").unwrap();

        let tmpl = D::write_params(job_num, params, molecule.template.clone());

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
                let (moles, taylor, atomic_numbers) =
                    sic.generate_pts(".", w, &mol, &pg, STEP_SIZE)?;

                // dir created in generate_pts but unused here
                let _ = std::fs::remove_dir_all("pts");

                // call build_jobs like before
                let jobs = <D as Driver>::build_jobs(
                    moles,
                    "inp",
                    start_index,
                    1.0,
                    job_num,
                    molecule.charge,
                    tmpl,
                );
                Ok((FreqParts::sic(sic.intder, taylor, atomic_numbers), jobs))
            }
            CoordType::Cart => {
                let (inner, jobs) = build_findiff(
                    3 * mol.atoms.len(),
                    &mol,
                    pg,
                    energy,
                    start_index,
                    molecule,
                    &tmpl,
                    &Cart,
                    false,
                );
                Ok((FreqParts::Cart { inner, mol }, jobs))
            }
            CoordType::Normal(..) => {
                // there must be a way to tie these types together more smoothly
                let Builder::Norm {
                    mut norm,
                    s,
                    o,
                    pg,
                    ref_energy,
                } = builder
                else {
                    unreachable!()
                };
                norm.prep_qff(w, &o, pg);
                if !norm.findiff {
                    // adapted from Normal::run_fitted in pbqff
                    let pg = if pg.is_d2h() {
                        pg.subgroup(symm::Pg::C2v).unwrap()
                    } else {
                        pg
                    };
                    let dir = "inp";
                    let (geoms, taylor, _atomic_numbers) = norm
                        .generate_pts(dir, w, &o.geom, &pg, STEP_SIZE)
                        .unwrap();
                    let jobs = <D as Driver>::build_jobs(
                        geoms,
                        dir,
                        start_index,
                        1.0,
                        job_num,
                        molecule.charge,
                        tmpl,
                    );
                    Ok((FreqParts::norm(norm, taylor, s, o), jobs))
                } else {
                    let (inner, jobs) = build_findiff(
                        norm.ncoords,
                        &o.geom,
                        pg,
                        ref_energy.unwrap(),
                        start_index,
                        molecule,
                        &tmpl,
                        &norm,
                        false,
                    );
                    Ok((
                        FreqParts::FinNorm {
                            normal: norm,
                            spectro: s,
                            output: o,
                            inner,
                        },
                        jobs,
                    ))
                }
            }
            CoordType::NormalHarm => {
                let (inner, jobs) = build_findiff(
                    3 * mol.atoms.len(),
                    &mol,
                    pg,
                    energy,
                    start_index,
                    molecule,
                    &tmpl,
                    &Cart,
                    true,
                );
                Ok((FreqParts::NormHarm { inner, mol, pg }, jobs))
            }
        }
    }

    /// call `pbqff::coord_type::freqs`, but sort the frequencies by irrep
    /// and then frequency and return the result as a DVector
    fn freqs<W: std::io::Write>(
        &self,
        w: &mut W,
        energies: &mut [f64],
        freq: FreqParts,
    ) -> DVector<f64> {
        let (_, summary) = match freq {
            FreqParts::Sic {
                intder,
                taylor,
                atomic_numbers,
            } => Sic::new(intder)
                .freqs(
                    w,
                    None::<&Path>,
                    energies,
                    &taylor,
                    &atomic_numbers,
                    STEP_SIZE,
                )
                .unwrap_or_else(|_| {
                    eprintln!("SIC freqs failed");
                    Default::default()
                }),
            FreqParts::Cart {
                mol,
                inner:
                    FinDiff {
                        mut fcs,
                        targets,
                        n,
                        nfc2,
                        nfc3,
                    },
            } => {
                let (fc2, f3, f4) = Cart.make_fcs(
                    targets,
                    energies,
                    &mut fcs,
                    n,
                    Derivative::Quartic(nfc2, nfc3, 0),
                    None::<&Path>,
                );
                cart::freqs(None::<&Path>, &mol, fc2, f3, f4)
            }
            FreqParts::FitNorm {
                normal,
                taylor,
                spectro,
                output,
            } => {
                let (f3, f4) = normal.fit_freqs(
                    None::<&Path>,
                    energies,
                    taylor,
                    STEP_SIZE,
                    w,
                    &output,
                );
                let o = spectro.finish(
                    DVector::from(output.harms.clone()),
                    F3qcm::new(f3),
                    F4qcm::new(f4),
                    output.irreps,
                    normal.lxm.as_ref().unwrap().clone(),
                    normal.lx.as_ref().unwrap().clone(),
                );
                (spectro, o)
            }
            FreqParts::FinNorm {
                normal,
                spectro,
                output,
                inner:
                    FinDiff {
                        mut fcs,
                        targets,
                        n,
                        nfc2,
                        nfc3,
                    },
            } => {
                normal.map_energies(targets, energies, &mut fcs);
                let cubs = &fcs[nfc2..nfc2 + nfc3];
                let quarts = &fcs[nfc2 + nfc3..];
                let (f3, f4) =
                    to_qcm(&output.harms, n, cubs, quarts, intder::HART);
                let o = spectro.finish(
                    DVector::from(output.harms.clone()),
                    F3qcm::new(f3),
                    F4qcm::new(f4),
                    output.irreps,
                    normal.lxm.as_ref().unwrap().clone(),
                    normal.lx.as_ref().unwrap().clone(),
                );
                (spectro, o)
            }
            FreqParts::NormHarm { .. } => unimplemented!(),
        };
        if *DEBUG {
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
        let freqs = (self.sort_fn)(&freqs, &summary.irreps);
        DVector::from(freqs)
    }

    /// generate all of the HFF energies for finding normal coordinates
    fn harm_jac_energies<Q: Queue<D> + Sync, D: Driver>(
        &self,
        rows: usize,
        params: &D::Params,
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
                pf.incr_value(row, self.delta);
                let mut pb = params.clone();
                pb.decr_value(row, self.delta);

                for (p, pf) in [pf, pb].iter().enumerate() {
                    // idx = job_num so use it twice
                    let (freq, fwd_jobs) = self
                        .build_jobs(
                            &mut output_stream(),
                            geoms[index + p].clone(),
                            pf,
                            idx,
                            idx,
                            molecule,
                            &CoordType::NormalHarm,
                            Builder::None,
                        )
                        .unwrap_or_else(|e| {
                            // we exit here because it's not entirely clear how
                            // to recover from one half of one column going off
                            // the rails.
                            eprintln!(
				"failed to build jobs with {e} at index {index}"
			    );
                            eprintln!("params:\n{pf}");
                            std::process::exit(1);
                        });
                    idx += fwd_jobs.len();
                    indices[i].push(idx);
                    jobs.extend(fwd_jobs);
                    freqs[i].push(freq);
                }
            }
        }

        // at this point, all of the jobs are built for running all of the HFFs,
        // but only the builders for non-normal molecules are built.
        let mut energies = vec![0.0; jobs.len()];
        submitter
            .drain("inp", jobs, &mut energies, Check::None)
            .expect("numjac optimizations failed");
        // reset for rest of num_jac
        setup();

        use rayon::prelude::*;
        freqs.reverse(); // for popping
        let mut builders = Vec::new();
        for m in 0..molecules.len() {
            let CoordType::Normal(b) = molecules[m].coord_type else {
                builders.extend(vec![Builder::None; 2 * params.len()]);
                continue;
            };

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
                .into_iter()
                .zip(freqs.pop().unwrap())
                .collect();
            // iterate over the energy/freqparts pairs in parallel
            let bs: Vec<_> = pairs
                .par_iter()
                .enumerate()
                .map(|(i, (energy, freq))| {
                    // there are 2 entries per molecule per row + i, which is
                    // the current entry
                    let index = 2 * m * rows + i;
                    let freq = freq.clone();
                    let norm = Normal::findiff(b);
                    let FreqParts::NormHarm {
                        inner:
                            FinDiff {
                                mut fcs,
                                targets,
                                n,
                                nfc2,
                                nfc3,
                            },
                        mol,
                        pg,
                    } = freq
                    else {
                        unimplemented!()
                    };
                    let (fc2, _, _) = norm.make_fcs(
                        targets,
                        energy,
                        &mut fcs,
                        n,
                        Derivative::Quartic(nfc2, nfc3, 0),
                        None::<&Path>,
                    );
                    let (s, o) = norm.harm_freqs(None, &mol, fc2);
                    Builder::Norm {
                        norm,
                        s,
                        o,
                        pg,
                        ref_energy: Some(geoms[index].energy),
                    }
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
    fn jac_energies<D: Driver>(
        &self,
        rows: usize,
        params: &D::Params,
        geoms: &[ProgramResult],
        molecules: &[config::Molecule],
        mut builders: Vec<Builder>,
    ) -> (Vec<Job<D>>, Vec<Vec<FreqParts>>, Vec<Vec<usize>>) {
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
                pf.incr_value(row, self.delta);
                let mut pb = params.clone();
                pb.decr_value(row, self.delta);

                for (p, pf) in [pf, pb].iter().enumerate() {
                    // idx = job_num so use it twice
                    let (freq, fwd_jobs) = self
                        .build_jobs(
                            &mut output_stream(),
                            geoms[index + p].clone(),
                            pf,
                            idx,
                            idx,
                            molecule,
                            &molecule.coord_type,
                            builders.pop().unwrap(),
                        )
                        .unwrap();
                    // we unwrap here because it's not entirely clear how to
                    // recover from one half of one column going off the rails.
                    idx += fwd_jobs.len();
                    indices[i].push(idx);
                    jobs.extend(fwd_jobs);
                    freqs[i].push(freq);
                }
            }
        }
        (jobs, freqs, indices)
    }
}

#[allow(clippy::too_many_arguments)]
fn build_findiff<F: FiniteDifference, D: Driver>(
    ncoords: usize,
    mol: &Molecule,
    pg: PointGroup,
    energy: f64,
    start_index: usize,
    molecule: &config::Molecule,
    tmpl: &Template,
    coord: &F,
    harmonic: bool,
) -> (FinDiff, Vec<Job<D>>) {
    let n = ncoords;
    let nfc2 = n * n;
    let nfc3 = n * (n + 1) * (n + 2) / 6;
    let nfc4 = n * (n + 1) * (n + 2) * (n + 3) / 24;
    let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];
    let mut target_map = BigHash::new(mol.clone(), pg);
    let deriv = if harmonic {
        Derivative::Harmonic(nfc2)
    } else {
        Derivative::Quartic(nfc2, nfc3, nfc4)
    };
    let geoms = coord.build_points(
        Geom::Xyz(mol.atoms.clone()),
        STEP_SIZE,
        energy,
        deriv,
        &mut fcs,
        &mut target_map,
        n,
    );
    let targets = target_map.values();
    let dir = "inp";
    let mut job_num = start_index;
    let mut jobs = Vec::new();
    for mol in geoms {
        let filename = format!("{dir}/job.{job_num:08}");
        job_num += 1;
        jobs.push(Job::new(
            D::new(filename, tmpl.clone(), molecule.charge, mol.geom.clone()),
            mol.index + start_index,
        ));
    }
    (
        FinDiff {
            n,
            nfc2,
            nfc3,
            fcs,
            targets,
        },
        jobs,
    )
}

impl Default for Frequency {
    fn default() -> Self {
        Self::new(1e-4, false)
    }
}

impl<D> Optimize<D> for Frequency
where
    D: Driver,
{
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q>(
        &self,
        params: &D::Params,
        submitter: &Q,
        molecules: &[config::Molecule],
        ntrue: usize,
    ) -> Option<DVector<f64>>
    where
        Q: Queue<D> + Sync,
    {
        let mut w = output_stream();

        writeln!(w, "Params:\n{params}").unwrap();

        let mut ret = Vec::new();

        for molecule in molecules {
            setup();
            let geom = D::optimize_geometry(
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

            let builder = if let CoordType::Normal(b) = molecule.coord_type {
                // NOTE using fitted normal coordinates
                let norm = Normal::findiff(b);
                // safe to use 0 as job_num because we have to reset directories
                // after the harmonic part anyway
                let tmpl =
                    D::write_params(0, params, molecule.template.clone());
                let CartPart {
                    spectro: s,
                    output: o,
                    ref_energy,
                    pg,
                    ..
                } = norm
                    .cart_part(
                        &FirstPart {
                            template: tmpl.header,
                            optimize: false,
                            geometry: Geom::Xyz(
                                geom.cart_geom.as_ref().unwrap().clone(),
                            ),
                            charge: molecule.charge,
                            step_size: STEP_SIZE,
                            weights: None,
                            dummy_atoms: None,
                            check_int: 0,
                            norm_resume_hff: false,
                        },
                        submitter,
                        &mut w,
                        "inp",
                    )
                    .ok()?;
                setup();
                Builder::Norm {
                    norm,
                    s,
                    o,
                    pg,
                    ref_energy: Some(ref_energy),
                }
            } else {
                Builder::None
            };

            let Ok((freq, jobs)) = self.build_jobs(
                &mut w,
                geom,
                params,
                0,
                0,
                molecule,
                &molecule.coord_type,
                builder,
            ) else {
                return None;
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
            let status =
                submitter.drain("inp", jobs, &mut energies, Check::None);

            if status.is_err() {
                return None;
            }

            takedown();

            // convert energies to frequencies and return those
            let res = self.freqs(&mut w, &mut energies, freq);
            ret.extend(res.iter());
        }
        ret.resize(ntrue, 0.0);
        Some(DVector::from(ret))
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned. `ntrue` is the expected ("true") size of
    /// the frequency vectors in case something goes wrong and a different
    /// number of entries comes out of freqs
    fn num_jac<Q: Queue<D> + Sync>(
        &self,
        params: &D::Params,
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
            .drain("inp", jobs, &mut energies, Check::None)
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
                .into_iter()
                .zip(freqs.pop().unwrap())
                .collect();
            // iterate over the energy/freqparts pairs in parallel
            let freqs: Vec<_> = pairs
                .par_iter()
                .map(|(energy, freq)| {
                    let mut energy = energy.clone();
                    let freq = freq.clone();
                    self.freqs(&mut output_stream(), &mut energy, freq)
                })
                .collect();

            // should be a forward and backward step for each parameter
            assert_eq!(freqs.len(), 2 * params.len());

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
            if tmp.len() != ntrue {
                tmp.resize(ntrue, 0.0);
            }
            // if the whole column is zeros, make the first element 1 to avoid
            // linear algebra issues
            if tmp.iter().all(|s| *s == 0.0) {
                eprintln!("singular column {i} in jacobian, fixing");
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
        write!(logger, "{iter:5}").unwrap();
        for g in got {
            write!(logger, "{g:8.1}").unwrap();
        }
        writeln!(logger, "{:8.1}", mae(got, want)).unwrap();
    }
}

impl Frequency {
    /// helper function for generating the optimizations for the jacobian.
    /// returns a vec of Jobs and also a vec of indices separating each molecule
    fn jac_opt<D: Driver>(
        &self,
        rows: usize,
        params: &D::Params,
        molecules: &[config::Molecule],
    ) -> Vec<Job<D>> {
        let mut opts = Vec::new();
        // each row corresponds to one parameter
        for (i, molecule) in molecules.iter().enumerate() {
            for row in 0..rows {
                let index = 2 * i * rows + 2 * row;
                // forward
                {
                    let mut pf = params.clone();
                    pf.incr_value(row, self.delta);
                    opts.push(Job::new(
                        D::new_full(
                            format!("inp/opt{row}_fwd{i}"),
                            pf,
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
                    pb.decr_value(row, self.delta);
                    opts.push(Job::new(
                        D::new_full(
                            format!("inp/opt{row}_bwd{i}"),
                            pb,
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
