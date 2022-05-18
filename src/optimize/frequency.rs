#![allow(unused)]
use psqs::geom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::program::Job;
use psqs::queue::Queue;

use crate::{build_jobs, relative, setup, takedown, DEBUG, MOPAC_TMPL};
use nalgebra as na;
use std::rc::Rc;

use super::Optimize;

static DELTA: f64 = 1e-8;
static DELTA_FWD: f64 = 5e7; // 1 / 2Δ
static DELTA_BWD: f64 = -5e7; // -1 / 2Δ

pub struct Frequency {
    pub config: rust_pbqff::config::Config,
    pub intder: rust_pbqff::Intder,
    pub spectro: rust_pbqff::Spectro,
}

pub fn optimize_geometry<Q: Queue<Mopac>>(geom: Geom, queue: &Q) -> Geom {
    // TODO handle error
    let _ = std::fs::create_dir("opt");
    let opt = Job::new(
        Mopac::new("opt/opt".to_string(), None, Rc::new(geom), 0, &MOPAC_TMPL),
        0,
    );
    queue.optimize(opt)
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
        let mut intder = self.intder.clone();
        // optimize
        let geom = optimize_geometry(self.config.geometry, submitter);
        // generate pts == moles
        let (moles, taylor, taylor_disps, atomic_numbers) =
            rust_pbqff::generate_pts(geom, &mut intder);
        // call build_jobs like before
        let mut jobs = build_jobs(&moles, params, 0, 1.0, 0, charge);
        let mut energies = vec![0.0; jobs.len()];
        // drain to get energies
        setup();
        submitter.drain(&mut jobs, &mut energies);
        takedown();
        // convert energies to frequencies and return those
        println!(
            "{}",
            na::DVector::from(
                rust_pbqff::freqs(
                    energies,
                    &mut intder,
                    &taylor,
                    &taylor_disps,
                    &atomic_numbers,
                    &self.spectro,
                    &self.config.spectro,
                )
                .corr,
            )
        );
        todo!()
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned
    fn num_jac<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Geom>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DMatrix<f64> {
        todo!()
    }
}
