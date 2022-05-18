#![allow(unused)]
use psqs::geom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::program::{Job, Template};
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
        setup(); // setup because submitter tries to write in inp, not opt
        let geom =
            optimize_geometry(self.config.geometry.clone(), params, submitter);
        // generate pts == moles
        let (moles, taylor, taylor_disps, atomic_numbers) =
            rust_pbqff::generate_pts(geom, &mut intder);
        // dir created in generate_pts but unused here
        let _ = std::fs::remove_dir_all("pts");
        // call build_jobs like before
        let mut jobs = build_jobs(&moles, params, 0, 1.0, 0, charge);
        let mut energies = vec![0.0; jobs.len()];
        // drain to get energies
        submitter.drain(&mut jobs, &mut energies);
        takedown();
        // convert energies to frequencies and return those
        let res = na::DVector::from(
            rust_pbqff::freqs(
                energies,
                &mut intder,
                &taylor,
                &taylor_disps,
                &atomic_numbers,
                &self.spectro,
                &self.config.gspectro_cmd,
                &self.config.spectro_cmd,
                &self.config.summary_cmd,
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
        moles: &Vec<Rc<Geom>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DMatrix<f64> {
        // TODO do this in parallel, but as a first pass just use semi_empirical
        // serially
        let rows = params.values.len();
        let mut jac_t = Vec::new();
        for row in 0..rows {
            let mut pf = params.clone();
            pf.values[row] += DELTA;
            let fwd_freqs = self.semi_empirical(moles, &pf, submitter, charge);

            let mut pb = params.clone();
            pb.values[row] -= DELTA;
            let bwd_freqs = self.semi_empirical(moles, &pb, submitter, charge);

            // (fwd_freqs - bwd_freqs) / 2DELTA gives a column of jac_t
            jac_t.extend_from_slice(
                ((fwd_freqs - bwd_freqs) / (2. * DELTA)).as_slice(),
            );
        }
        let cols = jac_t.len() / rows;
        na::DMatrix::from_vec(cols, rows, jac_t)
    }

    fn stat_multiplier(&self) -> f64 {
        1.0
    }
}
