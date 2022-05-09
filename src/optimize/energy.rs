use crate::{
    build_jobs,
    mopac::{Mopac, Params},
    queue::Queue,
    relative, Atom, DEBUG, setup, takedown,
};
use nalgebra as na;
use std::rc::Rc;

use super::Optimize;

static DELTA: f64 = 1e-8;
static DELTA_FWD: f64 = 5e7; // 1 / 2Δ
static DELTA_BWD: f64 = -5e7; // -1 / 2Δ

pub struct Energy;

impl Optimize for Energy {
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Vec<Atom>>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DVector<f64> {
        let mut jobs = build_jobs(moles, params, 0, 1.0, 0, charge);
        let mut got = vec![0.0; jobs.len()];
	setup();
        submitter.drain(&mut jobs, &mut got);
	takedown();
        relative(&na::DVector::from(got))
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned
    fn num_jac<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Vec<Atom>>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DMatrix<f64> {
        let rows = params.values.len();
        let cols = moles.len();
        let mut jac_t = vec![0.; rows * cols];
        // front and back for each row and col
        let mut jobs = Vec::with_capacity(2 * rows * cols);
        let mut job_num = 0;
        for row in 0..rows {
            // loop over params
            let mut pf = params.clone();
            let mut pb = params.clone();
            let idx = row * cols;
            pf.values[row] += DELTA;
            let mut fwd_jobs =
                build_jobs(&moles, &pf, idx, DELTA_FWD, job_num, charge);
            job_num += fwd_jobs.len();
            jobs.append(&mut fwd_jobs);

            pb.values[row] -= DELTA;
            let mut bwd_jobs =
                build_jobs(&moles, &pb, idx, DELTA_BWD, job_num, charge);
            job_num += bwd_jobs.len();
            jobs.append(&mut bwd_jobs);
        }
        if DEBUG {
            eprintln!("num_jac: running {} jobs", jobs.len());
        }
	setup();
        submitter.drain(&mut jobs, &mut jac_t);
	takedown();
        // nalgebra does from_vec in col-major order, so lead with cols and I get
        // jac_t_t or jac back
        let jac = na::DMatrix::from_vec(cols, rows, jac_t);
        jac
    }
}
