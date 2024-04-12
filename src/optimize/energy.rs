use psqs::geom::Geom;
use psqs::program::mopac::Mopac;
use psqs::queue::{Check, Queue};

use crate::config::Molecule;
use crate::params::mopac::{self, MopacParams};
use crate::utils::relative;
use crate::Dvec;
use crate::{utils::setup, utils::takedown, DEBUG};
use nalgebra as na;

use super::Optimize;

static DELTA: f64 = 1e-8;
static DELTA_FWD: f64 = 5e7; // 1 / 2Δ
static DELTA_BWD: f64 = -5e7; // -1 / 2Δ

pub struct Energy {
    pub moles: Vec<Geom>,
}

impl Optimize<Mopac> for Energy {
    /// compute the semi-empirical energies of `moles` for the given `params`
    fn semi_empirical<Q>(
        &self,
        params: &MopacParams,
        submitter: &Q,
        molecules: &[Molecule],
        _ntrue: usize,
    ) -> Option<Dvec>
    where
        Q: Queue<Mopac> + Sync,
    {
        // NOTE: still no loop over molecules here
        let jobs = Mopac::build_jobs(
            self.moles.clone(),
            Some(&params.0),
            "inp",
            0,
            1.0,
            0,
            molecules[0].charge,
            molecules[0].template.clone(),
        );
        let mut got = vec![0.0; jobs.len()];
        setup();
        submitter
            .drain("inp", jobs, &mut got, Check::None)
            .expect("energies failed");
        takedown();
        Some(relative(&na::DVector::from(got)))
    }

    /// Compute the numerical Jacobian for the geomeries in `moles` and the
    /// parameters in `params`. For convenience of indexing, the transpose is
    /// actually computed and returned
    fn num_jac<Q: Queue<Mopac> + std::marker::Sync>(
        &self,
        params: &mopac::MopacParams,
        submitter: &Q,
        molecules: &[Molecule],
        _ntrue: usize,
    ) -> na::DMatrix<f64> {
        // NOTE: no loop over molecules here
        let rows = params.0.values.len();
        let cols = self.moles.len();
        let mut jac_t = vec![0.; rows * cols];
        // front and back for each row and col
        let mut jobs = Vec::with_capacity(2 * rows * cols);
        let mut job_num = 0;
        for row in 0..rows {
            // loop over params
            let mut pf = params.clone();
            let mut pb = params.clone();
            let idx = row * cols;
            pf.0.values[row] += DELTA;
            let mut fwd_jobs = Mopac::build_jobs(
                self.moles.clone(),
                Some(&pf.0),
                "inp",
                idx,
                DELTA_FWD,
                job_num,
                molecules[0].charge,
                molecules[0].template.clone(),
            );
            job_num += fwd_jobs.len();
            jobs.append(&mut fwd_jobs);

            pb.0.values[row] -= DELTA;
            let mut bwd_jobs = Mopac::build_jobs(
                self.moles.clone(),
                Some(&pb.0),
                "inp",
                idx,
                DELTA_BWD,
                job_num,
                molecules[0].charge,
                molecules[0].template.clone(),
            );
            job_num += bwd_jobs.len();
            jobs.append(&mut bwd_jobs);
        }
        if *DEBUG == "jobs" {
            eprintln!("num_jac: running {} jobs", jobs.len());
        }
        setup();
        submitter
            .drain("inp", jobs, &mut jac_t, Check::None)
            .expect("numjac energies failed");
        takedown();
        // nalgebra does from_vec in col-major order, so lead with cols and I get
        // jac_t_t or jac back

        na::DMatrix::from_vec(cols, rows, jac_t)
    }

    fn stat_multiplier(&self) -> f64 {
        static HT_TO_CM: f64 = 219_474.5459784;
        HT_TO_CM
    }

    fn log(&self, _iter: usize, _got: &Dvec, _want: &Dvec) {}
}
