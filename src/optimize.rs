use nalgebra as na;

use psqs::program::mopac::{Mopac, Params};
use psqs::queue::Queue;

use crate::config::Molecule;

pub mod energy;
pub mod frequency;

pub trait Optimize {
    fn semi_empirical<Q: Queue<Mopac> + Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[Molecule],
    ) -> Option<na::DVector<f64>>;

    fn num_jac<Q: Queue<Mopac> + Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[Molecule],
        ntrue: usize,
    ) -> na::DMatrix<f64>;

    fn stat_multiplier(&self) -> f64;

    /// optional method for logging the actual values on each iteration
    fn log(
        &self,
        _iter: usize,
        _got: &na::DVector<f64>,
        _want: &na::DVector<f64>,
    ) {
    }
}
