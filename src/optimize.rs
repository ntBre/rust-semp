use psqs::queue::Queue;

use crate::config::Molecule;
use crate::driver::Driver;
use crate::{Dmat, Dvec};

pub mod energy;
pub mod frequency;

pub trait Optimize<D: Driver> {
    fn semi_empirical<Q: Queue<D> + Sync>(
        &self,
        params: &D::Params,
        submitter: &Q,
        molecules: &[Molecule],
        ntrue: usize,
    ) -> Option<Dvec>;

    fn num_jac<Q: Queue<D> + Sync>(
        &self,
        params: &D::Params,
        submitter: &Q,
        molecules: &[Molecule],
        ntrue: usize,
    ) -> Dmat;

    fn stat_multiplier(&self) -> f64;

    /// optional method for logging the actual values on each iteration
    fn log(&self, _iter: usize, _got: &Dvec, _want: &Dvec) {}
}
