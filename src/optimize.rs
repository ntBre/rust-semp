use psqs::program::mopac::{Mopac, Params};
use psqs::queue::Queue;

use crate::config::Molecule;
use crate::{Dmat, Dvec};

pub mod energy;
pub mod frequency;

pub trait Optimize {
    fn semi_empirical<Q: Queue<Mopac> + Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[Molecule],
        ntrue: usize,
    ) -> Option<Dvec>;

    fn num_jac<Q: Queue<Mopac> + Sync>(
        &self,
        params: &Params,
        submitter: &Q,
        molecules: &[Molecule],
        ntrue: usize,
    ) -> Dmat;

    fn stat_multiplier(&self) -> f64;

    /// optional method for logging the actual values on each iteration
    fn log(&self, _iter: usize, _got: &Dvec, _want: &Dvec) {}
}
