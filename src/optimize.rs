use std::rc::Rc;

use nalgebra as na;

use psqs::atom::Geom;
use psqs::program::mopac::{Mopac, Params};
use psqs::queue::Queue;

pub mod energy;
pub mod frequency;

pub trait Optimize {
    fn semi_empirical<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Geom>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DVector<f64>;

    fn num_jac<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Geom>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DMatrix<f64>;
}
