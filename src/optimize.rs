use std::rc::Rc;

use nalgebra as na;

use crate::{
    mopac::{Mopac, Params},
    queue::Queue,
    Atom,
};

pub mod energy;
pub mod frequency;

pub trait Optimize {
    fn semi_empirical<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Vec<Atom>>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DVector<f64>;

    fn num_jac<Q: Queue<Mopac>>(
        &self,
        moles: &Vec<Rc<Vec<Atom>>>,
        params: &Params,
        submitter: &Q,
        charge: isize,
    ) -> na::DMatrix<f64>;
}
