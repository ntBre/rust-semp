use psqs::program::{mopac::Mopac, Program};
use serde::{Deserialize, Serialize};

use crate::Dvec;

pub trait Params {
    fn new(names: Vec<String>, atoms: Vec<String>, values: Dvec) -> Self;
    fn names(&self) -> &Vec<String>;
    fn atoms(&self) -> &Vec<String>;
    fn values(&self) -> &Dvec;
}

pub trait Driver:
    Program + Clone + Sync + Send + Serialize + for<'a> Deserialize<'a>
{
    type Params: Params;
}

impl Params for psqs::program::mopac::Params {
    fn names(&self) -> &Vec<String> {
        &self.names
    }

    fn atoms(&self) -> &Vec<String> {
        &self.atoms
    }

    fn values(&self) -> &Dvec {
        &self.values
    }

    fn new(names: Vec<String>, atoms: Vec<String>, values: Dvec) -> Self {
        psqs::program::mopac::params::Params::new(names, atoms, values)
    }
}

impl Driver for Mopac {
    type Params = psqs::program::mopac::Params;
}
