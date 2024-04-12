use std::{fmt::Display, ops::Add, str::FromStr};

use psqs::program::mopac::Params;

use crate::{driver, Dvec};

#[derive(Clone)]
pub struct MopacParams(pub psqs::program::mopac::Params);

impl Display for MopacParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for MopacParams {
    type Err = std::string::ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Params::from_str(s).map(Self)
    }
}

impl driver::Params for MopacParams {
    fn names(&self) -> &Vec<String> {
        self.0.names()
    }

    fn atoms(&self) -> &Vec<String> {
        self.0.atoms()
    }

    fn values(&self) -> &Dvec {
        self.0.values()
    }

    fn incr_value(&mut self, idx: usize, delta: f64) {
        self.0.incr_value(idx, delta)
    }
}

impl Add<&Dvec> for MopacParams {
    type Output = Self;

    fn add(self, rhs: &Dvec) -> Self::Output {
        let MopacParams(Params {
            names,
            atoms,
            values,
        }) = self;
        Self(Params {
            names,
            atoms,
            values: values + rhs,
        })
    }
}
