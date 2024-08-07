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
    fn incr_value(&mut self, idx: usize, delta: f64) {
        self.0.incr_value(idx, delta)
    }

    fn len(&self) -> usize {
        self.0.names.len()
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
            values: (Dvec::from(values) + rhs).data.as_vec().clone(),
        })
    }
}
