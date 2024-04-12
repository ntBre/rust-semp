//! Implementations of the [Params] interface

pub mod mopac;

pub mod molpro {
    use std::{fmt::Display, ops::Add, str::FromStr};

    use crate::{driver::Params, Dvec};

    /// A basis set
    #[derive(Clone)]
    pub struct MolproParams {}

    impl FromStr for MolproParams {
        type Err = std::string::ParseError;

        fn from_str(_s: &str) -> Result<Self, Self::Err> {
            todo!()
        }
    }

    impl Display for MolproParams {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}", std::any::type_name::<Self>())
        }
    }

    impl Params for MolproParams {
        #[allow(unused)]
        fn incr_value(&mut self, idx: usize, delta: f64) {
            todo!()
        }

        fn len(&self) -> usize {
            todo!()
        }
    }

    impl Add<&Dvec> for MolproParams {
        type Output = Self;

        fn add(self, _rhs: &Dvec) -> Self::Output {
            todo!()
        }
    }
}
