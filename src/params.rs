//! Implementations of the [Params] interface

pub mod mopac;

pub mod molpro {
    use std::fmt::Display;

    use crate::driver::Params;

    /// A basis set
    #[derive(Clone)]
    pub struct MolproParams {}

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
}
