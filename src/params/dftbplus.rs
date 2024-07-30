use std::{fmt::Display, ops::Add, str::FromStr};

use crate::{
    driver::{self, Params},
    Dvec,
};

#[derive(Clone)]
pub struct DFTBPlusParams;

impl Display for DFTBPlusParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

impl Params for DFTBPlusParams {
    fn incr_value(&mut self, idx: usize, delta: f64) {
        todo!()
    }

    fn len(&self) -> usize {
        todo!()
    }
}
