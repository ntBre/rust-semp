use std::str::FromStr;

use psqs::{program::dftbplus::DFTBPlus, queue::local::Local};
use rust_semp::{
    config::Config,
    optimize::frequency::Frequency,
    params::dftbplus::DFTBPlusParams,
    run_algo,
    utils::{setup, takedown},
};

fn main() {
    env_logger::init();
    setup();
    let config = Config::load("test_files/dftb.toml");
    let queue = Local {
        dir: "inp".to_owned(),
        chunk_size: 128,
        template: None,
    };
    run_algo::<DFTBPlus, _, _, _>(
        &mut std::io::sink(),
        &config.molecules,
        DFTBPlusParams::from_str(&config.params).unwrap(),
        nalgebra::DVector::from(rust_semp::utils::sort_freqs(&config)),
        5,
        true,
        5,
        queue,
        false,
        Frequency::new(config.delta, false),
    );
    takedown();
}
