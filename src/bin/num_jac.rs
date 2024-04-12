use psqs::queue::local::Local;
use rust_semp::{
    config::Config,
    optimize::{energy::Energy, Optimize},
    params::MopacParams,
    utils::{load_geoms, load_params},
    *,
};

fn main() {
    let config = Config::load("test_files/test.toml");
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/small07", &names);
    let params = load_params("test_files/small.params");
    Energy { moles }.num_jac(
        &MopacParams(params),
        &Local::new(128, 128, 2, "inp", true, None),
        &config.molecules,
        9,
    );
}
