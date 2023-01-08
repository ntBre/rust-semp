use psqs::queue::local::Local;
use rust_semp::{
    config::Config,
    optimize::{energy::Energy, Optimize},
    utils::{load_geoms, load_params},
    *,
};

fn main() {
    let config = Config::load("test_files/test.toml");
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/small07", &names);
    let params = load_params("test_files/small.params");
    Energy { moles }.num_jac(
        &params,
        &Local::new("inp", 128, "/opt/mopac/mopac"),
        &config.molecules,
	9
    );
}
