use psqs::queue::local::LocalQueue;
use rust_semp::{
    config::Config,
    optimize::{energy::Energy, Optimize},
    *,
};

fn main() {
    let config = Config::load("test_files/test.toml");
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/small07", &names);
    let params = load_params("test_files/small.params");
    Energy {
        moles: moles.clone(),
    }
    .num_jac(&params, &LocalQueue::new("inp", 128), &config.molecules);
}
