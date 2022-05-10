use psqs::queue::local::LocalQueue;
use rust_semp::{
    optimize::{energy::Energy, Optimize},
    *,
};

fn main() {
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/small07", &names);
    let params = load_params("test_files/small.params");
    Energy.num_jac(&moles, &params, &LocalQueue{ dir: "inp".to_string() }, 0);
}
