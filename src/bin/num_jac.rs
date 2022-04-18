use crate::local::LocalQueue;
use rust_semp::*;

fn main() {
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/small07", &names);
    let params = load_params("test_files/small.params");
    num_jac(&moles, &params, &LocalQueue, 0);
}
