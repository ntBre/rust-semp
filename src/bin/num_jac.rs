use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("small07", names);
    let params = load_params("small.params");
    setup();
    num_jac(&moles, &params, LocalQueue {});
    takedown();
}
