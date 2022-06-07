use psqs::queue::local::LocalQueue;
use rust_semp::optimize::energy::Energy;
use rust_semp::*;

fn main() {
    run_algo(
        &mut std::io::stdout(),
        string!["C", "C", "C", "H", "H"],
        "test_files/small07",
        load_params("test_files/small.params"),
        "test_files/25.dat",
        10,
        true,
        5,
        LocalQueue::new("inp", 128),
        0,
        Energy,
    );
}
