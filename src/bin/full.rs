use psqs::queue::local::LocalQueue;
use rust_semp::config::Molecule;
use rust_semp::optimize::energy::Energy;
use rust_semp::*;

fn main() {
    run_algo(
        &mut std::io::stdout(),
        vec![Molecule {
            atom_names: string!["C", "C", "C", "H", "H"],
            charge: 0,
            dummies: vec![],
        }],
        "test_files/small07",
        load_params("test_files/small.params"),
        load_energies("test_files/25.dat"),
        10,
        true,
        5,
        LocalQueue::new("inp", 128),
        false,
        Energy,
    );
}
