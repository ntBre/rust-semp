use psqs::queue::local::LocalQueue;
use rust_semp::config::Molecule;
use rust_semp::optimize::energy::Energy;
use rust_semp::*;

fn main() {
    let names = string!["C", "C", "C", "H", "H"];
    run_algo(
        &mut std::io::stdout(),
        &[Molecule {
            atom_names: names.clone(),
            charge: 0,
            dummies: vec![],
        }],
        load_params("test_files/small.params"),
        load_energies("test_files/25.dat"),
        10,
        true,
        5,
        LocalQueue::new("inp", 128),
        false,
        Energy {
            moles: load_geoms("test_files/small07", &names),
        },
    );
}
