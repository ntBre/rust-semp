use rust_semp::slurm::Slurm;
use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    run_algo(names, "file07", "params.dat", "rel.dat", Slurm::default());
}
