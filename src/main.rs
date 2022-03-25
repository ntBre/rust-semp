use std::fs::File;

use rust_semp::slurm::Slurm;
use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let mut outfile =
        File::create("rust.out").expect("failed to create output file");
    let mut logfile =
        File::create("rust.log").expect("failed to create log file");
    run_algo(
        names,
        "file07",
        "params.dat",
        "rel.dat",
        Slurm::default(),
        &mut outfile,
        &mut logfile,
    );
}
