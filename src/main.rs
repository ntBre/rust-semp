use std::fs::File;
use std::os::unix::io::AsRawFd;

use rust_semp::slurm::Slurm;
use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let outfile = File::create("rust.out").expect("failed to create outfile");
    let logfile = File::create("rust.log").expect("failed to create log file");
    let out_fd = outfile.as_raw_fd();
    let log_fd = logfile.as_raw_fd();
    // redirect stdout to outfile and stderr to logfile
    unsafe {
        libc::dup2(out_fd, 1);
        libc::dup2(log_fd, 2);
    }
    run_algo(names, "file07", "params.dat", "rel.dat", Slurm::default());
}
