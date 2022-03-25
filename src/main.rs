use std::fs::File;
use std::os::unix::io::AsRawFd;

use rust_semp::config::Config;
use rust_semp::slurm::Slurm;
use rust_semp::*;

fn main() {
    // load first arg or default to `semp.toml`
    let args: Vec<String> = std::env::args().collect();
    let mut conf_name = "semp.toml".to_string();
    if args.len() == 2 {
        conf_name = args[1].clone();
    }
    // TODO take these names from conf_name
    let outfile = File::create("rust.out").expect("failed to create outfile");
    let logfile = File::create("rust.log").expect("failed to create log file");
    let out_fd = outfile.as_raw_fd();
    let log_fd = logfile.as_raw_fd();
    // redirect stdout to outfile and stderr to logfile
    unsafe {
        libc::dup2(out_fd, 1);
        libc::dup2(log_fd, 2);
    }
    let conf = Config::load(&conf_name);
    let queue = Slurm::new(conf.chunk_size, conf.job_limit, conf.sleep_int);
    let mut param_log = File::create("params.log")
        .expect("failed to create parameter log file");
    run_algo(
        &mut param_log,
        conf.atom_names,
        "file07",
        "params.dat",
        "rel.dat",
        conf.max_iter,
        queue,
    );
}
