use std::fs::File;
use std::os::unix::io::AsRawFd;

use nalgebra as na;

use psqs::queue::slurm::Slurm;
use rust_semp::config::Config;
use rust_semp::optimize::energy::Energy;
use rust_semp::optimize::frequency::Frequency;
use rust_semp::*;

fn main() {
    // load first arg or default to `semp.toml`
    let args: Vec<String> = std::env::args().collect();
    let mut conf_name = "semp.toml".to_string();
    if args.len() == 2 {
        conf_name = args[1].clone();
    }
    // TODO take these names from conf_name
    let outfile = File::create("semp.out").expect("failed to create outfile");
    let logfile = File::create("semp.log").expect("failed to create log file");
    let out_fd = outfile.as_raw_fd();
    let log_fd = logfile.as_raw_fd();
    // redirect stdout to outfile and stderr to logfile
    unsafe {
        libc::dup2(out_fd, 1);
        libc::dup2(log_fd, 2);
    }
    let conf = Config::load(&conf_name);
    let queue =
        Slurm::new(conf.chunk_size, conf.job_limit, conf.sleep_int, "inp");
    let mut param_log = File::create("params.log")
        .expect("failed to create parameter log file");
    match conf.optimize {
        config::Protocol::Energy => {
            let ai = load_energies("rel.dat");
            run_algo(
                &mut param_log,
                &conf.molecules,
                parse_params(&conf.params),
                ai,
                conf.max_iter,
                conf.broyden,
                conf.broyd_int,
                queue,
                conf.reset_lambda,
                Energy {
                    moles: load_geoms("file07", &conf.molecules[0].atom_names),
                },
            );
        }
        // requirements for a Frequency calculation:
        // 1. pbqff input file called `pbqff.toml`
        // 2. intder template called `intder.in`
        // 3. spectro template called `spectro.in`
        // 4. optimize = "frequency" in `semp.toml`
        // 5. "true" frequencies in rel.dat
        // 6. an empty `file07`
        // 7. irreps for each of the frequencies in rel.dat in `symm`
        config::Protocol::Frequency => {
            // TODO load freqs and irreps for each molecule and combine them
            let ai = Vec::from(load_energies("rel.dat").as_slice());
            let irreps = Frequency::load_irreps("symm");
            let ai = sort_irreps(&ai, &irreps);
            eprintln!("symmetry-sorted true frequencies:");
            for a in &ai {
                eprintln!("{a:8.1}");
            }
            run_algo(
                &mut param_log,
                &conf.molecules,
                parse_params(&conf.params),
                na::DVector::from(ai),
                conf.max_iter,
                conf.broyden,
                conf.broyd_int,
                queue,
                conf.reset_lambda,
                Frequency::new(
                    conf.molecules[0].dummies.clone(),
                    conf.reorder,
                    conf.gspectro_cmd,
                    conf.spectro_cmd,
                ),
            );
        }
    }
}
