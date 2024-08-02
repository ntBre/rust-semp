use std::fs::File;
use std::os::unix::io::AsRawFd;
use std::str::FromStr;

use nalgebra as na;

use params::dftbplus::DFTBPlusParams;
use pbqff::config::Queue;
use psqs::program::dftbplus::DFTBPlus;
use psqs::program::molpro::Molpro;
use psqs::program::mopac::Mopac;
use psqs::queue::local::Local;
use psqs::queue::pbs::Pbs;
use psqs::queue::slurm::Slurm;
use rust_semp::config::Config;
use rust_semp::optimize::energy::Energy;
use rust_semp::optimize::frequency::Frequency;
use rust_semp::params::molpro::MolproParams;
use rust_semp::params::mopac::MopacParams;
use rust_semp::utils::{load_energies, load_geoms, parse_params};
use rust_semp::*;

fn main() {
    env_logger::init();
    // load first arg or default to `semp.toml`
    let args: Vec<String> = std::env::args().collect();
    let mut conf_name = "semp.toml".to_string();
    if args.len() == 2 {
        conf_name.clone_from(&args[1]);
    }
    // TODO take these names from conf_name
    let outfile = File::create("semp.out").expect("failed to create outfile");
    let logfile = File::create("semp.log").expect("failed to create log file");
    let out_fd = outfile.as_raw_fd();
    let log_fd = logfile.as_raw_fd();
    // redirect stdout to outfile and stderr to logfile
    unsafe {
        if libc::dup2(out_fd, 1) == -1 {
            panic!("failed to dup stdout");
        }
        if libc::dup2(log_fd, 2) == -1 {
            panic!("failed to dup stderr");
        }
    }
    let conf = Config::load(&conf_name);
    let mut param_log = File::create("params.log")
        .expect("failed to create parameter log file");
    match (conf.optimize, conf.program) {
        (config::Protocol::Energy, config::ProgramType::Mopac) => {
            let ai = load_energies("rel.dat");
            match conf.queue {
                Queue::Pbs => run_algo(
                    &mut param_log,
                    &conf.molecules,
                    MopacParams(parse_params(&conf.params)),
                    ai,
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Pbs::new(
                        conf.chunk_size,
                        conf.job_limit,
                        conf.sleep_int,
                        "inp",
                        false,
                        conf.queue_template,
                    ),
                    conf.reset_lambda,
                    Energy {
                        moles: load_geoms(
                            "file07",
                            &conf.molecules[0].atom_names,
                        ),
                    },
                ),
                Queue::Slurm => run_algo(
                    &mut param_log,
                    &conf.molecules,
                    MopacParams(parse_params(&conf.params)),
                    ai,
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Slurm::new(
                        conf.chunk_size,
                        conf.job_limit,
                        conf.sleep_int,
                        "inp",
                        false,
                        conf.queue_template,
                    ),
                    conf.reset_lambda,
                    Energy {
                        moles: load_geoms(
                            "file07",
                            &conf.molecules[0].atom_names,
                        ),
                    },
                ),
                _ => panic!("unsupported queue `{:#?}`", conf.queue),
            };
        }
        // requirements for a Frequency calculation:
        // 1. intder template called `intder.in` (for SICs)
        // 2. optimize = "frequency" in `semp.toml`
        (config::Protocol::Frequency, config::ProgramType::Mopac) => {
            let ai = crate::utils::sort_freqs(&conf);
            match conf.queue {
                Queue::Pbs => run_algo::<Mopac, _, _, _>(
                    &mut param_log,
                    &conf.molecules,
                    MopacParams(parse_params(&conf.params)),
                    na::DVector::from(ai),
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Pbs::new(
                        conf.chunk_size,
                        conf.job_limit,
                        conf.sleep_int,
                        "inp",
                        false,
                        conf.queue_template,
                    ),
                    conf.reset_lambda,
                    Frequency::new(conf.delta, conf.sort_ascending),
                ),
                Queue::Slurm => run_algo::<Mopac, _, _, _>(
                    &mut param_log,
                    &conf.molecules,
                    MopacParams(parse_params(&conf.params)),
                    na::DVector::from(ai),
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Slurm::new(
                        conf.chunk_size,
                        conf.job_limit,
                        conf.sleep_int,
                        "inp",
                        false,
                        conf.queue_template,
                    ),
                    conf.reset_lambda,
                    Frequency::new(conf.delta, conf.sort_ascending),
                ),
                Queue::Local => run_algo::<Mopac, _, _, _>(
                    &mut param_log,
                    &conf.molecules,
                    MopacParams(parse_params(&conf.params)),
                    na::DVector::from(ai),
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Local {
                        dir: "inp".to_owned(),
                        chunk_size: conf.chunk_size,
                        mopac: conf.mopac.expect(
                            "mopac field must be specified for local queue",
                        ),
                        template: None,
                    },
                    conf.reset_lambda,
                    Frequency::new(conf.delta, conf.sort_ascending),
                ),
            };
        }
        (config::Protocol::Energy, config::ProgramType::Molpro) => todo!(),
        (config::Protocol::Frequency, config::ProgramType::Molpro) => {
            let ai = crate::utils::sort_freqs(&conf);
            match conf.queue {
                Queue::Pbs => run_algo::<Molpro, _, _, _>(
                    &mut param_log,
                    &conf.molecules,
                    MolproParams::from_str(&conf.params).unwrap(),
                    na::DVector::from(ai),
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Pbs::new(
                        conf.chunk_size,
                        conf.job_limit,
                        conf.sleep_int,
                        "inp",
                        false,
                        conf.queue_template,
                    ),
                    conf.reset_lambda,
                    Frequency::new(conf.delta, conf.sort_ascending),
                ),
                _ => todo!(),
            };
        }
        (config::Protocol::Energy, config::ProgramType::Dftb) => todo!(),
        (config::Protocol::Frequency, config::ProgramType::Dftb) => {
            let ai = crate::utils::sort_freqs(&conf);
            match conf.queue {
                Queue::Pbs => run_algo::<DFTBPlus, _, _, _>(
                    &mut param_log,
                    &conf.molecules,
                    DFTBPlusParams::from_str(&conf.params).unwrap(),
                    na::DVector::from(ai),
                    conf.max_iter,
                    conf.broyden,
                    conf.broyd_int,
                    Pbs::new(
                        conf.chunk_size,
                        conf.job_limit,
                        conf.sleep_int,
                        "inp",
                        false,
                        conf.queue_template,
                    ),
                    conf.reset_lambda,
                    Frequency::new(conf.delta, conf.sort_ascending),
                ),
                _ => todo!(),
            };
        }
    }
}
