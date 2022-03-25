use std::fs::File;
use std::io::Write;

use crate::{mopac::Mopac, queue::Queue};

/// Slurm is a type for holding the information for submitting a slurm job.
/// `filename` is the name of the Slurm submission script
#[derive(Debug)]
pub struct Slurm {
    chunk_size: usize,
    job_limit: usize,
    sleep_int: usize,
}

impl Slurm {
    pub fn new(chunk_size: usize, job_limit: usize, sleep_int: usize) -> Self {
        Self {
            chunk_size,
            job_limit,
            sleep_int,
        }
    }

    pub fn default() -> Self {
        Self {
            chunk_size: 128,
            job_limit: 1600,
            sleep_int: 5,
        }
    }
}

impl Queue<Mopac> for Slurm {
    fn write_submit_script(&self, infiles: Vec<String>, filename: &str) {
        let mut body = format!(
            "#!/bin/bash
#SBATCH --job-name=semp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o {filename}.out
#SBATCH --no-requeue
#SBATCH --mem=1gb
export LD_LIBRARY_PATH=/home/qc/mopac2016/\n",
        );
        for f in infiles {
            body.push_str(&format!(
                "/home/qc/mopac2016/MOPAC2016.exe {f}.mop\n"
            ));
        }
        let mut file = match File::create(filename) {
            Ok(f) => f,
            Err(_) => {
                panic!("write_submit_script: failed to create {filename}");
            }
        };
        write!(file, "{}", body).expect("failed to write params file");
    }

    fn submit_command(&self) -> &str {
        "sbatch"
    }

    fn chunk_size(&self) -> usize {
        self.chunk_size
    }

    fn job_limit(&self) -> usize {
        self.job_limit
    }

    fn sleep_int(&self) -> usize {
        self.sleep_int
    }
}
