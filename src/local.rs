use std::fs::File;
use std::io::Write;

use crate::queue::{Program, Queue};

/// Minimal implementation for testing MOPAC locally
#[derive(Debug)]
pub struct LocalQueue;

impl<P: Program + Clone> Queue<P> for LocalQueue {
    fn write_submit_script(&self, infiles: Vec<String>, filename: &str) {
        let mut body = String::from("export LD_LIBRARY_PATH=/opt/mopac/\n");
        for f in infiles {
            body.push_str(&format!("/opt/mopac/mopac {f}.mop\n"));
        }
        body.push_str(&format!("date +%s\n"));
        let mut file =
            File::create(filename).expect("failed to create params file");
        write!(file, "{}", body).expect("failed to write params file");
    }

    fn submit_command(&self) -> &str {
        "bash"
    }

    fn chunk_size(&self) -> usize {
        128
    }

    fn job_limit(&self) -> usize {
        1600
    }

    fn sleep_int(&self) -> usize {
        1
    }
}
