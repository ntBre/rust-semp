use std::fs;

use rust_semp::{mopac::*, queue::Submit, *};

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = rust_semp::load_geoms("three07", names);
    let params = load_params("params.dat");
    // build jobs (in memory) -> write jobs (to disk) -> run jobs
    let mut count: usize = 0;
    let mut jobs = Vec::new();
    for mol in moles {
        let filename = format!("inp/job.{count:08}");
        count += 1;
        jobs.push(Mopac {
            filename,
            params: params.clone(),
            geom: mol,
        })
    }
    // write jobs - know all params are the same right now
    match fs::create_dir("inp") {
        Ok(_) => (),
        Err(_) => (), // ideally I'd only accept not exist error but idk
    }
    let mut chunk = Vec::new();
    jobs[0].write_params("inp/params.dat");
    for job in jobs {
        job.write_input("inp/params.dat");
        chunk.push(job.filename);
    }
    let slurm = LocalQueue::new("inp/main.slurm");
    slurm.write_submit_script(chunk);
    // run jobs
    println!("{}", slurm.submit());
}

// really params doesn't need to be a field of Mopac, each Mopac doesn't really
// have its own params, a group of Mopacs will always share a params. I guess I
// could have chunks contained within groups as well

struct Group {
    jobs: Vec<Mopac>,
    params: Vec<Param>,
}
