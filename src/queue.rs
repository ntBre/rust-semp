use core::time;
// this is going to make my psandqs (Ps and Qs) crate - programs and queues
use std::{collections::HashMap, process::Command, str, thread};

use crate::{dump::Dump, DEBUG};

pub trait Program {
    fn filename(&self) -> String;

    fn write_input(&mut self);

    fn read_output(&self) -> Option<f64>;

    /// Return all the filenames associated with the Program
    fn associated_files(&self) -> Vec<String>;
}

// TODO move this into a separate file
#[derive(Debug, Clone)]
pub struct Job<P: Program> {
    pub program: P,
    pub pbs_file: String,
    pub job_id: String,
    pub index: usize,
    pub coeff: f64,
}

impl<P: Program> Job<P> {
    pub fn new(program: P, index: usize) -> Self {
        Self {
            program,
            pbs_file: String::new(),
            job_id: String::new(),
            index,
            coeff: 1.0,
        }
    }
}

pub trait Queue<P>
where
    P: Program + Clone,
{
    fn write_submit_script(&self, infiles: Vec<String>, filename: &str);

    fn submit_command(&self) -> &str;

    fn chunk_size(&self) -> usize;

    fn job_limit(&self) -> usize;

    fn sleep_int(&self) -> usize;

    fn submit(&self, filename: &str) -> String {
        match Command::new(self.submit_command()).arg(filename).output() {
            Ok(s) => {
                return str::from_utf8(&s.stdout).unwrap().trim().to_string()
            }
            Err(_) => todo!(),
        };
    }

    /// Build a chunk of jobs by writing the Program input file and the
    /// corresponding submission script and then submitting the script
    fn build_chunk<'a>(
        &self,
        jobs: &'a [Job<P>],
        chunk_num: usize,
        slurm_jobs: &'a mut HashMap<String, usize>,
    ) -> Vec<Job<P>> {
        let queue_file = format!("inp/main{}.slurm", chunk_num);
        let mut chunk_jobs = Vec::new();
        for job in jobs {
            let mut job = (*job).clone();
            job.program.write_input();
            job.pbs_file = queue_file.to_string();
            chunk_jobs.push(job);
        }
        slurm_jobs.insert(queue_file.clone(), chunk_jobs.len());
        self.write_submit_script(
            chunk_jobs.iter().map(|j| j.program.filename()).collect(),
            &queue_file,
        );
        // run jobs
        let job_id = self.submit(&queue_file);
        for job in &mut chunk_jobs {
            job.job_id = job_id.clone();
        }
        chunk_jobs
    }

    fn drain(&self, jobs: &mut Vec<Job<P>>, dst: &mut [f64]) {
        let mut chunk_num: usize = 0;
        let mut cur_jobs = Vec::new();
        let mut slurm_jobs = HashMap::new();
        let mut cur = 0; // current index into jobs
        let tot_jobs = jobs.len();
        let mut dump = Dump::new(self.chunk_size() * 5);
        loop {
            let mut finished = 0;
            while cur_jobs.len() < self.job_limit() && cur < tot_jobs {
                let new_chunk = self.build_chunk(
                    &mut jobs
                        [cur..std::cmp::min(cur + self.chunk_size(), tot_jobs)],
                    chunk_num,
                    &mut slurm_jobs,
                );
                if DEBUG {
                    eprintln!("submitted chunk {}", chunk_num);
                }
                chunk_num += 1;
                cur += new_chunk.len();
                cur_jobs.extend(new_chunk);
            }
            // collect output
            let mut to_remove = Vec::new();
            for (i, job) in cur_jobs.iter_mut().enumerate() {
                match job.program.read_output() {
                    Some(val) => {
                        to_remove.push(i);
                        dst[job.index] += job.coeff * val;
                        dump.add(job.program.associated_files());
                        finished += 1;
                        let job_name = job.pbs_file.as_str();
                        let count = slurm_jobs.get_mut(job_name).unwrap();
                        *count -= 1;
                        if *count == 0 {
                            // delete the submit script
                            dump.add(vec![
                                job_name.to_string(),
                                format!("{}.out", job_name),
                            ]);
                        }
                    }
                    None => (),
                }
            }
            // have to remove the highest index first so sort and reverse
            to_remove.sort();
            to_remove.reverse();
            for i in to_remove {
                cur_jobs.swap_remove(i);
            }
            if cur_jobs.len() == 0 && cur == tot_jobs {
                return;
            }
            if finished == 0 {
                eprintln!("none finished, sleeping");
                thread::sleep(time::Duration::from_secs(
                    self.sleep_int() as u64
                ));
            } else {
                eprintln!("finished {}", finished);
            }
        }
    }
}
