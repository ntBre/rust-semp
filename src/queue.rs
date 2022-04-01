use core::time;
// this is going to make my psandqs (Ps and Qs) crate - programs and queues
use std::{
    collections::{HashMap, HashSet},
    fs,
    path::Path,
    process::Command,
    str, thread,
};

use crate::{dump::Dump, DEBUG};

#[derive(Debug, PartialEq)]
pub enum ProgramStatus {
    Success(f64),
    FileNotFound,
    ErrorInOutput,
    EnergyNotFound,
    EnergyParseError,
}

pub trait Program {
    fn filename(&self) -> String;

    fn write_input(&mut self);

    fn read_output(&self) -> ProgramStatus;

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

// TODO this really doesn't belong here
static DIRS: &'static [&str] = &["inp", "tmparam"];

/// set up the directories needed for the program after deleting existing ones
pub fn setup() {
    takedown();
    for dir in DIRS {
        match fs::create_dir(dir) {
            Ok(_) => (),
            Err(_) => {
                eprintln!("can't create '{}'", dir);
                std::process::exit(1);
            }
        }
    }
}

// TODO don't do this if -nodel flag
pub fn takedown() {
    for dir in DIRS {
        let path = Path::new(dir);
        if path.is_dir() {
            match fs::remove_dir_all(dir) {
                Ok(_) => (),
                Err(_) => {
                    eprintln!("can't remove '{}'", dir);
                    std::process::exit(1);
                }
            }
        }
    }
}

pub trait Queue<P>
where
    P: Program + Clone,
{
    /// the extension to append to submit scripts for this type of Queue
    const SCRIPT_EXT: &'static str;

    const DIR: &'static str = "inp";

    fn write_submit_script(&self, infiles: Vec<String>, filename: &str);

    fn submit_command(&self) -> &str;

    fn chunk_size(&self) -> usize;

    fn job_limit(&self) -> usize;

    fn sleep_int(&self) -> usize;

    fn submit(&self, filename: &str) -> String {
        match Command::new(self.submit_command()).arg(filename).output() {
            Ok(s) => {
                let raw = str::from_utf8(&s.stdout).unwrap().trim().to_string();
                return raw.split_whitespace().last().unwrap().to_string();
            }
            Err(_) => todo!(),
        };
    }

    /// the command to check the status of jobs in the queue
    fn stat_cmd(&self) -> String;

    /// take a name of a Program input file with the extension attached, replace
    /// the extension (ext) with _redo.ext and write _redo.SCRIPT_EXT, then
    /// submit the redo script
    fn resubmit(&self, filename: &str) -> String {
        let path = Path::new(filename);
        let dir = path.parent().unwrap().to_str().unwrap();
        let ext = path.extension().unwrap().to_str().unwrap();
        let base = path.file_stem().unwrap().to_str().unwrap();
        let redo_inp = format!("{}/{}_redo.{}", dir, base, ext);
        std::fs::copy(filename, &redo_inp).unwrap();
        let redo_pbs = format!("{}/{}_redo.{}", dir, base, Self::SCRIPT_EXT);
        self.write_submit_script(vec![String::from(redo_inp)], &redo_pbs);
        self.submit(&redo_pbs)
    }

    /// return a HashSet of jobs found in the queue based on the output of
    /// `stat_cmd`
    fn status(&self) -> HashSet<String> {
        let mut ret = HashSet::new();
        let lines = self.stat_cmd();
        let lines = lines.lines();
        for line in lines {
            if !line.contains("JOBID") {
                ret.insert(line.split_whitespace().next().unwrap().to_string());
            }
        }
        ret
    }

    /// Build a chunk of jobs by writing the Program input file and the
    /// corresponding submission script and then submitting the script
    fn build_chunk<'a>(
        &self,
        jobs: &'a [Job<P>],
        chunk_num: usize,
        slurm_jobs: &'a mut HashMap<String, usize>,
    ) -> Vec<Job<P>> {
        let queue_file =
            format!("{}/main{}.{}", Self::DIR, chunk_num, Self::SCRIPT_EXT);
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
        let mut remaining = tot_jobs;
        let mut dump = Dump::new(self.chunk_size() * 5);
        let mut qstat = HashSet::<String>::new();
        setup();
        loop {
            let mut finished = 0;
            while cur_jobs.len() < self.job_limit() && cur < tot_jobs {
                let new_chunk = self.build_chunk(
                    &mut jobs
                        [cur..std::cmp::min(cur + self.chunk_size(), tot_jobs)],
                    chunk_num,
                    &mut slurm_jobs,
                );
                let job_id = new_chunk[0].job_id.clone();
                qstat.insert(job_id);
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
                    ProgramStatus::Success(val) => {
                        to_remove.push(i);
                        dst[job.index] += job.coeff * val;
                        dump.add(job.program.associated_files());
                        finished += 1;
                        remaining -= 1;
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
                    e => {
                        // check if job was in the queue last time we checked,
                        // if not (maybe depending on the type of issue),
                        // resubmit
                        if !qstat.contains(&job.job_id) {
                            println!(
                                "{} not in there with {:?}",
                                job.job_id, e
                            );
                            // TODO actually figure out how to resub - have to
                            // add a field on job to hold the resub and add a
                            // check somewhere in here. I guess I could just
                            // while resub != nil in this arm to traverse resub
                            // linked list
                            todo!()
                        }
                    }
                }
            }
            // have to remove the highest index first so sort and reverse
            to_remove.sort();
            to_remove.reverse();
            for i in to_remove {
                cur_jobs.swap_remove(i);
            }
            if cur_jobs.len() == 0 && cur == tot_jobs {
                takedown();
                return;
            }
            if finished == 0 {
                eprintln!("{} jobs remaining", remaining);
                qstat = self.status();
                thread::sleep(time::Duration::from_secs(
                    self.sleep_int() as u64
                ));
            }
        }
    }
}
