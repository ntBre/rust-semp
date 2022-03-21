use std::{
    collections::HashMap,
    fs::{self, File},
    io::BufRead,
    io::BufReader,
    io::Write,
    thread, time,
};

use mopac::{Mopac, Param};
use queue::Submit;

pub mod mopac;
pub mod queue;

static DIRS: &'static [&str] = &["inp", "tmparam"];
static JOB_LIMIT: usize = 1024;
static CHUNK_SIZE: usize = 128;
static SLEEP_INT: usize = 1;

/// set up the directories needed for the program
pub fn setup() {
    for dir in DIRS {
        match fs::create_dir(dir) {
            Ok(_) => (),
            Err(_) => (),
        }
    }
}

// TODO don't do this if -nodel flag
pub fn takedown() {
    for dir in DIRS {
        match fs::remove_dir_all(dir) {
            Ok(_) => (),
            Err(_) => (),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub label: String,
    pub coord: Vec<f64>,
}

impl Atom {
    pub fn new(label: &str, coord: Vec<f64>) -> Self {
        Self {
            label: label.to_string(),
            coord,
        }
    }
}

impl ToString for Atom {
    fn to_string(&self) -> String {
        format!(
            "{:2} {:15.10} {:15.10} {:15.10}",
            self.label, self.coord[0], self.coord[1], self.coord[2]
        )
    }
}

pub fn geom_string(geom: &Vec<Atom>) -> String {
    use std::fmt::Write;
    let mut ret = String::new();
    for g in geom {
        write!(ret, "{}\n", g.to_string()).unwrap();
    }
    ret
}

/// Take an INTDER-style `file07` file and parse it into a Vec of geometries
pub fn load_geoms(filename: &str, atom_names: Vec<&str>) -> Vec<Vec<Atom>> {
    let f = match File::open(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("failed to open {filename} with {e}");
            std::process::exit(1);
        }
    };
    let lines = BufReader::new(f).lines();
    let mut ret = Vec::new();
    let mut buf = Vec::new();
    let mut idx: usize = 0;
    for (i, line) in lines.map(|x| x.unwrap()).enumerate() {
        if line.contains("# GEOM") {
            if i > 0 {
                ret.push(buf);
                buf = Vec::new();
                idx = 0;
            }
            continue;
        }
        let fields = line
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect();
        buf.push(Atom {
            label: atom_names[idx].to_string(),
            coord: fields,
        });
        idx += 1;
    }
    ret.push(buf);
    ret
}

/// parse a file containing lines like:
///
///   USS            H    -11.246958000000
///   ZS             H      1.268641000000
///
/// into a vec of Params
pub fn load_params(filename: &str) -> Vec<mopac::Param> {
    let mut ret = Vec::new();
    let f = match File::open(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "failed to open {filename} for reading parameters with {e}"
            );
            std::process::exit(1);
        }
    };
    let lines = BufReader::new(f).lines();
    for line in lines.map(|x| x.unwrap()) {
        let fields: Vec<&str> = line.split_whitespace().collect();
        ret.push(mopac::Param::new(
            fields[0],
            fields[1],
            fields[2].parse().unwrap(),
        ))
    }
    ret
}

/// Slurm is a type for holding the information for submitting a slurm job.
/// `filename` is the name of the Slurm submission script
pub struct Slurm;

impl queue::Submit for Slurm {
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

    fn new() -> Self {
        Slurm {}
    }
}

/// Minimal implementation for testing MOPAC locally
pub struct LocalQueue;

impl queue::Submit for LocalQueue {
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

    fn new() -> Self {
        Self {}
    }
}

#[derive(Debug, Clone)]
pub struct Job {
    pub mopac: mopac::Mopac,
    pub pbs_file: String,
    pub job_id: String,
    pub index: usize,
}

impl Job {
    pub fn new(mopac: mopac::Mopac, index: usize) -> Self {
        Self {
            mopac,
            pbs_file: String::new(),
            job_id: String::new(),
            index,
        }
    }
}

/// Build a chunk of jobs by writing the MOPAC input file and the corresponding
/// submission script and submitting the script
pub fn build_chunk<'a, S: Submit>(
    jobs: &'a [Job],
    chunk_num: usize,
    slurm_jobs: &'a mut HashMap<String, usize>,
    queue: S,
) -> Vec<Job> {
    let queue_file = format!("inp/main{}.slurm", chunk_num);
    let mut chunk_jobs = Vec::new();
    for job in jobs {
        let mut job = (*job).clone();
        job.mopac.write_input();
        job.pbs_file = queue_file.to_string();
        chunk_jobs.push(job);
    }
    slurm_jobs.insert(queue_file.clone(), chunk_jobs.len());
    queue.write_submit_script(
        chunk_jobs
            .iter()
            .map(|j| j.mopac.filename.clone())
            .collect(),
        &queue_file,
    );
    // run jobs
    let job_id = queue.submit(&queue_file);
    for job in &mut chunk_jobs {
        job.job_id = job_id.clone();
    }
    chunk_jobs
}

/// Build the jobs described by `moles` in memory, but don't write any of their
/// files yet
pub fn build_jobs(moles: Vec<Vec<Atom>>, params: Vec<Param>) -> Jobs {
    let mut count: usize = 0;
    let mut jobs = Jobs::new();
    for mol in moles {
        let filename = format!("inp/job.{count:08}");
        jobs.jobs
            .push(Job::new(Mopac::new(filename, params.clone(), mol), count));
        count += 1;
    }
    jobs.dst = vec![0.0; count];
    jobs
}

pub fn drain<S: Submit>(jobs: &mut Jobs, _submitter: S) {
    let mut chunk_num: usize = 0;
    let mut cur_jobs = Vec::new();
    let mut slurm_jobs = HashMap::new();
    let mut cur = 0; // current index into jobs
    let tot_jobs = jobs.jobs.len();
    loop {
        if cur_jobs.len() < JOB_LIMIT {
            let new_chunk = build_chunk(
                &mut jobs.jobs[cur..std::cmp::min(cur + CHUNK_SIZE, tot_jobs)],
                chunk_num,
                &mut slurm_jobs,
                S::new(),
            );
            chunk_num += 1;
            cur += new_chunk.len();
            cur_jobs.extend(new_chunk);
        }
        // collect output
        let mut to_remove = Vec::new();
        for (i, job) in cur_jobs.iter_mut().enumerate() {
            match job.mopac.read_output() {
                Some(val) => {
                    to_remove.push(i);
                    jobs.dst[job.index] = val;
                }
                None => (),
            }
            let count = slurm_jobs.get_mut(job.pbs_file.as_str()).unwrap();
            *count -= 1;
        }
        // remove finished jobs. TODO also delete the files

        // have to remove the highest index first so sort and reverse
        to_remove.sort();
        to_remove.reverse();
        for i in to_remove {
            cur_jobs.swap_remove(i);
        }
        if cur_jobs.len() == 0 && cur == tot_jobs {
            return;
        }
        thread::sleep(time::Duration::from_secs(SLEEP_INT as u64));
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use crate::mopac::*;

    #[test]
    fn test_load_geoms() {
        let got = load_geoms("three07", vec!["C", "C", "C", "H", "H"]);
        let want = vec![
            vec![
                Atom::new("C", vec![0.0000000000, 0.0000000000, -1.6794733900]),
                Atom::new("C", vec![0.0000000000, 1.2524327590, 0.6959098120]),
                Atom::new("C", vec![0.0000000000, -1.2524327590, 0.6959098120]),
                Atom::new("H", vec![0.0000000000, 3.0146272390, 1.7138963510]),
                Atom::new("H", vec![0.0000000000, -3.0146272390, 1.7138963510]),
            ],
            vec![
                Atom::new("C", vec![0.0000000000, 0.0000000000, -1.6929508795]),
                Atom::new("C", vec![0.0000000000, 1.2335354991, 0.6923003326]),
                Atom::new("C", vec![0.0000000000, -1.2335354991, 0.6923003326]),
                Atom::new("H", vec![0.0000000000, 2.9875928126, 1.7242445752]),
                Atom::new("H", vec![0.0000000000, -2.9875928126, 1.7242445752]),
            ],
            vec![
                Atom::new("C", vec![0.0000000000, 0.0000000000, -1.6826636598]),
                Atom::new("C", vec![0.0000000000, 1.2382598141, 0.6926064201]),
                Atom::new("C", vec![0.0000000000, -1.2382598141, 0.6926064201]),
                Atom::new("H", vec![0.0000000000, 2.9956906742, 1.7187948778]),
                Atom::new("H", vec![0.0000000000, -2.9956906742, 1.7187948778]),
            ],
        ];
        assert!(comp_geoms(got, want, 1e-10));
    }

    fn comp_geoms(got: Vec<Vec<Atom>>, want: Vec<Vec<Atom>>, eps: f64) -> bool {
        for (i, mol) in got.iter().enumerate() {
            for (j, atom) in mol.iter().enumerate() {
                let btom = &want[i][j];
                if atom.label != btom.label {
                    return false;
                }
                for k in 0..atom.coord.len() {
                    if (atom.coord[k] - btom.coord[k]).abs() > eps {
                        return false;
                    }
                }
            }
        }
        true
    }

    #[test]
    fn test_load_params() {
        let got = load_params("params.dat");
        let want = vec![
            Param::new("USS", "H", -11.246958000000),
            Param::new("ZS", "H", 1.268641000000),
            Param::new("BETAS", "H", -8.352984000000),
            Param::new("GSS", "H", 14.448686000000),
            Param::new("USS", "C", -51.089653000000),
            Param::new("UPP", "C", -39.937920000000),
            Param::new("ZS", "C", 2.047558000000),
            Param::new("ZP", "C", 1.702841000000),
            Param::new("BETAS", "C", -15.385236000000),
            Param::new("BETAP", "C", -7.471929000000),
            Param::new("GSS", "C", 13.335519000000),
            Param::new("GPP", "C", 10.778326000000),
            Param::new("GSP", "C", 11.528134000000),
            Param::new("GP2", "C", 9.486212000000),
            Param::new("HSP", "C", 0.717322000000),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn test_write_submit_script() {
        Slurm {}.write_submit_script(
            vec![
                String::from("input1"),
                String::from("input2"),
                String::from("input3"),
            ],
            "/tmp/submit.slurm",
        );
        let got = fs::read_to_string("/tmp/submit.slurm")
            .expect("failed to read /tmp/submit.slurm");
        let want = "#!/bin/bash
#SBATCH --job-name=semp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /tmp/submit.slurm.out
#SBATCH --no-requeue
#SBATCH --mem=1gb
export LD_LIBRARY_PATH=/home/qc/mopac2016/
/home/qc/mopac2016/MOPAC2016.exe input1.mop
/home/qc/mopac2016/MOPAC2016.exe input2.mop
/home/qc/mopac2016/MOPAC2016.exe input3.mop
";
        assert_eq!(got, want);
        fs::remove_file("/tmp/submit.slurm").unwrap();
    }

    #[test]
    fn test_full() {
        let names = vec!["C", "C", "C", "H", "H"];
        let moles = load_geoms("three07", names);
        let params = load_params("params.dat");
        // build jobs (in memory) -> write jobs (to disk) -> run jobs
        let mut jobs = build_jobs(moles, params);
        setup();
        drain(&mut jobs, LocalQueue {});
        let want = vec![
            0.20374485388911504,
            0.20541305733965845,
            0.20511337069030972,
        ];
        let got = jobs.dst;
        let eps = 1e-20;
        for i in 0..want.len() {
            assert!((got[i] - want[i]).abs() < eps);
        }
        takedown();
    }
}

pub struct Jobs {
    pub jobs: Vec<Job>,
    pub dst: Vec<f64>,
}

impl Jobs {
    pub fn new() -> Self {
        Jobs {
            jobs: Vec::new(),
            dst: Vec::new(),
        }
    }
}
