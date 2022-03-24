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
static JOB_LIMIT: usize = 1600;
static CHUNK_SIZE: usize = 128;
static SLEEP_INT: usize = 1;
static DELTA: f64 = 1e-8;
static DELTA_FWD: f64 = 5e7; // 1 / 2Δ
static DELTA_BWD: f64 = -5e7; // -1 / 2Δ
static DEBUG: bool = false;
static HT_TO_CM: f64 = 219_474.5459784;

/// set up the directories needed for the program after deleting existing ones
pub fn setup() {
    takedown();
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

pub fn load_energies(filename: &str) -> Vec<f64> {
    let mut ret = Vec::new();
    let f = match File::open(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "failed to open {filename} for reading energies with {e}"
            );
            std::process::exit(1);
        }
    };
    let lines = BufReader::new(f).lines();
    for line in lines.map(|x| x.unwrap()) {
        ret.push(line.trim().parse().unwrap());
    }
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
    pub coeff: f64,
}

impl Job {
    pub fn new(mopac: mopac::Mopac, index: usize) -> Self {
        Self {
            mopac,
            pbs_file: String::new(),
            job_id: String::new(),
            index,
            coeff: 1.0,
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
pub fn build_jobs(
    moles: &Vec<Vec<Atom>>,
    params: &Vec<Param>,
    start_index: usize,
    coeff: f64,
    job_num: usize,
) -> Vec<Job> {
    let mut count: usize = start_index;
    let mut job_num = job_num;
    let mut jobs = Vec::new();
    for mol in moles {
        let filename = format!("inp/job.{:08}", job_num);
        job_num += 1;
        let mut job =
            Job::new(Mopac::new(filename, params.clone(), mol.clone()), count);
        job.coeff = coeff;
        jobs.push(job);
        count += 1;
    }
    jobs
}

struct Dump {
    buf: Vec<String>,
    ptr: usize,
    max: usize,
}

impl Dump {
    fn new(size: usize) -> Self {
        Self {
            buf: vec![String::new(); size],
            ptr: 0,
            max: size,
        }
    }
    fn add(&mut self, files: Vec<String>) {
        for file in files {
            if self.ptr == self.max - 1 {
                self.dump();
            }
            self.buf[self.ptr] = file;
            self.ptr += 1;
        }
    }

    fn dump(&mut self) {
        for file in &self.buf {
            let _ = std::fs::remove_file(file);
        }
        self.ptr = 0;
    }
}

pub fn drain<S: Submit>(jobs: &mut Vec<Job>, dst: &mut [f64], _submitter: S) {
    let mut chunk_num: usize = 0;
    let mut cur_jobs = Vec::new();
    let mut slurm_jobs = HashMap::new();
    let mut cur = 0; // current index into jobs
    let tot_jobs = jobs.len();
    let mut dump = Dump::new(CHUNK_SIZE * 5);
    loop {
        let mut finished = 0;
        while cur_jobs.len() < JOB_LIMIT && cur < tot_jobs {
            let new_chunk = build_chunk(
                &mut jobs[cur..std::cmp::min(cur + CHUNK_SIZE, tot_jobs)],
                chunk_num,
                &mut slurm_jobs,
                S::new(),
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
            match job.mopac.read_output() {
                Some(val) => {
                    to_remove.push(i);
                    dst[job.index] += job.coeff * val;
                    dump.add(vec![
                        format!("{}.mop", job.mopac.filename),
                        format!("{}.out", job.mopac.filename),
                        format!("{}.arc", job.mopac.filename),
                        format!("{}.aux", job.mopac.filename),
                        job.mopac.param_file.clone(),
                    ]);
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
            dbg!(
                cur_jobs.len(),
                &cur_jobs[0].mopac.filename,
                &cur_jobs[0].pbs_file
            );
            thread::sleep(time::Duration::from_secs(SLEEP_INT as u64));
        } else {
            eprintln!("finished {}", finished);
        }
    }
}

/// compute the semi-empirical energies of `moles` for the given `params`
pub fn semi_empirical(moles: &Vec<Vec<Atom>>, params: &Vec<Param>) -> Vec<f64> {
    let mut jobs = build_jobs(moles, params, 0, 1.0, 0);
    let mut got = vec![0.0; jobs.len()];
    drain(&mut jobs, &mut got, LocalQueue {});
    got
}

/// Compute the numerical Jacobian for the geomeries in `moles` and the
/// parameters in `params`. For convenience of indexing, the transpose is
/// actually computed and returned
pub fn num_jac<S: Submit>(
    moles: &Vec<Vec<Atom>>,
    params: &Vec<Param>,
    submitter: S,
) -> Vec<f64> {
    let rows = params.len();
    let cols = moles.len();
    let mut jac_t = vec![0.; rows * cols];
    // front and back for each row and col
    let mut jobs = Vec::with_capacity(2 * rows * cols);
    let mut job_num = 0;
    for row in 0..rows {
        // loop over params
        let mut pf = params.clone();
        let mut pb = params.clone();
        let idx = row * cols;
        // TODO do I really have to clone here?
        pf[row].value += DELTA;
        let fwd_jobs = build_jobs(&moles.clone(), &pf, idx, DELTA_FWD, job_num);
        job_num += fwd_jobs.len();
        jobs.extend_from_slice(&fwd_jobs);

        pb[row].value -= DELTA;
        let bwd_jobs = build_jobs(&moles.clone(), &pb, idx, DELTA_BWD, job_num);
        job_num += bwd_jobs.len();
        jobs.extend_from_slice(&bwd_jobs);
    }
    if DEBUG {
        eprintln!("num_jac: running {} jobs", jobs.len());
    }
    drain(&mut jobs, &mut jac_t, submitter);
    jac_t
}

#[derive(Debug)]
pub struct Stats {
    pub norm: f64,
    pub rmsd: f64,
    pub max: f64,
}

impl Stats {
    /// compute the Stats between `v` and `w` in cm-1 under the assumption they
    /// started out in Ht
    pub fn new(v: &Vec<f64>, w: &Vec<f64>) -> Self {
        let count = v.len();
        assert_eq!(count, w.len());
        let mut sq_diffs = 0.0;
        let mut max = v[0] - w[0];
        for i in 0..count {
            let diff = v[i] - w[i];
            sq_diffs += diff * diff;
            if diff.abs() > max {
                max = diff.abs();
            }
        }
        Self {
            norm: sq_diffs.sqrt() * HT_TO_CM,
            rmsd: (sq_diffs / count as f64).sqrt() * HT_TO_CM,
            max: max * HT_TO_CM,
        }
    }
    pub fn print_header() {
        println!(
            "{:>17}{:>12}{:>12}{:>12}{:>12}{:>12}",
            "cm-1", "cm-1", "cm-1", "cm-1", "cm-1", "s"
        );
        println!(
            "{:>5}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}",
            "Iter", "Norm", "ΔNorm", "RMSD", "ΔRMSD", "Max", "Time"
        );
    }

    // TODO bring time back in last column
    pub fn print_step(&self, iter: usize, last: &Self) {
        print!(
            "{:5}{:12.4}{:12.4}{:12.4}{:12.4}{:12.4}{:12.1}\n",
            iter,
            self.norm,
            self.norm - last.norm,
            self.rmsd,
            self.rmsd - last.rmsd,
            self.max,
            0.0
        );
    }
}

impl Default for Stats {
    fn default() -> Self {
        Self {
            norm: 0.0,
            rmsd: 0.0,
            max: 0.0,
        }
    }
}

/// return `energies` relative to its minimum element
pub fn relative(energies: &Vec<f64>) -> Vec<f64> {
    let mut ret = Vec::with_capacity(energies.len());
    let mut min = energies[0];
    // compute the min
    for e in energies {
        if e < &min {
            min = *e
        }
    }

    for e in energies {
        ret.push(e - min);
    }
    ret
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

    fn comp_vec(got: &[f64], want: &[f64], eps: f64) -> bool {
        if got.len() != want.len() {
            return false;
        }
        for (i, g) in got.iter().enumerate() {
            let diff = (g - want[i]).abs();
            if diff > eps {
                eprintln!(
                    "{:5} got: {}, wanted {}: diff {:e}",
                    i, g, want[i], diff
                );
                return false;
            }
        }
        true
    }

    fn comp_geoms(got: Vec<Vec<Atom>>, want: Vec<Vec<Atom>>, eps: f64) -> bool {
        for (i, mol) in got.iter().enumerate() {
            for (j, atom) in mol.iter().enumerate() {
                let btom = &want[i][j];
                if atom.label != btom.label {
                    return false;
                }
                if !comp_vec(&atom.coord, &btom.coord, eps) {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn test_load_energies() {
        let got = load_energies("small.dat");
        let want = vec![
            0.000000000000,
            0.000453458157,
            0.000258906149,
            0.000268926371,
            0.000247426244,
            0.000251632099,
            0.000259236055,
            0.000267395682,
            0.000273547427,
            0.000159696616,
            0.000136655279,
            0.000118039859,
            0.000119823282,
            0.000124994640,
            0.000136055791,
            0.000177574442,
            0.000124620943,
            0.000126842659,
            0.000132448769,
            0.000108785467,
            0.000107629675,
            0.000176541210,
            0.000144535949,
            0.000126349619,
            0.000129459053,
            0.000141660093,
            0.000176137923,
            0.000127547307,
            0.000128059983,
        ];
        assert!(comp_vec(&got, &want, 1e-12));
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
    fn test_one_iter() {
        let names = vec!["C", "C", "C", "H", "H"];
        let moles = load_geoms("three07", names);
        let params = load_params("params.dat");
        let want = vec![
            0.20374485388911504,
            0.20541305733965845,
            0.20511337069030972,
        ];
        setup();
        let got = semi_empirical(&moles, &params);
        let eps = 1e-20;
        for i in 0..want.len() {
            assert!((got[i] - want[i]).abs() < eps);
        }
        takedown();
    }

    fn load_mat(filename: &str) -> Vec<f64> {
        let f = match File::open(filename) {
            Ok(f) => f,
            Err(e) => {
                panic!("failed to open {} with {}", filename, e);
            }
        };
        let lines = BufReader::new(f).lines().map(|x| x.unwrap());
        let mut ret = Vec::new();
        for line in lines {
            let sp: Vec<f64> = line
                .trim()
                .split_whitespace()
                .skip(1)
                .map(|x| x.parse().unwrap())
                .collect();
            ret.extend_from_slice(&sp);
        }
        ret
    }

    fn transpose(orig: &[f64], rows: usize, cols: usize) -> Vec<f64> {
        let mut ret = Vec::with_capacity(rows * cols);
        for i in 0..cols {
            for j in 0..rows {
                ret.push(orig[i + j * cols])
            }
        }
        ret
    }

    #[test]
    fn test_num_jac() {
        let names = vec!["C", "C", "C", "H", "H"];
        let moles = load_geoms("small07", names);
        let params = load_params("small.params");
        let want = transpose(&load_mat("small.jac"), moles.len(), params.len());
        setup();
        let got = num_jac(&moles, &params, LocalQueue {});
        assert!(comp_vec(&got, &want, 2e-6));
        takedown();
    }

    #[test]
    fn test_lev_mar() {
        let names = vec!["C", "C", "C", "H", "H"];
        let moles = load_geoms("small07", names);
        let params = load_params("small.params");
        let ai = load_energies("25.dat");
        setup();
        let se = semi_empirical(&moles, &params);
        let rel = relative(&se);
        let mut stats = Stats::new(&ai, &rel);
        let mut last_stats = Stats::default();
        Stats::print_header();
        stats.print_step(0, &last_stats);
        let mut iter = 1;
        while iter < 5 {
	    // let jac_t = num_jac(&moles, &params, LocalQueue{});
	    // TODO write lev_mar and call it with jac_t
            stats = Stats::new(&ai, &rel);
            stats.print_step(iter, &last_stats);
            last_stats = stats;
            iter += 1;
        }
        takedown();
    }
}
