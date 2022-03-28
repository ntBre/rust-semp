use std::{
    clone::Clone,
    fs::File,
    io::BufReader,
    io::{BufRead, Write},
};

use nalgebra as na;

use mopac::{Mopac, Params};
use queue::{Job, Queue};
use stats::Stats;

pub mod config;
pub mod dump;
pub mod local;
pub mod mopac;
pub mod queue;
pub mod slurm;
pub mod stats;

static DELTA: f64 = 1e-8;
static DELTA_FWD: f64 = 5e7; // 1 / 2Δ
static DELTA_BWD: f64 = -5e7; // -1 / 2Δ
static DEBUG: bool = false;

pub static LAMBDA0: f64 = 1e-8;
pub static NU: f64 = 2.0;
pub static MAX_TRIES: usize = 5;

/// from [StackOverflow](https://stackoverflow.com/a/45145246)
#[macro_export]
macro_rules! string {
    // match a list of expressions separated by comma:
    ($($str:expr),*) => ({
        // create a Vec with this list of expressions,
        // calling String::from on each:
        vec![$(String::from($str),)*] as Vec<String>
    });
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
pub fn load_geoms(filename: &str, atom_names: Vec<String>) -> Vec<Vec<Atom>> {
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

pub fn load_energies(filename: &str) -> na::DVector<f64> {
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
    na::DVector::from(ret)
}

/// Read `filename` into a string and then call [`parse_params`]
pub fn load_params(filename: &str) -> Params {
    let params = match std::fs::read_to_string(filename) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("failed to read {filename} with {e}");
            std::process::exit(1);
        }
    };
    parse_params(&params)
}

/// parse a string containing lines like:
///
/// ```text
///   USS            H    -11.246958000000
///   ZS             H      1.268641000000
/// ```
/// into a vec of Params
pub fn parse_params(params: &str) -> Params {
    let lines = params.split('\n');
    let mut names = Vec::new();
    let mut atoms = Vec::new();
    let mut values = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() == 3 {
            names.push(fields[0].to_string());
            atoms.push(fields[1].to_string());
            values.push(fields[2].parse().unwrap());
        }
    }
    Params::from(names, atoms, values)
}

/// Build the jobs described by `moles` in memory, but don't write any of their
/// files yet
pub fn build_jobs(
    moles: &Vec<Vec<Atom>>,
    params: &Params,
    start_index: usize,
    coeff: f64,
    job_num: usize,
) -> Vec<Job<Mopac>> {
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

/// compute the semi-empirical energies of `moles` for the given `params`
pub fn semi_empirical<Q: Queue<Mopac>>(
    moles: &Vec<Vec<Atom>>,
    params: &Params,
    submitter: &Q,
) -> na::DVector<f64> {
    let mut jobs = build_jobs(moles, params, 0, 1.0, 0);
    let mut got = vec![0.0; jobs.len()];
    submitter.drain(&mut jobs, &mut got);
    na::DVector::from(got)
}

/// Compute the numerical Jacobian for the geomeries in `moles` and the
/// parameters in `params`. For convenience of indexing, the transpose is
/// actually computed and returned
pub fn num_jac<Q: Queue<Mopac>>(
    moles: &Vec<Vec<Atom>>,
    params: &Params,
    submitter: &Q,
) -> na::DMatrix<f64> {
    let rows = params.values.len();
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
        pf.values[row] += DELTA;
        let fwd_jobs = build_jobs(&moles.clone(), &pf, idx, DELTA_FWD, job_num);
        job_num += fwd_jobs.len();
        jobs.extend_from_slice(&fwd_jobs);

        pb.values[row] -= DELTA;
        let bwd_jobs = build_jobs(&moles.clone(), &pb, idx, DELTA_BWD, job_num);
        job_num += bwd_jobs.len();
        jobs.extend_from_slice(&bwd_jobs);
    }
    if DEBUG {
        eprintln!("num_jac: running {} jobs", jobs.len());
    }
    submitter.drain(&mut jobs, &mut jac_t);
    // nalgebra does from_vec in col-major order, so lead with cols and I get
    // jac_t_t or jac back
    let jac = na::DMatrix::from_vec(cols, rows, jac_t);
    jac
}

/// Solve (JᵀJ + λI)δ = Jᵀ[y - f(β)] for δ. y is the vector of "true" training
/// energies, and f(β) represents the current semi-empirical energies.
pub fn lev_mar(
    jac: &na::DMatrix<f64>,
    ai: &na::DVector<f64>,
    se: &na::DVector<f64>,
    lambda: f64,
) -> na::DVector<f64> {
    let jac_t = jac.transpose();
    let (rows, _) = jac_t.shape();
    // compute scaled matrix, A*, from JᵀJ
    let a = &jac_t * jac;
    let mut a_star = na::DMatrix::from_vec(rows, rows, vec![0.0; rows * rows]);
    for i in 0..rows {
        for j in 0..rows {
            a_star[(i, j)] = a[(i, j)] / (a[(i, i)].sqrt() * a[(j, j)].sqrt())
        }
    }
    // compute scaled right side, g*, from b
    let b = {
        let b = &jac_t * (ai - se);
        let mut g = na::DVector::from(vec![0.0; rows]);
        for j in 0..rows {
            g[j] = b[j] / a[(j, j)].sqrt()
        }
        g
    };
    // compute λI
    let li = {
        let mut i = na::DMatrix::<f64>::identity(rows, rows);
        i.scale_mut(lambda);
        i
    };
    // add A* to λI (left-hand side) and compute the Cholesky decomposition
    let lhs = a_star + li;
    let lhs = match na::linalg::Cholesky::new(lhs) {
        Some(a) => a,
        None => {
            eprintln!("cholesky decomposition failed");
            std::process::exit(1);
        }
    };
    let mut d = lhs.solve(&b);
    // convert back from δ* to δ
    for j in 0..rows {
        d[j] /= a[(j, j)].sqrt()
    }
    d
}

/// Approximately update the Jacobian matrix using Broyden's method:
///
/// Jₙ = Jₙ₋₁ + /// (Δfₙ - Jₙ₋₁Δxₙ)Δxₙᵀ / ||Δxₙ||², where Jₙ is the updated
/// Jacobian, Jₙ₋₁ is the Jacobian from the previous iteration, Δf is the change
/// in the function value between iterations (the new semi-empirical energies
/// minus the old), and Δx is the step (δ) determined by the last lev_mar
/// iteration
pub fn broyden_update(
    jac: &na::DMatrix<f64>,
    se_old: &na::DVector<f64>,
    se_new: &na::DVector<f64>,
    step: &na::DVector<f64>,
) -> na::DMatrix<f64> {
    let df = se_new - se_old;
    let numerator = (df - jac * step) * step.transpose();
    let denominator = 1.0 / step.norm();
    let update = numerator.scale(denominator);
    jac + update
}

/// return `energies` relative to its minimum element
pub fn relative(energies: &na::DVector<f64>) -> na::DVector<f64> {
    let min = energies.min();
    let min = na::DVector::from(vec![min; energies.len()]);
    let ret = energies.clone();
    ret - min
}

fn log_params<W: Write>(w: &mut W, iter: usize, params: &Params) {
    let _ = writeln!(w, "Iter {}", iter);
    let _ = writeln!(w, "{}", params.to_string());
}

pub fn run_algo<Q: Queue<Mopac>, W: Write>(
    param_log: &mut W,
    atom_names: Vec<String>,
    geom_file: &str,
    params: Params,
    energy_file: &str,
    max_iter: usize,
    broyden: bool,
    broyd_int: usize,
    queue: Q,
) -> Stats {
    let moles = load_geoms(geom_file, atom_names);
    let mut params = params.clone();
    log_params(param_log, 0, &params);
    let ai = load_energies(energy_file);
    // initial semi-empirical energies and stats
    let mut se = semi_empirical(&moles, &params, &queue);
    let mut old_se = se.clone();
    let rel = relative(&se);
    let mut stats = Stats::new(&ai, &rel);
    let mut last_stats = Stats::default();
    Stats::print_header();
    stats.print_step(0, &last_stats, 0);
    last_stats = stats;
    // start looping
    let mut iter = 1;
    let mut lambda = LAMBDA0;
    let mut del_norm: f64 = 1.0;
    let mut start = std::time::SystemTime::now();
    let mut jac = num_jac(&moles, &params, &queue);
    // have to "initialize" this to satisfy compiler, but any use should panic
    // since it has zero length
    let mut step = na::DVector::from(vec![]);
    while iter <= max_iter && del_norm.abs() > 1e-4 {
        if broyden && iter > 1 && iter % broyd_int != 1 {
            eprintln!("broyden on iter {}", iter);
            start = std::time::SystemTime::now();
            jac = broyden_update(&jac, &old_se, &se, &step);
        } else if iter > 1 {
            start = std::time::SystemTime::now();
            jac = num_jac(&moles, &params, &queue);
        } // else (first iteration) use jac from outside loop
        lambda /= NU;
        // BEGIN copy-paste
        step = lev_mar(&jac, &ai, &se, lambda);
        let try_params = Params::new(
            params.names.clone(),
            params.atoms.clone(),
            &params.values + &step,
        );
        let new_se = semi_empirical(&moles, &try_params, &queue);
        let rel = relative(&new_se);
        stats = Stats::new(&ai, &rel);
        // END copy-paste

        // cases ii. and iii. from Marquardt63; first iteration is case ii.
        let mut i = 0;
        let mut bad = false;
        while stats.norm > last_stats.norm {
            lambda *= NU;
            let dnorm = stats.norm - last_stats.norm;

            // BCP
            step = lev_mar(&jac, &ai, &se, lambda);
            let try_params = Params::new(
                params.names.clone(),
                params.atoms.clone(),
                &params.values + &step,
            );
            let new_se = semi_empirical(&moles, &try_params, &queue);
            let rel = relative(&new_se);
            stats = Stats::new(&ai, &rel);
            // ECP

            eprintln!(
                "\tλ_{} to {:e} with ΔNorm = {}",
                i,
                lambda,
                stats.norm - last_stats.norm
            );
            i += 1;
            if stats.norm - last_stats.norm > dnorm || i > MAX_TRIES {
                bad = true;
                break;
            }
        }
        let mut k = 1.0;
        let mut i = 2;
        while bad && stats.norm > last_stats.norm && k > 1e-14 {
            k = 1.0 / 10.0_f64.powf(i as f64);

            // BCP
            step = lev_mar(&jac, &ai, &se, lambda);
            let try_params = Params::new(
                params.names.clone(),
                params.atoms.clone(),
                &params.values + k * &step,
            );
            let new_se = semi_empirical(&moles, &try_params, &queue);
            let rel = relative(&new_se);
            stats = Stats::new(&ai, &rel);
            // ECP

            eprintln!(
                "\tk_{} to {:e} with ΔNorm = {}",
                i,
                k,
                stats.norm - last_stats.norm
            );
            i += 1;
        }
        let time = if let Ok(elapsed) = start.elapsed() {
            elapsed.as_millis()
        } else {
            0
        };
        stats.print_step(iter, &last_stats, time);
        // end of loop updates
        old_se = se;
        se = new_se;
        params = try_params;
        log_params(param_log, iter, &params);
        del_norm = stats.norm - last_stats.norm;
        last_stats = stats;
        iter += 1;
    }
    println!("\nconvergence reached after {} iterations", iter);
    last_stats
}

#[cfg(test)]
mod tests {
    use std::fs;

    use crate::{local::LocalQueue, slurm::Slurm, stats::Stats};

    use super::*;

    #[test]
    fn test_load_geoms() {
        let got =
            load_geoms("test_files/three07", string!["C", "C", "C", "H", "H"]);
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
        let got = load_energies("test_files/small.dat");
        let want = na::DVector::from(vec![
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
        ]);
        assert!(comp_dvec(got, want, 1e-12));
    }

    #[test]
    fn test_load_params() {
        let got = load_params("test_files/params.dat");
        #[rustfmt::skip]
	let want = Params::from_literal(
            vec![
                "USS", "ZS", "BETAS", "GSS", "USS", "UPP", "ZS", "ZP", "BETAS",
                "BETAP", "GSS", "GPP", "GSP", "GP2", "HSP",
            ],
            vec![
                "H", "H", "H", "H", "C", "C", "C", "C", "C", "C", "C", "C",
                "C", "C", "C",
            ],
            vec![
                -11.246958000000, 1.268641000000, -8.352984000000,
                14.448686000000, -51.089653000000, -39.937920000000,
                2.047558000000, 1.702841000000, -15.385236000000,
                -7.471929000000, 13.335519000000, 10.778326000000,
                11.528134000000, 9.486212000000, 0.717322000000,
            ],
        );
        assert_eq!(got, want);
    }

    #[test]
    fn test_write_submit_script() {
        Slurm::default().write_submit_script(
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
        let names = string!["C", "C", "C", "H", "H"];
        let moles = load_geoms("test_files/three07", names);
        let params = load_params("test_files/params.dat");
        let want = na::DVector::from(vec![
            0.20374485388911504,
            0.20541305733965845,
            0.20511337069030972,
        ]);
        let got = semi_empirical(&moles, &params, &LocalQueue);
        let eps = 1e-14;
        assert!(comp_dvec(got, want, eps));
    }

    fn load_mat(filename: &str) -> na::DMatrix<f64> {
        let f = match File::open(filename) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("failed to open {} with {}", filename, e);
                std::process::exit(1);
            }
        };
        let lines = BufReader::new(f).lines().map(|x| x.unwrap());
        let mut ret = Vec::new();
        let mut rows = 0;
        for line in lines {
            let sp: Vec<f64> = line
                .trim()
                .split_whitespace()
                .skip(1)
                .map(|x| x.parse().unwrap())
                .collect();
            ret.extend_from_slice(&sp);
            rows += 1;
        }
        let cols = ret.len() / rows;
        let ret = transpose(&ret, rows, cols);
        na::DMatrix::from_vec(rows, cols, ret)
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

    /// report whether or not the euclidean norm of the difference between `got`
    /// and `want` is less than epsilon
    fn comp_dvec(
        got: na::DVector<f64>,
        want: na::DVector<f64>,
        eps: f64,
    ) -> bool {
        let norm = (got - want).norm();
        if norm < eps {
            return true;
        }
        eprintln!("comp_vec: norm = {}", norm);
        false
    }

    fn comp_mat(
        got: na::DMatrix<f64>,
        want: na::DMatrix<f64>,
        eps: f64,
    ) -> bool {
        let norm = (got - want).norm();
        if norm < eps {
            return true;
        }
        eprintln!("comp_mat: norm = {}", norm);
        false
    }

    #[ignore]
    #[test]
    fn test_num_jac() {
        let names = string!["C", "C", "C", "H", "H"];
        let moles = load_geoms("test_files/small07", names);
        let params = load_params("test_files/small.params");
        let want = load_mat("test_files/small.jac");
        let got = num_jac(&moles, &params, &LocalQueue);
        assert!(comp_mat(got, want, 1e-5));
    }

    #[ignore]
    #[test]
    fn test_lev_mar() {
        // loading everything
        let names = string!["C", "C", "C", "H", "H"];
        let queue = LocalQueue;
        let geom_file = "test_files/small07";
        let param_file = "test_files/small.params";
        let energy_file = "test_files/25.dat";
        let got = run_algo(
            &mut std::io::sink(),
            names,
            geom_file,
            load_params(param_file),
            energy_file,
            5,
            true,
            5,
            queue,
        );
        let want = Stats {
            norm: 16.1347,
            rmsd: 3.2269,
            max: 9.2995,
        };
        assert_eq!(got, want);
    }
    /*
       current output:
    Iter        Norm       ΔNorm        RMSD       ΔRMSD         Max        Time
    0    828.6919    828.6919    165.7384    165.7384    266.6057         0.0
    1    325.2037   -503.4882     65.0407   -100.6976    126.9844        40.1

    after scaling:
    1    389.9481   -438.7438     77.9896    -87.7488    201.6500        34.2

    after fixing bounds and <= 5:
    1     70.7093   -757.9826     14.1419   -151.5965     27.7492        29.3
    2     28.7168    -41.9925      5.7434     -8.3985      8.0413        29.7
    3     22.7316     -5.9852      4.5463     -1.1970      7.0542        29.8
    4     22.8138      0.0823      4.5628      0.0165      6.2959        29.7
    5     21.3078     -1.5061      4.2616     -0.3012      6.1055        30.6

    this is weirdly better than the Go version, not sure it should be

    after rest of lev mar:
    1     70.5271   -758.1648     14.1054   -151.6330     27.1241        28.2
    2     27.8114    -42.7157      5.5623     -8.5431      9.7420        29.2
    3     21.9563     -5.8551      4.3913     -1.1710      7.0847        31.8
    4     19.8039     -2.1523      3.9608     -0.4305      5.4642        30.5
    5     19.1505     -0.6535      3.8301     -0.1307      5.3263        30.6

    after broyden:
    1     70.5271   -758.1648     14.1054   -151.6330     27.1241        26.7
    2     66.8706     -3.6565     13.3741     -0.7313     29.8608         1.5
    3     55.3589    -11.5117     11.0718     -2.3023     23.5758         4.3
    4    246.9214    191.5625     49.3843     38.3125     97.3837        12.5
    5     16.1347   -230.7867      3.2269    -46.1573      9.2995         0.9

    absolutely disastrous 4th iteration, but it recovers very nicely on the 5th,
    and it does all four of its iterations faster than one more numerical
    Jacobian so I guess it's worth it

    also converges on iteration 8 after another num_jac on 6, so I'd say it's
    working well

    6     15.9826     -0.1521      3.1965     -0.0304      9.2602        34.9
    7     14.2220     -1.7606      2.8444     -0.3521      5.4632         0.8
    8     14.2220     -0.0000      2.8444     -0.0000      5.4632         9.2

     */
}
