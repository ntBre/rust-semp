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
/// convergence threshold for the change in the norm, rmsd, and max
const DCONV_THRSH: f64 = 1e-5;

pub static LAMBDA0: f64 = 1e-8;
pub static NU: f64 = 2.0;
pub static MAX_TRIES: usize = 10;

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
pub fn load_geoms(filename: &str, atom_names: &[String]) -> Vec<Vec<Atom>> {
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
    charge: isize,
) -> Vec<Job<Mopac>> {
    let mut count: usize = start_index;
    let mut job_num = job_num;
    let mut jobs = Vec::new();
    for mol in moles {
        let filename = format!("inp/job.{:08}", job_num);
        job_num += 1;
        let mut job = Job::new(
            Mopac::new(filename, params.clone(), mol.clone(), charge),
            count,
        );
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
    charge: isize,
) -> na::DVector<f64> {
    let mut jobs = build_jobs(moles, params, 0, 1.0, 0, charge);
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
    charge: isize,
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
        let fwd_jobs =
            build_jobs(&moles.clone(), &pf, idx, DELTA_FWD, job_num, charge);
        job_num += fwd_jobs.len();
        jobs.extend_from_slice(&fwd_jobs);

        pb.values[row] -= DELTA;
        let bwd_jobs =
            build_jobs(&moles.clone(), &pb, idx, DELTA_BWD, job_num, charge);
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
    let g_star = {
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
            eprintln!("jacobian: \n{:12.8}", jac);
            std::process::exit(1);
        }
    };
    let mut d = lhs.solve(&g_star);
    // convert back from δ* to δ
    for j in 0..rows {
        d[j] /= a[(j, j)].sqrt()
    }
    d
}

/// Approximately update the Jacobian matrix using Broyden's method:
///
/// Jₙ = Jₙ₋₁ + (Δfₙ - Jₙ₋₁Δxₙ)Δxₙᵀ / ||Δxₙ||², where Jₙ is the updated Jacobian,
/// Jₙ₋₁ is the Jacobian from the previous iteration, Δf is the change in the
/// function value between iterations (the new semi-empirical energies minus the
/// old), and Δx is the step (δ) determined by the last lev_mar iteration
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

pub fn dump_vec<W: Write>(w: &mut W, vec: &na::DVector<f64>) {
    for (i, v) in vec.iter().enumerate() {
        writeln!(w, "{:>5}{:>20.12}", i, v).unwrap();
    }
}

pub fn dump_mat<W: Write>(w: &mut W, mat: &na::DMatrix<f64>) {
    let (rows, cols) = mat.shape();
    writeln!(w).unwrap();
    for i in 0..rows {
        write!(w, "{:>5}", i).unwrap();
        for j in 0..cols {
            write!(w, "{:>12.8}", mat[(i, j)]).unwrap();
        }
        writeln!(w).unwrap();
    }
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
    charge: isize,
) -> Stats {
    let moles = load_geoms(geom_file, &atom_names);
    let mut params = params.clone();
    log_params(param_log, 0, &params);
    let ai = load_energies(energy_file);
    // initial semi-empirical energies and stats
    let mut se = relative(&semi_empirical(&moles, &params, &queue, charge));
    let mut old_se = se.clone();
    let mut stats = Stats::new(&ai, &se);
    let mut last_stats = Stats::default();
    Stats::print_header();
    stats.print_step(0, &last_stats, 0);
    last_stats = stats;
    // start looping
    let mut iter = 1;
    let mut lambda = LAMBDA0;
    let mut del_norm: f64 = 1.0;
    let mut del_rmsd: f64 = 1.0;
    let mut del_max: f64 = 1.0;
    let mut in_broyden = false;
    let mut need_num_jac = false;
    let mut start = std::time::SystemTime::now();
    let mut jac = num_jac(&moles, &params, &queue, charge);
    // have to "initialize" this to satisfy compiler, but any use should panic
    // since it has zero length
    let mut step = na::DVector::from(vec![]);
    while iter <= max_iter
        && (del_norm.abs() > DCONV_THRSH
            || del_rmsd.abs() > DCONV_THRSH
            || del_max.abs() > DCONV_THRSH)
    {
        if broyden && !need_num_jac && iter > 1 && iter % broyd_int != 1 {
            eprintln!("broyden on iter {}", iter);
            in_broyden = true;
            start = std::time::SystemTime::now();
            jac = broyden_update(&jac, &old_se, &se, &step);
        } else if iter > 1 {
            in_broyden = false;
            need_num_jac = false;
            start = std::time::SystemTime::now();
            jac = num_jac(&moles, &params, &queue, charge);
        } // else (first iteration) use jac from outside loop
        let lambda_init = lambda;

        // first try with λ/ν
        lambda /= NU;
        step = lev_mar(&jac, &ai, &se, lambda);
        let mut try_params = Params::new(
            params.names.clone(),
            params.atoms.clone(),
            &params.values + &step,
        );
        let mut new_se =
            relative(&semi_empirical(&moles, &try_params, &queue, charge));
        stats = Stats::new(&ai, &new_se);

        // cases ii. and iii. from Marquardt63; first iteration is case ii.
        // increase λ*ν until the norm improves or we hit MAX_TRIES or Δnorm
        // increases
        let mut i = 0;
        let mut bad = false;
        while stats.norm > last_stats.norm {
            eprintln!(
                "\tλ_{} to {:e} with ΔNorm = {}",
                i,
                lambda,
                stats.norm - last_stats.norm
            );
            lambda *= NU;
            let dnorm = stats.norm - last_stats.norm;

            step = lev_mar(&jac, &ai, &se, lambda);
            try_params = Params::new(
                params.names.clone(),
                params.atoms.clone(),
                &params.values + &step,
            );
            new_se =
                relative(&semi_empirical(&moles, &try_params, &queue, charge));
            stats = Stats::new(&ai, &new_se);

            i += 1;
            if stats.norm - last_stats.norm > dnorm || i > MAX_TRIES {
                bad = true;
                break;
            }
        }

        // adjusting λ failed to decrease the norm. try decreasing the step size
        // K until the norm improves or k ≈ 0 or Δnorm increases
        let mut k = 1.0;
        let mut i = 2;
        while bad && stats.norm > last_stats.norm && k > 1e-14 {
            eprintln!(
                "\tk_{} to {:e} with ΔNorm = {}",
                i - 2,
                k,
                stats.norm - last_stats.norm
            );
            k = 1.0 / 10.0_f64.powf(i as f64);
            let dnorm = stats.norm - last_stats.norm;

            step = lev_mar(&jac, &ai, &se, lambda);
            try_params = Params::new(
                params.names.clone(),
                params.atoms.clone(),
                &params.values + k * &step,
            );
            new_se =
                relative(&semi_empirical(&moles, &try_params, &queue, charge));
            stats = Stats::new(&ai, &new_se);

            if stats.norm - last_stats.norm > dnorm {
                break;
            }
            i += 1;
        }

        // log the time
        let time = if let Ok(elapsed) = start.elapsed() {
            elapsed.as_millis()
        } else {
            0
        };

        // don't accept a bad step from broyden, do numjac
        if in_broyden && stats.norm > last_stats.norm {
            eprintln!(
                "bad broyden step, doing num_jac after {:.1} sec",
                time as f64 / 1000.
            );
            need_num_jac = true;
            lambda = lambda_init;
            continue;
        }

        // end of loop updates
        stats.print_step(iter, &last_stats, time);
        old_se = se;
        se = new_se;
        params = try_params;
        log_params(param_log, iter, &params);
        del_norm = stats.norm - last_stats.norm;
        del_rmsd = stats.rmsd - last_stats.rmsd;
        del_max = stats.max - last_stats.max;
        last_stats = stats;
        iter += 1;
    }
    println!("\nconvergence reached after {} iterations", iter - 1);
    println!("DNORM = {:.2e}", del_norm.abs());
    println!("DRMSD = {:.2e}", del_rmsd.abs());
    println!("DMAX = {:.2e}", del_max.abs());
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
            load_geoms("test_files/three07", &string!["C", "C", "C", "H", "H"]);
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
            &string!["input1", "input2", "input3"],
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
        let moles = load_geoms("test_files/three07", &names);
        let params = load_params("test_files/params.dat");
        let want = na::DVector::from(vec![
            0.20374485388911504,
            0.20541305733965845,
            0.20511337069030972,
        ]);
        let got = semi_empirical(&moles, &params, &LocalQueue, 0);
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
        eprintln!("comp_vec: norm = {:e}", norm);
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
        {
            let moles = load_geoms("test_files/small07", &names);
            let params = load_params("test_files/small.params");
            let want = load_mat("test_files/small.jac");
            let got = num_jac(&moles, &params, &LocalQueue, 0);
            assert!(comp_mat(got, want, 1e-5));
        }
        {
            // want jac straight from the Go version
            let moles = load_geoms("test_files/three07", &names);
            let params = load_params("test_files/three.params");
            let want = load_mat("test_files/three.jac");
            let got = num_jac(&moles, &params, &LocalQueue, 0);
            assert!(comp_mat(got, want, 1e-8));
        }
    }

    #[test]
    fn test_norm() {
        // making sure the nalgebra vector norm is what I think it is
        let v = na::DVector::<f64>::from(vec![1.0, 2.0, 3.0]);
        let got = v.norm();
        let want = ((1.0 + 4.0 + 9.0) as f64).sqrt();
        assert_eq!(got, want);
    }

    #[test]
    fn test_broyden() {
        let jac = load_mat("test_files/three.jac");
        let se_old = na::DVector::from(vec![
            0.000000000000,
            0.001668203451,
            0.001368516801,
        ]);
        let se_new = na::DVector::from(vec![
            0.000000000000,
            0.000375747088,
            0.000113449727,
        ]);
        let step = na::DVector::from(vec![
            0.879074952563,
            -0.053372779748,
            -0.386443051094,
            -1.723741422927,
            0.085500727064,
            -0.083919956047,
            -0.006440778691,
            -0.008037008147,
            -0.161085646344,
            -0.106676930451,
            0.083518030316,
            0.454558059306,
            0.599560620927,
            -0.065880893329,
            0.944895856082,
        ]);
        let got = broyden_update(&jac, &se_old, &se_new, &step);
        let want = load_mat("test_files/broyden.jac");
        assert!(comp_mat(got, want, 3e-8));
    }

    #[test]
    fn test_lev_mar() {
        // test input taken from the output of go TestLevMar
        #[rustfmt::skip]
        let jac = load_mat("test_files/levmar.jac");
        let got = lev_mar(
            &jac,
            &na::DVector::from(vec![
                0.000000000000,
                0.000453458157,
                0.000258906149,
            ]),
            &na::DVector::from(vec![
                0.203744853890,
                0.205413057340,
                0.205113370691,
            ]),
            1e-8,
        );
        let want = na::DVector::from(vec![
            0.8790749525625414,
            -0.053372779747915,
            -0.38644305109436294,
            -1.7237414229267645,
            0.08550072706391243,
            -0.08391995604728157,
            -0.006440778691258624,
            -0.008037008146972191,
            -0.161085646343529,
            -0.10667693045073838,
            0.08351803031611686,
            0.45455805930613674,
            0.599560620927178,
            -0.06588089332916275,
            0.9448958560822209,
        ]);
        assert!(comp_dvec(got, want, 3.2e-7));
    }

    #[test]
    fn test_solve() {
        // try nalgebra matrix equation solution with input from Go version
        #[rustfmt::skip]
        let lhs = load_mat("test_files/go.lhs");
        let rhs = na::dvector![
            0.35423900602931,
            -0.35423663730942,
            -0.35423845555583,
            -0.35423646488059,
            0.35423942795133,
            -0.35423939669513,
            -0.35423686878446,
            -0.35423575920428,
            -0.35423924460534,
            -0.35423965293655,
            0.35423629270772,
            -0.05025294054723,
            0.35423682777338,
            -0.35423906998699,
            0.35423447419595
        ];
        let want = na::dvector![
            0.022821308059,
            -0.019004630312,
            -0.030395950124,
            -0.022549225832,
            0.011349968251,
            -0.013318764149,
            -0.006731748989,
            -0.017577820182,
            -0.018559932276,
            -0.026805396091,
            0.009825237417,
            0.000558292732,
            0.052824051563,
            -0.029828987570,
            0.072728316382
        ];
        let lhs = match na::linalg::Cholesky::new(lhs) {
            Some(a) => a,
            None => {
                eprintln!("cholesky decomposition failed");
                std::process::exit(1);
            }
        };
        let got = lhs.solve(&rhs);
        assert!(comp_dvec(got, want, 1.04e-6));
    }

    #[ignore]
    #[test]
    fn test_algo() {
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
            0,
        );
        let want = Stats {
            norm: 7.1820,
            rmsd: 1.4364,
            max: 3.7770,
        };
        assert_eq!(got, want);
    }
}
