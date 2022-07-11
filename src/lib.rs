#[cfg(test)]
mod tests;

use std::{
    clone::Clone,
    fs::{self, File},
    io::BufReader,
    io::{BufRead, Write},
    path::Path,
    rc::Rc,
};

use config::Molecule;
use nalgebra as na;
use std::cmp::Ordering;
use std::iter::zip;
use symm::Irrep;

use psqs::program::mopac::{Mopac, Params};
use psqs::queue::Queue;
use psqs::{geom::Geom, program::Template};
use stats::Stats;
use symm::atom::Atom;

use crate::optimize::Optimize;

pub mod config;
pub mod optimize;
pub mod stats;

static DEBUG: bool = false;
/// convergence threshold for the change in the norm, rmsd, and max
const DCONV_THRSH: f64 = 1e-5;

pub static LAMBDA0: f64 = 1e-8;
pub static NU: f64 = 2.0;
pub static MAX_TRIES: usize = 10;

static MOPAC_TMPL: Template =
    Template::from("XYZ A0 scfcrt=1.D-21 aux(precision=14) PM6");

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

/// Take an INTDER-style `file07` file and parse it into a Vec of geometries
pub fn load_geoms(filename: &str, atom_names: &[String]) -> Vec<Rc<Geom>> {
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
                ret.push(Rc::new(Geom::Xyz(buf)));
                buf = Vec::new();
                idx = 0;
            }
            continue;
        }
        let fields: Vec<_> = line
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect();
        buf.push(Atom::new_from_label(
            &atom_names[idx],
            fields[0],
            fields[1],
            fields[2],
        ));
        idx += 1;
    }
    ret.push(Rc::new(Geom::Xyz(buf)));
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
    let mut d = match na::linalg::Cholesky::new(lhs.clone()) {
        Some(a) => {
            let lhs = a;
            lhs.solve(&g_star)
        }
        None => {
            eprintln!("cholesky decomposition failed");
            let lhs = na::linalg::LU::new(lhs);
            lhs.solve(&g_star).expect("LU decomposition also failed")
        }
    };
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

/// sort freqs by the irrep in the same position and then by frequency. panics
/// if `freqs` and `irreps` are not the same length.
pub fn sort_irreps(freqs: &[f64], irreps: &[Irrep]) -> Vec<f64> {
    assert!(freqs.len() == irreps.len());
    let mut pairs: Vec<_> = zip(irreps, freqs).collect();
    pairs.sort_by(|a, b| match a.0.cmp(b.0) {
        Ordering::Equal => a.1.partial_cmp(b.1).unwrap(),
        other => other,
    });
    pairs.iter().map(|x| x.1.clone()).collect()
}

pub fn run_algo<O: Optimize, Q: Queue<Mopac>, W: Write>(
    param_log: &mut W,
    molecules: Vec<Molecule>,
    geom_file: &str,
    params: Params,
    ai: na::DVector<f64>,
    max_iter: usize,
    broyden: bool,
    broyd_int: usize,
    queue: Q,
    reset_lambda: bool,
    optimizer: O,
) -> Stats {
    let conv = optimizer.stat_multiplier();
    // this is only needed for energies, can I take it out in general?
    let moles = load_geoms(geom_file, &molecules[0].atom_names);
    let mut params = params.clone();
    log_params(param_log, 0, &params);
    // initial semi-empirical energies and stats
    let mut se =
        optimizer.semi_empirical(&moles, &params, &queue, molecules[0].charge);
    optimizer.log(0, &se, &ai);
    let mut old_se = se.clone();
    let mut stats = Stats::new(&ai, &se, conv);
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
    let mut jac =
        optimizer.num_jac(&moles, &params, &queue, molecules[0].charge);
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
            jac =
                optimizer.num_jac(&moles, &params, &queue, molecules[0].charge);
        } // else (first iteration) use jac from outside loop

        if DEBUG {
            eprintln!("{:8.3e}", jac);
        }

        let lambda_init = lambda;

        // first try with λ/ν
        lambda /= NU;
        step = lev_mar(&jac, &ai, &se, lambda);
        let mut try_params = Params::new(
            params.names.clone(),
            params.atoms.clone(),
            &params.values + &step,
        );
        let mut new_se = optimizer.semi_empirical(
            &moles,
            &try_params,
            &queue,
            molecules[0].charge,
        );
        stats = Stats::new(&ai, &new_se, conv);

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
            new_se = optimizer.semi_empirical(
                &moles,
                &try_params,
                &queue,
                molecules[0].charge,
            );
            stats = Stats::new(&ai, &new_se, conv);

            i += 1;
            if stats.norm - last_stats.norm > dnorm || i > MAX_TRIES {
                bad = true;
                break;
            }
        }

        if reset_lambda {
            lambda /= NU.powi(i as i32);
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
            new_se = optimizer.semi_empirical(
                &moles,
                &try_params,
                &queue,
                molecules[0].charge,
            );
            stats = Stats::new(&ai, &new_se, conv);

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
        optimizer.log(iter, &se, &ai);
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
