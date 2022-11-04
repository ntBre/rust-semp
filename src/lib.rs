use crate::optimize::Optimize;
use crate::utils::log_params;
use config::Molecule;
use nalgebra as na;
use psqs::program::mopac::{Mopac, Params};
use psqs::queue::Queue;
use stats::Stats;
use std::{clone::Clone, io::Write};

pub mod config;
pub mod optimize;
pub mod stats;
#[cfg(test)]
mod tests;
pub mod utils;

static DEBUG: bool = false;
/// convergence threshold for the change in the norm, rmsd, and max
const DCONV_THRSH: f64 = 1e-5;

pub static LAMBDA0: f64 = 1e-8;
pub static NU: f64 = 2.0;
pub static MAX_TRIES: usize = 10;

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
    let lhs = a_star.clone() + li;
    let mut d = match na::linalg::Cholesky::new(lhs.clone()) {
        Some(a) => {
            let lhs = a;
            lhs.solve(&g_star)
        }
        None => {
            eprintln!("cholesky decomposition failed");
            eprintln!("jac\n{:.8}", jac);
            eprintln!("a\n{:.8}", a);
            eprintln!("a*\n{:.8}", a_star);
            eprintln!("lhs\n{:.8}", lhs);
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

#[allow(clippy::too_many_arguments)]
pub fn run_algo<O: Optimize, Q: Queue<Mopac> + Sync, W: Write>(
    param_log: &mut W,
    molecules: &[Molecule],
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
    let mut params = params;
    log_params(param_log, 0, &params);
    let mut start = std::time::Instant::now();
    // initial semi-empirical energies and stats
    let mut se = optimizer
        .semi_empirical(&params, &queue, molecules)
        .unwrap_or_else(|| na::DVector::zeros(ai.len()));
    optimizer.log(0, &se, &ai);
    let mut old_se = se.clone();
    let mut stats = Stats::new(&ai, &se, conv);
    let mut last_stats = Stats::default();
    Stats::print_header();
    stats.print_step(0, &last_stats, start.elapsed().as_millis(), LAMBDA0);
    last_stats = stats;
    // start looping
    let mut iter = 1;
    let mut lambda = LAMBDA0;
    let mut del_norm: f64 = 1.0;
    let mut del_rmsd: f64 = 1.0;
    let mut del_max: f64 = 1.0;
    let mut in_broyden = false;
    let mut need_num_jac = false;
    start = std::time::Instant::now();
    let mut jac = optimizer.num_jac(&params, &queue, molecules);

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
            start = std::time::Instant::now();
            jac = broyden_update(&jac, &old_se, &se, &step);
        } else if iter > 1 {
            in_broyden = false;
            need_num_jac = false;
            start = std::time::Instant::now();
            jac = optimizer.num_jac(&params, &queue, molecules);
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
        let mut new_se = optimizer
            .semi_empirical(&try_params, &queue, molecules)
            .unwrap_or_else(|| na::DVector::zeros(ai.len()));
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
            new_se = optimizer
                .semi_empirical(&try_params, &queue, molecules)
                .unwrap_or_else(|| na::DVector::zeros(ai.len()));
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
            new_se = optimizer
                .semi_empirical(&try_params, &queue, molecules)
                .unwrap_or_else(|| na::DVector::zeros(ai.len()));
            stats = Stats::new(&ai, &new_se, conv);

            if stats.norm - last_stats.norm > dnorm {
                break;
            }
            i += 1;
        }

        // log the time
        let time = start.elapsed().as_millis();

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
        stats.print_step(iter, &last_stats, time, lambda);
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
