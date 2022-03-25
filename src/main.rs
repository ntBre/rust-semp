use rust_semp::mopac::Params;
use rust_semp::*;
use rust_semp::slurm::Slurm;
use rust_semp::stats::Stats;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("file07", names);
    let mut params = load_params("params.dat");
    let ai = load_energies("rel.dat");
    // initial semi-empirical energies and stats
    let queue = Slurm::default(); // TODO use config
    let mut se = semi_empirical(&moles, &params, &queue);
    let rel = relative(&se);
    let mut stats = Stats::new(&ai, &rel);
    let mut last_stats = Stats::default();
    Stats::print_header();
    stats.print_step(0, &last_stats, 0);
    last_stats = stats;
    // start looping
    let mut iter = 1;
    let mut lambda = LAMBDA0;
    while iter <= 500 {
        let start = std::time::SystemTime::now();
        // TODO if Broyden do broyden
        let jac = num_jac(&moles, &params, &queue);

        // BEGIN copy-paste
        let step = lev_mar(&jac, &ai, &se, lambda / NU);
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
            let step = lev_mar(&jac, &ai, &se, lambda);
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
            let step = lev_mar(&jac, &ai, &se, lambda);
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
        se = new_se;
        params = try_params;
        last_stats = stats;
        iter += 1;
    }
}
