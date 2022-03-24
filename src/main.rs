use rust_semp::mopac::Params;
use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("file07", names);
    let mut params = load_params("params.dat");
    let ai = load_energies("rel.dat");
    setup();
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
    let start = std::time::SystemTime::now();
    while iter <= 5 {
        let jac_t = num_jac(&moles, &params, &queue);
        let step = lev_mar(&jac_t, &ai, &se, LAMBDA0);
        let try_params = Params::new(
            params.names.clone(),
            params.atoms.clone(),
            &params.values + &step,
        );
        let new_se = semi_empirical(&moles, &try_params, &queue);
        let rel = relative(&new_se);
        stats = Stats::new(&ai, &rel);
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
    takedown();
}
