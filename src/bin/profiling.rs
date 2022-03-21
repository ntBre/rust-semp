use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("quarter07", names);
    let params = load_params("params.dat");
    let mut jobs = build_jobs(moles, params);
    setup();
    drain(&mut jobs, LocalQueue {});
    for (i, g) in jobs.dst.iter().enumerate() {
        println!("{i:5}{g:20.12}");
    }
    takedown();
}
