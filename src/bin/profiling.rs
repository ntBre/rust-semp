use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("quarter07", names);
    let ml = moles.len();
    let params = load_params("params.dat");
    let mut jobs = build_jobs(moles, params, 0, 1.0);
    setup();
    let mut energies = vec![0.; ml];
    drain(&mut jobs, &mut energies, LocalQueue {});
    for (i, e) in energies.iter().enumerate() {
        println!("{i:5}{e:20.12}");
    }
    takedown();
}
