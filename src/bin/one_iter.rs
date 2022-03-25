use rust_semp::*;
use rust_semp::local::LocalQueue;
use rust_semp::queue::Queue;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/quarter07", names);
    let ml = moles.len();
    let params = load_params("test_files/params.dat");
    let mut jobs = build_jobs(&moles, &params, 0, 1.0, 0);
    setup();
    let mut energies = vec![0.; ml];
    LocalQueue.drain(&mut jobs, &mut energies);
    for (i, e) in energies.iter().enumerate() {
        println!("{i:5}{e:20.12}");
    }
    takedown();
}
