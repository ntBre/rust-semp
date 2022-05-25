use psqs::queue::local::LocalQueue;
use psqs::queue::Queue;
use rust_semp::*;

fn main() {
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/quarter07", &names);
    let ml = moles.len();
    let params = load_params("test_files/params.dat");
    let mut jobs = build_jobs(&moles, &params, 0, 1.0, 0, 0);
    let mut energies = vec![0.; ml];
    setup();
    LocalQueue {
        dir: "inp".to_string(),
    }
    .drain(&mut jobs, &mut energies);
    takedown();
    for (i, e) in energies.iter().enumerate() {
        println!("{i:5}{e:20.12}");
    }
}
