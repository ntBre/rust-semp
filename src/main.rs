use rust_semp::*;

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let moles = load_geoms("file07", names);
    let ml = moles.len();
    let params = load_params("params.dat");
    let mut jobs = build_jobs(moles, params, 0, 1.0);
    setup();
    let mut energies = vec![0.; ml];
    drain(&mut jobs, &mut energies, Slurm {});
    for (i, g) in energies.iter().enumerate() {
        println!("{i:5}{g:20.12}");
    }
}
