use psqs::program::mopac::Mopac;
use psqs::program::Template;
use psqs::queue::local::LocalQueue;
use psqs::queue::Queue;
use rust_semp::*;

static MOPAC_TMPL: Template =
    Template::from("XYZ A0 scfcrt=1.D-21 aux(precision=14) PM6");

fn main() {
    let names = string!["C", "C", "C", "H", "H"];
    let moles = load_geoms("test_files/quarter07", &names);
    let ml = moles.len();
    let params = load_params("test_files/params.dat");
    let mut jobs = Mopac::build_jobs(
        &moles,
        Some(&params),
        "inp",
        0,
        1.0,
        0,
        0,
        &MOPAC_TMPL,
    );
    let mut energies = vec![0.; ml];
    setup();
    LocalQueue::new("inp", 128).drain(&mut jobs, &mut energies);
    takedown();
    for (i, e) in energies.iter().enumerate() {
        println!("{i:5}{e:20.12}");
    }
}
