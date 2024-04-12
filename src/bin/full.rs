use psqs::queue::local::Local;
use rust_semp::config::Config;
use rust_semp::optimize::energy::Energy;
use rust_semp::params::mopac::MopacParams;
use rust_semp::utils::*;
use rust_semp::*;

fn main() {
    let config = Config::load("test_files/test.toml");
    let names = string!["C", "C", "C", "H", "H"];
    run_algo(
        &mut std::io::stdout(),
        &config.molecules,
        MopacParams(load_params("test_files/small.params")),
        load_energies("test_files/25.dat"),
        10,
        true,
        5,
        Local::new(128, 128, 2, "inp", true, None),
        false,
        Energy {
            moles: load_geoms("test_files/small07", &names),
        },
    );
}
