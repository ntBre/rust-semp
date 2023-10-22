use psqs::program::mopac::Mopac;
use psqs::queue::local::Local;
use rust_semp::optimize::{frequency::Frequency, Optimize};
use rust_semp::utils::setup;

fn main() {
    let config = rust_semp::config::Config::load("test_files/test.toml");
    let freq = Frequency::default();
    setup();
    let queue = Local::new(128, 128, 2, "inp", true, None);
    Optimize::<Mopac>::num_jac(
        &freq,
        &"USS            H    -11.246958000000
    ZS             H      1.268641000000
    BETAS          H     -8.352984000000
    GSS            H     14.448686000000
    USS            C    -51.089653000000
    UPP            C    -39.937920000000
    ZS             C      2.047558000000
    ZP             C      1.702841000000
    BETAS          C    -15.385236000000
    BETAP          C     -7.471929000000
    GSS            C     13.335519000000
    GPP            C     10.778326000000
    GSP            C     11.528134000000
    GP2            C      9.486212000000
    HSP            C      0.717322000000
    FN11           C      0.046302000000"
            .parse()
            .unwrap(),
        &queue,
        &config.molecules,
        9,
    );
}
