use super::*;
use crate::{config::Config, tests::load_mat, utils::setup};
use psqs::queue::local::Local;

#[test]
#[ignore]
fn freq_semi_empirical() {
    let config = Config::load("test_files/test.toml");
    let freq = Frequency::default();
    setup();
    let queue = Local {
        dir: "inp".to_owned(),
        chunk_size: 128,
        mopac: "/opt/mopac/mopac".to_owned(),
    };
    let mut got = freq
        .semi_empirical(
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
        )
        .unwrap();
    let got = got.as_mut_slice();
    got.sort_by(|a, b| b.partial_cmp(a).unwrap());
    approx::assert_abs_diff_eq!(
        na::DVector::from(Vec::from(got)),
        na::DVector::from(vec![
            2783.961357017977,
            2764.3281931305864,
            1775.6737582101287,
            1177.145694605022,
            1040.6327552461862,
            960.1588882411834,
            927.0954879807678,
            919.9554364598707,
            905.3386838377455
        ]),
        epsilon = 0.1
    );
    takedown();
}

#[test]
#[ignore]
fn freq_num_jac() {
    // this test takes 23 minutes with the current implementation at work
    let config = Config::load("test_files/test.toml");
    let freq = Frequency::default();
    setup();
    let queue = Local {
        dir: "inp".to_owned(),
        chunk_size: 128,
        mopac: "/opt/mopac/mopac".to_owned(),
    };
    let got = freq.num_jac(
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
    let want = load_mat("test_files/freq.jac");
    // this should agree to 1e-12 depending on the computer I guess, I printed
    // `got` to 12 decimal places when obtaining the `want` value
    approx::assert_abs_diff_eq!(got, want, epsilon = 1e-5,);
    takedown();
}
