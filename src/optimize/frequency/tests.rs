use super::*;
use crate::{
    config::Config,
    tests::{hostname, load_mat},
    utils::setup,
};

use nalgebra as na;
use psqs::queue::local::Local;

/// this is kinda lame because these are essentially pbqff tests, but they do
/// utilize a lot of the same machinery
fn se_from_config(config: Config, want: na::DVector<f64>) {
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
        want,
        epsilon = 0.1
    );
    takedown();
}

#[test]
#[ignore]
fn sic_semi_empirical() {
    let config = Config::load("test_files/test.toml");
    let want = na::dvector![
        2783.961357017977,
        2764.3281931305864,
        1775.6737582101287,
        1177.145694605022,
        1040.6327552461862,
        960.1588882411834,
        927.0954879807678,
        919.9554364598707,
        905.3386838377455
    ];
    se_from_config(config, want);
}

#[test]
#[ignore]
fn norm_semi_empirical() {
    let config = Config::load("test_files/norm.toml");
    let want = na::dvector![
        2783.0355951341962,
        2763.275377957968,
        1776.3820512917491,
        1177.824074229593,
        1041.2641575043508,
        960.0159356063677,
        927.3499180637409,
        920.5868937867253,
        906.0837431489708
    ];
    se_from_config(config, want);
}

/// this tests:
/// 1) HFF part of normal coordinates
/// 2) fitted normal coordinates
/// 3) all parts of a numerical Jacobian on frequencies
#[test]
#[ignore]
fn freq_num_jac() {
    let config = Config::load("test_files/jac.toml");
    let freq = Frequency::default();
    setup();
    let queue = Local {
        dir: "inp".to_owned(),
        chunk_size: 128,
        mopac: "/opt/mopac/mopac".to_owned(),
    };
    let got = freq.num_jac(
        &"
USS            H    -11.246958000000
ZS             H      1.268641000000
USS            C    -51.089653000000
UPP            C    -39.937920000000
ZS             C      2.047558000000
"
        .parse()
        .unwrap(),
        &queue,
        &config.molecules,
        9,
    );
    let want = if hostname() == "keystone" {
        load_mat("test_files/freq.jac")
    } else {
        load_mat("test_files/github.jac")
    };
    if !approx::abs_diff_eq!(got, want, epsilon = 1e-12) {
        panic!("got\n{got}, wanted\n{want}");
    }
    takedown();
}
