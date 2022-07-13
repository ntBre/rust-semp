use serde::Deserialize;

use super::Protocol;

#[derive(Deserialize, Debug, PartialEq)]
pub(super) struct RawConfig {
    /// The maximum number of jobs that should be written/submitted to the Queue
    /// at one time.
    pub job_limit: usize,

    /// The number of jobs to group into a single Queue submission script. The
    /// higher this number, the less strain on the submission system. The lower
    /// this number, the more jobs that can theoretically run at one time.
    pub chunk_size: usize,

    /// The time in seconds to sleep between iterations of polling running jobs
    /// when no jobs finished on the previous iteration
    pub sleep_int: usize,

    /// The maximum number of iterations to run the Levenberg-Marquardt
    /// algorithm for
    pub max_iter: usize,

    /// String containing the initial parameters to be optimized
    pub params: String,

    /// Whether or not to use [Broyden's approximate update
    /// method](https://en.wikipedia.org/wiki/Broyden%27s_method) to update the
    /// Jacobian matrix
    pub broyden: bool,

    /// If `broyden` is true, the interval at which to calculate a numerical
    /// Jacobian instead of using Broyden's method
    pub broyd_int: usize,

    /// the values to base the optimization on. energy or frequency
    pub optimize: Protocol,

    /// whether or not to try re-aligning the coordinates onto the z-axis
    pub reorder: bool,

    /// whether or not to reset lambda after failing to improve it by successive
    /// multiplications by nu
    pub reset_lambda: bool,

    /// individual molecules
    pub molecule: Vec<RawMolecule>,

    /// path to the actual spectro program to run in gspectro
    pub spectro_cmd: String,

    /// path to gspectro
    pub gspectro_cmd: String,
}

/// deserialize into this for the String geometry before converting to a real
/// Geom with FromStr
#[derive(Clone, Deserialize, Debug, PartialEq)]
pub(super) struct RawMolecule {
    /// Array of string atomic symbols like ["C", "C", "C", "H", "H"]
    pub atom_names: Vec<String>,

    /// charge on the molecule. 0 for neutral, +1 for cation, -1 for anion, and
    /// so on
    pub charge: isize,

    /// dummy atoms of the form (axis, atom)
    pub dummies: Vec<(usize, usize)>,

    pub geometry: String,

    pub intder_file: String,
}

impl RawConfig {
    pub fn load(filename: &str) -> Result<Self, toml::de::Error> {
        let contents = std::fs::read_to_string(filename)
            .expect("failed to load config file");
        toml::from_str(&contents)
    }
}

#[cfg(test)]
mod tests {
    use crate::string;

    use super::*;

    #[test]
    fn test_load_full() {
        let got = RawConfig::load("test_files/test.toml").unwrap();
        let want = RawConfig {
            job_limit: 10000,
            chunk_size: 128,
            sleep_int: 5,
            max_iter: 5,
            params: String::from(
                "USS        H    -11.246958000000
ZS         H      1.268641000000
BETAS      H     -8.352984000000
GSS        H     14.448686000000
USS        C    -51.089653000000
UPP        C    -39.937920000000
ZS         C      2.047558000000
ZP         C      1.702841000000
BETAS      C    -15.385236000000
BETAP      C     -7.471929000000
GSS        C     13.335519000000
GPP        C     10.778326000000
GSP        C     11.528134000000
GP2        C      9.486212000000
HSP        C      0.717322000000
FN11       C      0.046302000000
",
            ),
            molecule: vec![RawMolecule {
                atom_names: string!["C", "C", "C", "H", "H"],
                charge: 0,
                dummies: vec![],
                geometry: "C
C 1 CC
C 1 CC 2 CCC
H 2 CH 1 HCC 3 180.0
H 3 CH 1 HCC 2 180.0

CC =                  1.42101898
CCC =                55.60133141
CH =                  1.07692776
HCC =               147.81488230
"
                .to_string(),
                intder_file: "test_files/intder.in".to_string(),
            }],
            broyden: false,
            broyd_int: 10,
            optimize: Protocol::Energy,
            reorder: false,
            reset_lambda: false,
            spectro_cmd: "/home/brent/Projects/pbqff/bin/spectro".to_string(),
            gspectro_cmd:
                "/home/brent/Projects/chemutils/spectro/spectro/spectro"
                    .to_string(),
        };
        assert_eq!(got, want);
    }
}
