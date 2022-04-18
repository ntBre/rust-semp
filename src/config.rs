use serde::Deserialize;

#[derive(Deserialize, Debug, PartialEq)]
pub struct Config {
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
    /// Array of string atomic symbols like ["C", "C", "C", "H", "H"]
    pub atom_names: Vec<String>,
    /// String containing the initial parameters to be optimized
    pub params: String,
    /// Whether or not to use [Broyden's approximate update
    /// method](https://en.wikipedia.org/wiki/Broyden%27s_method) to update the
    /// Jacobian matrix
    pub broyden: bool,
    /// If `broyden` is true, the interval at which to calculate a numerical
    /// Jacobian instead of using Broyden's method
    pub broyd_int: usize,
    /// charge on the molecule. 0 for neutral, +1 for cation, -1 for anion, and
    /// so on
    pub charge: isize,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let contents = std::fs::read_to_string(filename)
            .expect("failed to load config file");
        toml::from_str(&contents).expect("failed to deserialize config file")
    }
}

#[cfg(test)]
mod tests {
    use crate::string;

    use super::*;

    #[test]
    fn test_load_full() {
        let got = Config::load("test_files/test.toml");
        let want = Config {
            job_limit: 10000,
            chunk_size: 128,
            sleep_int: 5,
            max_iter: 5,
            atom_names: string!["C", "C", "C", "H", "H"],
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
            broyden: false,
            broyd_int: 10,
            charge: 0,
        };
        assert_eq!(got, want);
    }
}
