use serde::Deserialize;

#[derive(Deserialize, Debug, PartialEq)]
pub enum Protocol {
    #[serde(alias = "energy")]
    Energy,
    #[serde(alias = "frequency")]
    Frequency,
}

mod raw;
use raw::*;

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
    pub molecules: Vec<Molecule>,

    /// path to the actual spectro program to run in gspectro
    pub spectro_cmd: String,

    /// path to gspectro
    pub gspectro_cmd: String,
}

pub struct Molecule {
    /// Array of string atomic symbols like ["C", "C", "C", "H", "H"]
    pub atom_names: Vec<String>,

    /// charge on the molecule. 0 for neutral, +1 for cation, -1 for anion, and
    /// so on
    pub charge: isize,

    /// dummy atoms of the form (axis, atom)
    pub dummies: Vec<(usize, usize)>,

    pub geometry: psqs::geom::Geom,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let raw = match RawConfig::load(filename) {
            Ok(r) => r,
            Err(e) => panic!("failed to deserialize {} with {}", filename, e),
        };

        let mut molecules = Vec::new();
        for molecule in raw.molecule {
            molecules.push(Molecule {
                atom_names: molecule.atom_names,
                charge: molecule.charge,
                dummies: molecule.dummies,
                geometry: molecule.geometry.parse().unwrap(),
            });
        }
        Self {
            job_limit: raw.job_limit,
            chunk_size: raw.chunk_size,
            sleep_int: raw.sleep_int,
            max_iter: raw.max_iter,
            params: raw.params,
            broyden: raw.broyden,
            broyd_int: raw.broyd_int,
            optimize: raw.optimize,
            reorder: raw.reorder,
            reset_lambda: raw.reset_lambda,
            molecules,
            spectro_cmd: raw.spectro_cmd,
            gspectro_cmd: raw.gspectro_cmd,
        }
    }
}
