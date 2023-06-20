use psqs::program::Template;
use rust_pbqff::{config::Queue, Intder};
use serde::Deserialize;

#[derive(Deserialize, Debug, PartialEq, Eq)]
pub enum Protocol {
    #[serde(alias = "energy")]
    Energy,
    #[serde(alias = "frequency")]
    Frequency,
}

pub mod raw;
use raw::*;
use symm::Irrep;

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

    /// whether or not to reset lambda after failing to improve it by successive
    /// multiplications by nu
    pub reset_lambda: bool,

    /// individual molecules
    pub molecules: Vec<Molecule>,

    pub queue: Queue,

    /// the location of the mopac executable for a local queue
    pub mopac: Option<String>,

    /// the step size in the parameters
    pub delta: f64,

    pub queue_template: Option<String>,

    pub sort_ascending: bool,
}

#[allow(clippy::large_enum_variant)]
pub(crate) enum CoordType {
    Sic(Intder),
    Cart,
    Normal(Option<bool>),
    NormalHarm,
}

impl CoordType {
    /// Returns `true` if the coord type is [`Normal`].
    ///
    /// [`Normal`]: CoordType::Normal
    #[must_use]
    pub(crate) fn is_normal(&self) -> bool {
        matches!(self, Self::Normal(..))
    }
}

pub struct Molecule {
    /// Array of string atomic symbols like ["C", "C", "C", "H", "H"]
    pub atom_names: Vec<String>,

    /// charge on the molecule. 0 for neutral, +1 for cation, -1 for anion, and
    /// so on
    pub charge: isize,

    pub geometry: psqs::geom::Geom,

    /* from here down only needed for frequencies */
    pub(crate) coord_type: CoordType,

    pub true_freqs: Vec<f64>,

    pub irreps: Vec<Irrep>,

    pub template: Template,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let raw = match RawConfig::load(filename) {
            Ok(r) => r,
            Err(e) => panic!("failed to deserialize {filename} with {e}"),
        };

        let mut molecules = Vec::new();
        for molecule in raw.molecule {
            let coord_type = match molecule.coord_type {
                raw::CoordType::Sic(f) => CoordType::Sic(Intder::load_file(&f)),
                raw::CoordType::Cart => CoordType::Cart,
                raw::CoordType::Norm(b) => CoordType::Normal(b),
            };
            molecules.push(Molecule {
                atom_names: molecule.atom_names,
                charge: molecule.charge,
                geometry: molecule.geometry.parse().unwrap(),
                coord_type,
                true_freqs: molecule.true_freqs,
                irreps: molecule.irreps,
                template: Template::from(&molecule.template),
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
            reset_lambda: raw.reset_lambda,
            molecules,
            queue: raw.queue,
            mopac: raw.mopac,
            delta: raw.delta,
            queue_template: raw.queue_template,
            sort_ascending: raw.sort_ascending.unwrap_or(false),
        }
    }
}
