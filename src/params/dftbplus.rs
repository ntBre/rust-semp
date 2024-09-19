#![allow(unused, unreachable_code)]

use std::{
    error::Error,
    ffi::OsString,
    fmt::Display,
    fs::{read_to_string, File},
    io,
    ops::Add,
    path::{Path, PathBuf},
    str::FromStr,
};

use crate::{
    driver::{self, Params},
    Dvec,
};

mod utils {
    use std::{io::Error, num::ParseFloatError};

    /// It looks like SKF file entries can be separated by either commas or
    /// spaces, so this function will try to split on both and make sure we get
    /// all the entries we expect. It also converts everything to f64
    pub(super) fn split_somehow(
        line: &str,
    ) -> Result<Vec<f64>, ParseFloatError> {
        line.split(',')
            .flat_map(|s| s.split_ascii_whitespace())
            .map(|s| s.parse::<f64>())
            .collect()
    }
}

/// Representation of a single SKF.
#[derive(Clone, Debug)]
pub(crate) struct Skf {
    /// Base filename for writing back to disk
    pub(crate) basename: PathBuf,
    /// Header junk
    header: String,
    /// For now, the parameters to optimize. If empty, we can just write the
    /// file without this line. For my current understanding, that will imply a
    /// hetero-atomic SKF
    line2: Vec<f64>,
    /// Trailing junk
    footer: String,
}

impl Skf {
    /// Construct a [Self] from the Slater-Koster file located at `path`.
    ///
    /// The second line of the homo-atomic Slater-Koster files are the only ones
    /// that are easy to adjust, as far as I can tell. According to
    /// <https://dftb.org/fileadmin/DFTB/public/misc/slakoformat.pdf>, which is
    /// the only reference for these files I can find anywhere, this line
    /// contains the values: `Ed Ep Es SPE Ud Up Us fd fp fs`, but SPE is
    /// ignored by DFTB+. Quoting from the link above:
    /// * Ed, Ep and Es are the on-site energies for the angular momenta d, p
    ///   and s for the given atom
    /// * Ud, Up, Us are the Hubbard U values for the appropriate angular
    ///   momenta
    /// * fd, fp, and fs are the orbital occupations for the neutral atom in the
    ///   ground state
    ///
    /// So I don't think it makes much sense to optimize the f values either.
    /// That leaves us with only six parameters per atom type, and Ed is zero
    /// for both C-C and H-H in mio-1-1 anyway.
    ///
    /// The lines after 2 are the elements of the DFTB Hamiltonian, which I
    /// would generally like to include in the optimization, but there are ~500
    /// lines in each of them, and each one contains 20 entries. We sadly cannot
    /// optimize 10,000 parameters per atom pair. Similarly, the hetero-atom
    /// Slater-Koster files only have the big Hamiltonian blocks, so there's
    /// nothing we can do there either, as far as I can tell.
    ///
    /// The manual above also indicates that the `Spline` section is optional,
    /// so I think we can actually return as soon as we see it without breaking
    /// the calculations.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Box<dyn Error>> {
        let path = path.as_ref();
        let name = path.file_name();
        let s = read_to_string(path)?;
        let mut lines = s.lines();
        let header = lines
            .next()
            .expect("expecting at least 1 line in skf")
            .to_owned();
        let line = lines.next().expect("expecting at least 2 lines in skf");
        let mut rest = Vec::new();
        let mut line2 = Vec::new();
        if let Ok(fields) = utils::split_somehow(line) {
            // assume homo-atom skf here
            assert_eq!(fields.len(), 10, "expecting 10 fields for line2");
            line2 = fields;
        } else {
            log::warn!("assuming hetero-atom skf, enable tracing to see line");
            log::trace!("{name:?}:2: {line}");
            rest.push(line);
        };
        for line in lines {
            rest.push(line);
        }
        Ok(Self {
            basename: name.unwrap().to_owned().into(),
            header,
            line2,
            footer: rest.join("\n"),
        })
    }

    pub(crate) fn to_file(&self, param_file: PathBuf) -> std::io::Result<()> {
        use std::io::Write;
        let mut f = File::create(param_file)?;
        write!(f, "{self}")
    }
}

#[derive(Clone, Debug)]
pub struct DFTBPlusParams {
    pub(crate) files: Vec<Skf>,
}

impl DFTBPlusParams {
    /// take the first 7 entries from line2, skipping the orbital occupations
    const LINE2_SIZE: usize = 7;

    /// Returns a flattened vector of references into the first 7 entries of
    /// each line2 (the optimizable parameters).
    fn params_mut(&mut self) -> Vec<&mut f64> {
        self.files
            .iter_mut()
            .flat_map(|mut s| {
                let end = usize::min(s.line2.len(), Self::LINE2_SIZE);
                &mut s.line2[..end]
            })
            .collect()
    }

    fn params(&self) -> impl Iterator<Item = &f64> {
        self.files.iter().flat_map(|s| {
            &s.line2[..usize::min(s.line2.len(), Self::LINE2_SIZE)]
        })
    }
}

impl FromStr for DFTBPlusParams {
    type Err = Box<dyn Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self {
            files: s
                .lines()
                .filter_map(|s| {
                    let s = s.trim();
                    if s.is_empty() {
                        None
                    } else {
                        Some(s)
                    }
                })
                .map(Skf::from_file)
                .collect::<Result<Vec<Skf>, _>>()?,
        })
    }
}

impl Display for Skf {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.header)?; // should have trailing newline?
        if !self.line2.is_empty() {
            writeln!(f)?;
        }
        for v in self.line2.iter() {
            write!(f, " {v:20.12}")?;
        }
        writeln!(f, "\n{}", self.footer)?;
        Ok(())
    }
}

impl Display for DFTBPlusParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, file) in self.files.iter().enumerate() {
            writeln!(f, "=== START OF SKF {} ===", i + 1)?;
            write!(f, "{}", file)?;
            writeln!(f, "=== END OF SKF {} ===", i + 1)?;
        }
        Ok(())
    }
}

impl Params for DFTBPlusParams {
    fn incr_value(&mut self, idx: usize, delta: f64) {
        let mut v = self.params_mut();
        *v[idx] += delta;
    }

    fn len(&self) -> usize {
        self.params().count()
    }
}

impl Add<&Dvec> for DFTBPlusParams {
    type Output = Self;

    fn add(mut self, rhs: &Dvec) -> Self::Output {
        for (idx, p) in self.params_mut().into_iter().enumerate() {
            *p += rhs[idx];
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use insta::{assert_debug_snapshot, assert_snapshot};
    use psqs::{program::dftbplus::DFTBPlus, queue::local::Local};

    use crate::{
        config::Config,
        optimize::{energy::Energy, frequency::Frequency},
        run_algo,
        utils::{setup, takedown},
    };

    use super::*;

    #[test]
    fn test_from_file() {
        assert_debug_snapshot!(DFTBPlusParams::from_str(
            "test_files/H-O.skf
             test_files/H-H.skf
            ",
        ));
    }

    #[test]
    fn test_config() {
        assert_debug_snapshot!(Config::load("test_files/dftb.toml"));
    }

    #[test]
    fn test_display() {
        let params = DFTBPlusParams::from_str(
            "test_files/H-O.skf
             test_files/H-H.skf
            ",
        )
        .unwrap();
        assert_eq!(params.len(), 7);
        assert_snapshot!(params);
    }
}
