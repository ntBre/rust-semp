#![allow(unused, unreachable_code)]

use std::{
    error::Error, ffi::OsString, fmt::Display, fs::read_to_string, io,
    ops::Add, path::Path, str::FromStr,
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
struct SKF {
    /// Header junk
    header: String,
    /// For now, the parameters to optimize. If empty, we can just write the
    /// file without this line. For my current understanding, that will imply a
    /// hetero-atomic SKF
    line2: Vec<f64>,
    /// Trailing junk
    footer: String,
}

impl SKF {
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
            if line.starts_with("Spline") {
                break;
            }
            rest.push(line);
        }
        Ok(Self {
            header,
            line2,
            footer: rest.join("\n"),
        })
    }
}

#[derive(Clone, Debug)]
pub struct DFTBPlusParams {
    files: Vec<SKF>,
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
                .map(SKF::from_file)
                .collect::<Result<Vec<SKF>, _>>()?,
        })
    }
}

impl Display for DFTBPlusParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

impl Params for DFTBPlusParams {
    fn incr_value(&mut self, idx: usize, delta: f64) {
        todo!()
    }

    fn len(&self) -> usize {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_debug_snapshot;

    use super::*;

    #[test]
    fn test_from_file() {
        assert_debug_snapshot!(DFTBPlusParams::from_str(
            "test_files/H-O.skf
             test_files/H-H.skf
            ",
        ));
    }
}
