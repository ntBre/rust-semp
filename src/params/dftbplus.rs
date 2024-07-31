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
#[derive(Clone)]
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

#[derive(Clone)]
pub struct DFTBPlusParams {
    files: Vec<SKF>,
}

impl DFTBPlusParams {
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
        // TODO this whole function is actually for a single SKF, not for a
        // whole DFTBPlusParams
        for (i, line) in s.lines().enumerate() {
            if line.starts_with("Spline") {
                break;
            }
            if i == 0 {
                // I don't think we actually need these separately for our
                // purposes, but following the dftb+ parser for now.
                let mut sp = line.split(',').map(str::trim);
                let dist = sp.next().unwrap();
                let ngrid = sp.next().unwrap();
            } else if i == 2 {
                // the only line we care about
                if let Ok(fields) = utils::split_somehow(line) {
                    // assume homo-atom skf here
                    let &[ed, ep, es, spe, ud, up, us, fd, fp, fs] =
                        &fields[..]
                    else {
                        panic!("wrong number of skf fields for homo atoms")
                    };
                } else {
                    log::warn!(
                        "assuming hetero-atom skf, enable tracing to see line"
                    );
                    log::trace!("{name:?}:{i}: {line}");
                };
            }
            println!("{:?}", line);
        }
        Ok(Self)
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
    use super::*;

    #[test]
    fn test_from_file() {
        DFTBPlusParams::from_file("test_files/H-O.skf").unwrap();
        DFTBPlusParams::from_file("test_files/H-H.skf").unwrap();
    }
}
