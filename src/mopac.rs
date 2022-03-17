#![allow(dead_code)]

use super::*;
use std::fs::File;
use std::io::Write;

#[derive(Debug)]
pub struct Param {
    name: String,
    atom: String,
    value: f64,
}

impl ToString for Param {
    fn to_string(&self) -> String {
        format!("{:<8}{:>8}{:20.12}\n", self.name, self.atom, self.value)
    }
}

/// Mopac holds the information needed to write a MOPAC input file
pub struct Mopac {
    params: Vec<Param>,
    geom: Vec<Atom>,
}

impl Mopac {
    pub fn write_params(&self, filename: &str) {
        let mut body = String::new();
        for p in &self.params {
            body.push_str(&p.to_string());
        }
        let mut file = File::create(filename).expect("failed to create params file");
        write!(file, "{}", body).expect("failed to write params file");
    }

    /// Caller should ensure that paramfile is written before calling. This
    /// allows common parameter files to be shared between jobs
    pub fn write_input(&self, filename: &str, paramfile: &str) {
        let geom = geom_string(&self.geom);
        let mut file = File::create(filename).expect("failed to create input file");
        write!(
            file,
            "scfcrt=1.D-21 aux(precision=14) PM6 external={paramfile}
Comment line 1
Comment line 2
{geom}
"
        )
        .expect("failed to write input file");
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    fn test_mopac() -> Mopac {
        Mopac {
            params: vec![
                Param {
                    name: String::from("USS"),
                    atom: String::from("H"),
                    value: -11.246958000000,
                },
                Param {
                    name: String::from("ZS"),
                    atom: String::from("H"),
                    value: 1.268641000000,
                },
                Param {
                    name: String::from("BETAS"),
                    atom: String::from("H"),
                    value: -8.352984000000,
                },
                Param {
                    name: String::from("GSS"),
                    atom: String::from("H"),
                    value: 14.448686000000,
                },
                Param {
                    name: String::from("USS"),
                    atom: String::from("C"),
                    value: -51.089653000000,
                },
                Param {
                    name: String::from("UPP"),
                    atom: String::from("C"),
                    value: -39.937920000000,
                },
                Param {
                    name: String::from("ZS"),
                    atom: String::from("C"),
                    value: 2.047558000000,
                },
                Param {
                    name: String::from("ZP"),
                    atom: String::from("C"),
                    value: 1.702841000000,
                },
                Param {
                    name: String::from("BETAS"),
                    atom: String::from("C"),
                    value: -15.385236000000,
                },
                Param {
                    name: String::from("BETAP"),
                    atom: String::from("C"),
                    value: -7.471929000000,
                },
                Param {
                    name: String::from("GSS"),
                    atom: String::from("C"),
                    value: 13.335519000000,
                },
                Param {
                    name: String::from("GPP"),
                    atom: String::from("C"),
                    value: 10.778326000000,
                },
                Param {
                    name: String::from("GSP"),
                    atom: String::from("C"),
                    value: 11.528134000000,
                },
                Param {
                    name: String::from("GP2"),
                    atom: String::from("C"),
                    value: 9.486212000000,
                },
                Param {
                    name: String::from("HSP"),
                    atom: String::from("C"),
                    value: 0.717322000000,
                },
            ],
            geom: Vec::new(),
        }
    }

    #[test]
    fn test_write_input() {
        test_mopac().write_input("/tmp/test.mop", "params.dat");
        let got = fs::read_to_string("/tmp/test.mop").expect("file not found");
        let want = "scfcrt=1.D-21 aux(precision=14) PM6 external=params.dat
Comment line 1
Comment line 2

";
        assert_eq!(got, want);
    }

    #[test]
    fn test_write_params() {
        test_mopac().write_params("/tmp/params.dat");
        let got = fs::read_to_string("/tmp/params.dat").expect("file not found");
        let want = "USS            H    -11.246958000000
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
";
        assert_eq!(got, want);
    }
}
