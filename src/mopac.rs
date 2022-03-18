use super::*;
use std::fs::File;
use std::io::{BufReader, Write};

/// kcal/mol per hartree
const KCALHT: f64 = 627.5091809;

#[derive(Debug, Clone)]
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

impl PartialEq for Param {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.atom == other.atom
            && (self.value - other.value).abs() < 1e-12
    }
}

impl Param {
    pub fn new(name: &str, atom: &str, value: f64) -> Self {
        Self {
            name: name.to_string(),
            atom: atom.to_string(),
            value,
        }
    }
}

/// Mopac holds the information needed to write a MOPAC input file. `filename`
/// should not include an extension. `.mop` will be appended for input files,
/// and `.out` and `.aux` will be appended for output files.
#[derive(Debug)]
pub struct Mopac {
    pub filename: String,
    pub params: Vec<Param>,
    pub geom: Vec<Atom>,
}

pub enum OutputError {
    FileNotFound,
}

impl Mopac {
    pub fn write_params(&self, filename: &str) {
        let mut body = String::new();
        for p in &self.params {
            body.push_str(&p.to_string());
        }
        let mut file =
            File::create(filename).expect("failed to create params file");
        write!(file, "{}", body).expect("failed to write params file");
    }

    /// Caller should ensure that paramfile is written before calling. This
    /// allows common parameter files to be shared between jobs
    pub fn write_input(&self, paramfile: &str) {
        let geom = geom_string(&self.geom);
        let mut file = File::create(format!("{}.mop", self.filename))
            .expect("failed to create input file");
        write!(
            file,
            "XYZ 1SCF A0 scfcrt=1.D-21 aux(precision=14) PM6 external={paramfile}
Comment line 1
Comment line 2
{geom}
"
        )
        .expect("failed to write input file");
    }

    /// return filename with its extension replaced by ext
    fn with_ext(filename: &str, ext: &str) -> String {
        let base = Path::new(filename).file_stem().unwrap();
        let mut ret = String::from(base.to_str().unwrap());
        ret.push_str(ext);
        ret
    }

    /// Reads a MOPAC output file. If normal termination occurs, also try
    /// reading the `.aux` file to extract the energy from there. This function
    /// panics if an error is found in the output file. If a non-fatal error
    /// occurs (file not found, not written to yet, etc) None is returned.
    pub fn read_output(&self) -> Option<f64> {
        let outfile = Self::with_ext(&self.filename, ".out");
        let f = match File::open(outfile) {
            Ok(file) => file,
            Err(_) => {
                eprintln!("file {} not found", self.filename);
                return None;
            } // file not found
        };
        let f = BufReader::new(f);
        for line in f.lines().map(|x| x.unwrap()) {
            if line.to_uppercase().contains("PANIC") {
                panic!("panic requested in read_output");
            } else if line.to_uppercase().contains("ERROR") {
                panic!("error found in {}, exiting", self.filename);
            } else if line.contains(" == MOPAC DONE ==") {
                return self.read_aux();
            }
        }
        None
    }

    /// return the heat of formation from a MOPAC aux file in Hartrees
    fn read_aux(&self) -> Option<f64> {
        let auxfile = Self::with_ext(&self.filename, ".aux");
        let f = if let Ok(file) = File::open(auxfile) {
            file
        } else {
            return None;
        };
        let f = BufReader::new(f);
        for line in f.lines().map(|x| x.unwrap()) {
            // line like HEAT_OF_FORMATION:KCAL/MOL=+0.97127947459164715838D+02
            if line.contains("HEAT_OF_FORMATION") {
                let fields: Vec<&str> = line.split("=").collect();
                match fields[1].replace("D", "E").parse::<f64>() {
                    Ok(f) => return Some(f / KCALHT),
                    Err(_) => {
                        eprintln!("failed to parse {}", fields[1]);
                        return None;
                    }
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use crate::queue::*;

    fn test_mopac() -> Mopac {
        Mopac {
            filename: String::from("/tmp/test"),
            params: vec![
                Param::new("USS", "H", -11.246958000000),
                Param::new("ZS", "H", 1.268641000000),
                Param::new("BETAS", "H", -8.352984000000),
                Param::new("GSS", "H", 14.448686000000),
                Param::new("USS", "C", -51.089653000000),
                Param::new("UPP", "C", -39.937920000000),
                Param::new("ZS", "C", 2.047558000000),
                Param::new("ZP", "C", 1.702841000000),
                Param::new("BETAS", "C", -15.385236000000),
                Param::new("BETAP", "C", -7.471929000000),
                Param::new("GSS", "C", 13.335519000000),
                Param::new("GPP", "C", 10.778326000000),
                Param::new("GSP", "C", 11.528134000000),
                Param::new("GP2", "C", 9.486212000000),
                Param::new("HSP", "C", 0.717322000000),
            ],
            geom: Vec::new(),
        }
    }

    #[test]
    fn test_write_input() {
        test_mopac().write_input("params.dat");
        let got = fs::read_to_string("/tmp/test.mop").expect("file not found");
        let want = "XYZ 1SCF A0 scfcrt=1.D-21 aux(precision=14) PM6 external=params.dat
Comment line 1
Comment line 2

";
        assert_eq!(got, want);
        fs::remove_file("/tmp/test.mop").unwrap();
    }

    #[test]
    fn test_write_params() {
        test_mopac().write_params("/tmp/params.dat");
        let got =
            fs::read_to_string("/tmp/params.dat").expect("file not found");
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
        fs::remove_file("/tmp/params.dat").unwrap();
    }

    #[test]
    fn test_read_output() {
        let mp = Mopac {
            filename: String::from("job"),
            geom: Vec::new(),
            params: Vec::new(),
        };
        let got = mp.read_output().expect("expected a value");
        let want = 0.97127947459164715838e+02 / KCALHT;
        assert!((dbg!(got) - want).abs() < 1e-20);
    }

    /// minimal queue for testing general submission
    struct TestQueue;

    impl Submit for TestQueue {
        fn write_submit_script(&self, infiles: Vec<String>) {
            let mut body = String::new();
            for f in infiles {
                body.push_str(&format!("echo {f}\n"));
            }
            let mut file = File::create(self.filename())
                .expect("failed to create params file");
            write!(file, "{}", body).expect("failed to write params file");
        }

        fn filename(&self) -> &str {
            "/tmp/main.pbs"
        }

        fn submit_command(&self) -> &str {
            "bash"
        }
    }

    #[test]
    fn test_submit() {
        let tq = TestQueue;
        tq.write_submit_script(
            vec!["input1.mop", "input2.mop", "input3.mop"]
                .iter()
                .map(|x| x.to_string())
                .collect(),
        );
        let got = tq.submit();
        let want = "input1.mop\ninput2.mop\ninput3.mop\n";
        assert_eq!(got, want);
    }
}
