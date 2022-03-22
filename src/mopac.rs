use super::*;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
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
#[derive(Debug, Clone)]
pub struct Mopac {
    pub filename: String,
    pub params: Vec<Param>,
    pub geom: Vec<Atom>,
    pub param_file: String,
    pub param_dir: String,
}

impl Mopac {
    pub fn new(filename: String, params: Vec<Param>, geom: Vec<Atom>) -> Self {
        Self {
            filename,
            params,
            geom,
            param_file: String::new(),
            param_dir: "tmparam".to_string(),
        }
    }

    fn write_params(&self) {
        let mut body = String::new();
        for p in &self.params {
            body.push_str(&p.to_string());
        }
        let mut file = match File::create(&self.param_file) {
            Ok(f) => f,
            Err(e) => panic!("failed to create {} with {}", self.param_file, e),
        };
        write!(file, "{}", body).expect("failed to write params file");
    }

    /// Writes the parameters of self to a parameter file, then writes the MOPAC
    /// input file with external=paramfile. Also update self.paramfile to point
    /// to the generated name for the parameter file
    pub fn write_input(&mut self) {
        let mut s = DefaultHasher::new();
        self.filename.hash(&mut s);
        let paramfile = format!("{}/{}", self.param_dir, s.finish());
        self.param_file = paramfile;
        self.write_params();
        let geom = geom_string(&self.geom);
        let mut file = File::create(format!("{}.mop", self.filename))
            .expect("failed to create input file");
        write!(
            file,
            "XYZ 1SCF A0 scfcrt=1.D-21 aux(precision=14) PM6 external={paramfile}
Comment line 1
Comment line 2
{geom}
", paramfile=&self.param_file)
        .expect("failed to write input file");
    }

    /// Reads a MOPAC output file. If normal termination occurs, also try
    /// reading the `.aux` file to extract the energy from there. This function
    /// panics if an error is found in the output file. If a non-fatal error
    /// occurs (file not found, not written to yet, etc) None is returned.
    pub fn read_output(&self) -> Option<f64> {
        let outfile = format!("{}.out", &self.filename);
        let f = match File::open(&outfile) {
            Ok(file) => file,
            Err(_) => {
                return None;
            } // file not found
        };
        let mut f = BufReader::new(f);
        let mut line = String::new();
        while let Ok(b) = f.read_line(&mut line) {
            if b == 0 {
                break;
            }
            line.make_ascii_uppercase();
            if let Some(_) = line.find("PANIC") {
                panic!("panic requested in read_output");
            } else if let Some(_) = line.find("ERROR") {
                panic!("error found in {}, exiting", self.filename);
            } else if let Some(_) = line.find(" == MOPAC DONE ==") {
                return self.read_aux();
            }
            line.clear();
        }
        None
    }

    /// return the heat of formation from a MOPAC aux file in Hartrees
    fn read_aux(&self) -> Option<f64> {
        let auxfile = format!("{}.aux", &self.filename);
        let f = if let Ok(file) = File::open(auxfile) {
            file
        } else {
            return None;
        };
        let mut f = BufReader::new(f);
        let mut line = String::new();
        while let Ok(b) = f.read_line(&mut line) {
            if b == 0 {
                break;
            }
            // line like HEAT_OF_FORMATION:KCAL/MOL=+0.97127947459164715838D+02
            if line.contains("HEAT_OF_FORMATION") {
                let fields: Vec<&str> = line.trim().split("=").collect();
                match fields[1].replace("D", "E").parse::<f64>() {
                    Ok(f) => return Some(f / KCALHT),
                    Err(_) => {
                        eprintln!("failed to parse '{}'", fields[1]);
                        return None;
                    }
                }
            }
            line.clear();
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    fn test_mopac() -> Mopac {
        Mopac::new(
            String::from("/tmp/test"),
            vec![
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
            Vec::new(),
        )
    }

    #[test]
    fn test_write_input() {
        let mut tm = test_mopac();
        tm.param_dir = "/tmp".to_string();
        tm.write_input();
        let got = fs::read_to_string("/tmp/test.mop").expect("file not found");
        let want = format!("XYZ 1SCF A0 scfcrt=1.D-21 aux(precision=14) PM6 external={paramfile}
Comment line 1
Comment line 2

", paramfile=tm.param_file);
        assert_eq!(got, want);
        fs::remove_file("/tmp/test.mop").unwrap();
    }

    #[test]
    fn test_write_params() {
        let mut tm = test_mopac();
        tm.param_file = String::from("/tmp/params.dat");
        tm.write_params();
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
        // success
        let mp = Mopac::new(String::from("job"), Vec::new(), Vec::new());
        let got = mp.read_output().expect("expected a value");
        let want = 0.97127947459164715838e+02 / KCALHT;
        assert!((got - want).abs() < 1e-20);

        // failure in output
        let mp = Mopac::new(String::from("nojob"), Vec::new(), Vec::new());
        let got = mp.read_output();
        assert_eq!(got, None);

        // failure in aux
        let mp = Mopac::new(String::from("noaux"), Vec::new(), Vec::new());
        let got = mp.read_output();
        assert_eq!(got, None);
    }

    /// minimal queue for testing general submission
    struct TestQueue;

    impl Submit for TestQueue {
        fn write_submit_script(&self, infiles: Vec<String>, filename: &str) {
            let mut body = String::new();
            for f in infiles {
                body.push_str(&format!("echo {f}\n"));
            }
            let mut file =
                File::create(filename).expect("failed to create params file");
            write!(file, "{}", body).expect("failed to write params file");
        }

        fn submit_command(&self) -> &str {
            "bash"
        }

        fn new() -> Self {
            Self {}
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
            "/tmp/main.pbs",
        );
        let got = tq.submit("/tmp/main.pbs");
        let want = "input1.mop\ninput2.mop\ninput3.mop";
        assert_eq!(got, want);
    }
}
