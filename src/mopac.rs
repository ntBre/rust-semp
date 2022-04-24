use crate::queue::{Program, ProgramStatus};

use super::*;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufReader, Write};

/// kcal/mol per hartree
const KCALHT: f64 = 627.5091809;

#[derive(Debug, Clone)]
pub struct Params {
    pub names: Vec<String>,
    pub atoms: Vec<String>,
    pub values: na::DVector<f64>,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            names: Default::default(),
            atoms: Default::default(),
            values: na::DVector::from(vec![0.; 0]),
        }
    }
}

impl ToString for Params {
    fn to_string(&self) -> String {
        let mut ret = String::new();
        for (i, n) in self.names.iter().enumerate() {
            ret.push_str(
                &format!(
                    "{:<8}{:>8}{:20.12}\n",
                    n, self.atoms[i], self.values[i]
                )
                .to_string(),
            );
        }
        ret
    }
}

impl PartialEq for Params {
    fn eq(&self, other: &Self) -> bool {
        for (i, n) in self.names.iter().enumerate() {
            if *n != other.names[i] {
                #[cfg(test)]
                eprintln!("{}: {} != {}", i, *n, other.names[i]);
                return false;
            }
            if self.atoms[i] != other.atoms[i] {
                #[cfg(test)]
                eprintln!("{}: {} != {}", i, self.atoms[i], other.atoms[i]);
                return false;
            }
            let diff = (self.values[i] - other.values[i]).abs();
            if diff >= 1e-12 {
                #[cfg(test)]
                eprintln!(
                    "{}: {} != {}, diff = {}",
                    i, self.values[i], other.values[i], diff
                );
                return false;
            }
        }
        true
    }
}

impl Params {
    pub fn new(
        names: Vec<String>,
        atoms: Vec<String>,
        values: na::DVector<f64>,
    ) -> Self {
        Self {
            names,
            atoms,
            values,
        }
    }
    pub fn from(
        names: Vec<String>,
        atoms: Vec<String>,
        values: Vec<f64>,
    ) -> Self {
        Self {
            names,
            atoms,
            values: na::DVector::from(values),
        }
    }

    pub fn from_literal(
        names: Vec<&str>,
        atoms: Vec<&str>,
        values: Vec<f64>,
    ) -> Self {
        Self {
            names: names.iter().map(|s| s.to_string()).collect(),
            atoms: atoms.iter().map(|s| s.to_string()).collect(),
            values: na::DVector::from(values),
        }
    }

    pub fn len(&self) -> usize {
        assert_eq!(self.names.len(), self.atoms.len());
        assert_eq!(self.names.len(), self.values.len());
        self.names.len()
    }
}

/// Mopac holds the information needed to write a MOPAC input file. `filename`
/// should not include an extension. `.mop` will be appended for input files,
/// and `.out` and `.aux` will be appended for output files.
#[derive(Debug, Clone)]
pub struct Mopac {
    pub filename: String,
    /// The semi-empirical parameters to use in the calculation via the EXTERNAL
    /// keyword. These are wrapped in an Rc to allow the same set of parameters
    /// to be shared between calculations without an expensive `clone`
    /// operation.
    pub params: Rc<Params>,
    /// The initial geometry for the calculation. These are also wrapped in an
    /// Rc to avoid allocating multiple copies for calculations with the same
    /// geometry.
    pub geom: Rc<Vec<Atom>>,
    pub param_file: String,
    pub param_dir: String,
    pub charge: isize,
    // TODO add field for resub - linked list of jobs
}

impl Program for Mopac {
    fn filename(&self) -> String {
        self.filename.clone()
    }

    fn set_filename(&mut self, filename: &str) {
        self.filename = String::from(filename);
    }

    fn extension(&self) -> String {
        String::from("mop")
    }

    /// Writes the parameters of self to a parameter file, then writes the MOPAC
    /// input file with external=paramfile. Also update self.paramfile to point
    /// to the generated name for the parameter file
    fn write_input(&mut self) {
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
            "XYZ 1SCF A0 scfcrt=1.D-21 aux(precision=14) PM6 \
	     external={paramfile} charge={charge}
Comment line 1
Comment line 2
{geom}
",
            paramfile = self.param_file,
            charge = self.charge,
        )
        .expect("failed to write input file");
    }

    /// Reads a MOPAC output file. If normal termination occurs, also try
    /// reading the `.aux` file to extract the energy from there. This function
    /// panics if an error is found in the output file. If a non-fatal error
    /// occurs (file not found, not written to yet, etc) None is returned.
    fn read_output(&self) -> ProgramStatus {
        let outfile = format!("{}.out", &self.filename);
        let f = match File::open(&outfile) {
            Ok(file) => file,
            Err(_) => {
                return ProgramStatus::FileNotFound;
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
                eprintln!("panic requested in read_output");
                std::process::exit(1)
            } else if let Some(_) = line.find("ERROR") {
                return ProgramStatus::ErrorInOutput;
            } else if let Some(_) = line.find(" == MOPAC DONE ==") {
                return self.read_aux();
            }
            line.clear();
        }
        ProgramStatus::EnergyNotFound
    }

    fn associated_files(&self) -> Vec<String> {
        let fname = self.filename();
        vec![
            format!("{}.mop", fname),
            format!("{}.out", fname),
            format!("{}.arc", fname),
            format!("{}.aux", fname),
            self.param_file.clone(),
        ]
    }

    fn charge(&self) -> isize {
        self.charge
    }
}

impl Mopac {
    pub fn new(
        filename: String,
        params: Rc<Params>,
        geom: Rc<Vec<Atom>>,
        charge: isize,
    ) -> Self {
        Self {
            filename,
            params,
            geom,
            param_file: String::new(),
            param_dir: "tmparam".to_string(),
            charge,
        }
    }

    fn write_params(&self) {
        let mut body = String::new();
        body.push_str(&self.params.to_string());
        let mut file = match File::create(&self.param_file) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("failed to create {} with {}", self.param_file, e);
                std::process::exit(1);
            }
        };
        write!(file, "{}", body).expect("failed to write params file");
    }

    /// return the heat of formation from a MOPAC aux file in Hartrees
    fn read_aux(&self) -> ProgramStatus {
        let auxfile = format!("{}.aux", &self.filename);
        let f = if let Ok(file) = File::open(auxfile) {
            file
        } else {
            return ProgramStatus::FileNotFound;
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
                    Ok(f) => return ProgramStatus::Success(f / KCALHT),
                    Err(_) => {
                        return ProgramStatus::EnergyParseError;
                    }
                }
            }
            line.clear();
        }
        ProgramStatus::EnergyNotFound
    }
}

#[cfg(test)]
mod tests {
    use std::fs::{self, read_to_string};

    use super::*;

    fn test_mopac() -> Mopac {
        let names = vec![
            "USS", "ZS", "BETAS", "GSS", "USS", "UPP", "ZS", "ZP", "BETAS",
            "BETAP", "GSS", "GPP", "GSP", "GP2", "HSP",
        ];
        let atoms = vec![
            "H", "H", "H", "H", "C", "C", "C", "C", "C", "C", "C", "C", "C",
            "C", "C",
        ];
        let values = vec![
            -11.246958000000,
            1.268641000000,
            -8.352984000000,
            14.448686000000,
            -51.089653000000,
            -39.937920000000,
            2.047558000000,
            1.702841000000,
            -15.385236000000,
            -7.471929000000,
            13.335519000000,
            10.778326000000,
            11.528134000000,
            9.486212000000,
            0.717322000000,
        ];
        Mopac::new(
            String::from("/tmp/test"),
            Rc::new(Params::from(
                names.iter().map(|s| s.to_string()).collect(),
                atoms.iter().map(|s| s.to_string()).collect(),
                values,
            )),
            Rc::new(Vec::new()),
            0,
        )
    }

    #[test]
    fn test_write_input() {
        let mut tm = test_mopac();
        tm.param_dir = "/tmp".to_string();
        tm.write_input();
        let got = fs::read_to_string("/tmp/test.mop").expect("file not found");
        let want = format!("XYZ 1SCF A0 scfcrt=1.D-21 aux(precision=14) PM6 external={paramfile} charge=0
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
        let mp = Mopac::new(
            String::from("test_files/job"),
            Rc::new(Params::default()),
            Rc::new(Vec::new()),
            0,
        );
        let got = if let ProgramStatus::Success(v) = mp.read_output() {
            v
        } else {
            panic!()
        };
        let want = 0.97127947459164715838e+02 / KCALHT;
        assert!((got - want).abs() < 1e-20);

        // failure in output
        let mp = Mopac::new(
            String::from("test_files/nojob"),
            Rc::new(Params::default()),
            Rc::new(Vec::new()),
            0,
        );
        let got = mp.read_output();
        assert_eq!(got, ProgramStatus::EnergyNotFound);

        // failure in aux
        let mp = Mopac::new(
            String::from("test_files/noaux"),
            Rc::new(Params::default()),
            Rc::new(Vec::new()),
            0,
        );
        let got = mp.read_output();
        assert_eq!(got, ProgramStatus::EnergyNotFound);
    }

    /// minimal queue for testing general submission
    struct TestQueue;

    impl Queue<Mopac> for TestQueue {
        fn write_submit_script(&self, infiles: &[String], filename: &str) {
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

        fn chunk_size(&self) -> usize {
            128
        }

        fn job_limit(&self) -> usize {
            1600
        }

        fn sleep_int(&self) -> usize {
            1
        }

        const SCRIPT_EXT: &'static str = "pbs";

        const DIR: &'static str = "inp";

        fn stat_cmd(&self) -> String {
            todo!()
        }
    }

    #[test]
    fn test_submit() {
        let tq = TestQueue;
        tq.write_submit_script(
            &string!["input1.mop", "input2.mop", "input3.mop"],
            "/tmp/main.pbs",
        );
        let got = tq.submit("/tmp/main.pbs");
        let want = "input3.mop";
        assert_eq!(got, want);
    }

    #[test]
    fn test_resubmit() {
        use std::path::Path;
        let tq = TestQueue;
        std::fs::copy("test_files/job.mop", "/tmp/job.mop").unwrap();
        let got = tq.resubmit("/tmp/job.mop");
        assert!(Path::new("/tmp/job_redo.mop").exists());
        assert!(Path::new("/tmp/job_redo.pbs").exists());
        assert_eq!(
            read_to_string("/tmp/job.mop").unwrap(),
            read_to_string("/tmp/job_redo.mop").unwrap()
        );
        let want = queue::Resubmit {
            inp_file: String::from("/tmp/job_redo"),
            pbs_file: String::from("/tmp/job_redo.pbs"),
            job_id: String::from("/tmp/job_redo.mop"),
        };
        assert_eq!(got, want);

        for f in vec!["/tmp/job.mop", "/tmp/job_redo.mop", "/tmp/job_redo.pbs"]
        {
            std::fs::remove_file(f).unwrap();
        }
    }
}
