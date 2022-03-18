use std::{
    fs::{self, File},
    io::BufRead,
    io::BufReader,
    io::Write,
};

pub mod mopac;
pub mod queue;

static DIRS: &'static [&str] = &["inp", "tmparam"];

/// set up the directories needed for the program
pub fn setup() {
    for dir in DIRS {
        match fs::create_dir(dir) {
            Ok(_) => (),
            Err(_) => (),
        }
    }
}

// TODO don't do this if -nodel flag
pub fn takedown() {
    for dir in DIRS {
        match fs::remove_dir_all(dir) {
            Ok(_) => (),
            Err(_) => (),
        }
    }
}

#[derive(Debug)]
pub struct Atom {
    pub label: String,
    pub coord: Vec<f64>,
}

impl Atom {
    pub fn new(label: &str, coord: Vec<f64>) -> Self {
        Self {
            label: label.to_string(),
            coord,
        }
    }
}

impl ToString for Atom {
    fn to_string(&self) -> String {
        format!(
            "{:2} {:15.10} {:15.10} {:15.10}",
            self.label, self.coord[0], self.coord[1], self.coord[2]
        )
    }
}

pub fn geom_string(geom: &Vec<Atom>) -> String {
    use std::fmt::Write;
    let mut ret = String::new();
    for g in geom {
        write!(ret, "{}\n", g.to_string()).unwrap();
    }
    ret
}

/// Take a file containing geometries like:
///   # GEOMUP #################
///           0.0000000000        0.0000000000       -1.6794733900
///           0.0000000000        1.2524327590        0.6959098120
///           0.0000000000       -1.2524327590        0.6959098120
///           0.0000000000        3.0146272390        1.7138963510
///           0.0000000000       -3.0146272390        1.7138963510
///   # GEOMUP #################
///           0.0000000000        0.0000000000       -1.6929508795
///           0.0000000000        1.2335354991        0.6923003326
///           0.0000000000       -1.2335354991        0.6923003326
///           0.0000000000        2.9875928126        1.7242445752
///           0.0000000000       -2.9875928126        1.7242445752
///   # GEOMUP #################
///           0.0000000000        0.0000000000       -1.6826636598
///           0.0000000000        1.2382598141        0.6926064201
///           0.0000000000       -1.2382598141        0.6926064201
///           0.0000000000        2.9956906742        1.7187948778
///           0.0000000000       -2.9956906742        1.7187948778
/// and split it into separate vectors for each geometry
pub fn load_geoms(filename: &str, atom_names: Vec<&str>) -> Vec<Vec<Atom>> {
    let f = File::open(filename).expect("failed to open three07");
    let lines = BufReader::new(f).lines();
    let mut ret = Vec::new();
    let mut buf = Vec::new();
    let mut idx: usize = 0;
    for (i, line) in lines.map(|x| x.unwrap()).enumerate() {
        if line.contains("# GEOM") {
            if i > 0 {
                ret.push(buf);
                buf = Vec::new();
                idx = 0;
            }
            continue;
        }
        let fields = line
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect();
        buf.push(Atom {
            label: atom_names[idx].to_string(),
            coord: fields,
        });
        idx += 1;
    }
    ret.push(buf);
    ret
}

/// parse a file containing lines like:
///
///   USS            H    -11.246958000000
///   ZS             H      1.268641000000
///
/// into a vec of Params
pub fn load_params(filename: &str) -> Vec<mopac::Param> {
    let mut ret = Vec::new();
    let f = File::open(filename).expect("failed to open three07");
    let lines = BufReader::new(f).lines();
    for line in lines.map(|x| x.unwrap()) {
        let fields: Vec<&str> = line.split_whitespace().collect();
        ret.push(mopac::Param::new(
            fields[0],
            fields[1],
            fields[2].parse().unwrap(),
        ))
    }
    ret
}

/// Slurm is a type for holding the information for submitting a slurm job.
/// `filename` is the name of the Slurm submission script
pub struct Slurm<'a> {
    filename: &'a str,
}

impl<'a> Slurm<'a> {
    pub fn new(filename: &'a str) -> Self {
        Slurm { filename }
    }
}

impl<'a> queue::Submit for Slurm<'a> {
    fn write_submit_script(&self, infiles: Vec<String>) {
        let mut body = String::from(
            "#!/bin/bash
#SBATCH --job-name=semp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o {filename}.out
#SBATCH --no-requeue
#SBATCH --mem=1gb
export LD_LIBRARY_PATH=/home/qc/mopac2016/\n",
        );
        for f in infiles {
            body.push_str(&format!(
                "/home/qc/mopac2016/MOPAC2016.exe {f}.mop\n"
            ));
        }
        let mut file =
            File::create(self.filename).expect("failed to create params file");
        write!(file, "{}", body).expect("failed to write params file");
    }

    fn filename(&self) -> &str {
        self.filename
    }

    fn submit_command(&self) -> &str {
        "sbatch"
    }
}

/// Minimal implementation for testing MOPAC locally
pub struct LocalQueue<'a> {
    filename: &'a str,
}

impl<'a> LocalQueue<'a> {
    pub fn new(filename: &'a str) -> Self {
        Self { filename }
    }
}

impl<'a> queue::Submit for LocalQueue<'a> {
    fn write_submit_script(&self, infiles: Vec<String>) {
        let mut body = String::from("export LD_LIBRARY_PATH=/opt/mopac/\n");
        for f in infiles {
            body.push_str(&format!("/opt/mopac/mopac {f}.mop\n"));
            body.push_str(&format!("date +%s\n"));
        }
        let mut file =
            File::create(self.filename).expect("failed to create params file");
        write!(file, "{}", body).expect("failed to write params file");
    }

    fn filename(&self) -> &str {
        self.filename
    }

    fn submit_command(&self) -> &str {
        "bash"
    }
}
#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::fs;

    use super::*;
    use crate::mopac::*;
    use crate::queue::*;

    #[test]
    fn test_load_geoms() {
        let got = load_geoms("three07", vec!["C", "C", "C", "H", "H"]);
        let want = vec![
            vec![
                Atom::new("C", vec![0.0000000000, 0.0000000000, -1.6794733900]),
                Atom::new("C", vec![0.0000000000, 1.2524327590, 0.6959098120]),
                Atom::new("C", vec![0.0000000000, -1.2524327590, 0.6959098120]),
                Atom::new("H", vec![0.0000000000, 3.0146272390, 1.7138963510]),
                Atom::new("H", vec![0.0000000000, -3.0146272390, 1.7138963510]),
            ],
            vec![
                Atom::new("C", vec![0.0000000000, 0.0000000000, -1.6929508795]),
                Atom::new("C", vec![0.0000000000, 1.2335354991, 0.6923003326]),
                Atom::new("C", vec![0.0000000000, -1.2335354991, 0.6923003326]),
                Atom::new("H", vec![0.0000000000, 2.9875928126, 1.7242445752]),
                Atom::new("H", vec![0.0000000000, -2.9875928126, 1.7242445752]),
            ],
            vec![
                Atom::new("C", vec![0.0000000000, 0.0000000000, -1.6826636598]),
                Atom::new("C", vec![0.0000000000, 1.2382598141, 0.6926064201]),
                Atom::new("C", vec![0.0000000000, -1.2382598141, 0.6926064201]),
                Atom::new("H", vec![0.0000000000, 2.9956906742, 1.7187948778]),
                Atom::new("H", vec![0.0000000000, -2.9956906742, 1.7187948778]),
            ],
        ];
        assert!(comp_geoms(got, want, 1e-10));
    }

    fn comp_geoms(got: Vec<Vec<Atom>>, want: Vec<Vec<Atom>>, eps: f64) -> bool {
        for (i, mol) in got.iter().enumerate() {
            for (j, atom) in mol.iter().enumerate() {
                let btom = &want[i][j];
                if atom.label != btom.label {
                    return false;
                }
                for k in 0..atom.coord.len() {
                    if (atom.coord[k] - btom.coord[k]).abs() > eps {
                        return false;
                    }
                }
            }
        }
        true
    }

    #[test]
    fn test_load_params() {
        let got = load_params("params.dat");
        let want = vec![
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
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn test_write_submit_script() {
        Slurm {
            filename: "/tmp/submit.slurm",
        }
        .write_submit_script(vec![
            String::from("input1"),
            String::from("input2"),
            String::from("input3"),
        ]);
        let got = fs::read_to_string("/tmp/submit.slurm")
            .expect("failed to read /tmp/submit.slurm");
        let want = "#!/bin/bash
#SBATCH --job-name=semp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o {filename}.out
#SBATCH --no-requeue
#SBATCH --mem=1gb
export LD_LIBRARY_PATH=/home/qc/mopac2016/
/home/qc/mopac2016/MOPAC2016.exe input1.mop
/home/qc/mopac2016/MOPAC2016.exe input2.mop
/home/qc/mopac2016/MOPAC2016.exe input3.mop
";
        assert_eq!(got, want);
        fs::remove_file("/tmp/submit.slurm").unwrap();
    }
    #[test]
    fn test_full() {
        let names = vec!["C", "C", "C", "H", "H"];
        let moles = load_geoms("three07", names);
        let params = load_params("params.dat");
        // build jobs (in memory) -> write jobs (to disk) -> run jobs
        let mut count: usize = 0;
        let mut jobs = Vec::new();
        for mol in moles {
            let filename = format!("inp/job.{count:08}");
            count += 1;
            jobs.push(Job::new(Mopac::new(filename, params.clone(), mol)))
        }
        setup();
        // slurm_file has to be set here so you can associate it with each job
        // in the chunk.
        let slurm_file = "inp/main.slurm";
        let mut chunk = Vec::new();
        for job in &mut jobs {
            job.mopac.write_input();
            job.pbs_file = slurm_file.to_string();
            chunk.push(job.mopac.filename.clone());
        }
        let mut slurm_jobs = HashMap::new();
        slurm_jobs.insert(slurm_file, chunk.len());
        let slurm = LocalQueue::new(slurm_file);
        slurm.write_submit_script(chunk);
        // run jobs. TODO loop over jobs in this chunk and set job.job_id to the
        // return value of slurm.submit
        slurm.submit();
        // collect output
        let mut got = Vec::new();
        for job in jobs {
            got.push(job.mopac.read_output().unwrap());
            let count = slurm_jobs.get_mut(job.pbs_file.as_str()).unwrap();
            *count -= 1;
        }
        let want = vec![
            0.20374485388911504,
            0.20541305733965845,
            0.20511337069030972,
        ];
        let eps = 1e-20;
        for i in 0..got.len() {
            assert!((got[i] - want[i]).abs() < eps);
        }
        takedown();
    }
}

pub struct Job {
    pub mopac: mopac::Mopac,
    pub pbs_file: String,
    pub job_id: String,
}

impl Job {
    pub fn new(mopac: mopac::Mopac) -> Self {
        Self {
            mopac,
            pbs_file: String::new(),
            job_id: String::new(),
        }
    }
}
