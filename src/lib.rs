use std::{fs::File, io::BufRead, io::BufReader, io::Write};

pub mod mopac;
pub mod queue;

#[derive(Debug)]
pub struct Atom {
    pub label: String,
    pub coord: Vec<f64>,
}

impl Atom {
    pub fn new(label: String, coord: Vec<f64>) -> Self {
        Self { label, coord }
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
    let ret = String::new();
    for g in geom {
        println!("{:2}", g.label);
        todo!()
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
pub fn load_geoms(filename: &str) -> Vec<Vec<f64>> {
    let f = File::open(filename).expect("failed to open three07");
    let lines = BufReader::new(f).lines();
    let mut ret = Vec::new();
    let mut buf = Vec::new();
    for (i, line) in lines.map(|x| x.unwrap()).enumerate() {
        if line.contains("# GEOM") {
            if i > 0 {
                ret.push(buf);
                buf = Vec::new();
            }
            continue;
        }
        let fields = line.split_whitespace();
        for f in fields {
            buf.push(f.parse::<f64>().unwrap());
        }
    }
    ret.push(buf);
    ret
}

/// `filename` is the name of the Slurm submission script
pub struct Slurm<'a> {
    filename: &'a str,
}

impl<'a> queue::Submit for Slurm<'a> {
    fn write_submit_script(&self, infiles: Vec<&str>) {
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
            body.push_str(&format!("/home/qc/mopac2016/MOPAC2016.exe {f}\n"));
        }
        let mut file =
            File::create(self.filename).expect("failed to create params file");
        write!(file, "{}", body).expect("failed to write params file");
    }

    fn submit(&self) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use crate::queue::*;

    #[test]
    fn test_load_geoms() {
        let got = load_geoms("three07");
        let want = vec![
            vec![
                0.0000000000,
                0.0000000000,
                -1.6794733900,
                0.0000000000,
                1.2524327590,
                0.6959098120,
                0.0000000000,
                -1.2524327590,
                0.6959098120,
                0.0000000000,
                3.0146272390,
                1.7138963510,
                0.0000000000,
                -3.0146272390,
                1.7138963510,
            ],
            vec![
                0.0000000000,
                0.0000000000,
                -1.6929508795,
                0.0000000000,
                1.2335354991,
                0.6923003326,
                0.0000000000,
                -1.2335354991,
                0.6923003326,
                0.0000000000,
                2.9875928126,
                1.7242445752,
                0.0000000000,
                -2.9875928126,
                1.7242445752,
            ],
            vec![
                0.0000000000,
                0.0000000000,
                -1.6826636598,
                0.0000000000,
                1.2382598141,
                0.6926064201,
                0.0000000000,
                -1.2382598141,
                0.6926064201,
                0.0000000000,
                2.9956906742,
                1.7187948778,
                0.0000000000,
                -2.9956906742,
                1.7187948778,
            ],
        ];
        assert!(comp_geoms(got, want, 1e-10));
    }

    fn comp_geoms(got: Vec<Vec<f64>>, want: Vec<Vec<f64>>, eps: f64) -> bool {
        for (i, geom) in got.iter().enumerate() {
            for (j, coord) in geom.iter().enumerate() {
                if (coord - want[i][j]).abs() > eps {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn test_write_submit_script() {
        Slurm {
            filename: "/tmp/submit.slurm",
        }
        .write_submit_script(vec![
            "input1.mop",
            "input2.mop",
            "input3.mop",
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

    struct TestQueue<'a> {
        filename: &'a str,
    }

    impl<'a> Submit for TestQueue<'a> {
        fn write_submit_script(&self, infiles: Vec<&str>) {
            Slurm {
                filename: self.filename,
            }
            .write_submit_script(infiles);
        }

        fn submit(&self) {
            todo!()
        }
    }
}
