use nalgebra as na;
use psqs::geom::Geom;
use psqs::program::mopac::Params;
use std;
use std::cmp::Ordering;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::iter::zip;
use std::path::Path;
use std::rc::Rc;
use symm::atom::Atom;
use symm::Irrep;

/// from [StackOverflow](https://stackoverflow.com/a/45145246)
#[macro_export]
macro_rules! string {
    // match a list of expressions separated by comma:
    ($($str:expr),*) => ({
        // create a Vec with this list of expressions,
        // calling String::from on each:
        vec![$(String::from($str),)*] as Vec<String>
    });
}

/// Take an INTDER-style `file07` file and parse it into a Vec of geometries
pub fn load_geoms(filename: &str, atom_names: &[String]) -> Vec<Rc<Geom>> {
    let f = match File::open(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("failed to open {filename} with {e}");
            std::process::exit(1);
        }
    };
    let lines = BufReader::new(f).lines();
    let mut ret = Vec::new();
    let mut buf = Vec::new();
    let mut idx: usize = 0;
    for (i, line) in lines.map(|x| x.unwrap()).enumerate() {
        if line.contains("# GEOM") {
            if i > 0 {
                ret.push(Rc::new(Geom::Xyz(buf)));
                buf = Vec::new();
                idx = 0;
            }
            continue;
        }
        let fields: Vec<_> = line
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect();
        buf.push(Atom::new_from_label(
            &atom_names[idx],
            fields[0],
            fields[1],
            fields[2],
        ));
        idx += 1;
    }
    ret.push(Rc::new(Geom::Xyz(buf)));
    ret
}

pub fn load_energies(filename: &str) -> na::DVector<f64> {
    let mut ret = Vec::new();
    let f = match File::open(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "failed to open {filename} for reading energies with {e}"
            );
            std::process::exit(1);
        }
    };
    let lines = BufReader::new(f).lines();
    for line in lines.map(|x| x.unwrap()) {
        ret.push(line.trim().parse().unwrap());
    }
    na::DVector::from(ret)
}

/// Read `filename` into a string and then call [`parse_params`]
pub fn load_params(filename: &str) -> Params {
    let params = match std::fs::read_to_string(filename) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("failed to read {filename} with {e}");
            std::process::exit(1);
        }
    };
    parse_params(&params)
}

pub(crate) static DIRS: &[&str] = &["inp", "tmparam"];

/// set up the directories needed for the program after deleting existing ones
pub fn setup() {
    takedown();
    for dir in DIRS {
        match fs::create_dir(dir) {
            Ok(_) => (),
            Err(e) => {
                panic!("can't create '{}' for '{}'", dir, e);
            }
        }
    }
}

pub fn takedown() {
    'outer: for dir in DIRS {
        let path = Path::new(dir);
        if path.is_dir() {
            while let Err(e) = fs::remove_dir_all(dir) {
                // assume it's okay if the directory already doesn't exist
                if !Path::new(dir).exists() {
                    continue 'outer;
                }
                eprintln!("trying to remove '{}' with '{}'", dir, e);
                std::thread::sleep(std::time::Duration::from_secs(1));
            }
        }
    }
}

/// parse a string containing lines like:
///
/// ```text
///   USS            H    -11.246958000000
///   ZS             H      1.268641000000
/// ```
/// into a vec of Params
pub fn parse_params(params: &str) -> Params {
    let lines = params.split('\n');
    let mut names = Vec::new();
    let mut atoms = Vec::new();
    let mut values = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() == 3 {
            names.push(fields[0].to_string());
            atoms.push(fields[1].to_string());
            values.push(fields[2].parse().unwrap());
        }
    }
    Params::from(names, atoms, values)
}

pub fn dump_vec<W: Write>(w: &mut W, vec: &na::DVector<f64>) {
    for (i, v) in vec.iter().enumerate() {
        writeln!(w, "{:>5}{:>20.12}", i, v).unwrap();
    }
}

pub fn dump_mat<W: Write>(w: &mut W, mat: &na::DMatrix<f64>) {
    let (rows, cols) = mat.shape();
    writeln!(w).unwrap();
    for i in 0..rows {
        write!(w, "{:>5}", i).unwrap();
        for j in 0..cols {
            write!(w, "{:>12.8}", mat[(i, j)]).unwrap();
        }
        writeln!(w).unwrap();
    }
}

/// return `energies` relative to its minimum element
pub fn relative(energies: &na::DVector<f64>) -> na::DVector<f64> {
    let min = energies.min();
    let min = na::DVector::from(vec![min; energies.len()]);
    let ret = energies.clone();
    ret - min
}

pub(crate) fn log_params<W: Write>(w: &mut W, iter: usize, params: &Params) {
    let _ = writeln!(w, "Iter {}", iter);
    let _ = writeln!(w, "{}", params.to_string());
}

/// sort freqs by the irrep in the same position and then by frequency. if
/// `freqs` and `irreps` are not the same length, just return the sorted
/// frequencies.
pub fn sort_irreps(freqs: &[f64], irreps: &[Irrep]) -> Vec<f64> {
    if freqs.len() != irreps.len() {
        let mut ret = Vec::from(freqs);
        ret.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return ret;
    }
    let mut pairs: Vec<_> = zip(irreps, freqs).collect();
    pairs.sort_by(|a, b| match a.0.cmp(b.0) {
        Ordering::Equal => a.1.partial_cmp(b.1).unwrap(),
        other => other,
    });
    pairs.iter().map(|x| *x.1).collect()
}
