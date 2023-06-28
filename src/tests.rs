use std::{
    fs::{self, File},
    io::{BufRead, BufReader},
};

use std::str::FromStr;

use crate::{
    config::Config,
    optimize::energy::Energy,
    stats::Stats,
    utils::{load_energies, load_params},
};

use psqs::{
    geom::Geom,
    queue::{local::Local, slurm::Slurm},
};
use symm::Molecule;

use super::*;

#[test]
fn test_load_geoms() {
    let got = utils::load_geoms(
        "test_files/three07",
        &string!["C", "C", "C", "H", "H"],
    );
    let moles = [
        Molecule::from_str(
            "
        C 0.0000000000        0.0000000000       -1.6794733900
        C 0.0000000000        1.2524327590        0.6959098120
        C 0.0000000000       -1.2524327590        0.6959098120
        H 0.0000000000        3.0146272390        1.7138963510
        H 0.0000000000       -3.0146272390        1.7138963510
",
        )
        .unwrap(),
        Molecule::from_str(
            "
        C 0.0000000000        0.0000000000       -1.6929508795
        C 0.0000000000        1.2335354991        0.6923003326
        C 0.0000000000       -1.2335354991        0.6923003326
        H 0.0000000000        2.9875928126        1.7242445752
        H 0.0000000000       -2.9875928126        1.7242445752
",
        )
        .unwrap(),
        Molecule::from_str(
            "
        C 0.0000000000        0.0000000000       -1.6826636598
        C 0.0000000000        1.2382598141        0.6926064201
        C 0.0000000000       -1.2382598141        0.6926064201
        H 0.0000000000        2.9956906742        1.7187948778
        H 0.0000000000       -2.9956906742        1.7187948778
",
        )
        .unwrap(),
    ];
    let want = moles
        .into_iter()
        .map(|mut m| {
            m.to_angstrom();
            Geom::Xyz(m.atoms)
        })
        .collect();
    assert!(comp_geoms(got, want, 1e-10));
}

fn comp_vec(got: &[f64], want: &[f64], eps: f64) -> bool {
    if got.len() != want.len() {
        return false;
    }
    for (i, g) in got.iter().enumerate() {
        let diff = (g - want[i]).abs();
        if diff > eps {
            eprintln!(
                "{:5} got: {}, wanted {}: diff {:e}",
                i, g, want[i], diff
            );
            return false;
        }
    }
    true
}

fn comp_geoms(got: Vec<Geom>, want: Vec<Geom>, eps: f64) -> bool {
    for (i, mol) in got.iter().enumerate() {
        let mol = mol.xyz().unwrap();
        for (j, atom) in mol.iter().enumerate() {
            let btom = &want[i].xyz().unwrap()[j];
            if atom.atomic_number != btom.atomic_number {
                return false;
            }
            if !comp_vec(&atom.coord(), &btom.coord(), eps) {
                return false;
            }
        }
    }
    true
}

#[test]
fn test_load_energies() {
    let got = load_energies("test_files/small.dat");
    let want = na::DVector::from(vec![
        0.000000000000,
        0.000453458157,
        0.000258906149,
        0.000268926371,
        0.000247426244,
        0.000251632099,
        0.000259236055,
        0.000267395682,
        0.000273547427,
        0.000159696616,
        0.000136655279,
        0.000118039859,
        0.000119823282,
        0.000124994640,
        0.000136055791,
        0.000177574442,
        0.000124620943,
        0.000126842659,
        0.000132448769,
        0.000108785467,
        0.000107629675,
        0.000176541210,
        0.000144535949,
        0.000126349619,
        0.000129459053,
        0.000141660093,
        0.000176137923,
        0.000127547307,
        0.000128059983,
    ]);
    assert!(comp_dvec(got, want, 1e-12));
}

#[test]
fn test_load_params() {
    let got = load_params("test_files/params.dat");
    #[rustfmt::skip]
	let want = Params::from_literal(
            vec![
                "USS", "ZS", "BETAS", "GSS", "USS", "UPP", "ZS", "ZP", "BETAS",
                "BETAP", "GSS", "GPP", "GSP", "GP2", "HSP",
            ],
            vec![
                "H", "H", "H", "H", "C", "C", "C", "C", "C", "C", "C", "C",
                "C", "C", "C",
            ],
            vec![
                -11.246958000000, 1.268641000000, -8.352984000000,
                14.448686000000, -51.089653000000, -39.937920000000,
                2.047558000000, 1.702841000000, -15.385236000000,
                -7.471929000000, 13.335519000000, 10.778326000000,
                11.528134000000, 9.486212000000, 0.717322000000,
            ],
        );
    assert_eq!(got, want);
}

#[test]
fn test_write_submit_script() {
    <Slurm as Queue<Mopac>>::write_submit_script(
        &Slurm::new(3, 3, 2, ".", false, None),
        &string!["input1", "input2", "input3"],
        "/tmp/submit.slurm",
    );
    let got = fs::read_to_string("/tmp/submit.slurm")
        .expect("failed to read /tmp/submit.slurm");
    let want = "#!/bin/bash
#SBATCH --job-name=semp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /tmp/submit.slurm.out
#SBATCH --no-requeue
#SBATCH --mem=1gb
export LD_LIBRARY_PATH=/home/qc/mopac2016/
echo $SLURM_JOB_ID
date
hostname
/home/qc/mopac2016/MOPAC2016.exe input1.mop
/home/qc/mopac2016/MOPAC2016.exe input2.mop
/home/qc/mopac2016/MOPAC2016.exe input3.mop
";

    assert_eq!(got, want);
    fs::remove_file("/tmp/submit.slurm").unwrap();
}

#[test]
fn test_one_iter() {
    let config = Config::load("test_files/test.toml");
    let names = string!["C", "C", "C", "H", "H"];
    let moles = utils::load_geoms("test_files/three07", &names);
    let params = utils::load_params("test_files/params.dat");
    let want = na::DVector::from(vec![
        0.0,
        0.0016682056863255301,
        0.0013685188198143683,
    ]);
    let got = Energy { moles }.semi_empirical(
        &params,
        &Local {
            chunk_size: 128,
            dir: "inp".to_string(),
            ..Default::default()
        },
        &config.molecules,
        want.len(),
    );
    let eps = 4e-8;
    assert!(comp_dvec(got.unwrap(), want, eps));
}

/// load a matrix from a file. each line becomes a row in the resulting matrix.
/// NOTE: assumes the rows are numbered, so it skips the first field of each
/// line
pub(crate) fn load_mat(filename: &str) -> na::DMatrix<f64> {
    let f = match File::open(filename) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("failed to open {filename} with {e}");
            std::process::exit(1);
        }
    };
    let lines = BufReader::new(f).lines().map(|x| x.unwrap());
    let mut ret = Vec::new();
    let mut rows = 0;
    for line in lines {
        let sp: Vec<f64> = line
            .split_whitespace()
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        ret.extend_from_slice(&sp);
        rows += 1;
    }
    let cols = ret.len() / rows;
    let ret = transpose(&ret, rows, cols);
    na::DMatrix::from_vec(rows, cols, ret)
}

fn transpose(orig: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut ret = Vec::with_capacity(rows * cols);
    for i in 0..cols {
        for j in 0..rows {
            ret.push(orig[i + j * cols])
        }
    }
    ret
}

/// report whether or not the euclidean norm of the difference between `got`
/// and `want` is less than epsilon
fn comp_dvec(got: na::DVector<f64>, want: na::DVector<f64>, eps: f64) -> bool {
    let norm = (got - want).norm();
    if norm < eps {
        return true;
    }
    eprintln!("comp_vec: norm = {norm:e}");
    false
}

fn comp_mat(got: na::DMatrix<f64>, want: na::DMatrix<f64>, eps: f64) -> bool {
    let norm = (got - want).norm();
    if norm < eps {
        return true;
    }
    eprintln!("comp_mat: norm = {norm}");
    false
}

#[ignore]
#[test]
fn test_num_jac() {
    let config = crate::config::Config::load("test_files/test.toml");
    let names = string!["C", "C", "C", "H", "H"];
    {
        let moles = utils::load_geoms("test_files/small07", &names);
        let params = utils::load_params("test_files/small.params");
        let want = load_mat("test_files/small.jac");
        let got = Energy { moles }.num_jac(
            &params,
            &Local {
                dir: "inp".to_owned(),
                chunk_size: 128,
                mopac: "/opt/mopac/mopac".to_owned(),
            },
            &config.molecules,
            9,
        );
        assert!(comp_mat(got, want, 2e-5));
    }
    {
        // want jac straight from the Go version
        let moles = utils::load_geoms("test_files/three07", &names);
        let params = utils::load_params("test_files/three.params");
        let want = load_mat("test_files/three.jac");
        let got = Energy { moles }.num_jac(
            &params,
            &Local {
                dir: "inp".to_owned(),
                chunk_size: 128,
                mopac: "/opt/mopac/mopac".to_owned(),
            },
            &config.molecules,
            9,
        );
        let tol = 4.2e-6;
        assert!(comp_mat(got, want, tol));
    }
}

#[test]
fn test_norm() {
    // making sure the nalgebra vector norm is what I think it is
    let v = na::DVector::<f64>::from(vec![1.0, 2.0, 3.0]);
    let got = v.norm();
    let want = f64::sqrt(1.0 + 4.0 + 9.0);
    assert_eq!(got, want);
}

#[test]
fn test_broyden() {
    let jac = load_mat("test_files/three.jac");
    let se_old =
        na::DVector::from(vec![0.000000000000, 0.001668203451, 0.001368516801]);
    let se_new =
        na::DVector::from(vec![0.000000000000, 0.000375747088, 0.000113449727]);
    let step = na::DVector::from(vec![
        0.879074952563,
        -0.053372779748,
        -0.386443051094,
        -1.723741422927,
        0.085500727064,
        -0.083919956047,
        -0.006440778691,
        -0.008037008147,
        -0.161085646344,
        -0.106676930451,
        0.083518030316,
        0.454558059306,
        0.599560620927,
        -0.065880893329,
        0.944895856082,
    ]);
    let got = broyden_update(&jac, &se_old, &se_new, &step);
    let want = load_mat("test_files/broyden.jac");
    assert!(comp_mat(got, want, 3e-8));
}

#[test]
fn test_lev_mar() {
    // test input taken from the output of go TestLevMar
    #[rustfmt::skip]
        let jac = load_mat("test_files/levmar.jac");
    let got = lev_mar(
        &jac,
        &na::DVector::from(vec![
            0.000000000000,
            0.000453458157,
            0.000258906149,
        ]),
        &na::DVector::from(vec![
            0.203744853890,
            0.205413057340,
            0.205113370691,
        ]),
        1e-8,
    );
    let want = na::DVector::from(vec![
        0.8790749525625414,
        -0.053372779747915,
        -0.38644305109436294,
        -1.7237414229267645,
        0.08550072706391243,
        -0.08391995604728157,
        -0.006440778691258624,
        -0.008037008146972191,
        -0.161085646343529,
        -0.10667693045073838,
        0.08351803031611686,
        0.45455805930613674,
        0.599560620927178,
        -0.06588089332916275,
        0.9448958560822209,
    ]);
    assert!(comp_dvec(got, want, 3.2e-7));
}

#[test]
fn test_solve() {
    // try nalgebra matrix equation solution with input from Go version
    #[rustfmt::skip]
        let lhs = load_mat("test_files/go.lhs");
    let rhs = na::dvector![
        0.35423900602931,
        -0.35423663730942,
        -0.35423845555583,
        -0.35423646488059,
        0.35423942795133,
        -0.35423939669513,
        -0.35423686878446,
        -0.35423575920428,
        -0.35423924460534,
        -0.35423965293655,
        0.35423629270772,
        -0.05025294054723,
        0.35423682777338,
        -0.35423906998699,
        0.35423447419595
    ];
    let want = na::dvector![
        0.022821308059,
        -0.019004630312,
        -0.030395950124,
        -0.022549225832,
        0.011349968251,
        -0.013318764149,
        -0.006731748989,
        -0.017577820182,
        -0.018559932276,
        -0.026805396091,
        0.009825237417,
        0.000558292732,
        0.052824051563,
        -0.029828987570,
        0.072728316382
    ];
    let got = match na::linalg::Cholesky::new(lhs.clone()) {
        Some(a) => {
            let lhs = a;
            lhs.solve(&rhs)
        }
        None => {
            eprintln!("cholesky decomposition failed");
            let lhs = na::linalg::LU::new(lhs);
            lhs.solve(&rhs).expect("LU decomposition also failed")
        }
    };
    assert!(comp_dvec(got, want, 1.04e-6));
}

pub(crate) fn hostname() -> String {
    String::from_utf8(
        std::process::Command::new("hostname")
            .output()
            .expect("failed to get hostname")
            .stdout,
    )
    .unwrap()
    .trim()
    .to_string()
}

#[ignore]
#[test]
fn test_algo() {
    // loading everything
    let want = match hostname().as_str() {
        "cactus" => Stats {
            norm: 17.68655029777878,
            rmsd: 3.537310059555756,
            max: 7.453759816083651,
        },
        "bonsai" | "lintel" => Stats {
            norm: 7.356730516402555,
            rmsd: 1.471346103280511,
            max: 3.4654645038349137,
        },
        "keystone" => Stats {
            norm: 9.253577484008394,
            rmsd: 1.850715496801679,
            max: 3.4100922731486816,
        },
        _ => {
            // github actions answer because hostname is hard there
            Stats {
                norm: 7.898535349571516,
                rmsd: 1.5797070699143032,
                max: 4.1090995758655415,
            }
        }
    };
    let config = Config::load("test_files/test.toml");
    let queue = Local {
        dir: "inp".to_owned(),
        chunk_size: 128,
        mopac: "/opt/mopac/mopac".to_owned(),
    };
    let geom_file = "test_files/small07";
    let param_file = "test_files/small.params";
    let energy_file = "test_files/25.dat";
    let got = run_algo(
        &mut std::io::sink(),
        &config.molecules,
        utils::load_params(param_file),
        utils::load_energies(energy_file),
        5,
        true,
        5,
        queue,
        false,
        Energy {
            moles: utils::load_geoms(
                geom_file,
                &config.molecules[0].atom_names,
            ),
        },
    );
    assert_eq!(got, want);
}
