use rust_semp::{self, Atom};

fn main() {
    let names = vec!["C", "C", "C", "H", "H"];
    let coords = rust_semp::load_geoms("three07");
    let _atoms = to_geoms(names, coords);
}

/// each coord is a vector of length 3N, need to split into pieces of length 3
/// and pair them with atom labels
fn to_geoms(names: Vec<&str>, coords: Vec<Vec<f64>>) -> Vec<Vec<Atom>> {
    let mut geoms = Vec::new();
    for coord in coords {
        let mut geom = Vec::new();
        for (c, _) in coord.iter().enumerate().step_by(3) {
            geom.push(Atom::new(
                String::from(names[c / 3]),
                coord[c..c + 3].to_vec(),
            ));
        }
        geoms.push(geom);
    }
    geoms
}
