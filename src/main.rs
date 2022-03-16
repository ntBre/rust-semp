use std::{fs::File, io::BufRead, io::BufReader};

fn main() {
    let f = File::open("../three07").expect("failed to open three07");
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
    dbg!(ret);
}
