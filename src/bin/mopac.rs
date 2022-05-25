use rayon::prelude::*;
use std::{
    env::args,
    fs::File,
    io::{BufRead, BufReader},
    process::Command,
};

fn main() {
    let infile = args().nth(1).unwrap();
    let f = File::open(infile).unwrap();
    let lines = BufReader::new(f).lines().flatten();
    let mut args = Vec::new();
    for line in lines {
        if line.contains("export") || line.contains("date") {
            continue;
        }
        args.push(line.split_whitespace().nth(1).unwrap().to_string());
    }
    args.iter().for_each(|arg| {
        Command::new("/opt/mopac/mopac")
            .env("LD_LIBRARY_PATH", "/opt/mopac")
            .arg(arg)
            .output()
            .unwrap();
    });
}
