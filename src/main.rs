use rust_semp::mopac::*;

fn main() {}

// really params doesn't need to be a field of Mopac, each Mopac doesn't really
// have its own params, a group of Mopacs will always share a params. I guess I
// could have chunks contained within groups as well

struct Group {
    jobs: Vec<Mopac>,
    params: Vec<Param>,
}
