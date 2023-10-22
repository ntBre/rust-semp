use std::fmt::Display;

use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Job, Program, ProgramResult, Template},
    queue::Queue,
};
use serde::{Deserialize, Serialize};

use crate::Dvec;

pub trait Params {
    fn new(names: Vec<String>, atoms: Vec<String>, values: Dvec) -> Self;
    fn names(&self) -> &Vec<String>;
    fn atoms(&self) -> &Vec<String>;
    fn values(&self) -> &Dvec;

    /// increment the `idx`th parameter value by `delta`
    fn incr_value(&mut self, idx: usize, delta: f64);

    /// decrement the `idx`th parameter value by `delta`. automatically
    /// implemented using `incr_value`
    fn decr_value(&mut self, idx: usize, delta: f64) {
        self.incr_value(idx, -delta);
    }

    /// returns the length of `self.names()`. assumes `names`, `atoms`, and
    /// `values` have the same length
    fn len(&self) -> usize {
        self.names().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Params for psqs::program::mopac::Params {
    fn names(&self) -> &Vec<String> {
        &self.names
    }

    fn atoms(&self) -> &Vec<String> {
        &self.atoms
    }

    fn values(&self) -> &Dvec {
        &self.values
    }

    fn new(names: Vec<String>, atoms: Vec<String>, values: Dvec) -> Self {
        psqs::program::mopac::params::Params::new(names, atoms, values)
    }

    fn incr_value(&mut self, idx: usize, delta: f64) {
        self.values[idx] += delta;
    }
}

/// [Driver] is an extension of [psqs::program::Program] with the additional
/// methods needed to call `run_algo`.
pub trait Driver:
    Program + Clone + Sync + Send + Serialize + for<'a> Deserialize<'a>
{
    type Params: Params + Display + Clone;

    fn optimize_geometry<Q: Queue<Self> + Sync>(
        geom: Geom,
        params: &Self::Params,
        queue: &Q,
        dir: &str,
        name: &str,
        charge: isize,
        template: Template,
    ) -> Option<ProgramResult>;

    fn write_params(
        job_num: usize,
        params: &Self::Params,
        template: Template,
    ) -> Template;

    #[allow(clippy::too_many_arguments)]
    fn build_jobs(
        moles: Vec<Geom>,
        params: Option<&Self::Params>,
        dir: &'static str,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>>;

    fn new_full(
        filename: String,
        params: Option<Self::Params>,
        geom: Geom,
        charge: isize,
        template: Template,
    ) -> Self;
}

impl Driver for Mopac {
    type Params = psqs::program::mopac::Params;

    fn optimize_geometry<Q: Queue<Self> + Sync>(
        geom: Geom,
        params: &Self::Params,
        queue: &Q,
        dir: &str,
        name: &str,
        charge: isize,
        template: Template,
    ) -> Option<ProgramResult> {
        let opt = Job::new(
            Mopac::new_full(
                format!("{dir}/{name}"),
                Some(params.clone()),
                geom,
                charge,
                template,
            ),
            0,
        );
        let mut res = vec![Default::default(); 1];
        let status = queue.energize(dir, vec![opt], &mut res);
        if status.is_err() {
            return None;
        }
        Some(res.pop().unwrap())
    }

    /// generate a parameter filename using `job_num`, write the parameters to
    /// that file, and return the template with the parameter filename added
    fn write_params(
        job_num: usize,
        params: &<Self as Driver>::Params,
        mut template: Template,
    ) -> Template {
        let param_file = format!("tmparam/{job_num}.dat");
        Mopac::write_params(params, &param_file);
        use std::fmt::Write;
        write!(template.header, " external={param_file}").unwrap();
        template
    }

    fn build_jobs(
        moles: Vec<Geom>,
        params: Option<&Self::Params>,
        dir: &'static str,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>> {
        Mopac::build_jobs(
            moles,
            params,
            dir,
            start_index,
            coeff,
            job_num,
            charge,
            tmpl,
        )
    }

    fn new_full(
        filename: String,
        params: Option<Self::Params>,
        geom: Geom,
        charge: isize,
        template: Template,
    ) -> Self {
        Mopac::new_full(filename, params, geom, charge, template)
    }
}
