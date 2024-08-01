use std::fmt::Display;

use psqs::{
    geom::Geom,
    program::{
        dftbplus::DFTBPlus, molpro::Molpro, mopac::Mopac, Job, Program,
        ProgramResult, Template,
    },
    queue::Queue,
};
use serde::{Deserialize, Serialize};

use crate::params::{
    dftbplus::DFTBPlusParams, molpro::MolproParams, mopac::MopacParams,
};

pub trait Params {
    /// increment the `idx`th parameter value by `delta`
    fn incr_value(&mut self, idx: usize, delta: f64);

    /// decrement the `idx`th parameter value by `delta`. automatically
    /// implemented using `incr_value`
    fn decr_value(&mut self, idx: usize, delta: f64) {
        self.incr_value(idx, -delta);
    }

    /// returns the length of `self.names()`. assumes `names`, `atoms`, and
    /// `values` have the same length
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Params for psqs::program::mopac::Params {
    fn incr_value(&mut self, idx: usize, delta: f64) {
        self.values[idx] += delta;
    }

    fn len(&self) -> usize {
        self.names.len()
    }
}

/// [Driver] is an extension of [psqs::program::Program] with the additional
/// methods needed to call `run_algo`.
pub trait Driver:
    Program + Clone + Sync + Send + Serialize + for<'a> Deserialize<'a>
{
    type Params: Params
        + Display
        + Clone
        + for<'a> std::ops::Add<&'a crate::Dvec>;

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
        dir: &'static str,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>>;

    fn new_full(
        filename: String,
        params: Self::Params,
        geom: Geom,
        charge: isize,
        template: Template,
    ) -> Self;
}

impl Driver for Mopac {
    type Params = MopacParams;

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
                Some(params.0.clone()),
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
        Mopac::write_params(&params.0, &param_file);
        use std::fmt::Write;
        write!(template.header, " external={param_file}").unwrap();
        template
    }

    fn build_jobs(
        moles: Vec<Geom>,
        dir: &'static str,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>> {
        Mopac::build_jobs(
            moles,
            None,
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
        params: Self::Params,
        geom: Geom,
        charge: isize,
        template: Template,
    ) -> Self {
        Mopac::new_full(filename, Some(params.0), geom, charge, template)
    }
}

impl Driver for Molpro {
    type Params = MolproParams;

    fn optimize_geometry<Q: Queue<Self> + Sync>(
        geom: Geom,
        params: &Self::Params,
        queue: &Q,
        dir: &str,
        name: &str,
        charge: isize,
        mut template: Template,
    ) -> Option<ProgramResult> {
        template.header =
            template.header.replace("{{.basis}}", &params.to_string());
        let opt = Job::new(
            Molpro::new(format!("{dir}/{name}"), template, charge, geom),
            0,
        );
        let mut res = vec![Default::default(); 1];
        let status = queue.energize(dir, vec![opt], &mut res);
        if status.is_err() {
            return None;
        }
        Some(res.pop().unwrap())
    }

    fn write_params(
        _job_num: usize,
        params: &Self::Params,
        mut template: Template,
    ) -> Template {
        template.header =
            template.header.replace("{{.basis}}", &params.to_string());
        template
    }

    fn build_jobs(
        moles: Vec<Geom>,
        dir: &'static str,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>> {
        <Molpro as Program>::build_jobs(
            moles,
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
        params: Self::Params,
        geom: Geom,
        charge: isize,
        mut template: Template,
    ) -> Self {
        // TODO do I need to do this?
        template.header =
            template.header.replace("{{.basis}}", &params.to_string());
        Self::new(filename, template, charge, geom)
    }
}

#[allow(unused)]
impl Driver for DFTBPlus {
    type Params = DFTBPlusParams;

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
            DFTBPlus::new_full(
                format!("{dir}/{name}"),
                params.clone(),
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

    fn write_params(
        job_num: usize,
        params: &Self::Params,
        template: Template,
    ) -> Template {
        todo!()
    }

    fn build_jobs(
        moles: Vec<Geom>,
        dir: &'static str,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>> {
        todo!()
    }

    fn new_full(
        filename: String,
        params: Self::Params,
        geom: Geom,
        charge: isize,
        template: Template,
    ) -> Self {
        // TODO put params into template somehow
        Self::new(filename, template, charge, geom)
    }
}
