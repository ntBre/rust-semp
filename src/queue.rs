pub trait Submit {
    fn write_submit_script(&self, infiles: Vec<&str>);
    fn submit(&self);
}
