use std::{process::Command, str};

pub trait Submit {
    fn new() -> Self;

    fn write_submit_script(&self, infiles: Vec<String>, filename: &str);

    fn submit_command(&self) -> &str;

    fn submit(&self, filename: &str) -> String {
        match Command::new(self.submit_command()).arg(filename).output() {
            Ok(s) => {
                return str::from_utf8(&s.stdout).unwrap().trim().to_string()
            }
            Err(_) => todo!(),
        };
    }
}
