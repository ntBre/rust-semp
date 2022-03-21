use std::{process::Command, str};

pub trait Submit {
    fn write_submit_script(&self, infiles: Vec<String>);

    fn filename(&self) -> &str;

    fn submit_command(&self) -> &str;

    fn submit(&self) -> String {
        match Command::new(self.submit_command())
            .arg(self.filename())
            .output()
        {
            Ok(s) => {
                return str::from_utf8(&s.stdout).unwrap().trim().to_string()
            }
            Err(_) => todo!(),
        };
    }
}
