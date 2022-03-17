use std::{path::Path, process::Command};

pub trait Submit {
    fn write_submit_script(&self, infiles: Vec<&str>);

    fn filename(&self) -> &str;

    fn submit_command(&self) -> &str;

    fn submit(&self) {
        let path = Path::new(self.filename());
        let dir = match path.parent() {
            Some(d) => d,
            None => Path::new("."),
        };
        let base = path.file_name().expect("expected filename");
        let output = match Command::new(self.submit_command())
            .current_dir(dir)
            .arg(base)
            .output()
        {
            Ok(s) => s,
            Err(_) => todo!(),
        };
        dbg!(path, dir, base, output);
    }
}
