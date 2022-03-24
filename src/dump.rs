pub(crate) struct Dump {
    pub(crate) buf: Vec<String>,
    pub(crate) ptr: usize,
    pub(crate) max: usize,
}

impl Dump {
    pub(crate) fn new(size: usize) -> Self {
        Self {
            buf: vec![String::new(); size],
            ptr: 0,
            max: size,
        }
    }
    pub(crate) fn add(&mut self, files: Vec<String>) {
        for file in files {
            if self.ptr == self.max - 1 {
                self.dump();
            }
            self.buf[self.ptr] = file;
            self.ptr += 1;
        }
    }

    fn dump(&mut self) {
        for file in &self.buf {
            let _ = std::fs::remove_file(file);
        }
        self.ptr = 0;
    }
}
