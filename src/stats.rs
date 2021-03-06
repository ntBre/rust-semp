use nalgebra as na;

#[derive(Debug)]
pub struct Stats {
    pub norm: f64,
    pub rmsd: f64,
    pub max: f64,
}

impl Stats {
    /// compute the Stats between `v` and `w` in cm-1 under the assumption they
    /// started out in Ht
    pub fn new(v: &na::DVector<f64>, w: &na::DVector<f64>, conv: f64) -> Self {
        let count = v.len();
        assert_eq!(count, w.len());
        let mut sq_diffs = 0.0;
        let mut max = v[0] - w[0];
        for i in 0..count {
            let diff = v[i] - w[i];
            sq_diffs += diff * diff;
            if diff.abs() > max {
                max = diff.abs();
            }
        }
        Self {
            norm: sq_diffs.sqrt() * conv,
            rmsd: (sq_diffs / count as f64).sqrt() * conv,
            max: max * conv,
        }
    }
    pub fn print_header() {
        println!(
            "{:>17}{:>12}{:>12}{:>12}{:>12}{:>12}",
            "cm-1", "cm-1", "cm-1", "cm-1", "cm-1", "s"
        );
        println!(
            "{:>5}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}",
            "Iter", "Norm", "ΔNorm", "RMSD", "ΔRMSD", "Max", "Time"
        );
    }

    pub fn print_step(&self, iter: usize, last: &Self, time_milli: u128) {
        print!(
            "{:5}{:12.4}{:12.4}{:12.4}{:12.4}{:12.4}{:12.1}\n",
            iter,
            self.norm,
            self.norm - last.norm,
            self.rmsd,
            self.rmsd - last.rmsd,
            self.max,
            time_milli as f64 / 1000.,
        );
    }
}

impl Default for Stats {
    fn default() -> Self {
        Self {
            norm: 0.0,
            rmsd: 0.0,
            max: 0.0,
        }
    }
}

impl PartialEq for Stats {
    fn eq(&self, other: &Self) -> bool {
        fn close(a: f64, b: f64, eps: f64) -> bool {
            if (a - b).abs() < eps {
                return true;
            }
            false
        }
        let eps = 1e-4;
        close(self.norm, other.norm, eps)
            && close(self.rmsd, other.rmsd, eps)
            && close(self.max, other.max, eps)
    }
}
