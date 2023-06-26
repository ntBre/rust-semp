use nalgebra as na;

#[derive(Debug)]
pub struct Stats {
    pub norm: f64,
    pub rmsd: f64,
    pub max: f64,
}

impl Stats {
    /// compute the Stats between `v` and `w` in cm-1 under the assumption they
    /// started out in Ht. `tru` is taken to be the "true" value, so length
    /// mismatches are resolved in its favor.
    pub fn new(
        tru: &na::DVector<f64>,
        w: &na::DVector<f64>,
        conv: f64,
    ) -> Self {
        let count = tru.len();
        let mut w = w.clone();
        if w.len() < count {
            w.extend(vec![0.0; count - w.len()]);
        }
        let mut sq_diffs = 0.0;
        let mut max = tru[0] - w[0];
        for i in 0..count {
            let wi = match w.get(i) {
                Some(i) => *i,
                None => {
                    eprintln!("failed to get w[{i}], replacing with 0.0");
                    0.0
                }
            };
            let diff = tru[i] - wi;
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
            "{:>17}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}",
            "cm-1", "cm-1", "cm-1", "cm-1", "cm-1", "arb", "s"
        );
        println!(
            "{:>5}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}{:>12}",
            "Iter", "Norm", "ΔNorm", "RMSD", "ΔRMSD", "Max", "λ", "Time"
        );
    }

    pub fn print_step(
        &self,
        iter: usize,
        last: &Self,
        time_milli: u128,
        lambda: f64,
    ) {
        println!(
            "{:5}{:12.4}{:12.4}{:12.4}{:12.4}{:12.4}{:12.3e}{:12.1}",
            iter,
            self.norm,
            self.norm - last.norm,
            self.rmsd,
            self.rmsd - last.rmsd,
            self.max,
            lambda,
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
        const EPS: f64 = 1e-4;
        fn close(a: f64, b: f64) -> bool {
            (a - b).abs() < EPS
        }
        close(self.norm, other.norm)
            && close(self.rmsd, other.rmsd)
            && close(self.max, other.max)
    }
}
