use std::{fmt::Display, ops::Add, str::FromStr};

use crate::{driver::Params, Dvec};

/// A basis set
#[derive(Clone, Debug, PartialEq)]
pub struct MolproParams {
    col1: Vec<String>,
    col2: Vec<String>,
    rest: Vec<Vec<f64>>,
}

impl FromStr for MolproParams {
    type Err = std::string::ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut col1 = Vec::new();
        let mut col2 = Vec::new();
        let mut rest = Vec::new();
        for line in s.lines() {
            let sp: Vec<_> = line
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .collect();
            if sp.is_empty() || sp[0].starts_with('!') {
                continue;
            }
            let mut sp = sp.into_iter();
            col1.push(sp.next().unwrap().to_string());
            col2.push(sp.next().unwrap().to_string());
            let mut tmp = Vec::new();
            for s in sp {
                tmp.push(s.parse().unwrap());
            }
            rest.push(tmp);
        }
        Ok(Self { col1, col2, rest })
    }
}

impl Display for MolproParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for ((c1, c2), rest) in self.col1.iter().zip(&self.col2).zip(&self.rest)
        {
            write!(f, "{c1}, {c2}, ")?;
            for (i, r) in rest.iter().enumerate() {
                write!(f, "{r:.10E}")?;
                if i < rest.len() - 1 {
                    write!(f, ", ")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl Params for MolproParams {
    fn incr_value(&mut self, idx: usize, delta: f64) {
        let mut cur = 0;
        for vs in self.rest.iter_mut() {
            for v in vs.iter_mut() {
                if cur == idx {
                    *v += delta;
                    return;
                }
                cur += 1;
            }
        }
    }

    /// Sum all of the basis set entries from each row
    fn len(&self) -> usize {
        self.rest.iter().map(|v| v.len()).sum()
    }
}

impl Add<&Dvec> for MolproParams {
    type Output = Self;

    fn add(self, _rhs: &Dvec) -> Self::Output {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use std::fs::read_to_string;

    use crate::string;

    use super::*;

    #[test]
    fn parse_params() {
        let got = MolproParams::from_str(
            &read_to_string("test_files/molpro.params").unwrap(),
        )
        .unwrap();
        let want = MolproParams {
            col1: string!["s", "c", "c", "s", "c", "c", "c", "p", "c", "c"],
            col2: string![
                "H", "1.2", "3.3", "C", "1.3", "4.5", "6.6", "C", "1.2", "3.3"
            ],
            rest: vec![
                vec![0.5447178000E+01, 0.8245472400E+00, 0.1831915800E+00],
                vec![0.1562849787E+00, 0.9046908767E+00],
                vec![1.0000000],
                vec![
                    0.1722560000E+03,
                    0.2591090000E+02,
                    0.5533350000E+01,
                    0.3664980000E+01,
                    0.7705450000E+00,
                    0.1958570000E+00,
                ],
                vec![0.6176690738E-01, 0.3587940429E+00, 0.7007130837E+00],
                vec![-0.3958951621E+00, 0.1215834356E+01],
                vec![0.1000000000E+01],
                vec![0.3664980000E+01, 0.7705450000E+00, 0.1958570000E+00],
                vec![0.2364599466E+00, 0.8606188057E+00],
                vec![0.1000000000E+01],
            ],
        };
        assert_eq!(got, want);
    }

    #[test]
    fn write_params() {
        let params = MolproParams::from_str(
            &read_to_string("test_files/molpro.params").unwrap(),
        )
        .unwrap();
        let got = params.to_string();
        let want = "s, H, 5.4471780000E0, 8.2454724000E-1, 1.8319158000E-1
c, 1.2, 1.5628497870E-1, 9.0469087670E-1
c, 3.3, 1.0000000000E0
s, C, 1.7225600000E2, 2.5910900000E1, 5.5333500000E0, 3.6649800000E0, 7.7054500000E-1, 1.9585700000E-1
c, 1.3, 6.1766907380E-2, 3.5879404290E-1, 7.0071308370E-1
c, 4.5, -3.9589516210E-1, 1.2158343560E0
c, 6.6, 1.0000000000E0
p, C, 3.6649800000E0, 7.7054500000E-1, 1.9585700000E-1
c, 1.2, 2.3645994660E-1, 8.6061880570E-1
c, 3.3, 1.0000000000E0
";
        assert_eq!(got, want);
    }
}
