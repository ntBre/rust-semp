use serde::Deserialize;

#[derive(Deserialize, Debug, PartialEq)]
pub struct Config {
    pub job_limit: usize,
    pub chunk_size: usize,
    pub sleep_int: usize,
    pub max_iter: usize,
    pub atom_names: Vec<String>,
    pub params: String,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let contents = std::fs::read_to_string(filename)
            .expect("failed to load config file");
        toml::from_str(&contents).expect("failed to deserialize config file")
    }
}

#[cfg(test)]
mod tests {
    use crate::string;

    use super::*;

    #[test]
    fn test_load_full() {
        let got = Config::load("test_files/test.toml");
        let want = Config {
            job_limit: 10000,
            chunk_size: 128,
            sleep_int: 5,
            max_iter: 5,
            atom_names: string!["C", "C", "C", "H", "H"],
            params: String::from(
                "USS        H    -11.246958000000
ZS         H      1.268641000000
BETAS      H     -8.352984000000
GSS        H     14.448686000000
USS        C    -51.089653000000
UPP        C    -39.937920000000
ZS         C      2.047558000000
ZP         C      1.702841000000
BETAS      C    -15.385236000000
BETAP      C     -7.471929000000
GSS        C     13.335519000000
GPP        C     10.778326000000
GSP        C     11.528134000000
GP2        C      9.486212000000
HSP        C      0.717322000000
FN11       C      0.046302000000
",
            ),
        };
        assert_eq!(got, want);
    }
}
