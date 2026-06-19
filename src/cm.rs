//! Port of `ibdseq.CmConverter`: bp -> cM conversion and segment-length bins.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

const BASES_PER_MB: f64 = 1_000_000.0;

pub struct CmConverter {
    positions: Option<Vec<f64>>,
    centimorgans: Option<Vec<f64>>,
    bins: Vec<f64>,
    bin_labels: Vec<String>,
    cm_per_mb: f64,
}

impl CmConverter {
    pub fn new(map_file: Option<&str>, bins_arg: &str, cm_per_mb: f64) -> Self {
        let bins = parse_bins(bins_arg);
        let bin_labels = bins.iter().map(|b| fmt_label(*b)).collect();
        match map_file {
            None => CmConverter {
                positions: None,
                centimorgans: None,
                bins,
                bin_labels,
                cm_per_mb,
            },
            Some(f) => {
                let (p, c) = read_map(f);
                CmConverter {
                    positions: Some(p),
                    centimorgans: Some(c),
                    bins,
                    bin_labels,
                    cm_per_mb,
                }
            }
        }
    }

    pub fn cm(&self, pos: i32) -> f64 {
        let positions = match &self.positions {
            None => return (pos as f64 / BASES_PER_MB) * self.cm_per_mb,
            Some(p) => p,
        };
        let cms = self.centimorgans.as_ref().unwrap();
        let x = pos as f64;
        if x <= positions[0] {
            return cms[0];
        }
        let last = positions.len() - 1;
        if x >= positions[last] {
            return cms[last];
        }
        match positions.binary_search_by(|p| p.partial_cmp(&x).unwrap()) {
            Ok(i) => cms[i],
            Err(hi) => {
                let lo = hi - 1;
                let scale = (x - positions[lo]) / (positions[hi] - positions[lo]);
                cms[lo] + scale * (cms[hi] - cms[lo])
            }
        }
    }

    pub fn bin(&self, length_cm: f64) -> String {
        let hi = match self.bins.binary_search_by(|b| b.partial_cmp(&length_cm).unwrap()) {
            Ok(i) => i + 1,
            Err(i) => i,
        };
        if hi == 0 || hi >= self.bins.len() {
            panic!(
                "Segment cM length outside bins: {}. Increase bins upper bound.",
                length_cm
            );
        }
        format!("{}-{}", self.bin_labels[hi - 1], self.bin_labels[hi])
    }
}

fn parse_bins(arg: &str) -> Vec<f64> {
    let v: Vec<f64> = arg
        .split(',')
        .map(|s| s.trim().parse::<f64>().expect("invalid bins value"))
        .collect();
    assert!(v.len() >= 2, "bins must contain at least two values");
    for j in 1..v.len() {
        assert!(v[j] > v[j - 1], "bins values must be strictly increasing");
    }
    v
}

/// Mimic Java Double.toString for the common integer-valued bin edges
/// (e.g. 0.0, 1.0, 3000.0); falls back to the shortest representation otherwise.
fn fmt_label(b: f64) -> String {
    if b.fract() == 0.0 && b.abs() < 1e15 {
        format!("{:.1}", b)
    } else {
        let s = format!("{}", b);
        if s.contains('.') {
            s
        } else {
            format!("{}.0", s)
        }
    }
}

fn read_map(path: &str) -> (Vec<f64>, Vec<f64>) {
    let file = File::open(path).expect("cannot open map file");
    let reader: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut lines = BufReader::new(reader).lines();
    let header = lines.next().expect("empty map file").unwrap();
    let hdr: Vec<&str> = header.split('\t').map(|s| s.trim()).collect();
    let pos_col = hdr.iter().position(|c| *c == "Position(bp)").expect("map missing Position(bp)");
    let cm_col = hdr.iter().position(|c| *c == "Map(cM)").expect("map missing Map(cM)");
    let mut positions = Vec::new();
    let mut cms = Vec::new();
    for line in lines {
        let line = line.unwrap();
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let p: f64 = f[pos_col].trim().parse().unwrap();
        let c: f64 = f[cm_col].trim().parse().unwrap();
        positions.push(p);
        cms.push(c);
    }
    assert!(!positions.is_empty(), "map file contains no rows");
    (positions, cms)
}
