//! Reader for the `.scorefreq` file (frozen minor allele + frequency + LD-pruned
//! flag from a reference full run), mirroring `VcfMarkerData.readScoreFile` /
//! `ScoreFreqTool`. Consuming a scorefreq file lets new individuals be scored
//! against the reference run's frozen allele frequencies and LD pruning, so a
//! focus run is equivalent to the full merged run for the scored pairs.
//!
//! Format (tab-delimited): `CHROM POS ID REF ALT ALLELE FREQ LD_PRUNED`.
//! `LD_PRUNED` is 0 (retained) or 1 (LD-pruned, exclusion-only). A 7-column file
//! (no `LD_PRUNED`) is read as all-retained, matching the Java legacy behavior.

use std::collections::HashSet;
use std::io::Read;

pub struct ScoreRecord {
    pub chrom: String,
    pub pos: i32,
    pub ref_a: String,
    pub alt_a: String,
    pub allele: String, // the scored (minor) allele; must equal REF or ALT in the VCF
    pub freq: f32,
    pub correlated: bool, // LD_PRUNED == 1
}

fn is_header(f: &[&str]) -> bool {
    f.len() >= 7
        && f[0] == "CHROM"
        && f[1] == "POS"
        && f[2] == "ID"
        && f[3] == "REF"
        && f[4] == "ALT"
        && f[5] == "ALLELE"
        && f[6] == "FREQ"
}

/// Reads and validates a scorefreq file (plain or gzipped, auto-detected by
/// magic bytes). Records are returned in file order; duplicate (CHROM,POS,REF,ALT)
/// keys are rejected.
pub fn read_score_file(path: &str) -> Vec<ScoreRecord> {
    let raw = std::fs::read(path).unwrap_or_else(|e| panic!("cannot read scorefreq {}: {}", path, e));
    let text = if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut s = String::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_string(&mut s)
            .expect("scorefreq gzip read error");
        s
    } else {
        String::from_utf8(raw).expect("scorefreq is not valid UTF-8")
    };

    let mut records: Vec<ScoreRecord> = Vec::new();
    let mut seen: HashSet<(String, i32, String, String)> = HashSet::new();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if is_header(&f) {
            continue;
        }
        if f.len() != 7 && f.len() != 8 {
            panic!(
                "scorefreq line {}: expected 7 or 8 tab-delimited fields, found {}",
                lineno + 1,
                f.len()
            );
        }
        if f[0].is_empty() || f[3].is_empty() || f[4].is_empty() || f[5].is_empty() {
            panic!("scorefreq line {}: missing CHROM/REF/ALT/ALLELE field", lineno + 1);
        }
        let pos: i32 = f[1]
            .parse()
            .unwrap_or_else(|_| panic!("scorefreq line {}: invalid POS '{}'", lineno + 1, f[1]));
        if pos < 1 {
            panic!("scorefreq line {}: invalid POS '{}'", lineno + 1, f[1]);
        }
        let freq: f32 = f[6]
            .parse()
            .unwrap_or_else(|_| panic!("scorefreq line {}: invalid FREQ '{}'", lineno + 1, f[6]));
        if !(freq > 0.0 && freq < 1.0) {
            panic!(
                "scorefreq line {}: FREQ {} out of range (must be >0 and <1)",
                lineno + 1,
                freq
            );
        }
        let correlated = if f.len() == 8 {
            match f[7] {
                "0" => false,
                "1" => true,
                other => panic!(
                    "scorefreq line {}: invalid LD_PRUNED '{}' (must be 0 or 1)",
                    lineno + 1,
                    other
                ),
            }
        } else {
            false
        };
        let key = (f[0].to_string(), pos, f[3].to_string(), f[4].to_string());
        if !seen.insert(key) {
            panic!(
                "scorefreq line {}: duplicate marker {}:{} REF={} ALT={}",
                lineno + 1,
                f[0],
                pos,
                f[3],
                f[4]
            );
        }
        records.push(ScoreRecord {
            chrom: f[0].to_string(),
            pos,
            ref_a: f[3].to_string(),
            alt_a: f[4].to_string(),
            allele: f[5].to_string(),
            freq,
            correlated,
        });
    }
    if records.is_empty() {
        panic!("empty scorefreq file: {}", path);
    }
    records
}
