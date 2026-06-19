//! Streaming reader for a single-chromosome, GT-only VCF (as produced by the
//! workflow `get_vcf` rule). Computes per-marker allele stats matching
//! `vcf.MarkerData`. Biallelic sites only (errors on multiallelic ALT).

use crate::scorefreq::ScoreRecord;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct Markers {
    pub n_samples: usize,
    pub sample_ids: Vec<String>,
    pub chrom: String,
    pub pos: Vec<i32>,
    pub id: Vec<String>,
    pub ref_a: Vec<String>,
    pub alt_a: Vec<String>,
    pub minor_is_alt: Vec<bool>,
    pub minor_freq: Vec<f32>,
    pub major_freq: Vec<f32>,   // matches MarkerData.alleleFrequency(majorAllele)
    pub alt_dose: Vec<Vec<u8>>, // [marker][sample]; 0,1,2 ALT-allele dose, 3 = missing
}

impl Markers {
    pub fn n_markers(&self) -> usize {
        self.pos.len()
    }
}

pub fn read_vcf(path: &str, min_alleles: i32) -> Markers {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
    let mut reader = BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file));

    let mut line: Vec<u8> = Vec::with_capacity(1 << 16);
    let mut sample_ids: Vec<String> = Vec::new();

    // Header
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).expect("read error");
        if n == 0 {
            panic!("VCF has no #CHROM header");
        }
        if line.starts_with(b"##") {
            continue;
        }
        if line.starts_with(b"#CHROM") {
            let trimmed = trim_newline(&line);
            let cols: Vec<&[u8]> = trimmed.split(|&b| b == b'\t').collect();
            for c in cols.iter().skip(9) {
                sample_ids.push(String::from_utf8_lossy(c).into_owned());
            }
            break;
        }
        panic!("unexpected line before #CHROM header");
    }
    let n_samples = sample_ids.len();
    assert!(n_samples > 0, "no samples in VCF");

    let mut m = Markers {
        n_samples,
        sample_ids,
        chrom: String::new(),
        pos: Vec::new(),
        id: Vec::new(),
        ref_a: Vec::new(),
        alt_a: Vec::new(),
        minor_is_alt: Vec::new(),
        minor_freq: Vec::new(),
        major_freq: Vec::new(),
        alt_dose: Vec::new(),
    };

    let denom_alleles = 2 * n_samples;
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).expect("read error");
        if n == 0 {
            break;
        }
        let rec = trim_newline(&line);
        if rec.is_empty() {
            continue;
        }
        parse_record(rec, n_samples, denom_alleles, min_alleles, &mut m);
    }
    m
}

fn parse_record(
    rec: &[u8],
    n_samples: usize,
    denom_alleles: usize,
    min_alleles: i32,
    m: &mut Markers,
) {
    let mut it = rec.split(|&b| b == b'\t');
    let chrom = it.next().expect("missing CHROM");
    let pos = it.next().expect("missing POS");
    let id = it.next().expect("missing ID");
    let ref_a = it.next().expect("missing REF");
    let alt_a = it.next().expect("missing ALT");
    let _qual = it.next();
    let _filter = it.next();
    let _info = it.next();
    let _format = it.next();

    if alt_a.contains(&b',') {
        panic!(
            "multiallelic marker not supported (ALT={}) at pos {}",
            String::from_utf8_lossy(alt_a),
            String::from_utf8_lossy(pos)
        );
    }

    let mut alt_dose = vec![0u8; n_samples];
    let mut ref_count: usize = 0;
    let mut alt_count: usize = 0;
    let mut missing_alleles: usize = 0;

    let mut s = 0usize;
    for gt in it {
        let (a1, a2) = parse_gt(gt);
        match a1 {
            -1 => missing_alleles += 1,
            0 => ref_count += 1,
            1 => alt_count += 1,
            _ => panic!("multiallelic allele index in GT"),
        }
        match a2 {
            -1 => missing_alleles += 1,
            0 => ref_count += 1,
            1 => alt_count += 1,
            _ => panic!("multiallelic allele index in GT"),
        }
        alt_dose[s] = if a1 < 0 || a2 < 0 {
            3
        } else {
            (a1 == 1) as u8 + (a2 == 1) as u8
        };
        s += 1;
    }
    assert_eq!(s, n_samples, "wrong number of sample columns");

    // major = allele with larger count (ties -> REF, the lower index), matching
    // MarkerData.largestIndex; minor is the other.
    let minor_is_alt = ref_count >= alt_count;
    let minor_count = if minor_is_alt { alt_count } else { ref_count };
    if (minor_count as i32) < min_alleles {
        return; // fails minor-allele-count filter
    }

    let non_missing = (denom_alleles - missing_alleles) as f32;
    let minor_freq = minor_count as f32 / non_missing;
    let major_count = (denom_alleles - missing_alleles) - minor_count;
    let major_freq = major_count as f32 / non_missing;

    if m.chrom.is_empty() {
        m.chrom = String::from_utf8_lossy(chrom).into_owned();
    } else if m.chrom.as_bytes() != chrom {
        // single-chromosome input expected; stop at a new chromosome
        return;
    }
    m.pos.push(parse_i32(pos));
    m.id.push(String::from_utf8_lossy(id).into_owned());
    m.ref_a.push(String::from_utf8_lossy(ref_a).into_owned());
    m.alt_a.push(String::from_utf8_lossy(alt_a).into_owned());
    m.minor_is_alt.push(minor_is_alt);
    m.minor_freq.push(minor_freq);
    m.major_freq.push(major_freq);
    m.alt_dose.push(alt_dose);
}

/// Reads the VCF using frozen per-marker stats from a `.scorefreq` file: the
/// scored (minor) allele, its frequency, and the LD-pruned flag are taken from
/// `records` rather than computed from these samples. Markers are assembled in
/// scorefreq-file order (mirroring `VcfMarkerData.readScoreData`); the
/// minor-allele filter and LD pruning are NOT applied. Returns the markers and a
/// per-marker `correlated` (LD-pruned) flag vector.
pub fn read_vcf_scorefreq(path: &str, records: &[ScoreRecord]) -> (Markers, Vec<bool>) {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
    let mut reader = BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file));
    let mut line: Vec<u8> = Vec::with_capacity(1 << 16);
    let mut sample_ids: Vec<String> = Vec::new();

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).expect("read error");
        if n == 0 {
            panic!("VCF has no #CHROM header");
        }
        if line.starts_with(b"##") {
            continue;
        }
        if line.starts_with(b"#CHROM") {
            let trimmed = trim_newline(&line);
            for c in trimmed.split(|&b| b == b'\t').skip(9) {
                sample_ids.push(String::from_utf8_lossy(c).into_owned());
            }
            break;
        }
        panic!("unexpected line before #CHROM header");
    }
    let n_samples = sample_ids.len();
    assert!(n_samples > 0, "no samples in VCF");

    // key (chrom,pos,ref,alt) -> scorefreq record index
    let mut map: HashMap<(&str, i32, &str, &str), usize> = HashMap::with_capacity(records.len() * 2);
    for (i, r) in records.iter().enumerate() {
        map.insert((r.chrom.as_str(), r.pos, r.ref_a.as_str(), r.alt_a.as_str()), i);
    }

    // matched[i] = (alt_dose row, id, minor_is_alt) for record i
    let mut matched: Vec<Option<(Vec<u8>, String, bool)>> = (0..records.len()).map(|_| None).collect();

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).expect("read error");
        if n == 0 {
            break;
        }
        let rec = trim_newline(&line);
        if rec.is_empty() {
            continue;
        }
        let mut it = rec.split(|&b| b == b'\t');
        let chrom = it.next().expect("missing CHROM");
        let pos_b = it.next().expect("missing POS");
        let id = it.next().expect("missing ID");
        let ref_b = it.next().expect("missing REF");
        let alt_b = it.next().expect("missing ALT");
        if alt_b.contains(&b',') {
            panic!(
                "multiallelic marker not supported (ALT={}) at pos {}",
                String::from_utf8_lossy(alt_b),
                String::from_utf8_lossy(pos_b)
            );
        }
        let pos = parse_i32(pos_b);
        let chrom_s = std::str::from_utf8(chrom).expect("CHROM not UTF-8");
        let ref_s = std::str::from_utf8(ref_b).expect("REF not UTF-8");
        let alt_s = std::str::from_utf8(alt_b).expect("ALT not UTF-8");
        let idx = match map.get(&(chrom_s, pos, ref_s, alt_s)) {
            Some(&i) => i,
            None => continue, // marker not in scorefreq -> dropped
        };
        if matched[idx].is_some() {
            panic!("duplicate VCF marker matches scorefreq marker at {}:{}", chrom_s, pos);
        }
        // scored allele must equal REF or ALT
        let allele = records[idx].allele.as_str();
        let minor_is_alt = if allele == alt_s {
            true
        } else if allele == ref_s {
            false
        } else {
            panic!(
                "scorefreq allele '{}' not found in VCF marker at {}:{} (REF={}, ALT={})",
                allele, chrom_s, pos, ref_s, alt_s
            );
        };
        // genotype columns -> ALT dose
        let _qual = it.next();
        let _filter = it.next();
        let _info = it.next();
        let _format = it.next();
        let mut alt_dose = vec![0u8; n_samples];
        let mut s = 0usize;
        for gt in it {
            let (a1, a2) = parse_gt(gt);
            if a1 > 1 || a2 > 1 {
                panic!("multiallelic allele index in GT at {}:{}", chrom_s, pos);
            }
            alt_dose[s] = if a1 < 0 || a2 < 0 {
                3
            } else {
                (a1 == 1) as u8 + (a2 == 1) as u8
            };
            s += 1;
        }
        assert_eq!(s, n_samples, "wrong number of sample columns");
        matched[idx] = Some((alt_dose, String::from_utf8_lossy(id).into_owned(), minor_is_alt));
    }

    // assemble markers in scorefreq order
    let mut m = Markers {
        n_samples,
        sample_ids,
        chrom: records[0].chrom.clone(),
        pos: Vec::with_capacity(records.len()),
        id: Vec::with_capacity(records.len()),
        ref_a: Vec::with_capacity(records.len()),
        alt_a: Vec::with_capacity(records.len()),
        minor_is_alt: Vec::with_capacity(records.len()),
        minor_freq: Vec::with_capacity(records.len()),
        major_freq: Vec::with_capacity(records.len()),
        alt_dose: Vec::with_capacity(records.len()),
    };
    let mut correlated = Vec::with_capacity(records.len());
    for (i, r) in records.iter().enumerate() {
        let (alt_dose, id, minor_is_alt) = matched[i].take().unwrap_or_else(|| {
            panic!(
                "scorefreq marker not found in VCF: {}:{} REF={} ALT={}",
                r.chrom, r.pos, r.ref_a, r.alt_a
            )
        });
        m.pos.push(r.pos);
        m.id.push(id);
        m.ref_a.push(r.ref_a.clone());
        m.alt_a.push(r.alt_a.clone());
        m.minor_is_alt.push(minor_is_alt);
        m.minor_freq.push(r.freq);
        m.major_freq.push(1.0 - r.freq); // unused (LD pruning skipped in scorefreq mode)
        m.alt_dose.push(alt_dose);
        correlated.push(r.correlated);
    }
    (m, correlated)
}

#[inline]
fn parse_gt(gt: &[u8]) -> (i8, i8) {
    // Expect diploid "a|b" or "a/b"; alleles single-char 0/1 or '.'.
    if gt.len() < 3 {
        // haploid or malformed; treat single allele, second missing
        let a1 = allele_byte(gt.first().copied().unwrap_or(b'.'));
        return (a1, -1);
    }
    let a1 = allele_byte(gt[0]);
    let a2 = allele_byte(gt[2]);
    (a1, a2)
}

#[inline]
fn allele_byte(b: u8) -> i8 {
    match b {
        b'.' => -1,
        b'0' => 0,
        b'1' => 1,
        d if d.is_ascii_digit() => (d - b'0') as i8, // >1 -> triggers multiallelic panic upstream
        _ => -1,
    }
}

#[inline]
fn parse_i32(b: &[u8]) -> i32 {
    let mut v: i32 = 0;
    for &c in b {
        v = v * 10 + (c - b'0') as i32;
    }
    v
}

#[inline]
fn trim_newline(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}
