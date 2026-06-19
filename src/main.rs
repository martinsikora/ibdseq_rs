//! Optimized Rust port of the modified ibdseq r1206 IBD/HBD detection pipeline
//! (full-run path: read VCF -> minor-allele filter -> LD thinning -> LOD scoring
//! -> parallel segment detection -> cM tables). Output matches the Java
//! `<out>.ibd.gz` / `<out>.hbd.gz` format.

mod cm;
mod detect;
mod ld;
mod scorer;
mod vcf;

use cm::CmConverter;
use detect::Segment;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use scorer::IbdScorer;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

struct Config {
    gt: String,
    out: String,
    min_alleles: i32,
    ibd_lod: f32,
    ibd_trim: f32,
    error_max: f64,
    error_prop: f64,
    r2_window: i32,
    r2_max: f32,
    nthreads: usize,
    focussamples: Option<String>,
    map: Option<String>,
    cmpermb: f64,
    bins: String,
    noout: bool,
}

fn parse_args() -> Config {
    let mut c = Config {
        gt: String::new(),
        out: String::new(),
        min_alleles: 2,
        ibd_lod: 3.0,
        ibd_trim: 0.0,
        error_max: 0.001,
        error_prop: 0.25,
        r2_window: 500,
        r2_max: 0.15,
        nthreads: 0,
        focussamples: None,
        map: None,
        cmpermb: 1.0,
        bins: "0,1,2,4,8,20,30,3000".to_string(),
        noout: false,
    };
    for arg in std::env::args().skip(1) {
        let (k, v) = arg.split_once('=').unwrap_or_else(|| {
            eprintln!("invalid argument (expected key=value): {}", arg);
            std::process::exit(2);
        });
        match k {
            "gt" => c.gt = v.to_string(),
            "out" => c.out = v.to_string(),
            "minalleles" => c.min_alleles = v.parse().unwrap(),
            "ibdlod" => c.ibd_lod = v.parse().unwrap(),
            "ibdtrim" => c.ibd_trim = v.parse().unwrap(),
            "errormax" => c.error_max = v.parse().unwrap(),
            "errorprop" => c.error_prop = v.parse().unwrap(),
            "r2window" => c.r2_window = v.parse().unwrap(),
            "r2max" => c.r2_max = v.parse().unwrap(),
            "nthreads" => c.nthreads = v.parse().unwrap(),
            "focussamples" => c.focussamples = Some(v.to_string()),
            "map" => c.map = Some(v.to_string()),
            "cmpermb" => c.cmpermb = v.parse().unwrap(),
            "bins" => c.bins = v.to_string(),
            "noout" => c.noout = v.parse().unwrap(),
            _ => {
                eprintln!("unknown argument: {}", k);
                std::process::exit(2);
            }
        }
    }
    if c.gt.is_empty() || c.out.is_empty() {
        eprintln!("usage: ibdseq_rs gt=<vcf> out=<prefix> [minalleles= ibdlod= ibdtrim= errormax= errorprop= r2window= r2max= nthreads= focussamples= map= cmpermb= bins= noout=]");
        std::process::exit(2);
    }
    c
}

fn main() {
    let c = parse_args();
    let nthreads = if c.nthreads == 0 {
        rayon::current_num_threads()
    } else {
        c.nthreads
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build_global()
        .unwrap();

    let t_start = Instant::now();

    // 1. Read VCF + per-marker stats.
    let t0 = Instant::now();
    let markers = vcf::read_vcf(&c.gt, c.min_alleles);
    let n_markers = markers.n_markers();
    let n_samples = markers.n_samples;
    let read_s = t0.elapsed().as_secs_f64();

    // 2. LD thinning -> correlated (LD-pruned) flags; markers kept (exclusion-only).
    let t0 = Instant::now();
    let correlated = ld::ld_prune(&markers, c.r2_window, c.r2_max);
    let n_correlated = correlated.iter().filter(|&&b| b).count();
    let prune_s = t0.elapsed().as_secs_f64();

    // 3. ALT-keyed score-cell tables (minor transform folded in).
    let t0 = Instant::now();
    let scorer = IbdScorer::new(c.error_max, c.error_prop);
    let tables = detect::build_tables(&markers, &correlated, &scorer);
    let prep_s = t0.elapsed().as_secs_f64();

    // focus mask
    let focus: Option<Vec<bool>> = c.focussamples.as_ref().map(|f| read_focus(f, &markers.sample_ids));

    // 4. Detection (the kernel metric).
    let t0 = Instant::now();
    let segments = detect::detect(
        &markers.alt_dose,
        &tables,
        n_markers,
        n_samples,
        focus.as_deref(),
        c.ibd_lod,
        c.ibd_trim,
    );
    let detect_s = t0.elapsed().as_secs_f64();

    // 5. Output.
    let t0 = Instant::now();
    let (n_ibd, n_hbd) = if c.noout {
        let mut ni = 0u64;
        let mut nh = 0u64;
        for s in &segments {
            if s.hbd {
                nh += 1;
            } else {
                ni += 1;
            }
        }
        (ni, nh)
    } else {
        write_segments(&c, &markers, &segments)
    };
    let out_s = t0.elapsed().as_secs_f64();

    let total_s = t_start.elapsed().as_secs_f64();
    let retained = n_markers - n_correlated;
    eprintln!("ibdseq_rs");
    eprintln!("  gt              : {}", c.gt);
    eprintln!("  out             : {}", c.out);
    eprintln!("  nthreads        : {}", nthreads);
    eprintln!("  samples         : {}", n_samples);
    eprintln!("  markers(>=maf)  : {}", n_markers);
    eprintln!("  LD-thinned      : {} ({} correlated kept exclusion-only)", retained, n_correlated);
    eprintln!("  IBD segments    : {}", n_ibd);
    eprintln!("  HBD segments    : {}", n_hbd);
    eprintln!("  time read       : {:.2} s", read_s);
    eprintln!("  time ld-prune   : {:.2} s", prune_s);
    eprintln!("  time prep       : {:.2} s", prep_s);
    eprintln!("  time DETECTION  : {:.2} s", detect_s);
    eprintln!("  time output     : {:.2} s", out_s);
    eprintln!("  time TOTAL      : {:.2} s", total_s);
}

fn read_focus(path: &str, sample_ids: &[String]) -> Vec<bool> {
    let bytes = std::fs::read(path).expect("cannot read focussamples");
    let text = String::from_utf8_lossy(&bytes);
    let wanted: HashSet<&str> = text.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).collect();
    let mask: Vec<bool> = sample_ids.iter().map(|id| wanted.contains(id.as_str())).collect();
    let found = mask.iter().filter(|&&b| b).count();
    if found != wanted.len() {
        eprintln!(
            "warning: {} focus IDs requested, {} matched in VCF",
            wanted.len(),
            found
        );
    }
    mask
}

const HEADER: &[u8] =
    b"sample1\tsample2\tchromosome\tpos_start\tpos_end\tlod\tpos_start_cm\tpos_end_cm\tl_cm\tl_cm_bin\n";

/// Append one segment line to `buf`, byte-for-byte identical to the serial
/// `format!` path (itoa for the integer columns, `write!` to preserve the
/// fixed-precision float rounding).
#[inline]
fn fmt_line(
    buf: &mut Vec<u8>,
    markers: &vcf::Markers,
    marker_cm: &[f64],
    conv: &CmConverter,
    seg: &Segment,
) {
    use std::io::Write as _;
    let mut ib = itoa::Buffer::new();
    let pos_start = markers.pos[seg.start as usize];
    let pos_end = markers.pos[seg.end as usize];
    let cm_start = marker_cm[seg.start as usize];
    let cm_end = marker_cm[seg.end as usize];
    let l_cm = cm_end - cm_start;
    buf.extend_from_slice(markers.sample_ids[seg.id1 as usize].as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(markers.sample_ids[seg.id2 as usize].as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(markers.chrom.as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(ib.format(pos_start).as_bytes());
    buf.push(b'\t');
    buf.extend_from_slice(ib.format(pos_end).as_bytes());
    write!(buf, "\t{:.2}\t{:.3}\t{:.3}\t{:.3}\t", seg.score, cm_start, cm_end, l_cm).unwrap();
    buf.extend_from_slice(conv.bin_str(l_cm).as_bytes());
    buf.push(b'\n');
}

#[inline]
fn gzip_block(bytes: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::with_capacity(bytes.len() / 3 + 64), Compression::default());
    enc.write_all(bytes).unwrap();
    enc.finish().unwrap()
}

fn write_segments(c: &Config, markers: &vcf::Markers, segments: &[Segment]) -> (u64, u64) {
    let conv = CmConverter::new(c.map.as_deref(), &c.bins, c.cmpermb);
    // Precompute cM per marker once (kills the per-line binary search / cm()).
    let marker_cm: Vec<f64> = markers.pos.iter().map(|&p| conv.cm(p)).collect();

    // Format + compress shards in parallel; gzip multi-member output concatenates.
    const SHARD: usize = 262_144;
    let blocks: Vec<(Vec<u8>, Vec<u8>, u64, u64)> = segments
        .par_chunks(SHARD)
        .map(|chunk| {
            let mut ibd_buf: Vec<u8> = Vec::with_capacity(SHARD * 48);
            let mut hbd_buf: Vec<u8> = Vec::new();
            let mut ni = 0u64;
            let mut nh = 0u64;
            for seg in chunk {
                if seg.hbd {
                    fmt_line(&mut hbd_buf, markers, &marker_cm, &conv, seg);
                    nh += 1;
                } else {
                    fmt_line(&mut ibd_buf, markers, &marker_cm, &conv, seg);
                    ni += 1;
                }
            }
            (gzip_block(&ibd_buf), gzip_block(&hbd_buf), ni, nh)
        })
        .collect();

    let mut ibd = BufWriter::with_capacity(1 << 20, File::create(format!("{}.ibd.gz", c.out)).unwrap());
    let mut hbd = BufWriter::with_capacity(1 << 20, File::create(format!("{}.hbd.gz", c.out)).unwrap());
    ibd.write_all(&gzip_block(HEADER)).unwrap();
    hbd.write_all(&gzip_block(HEADER)).unwrap();
    let mut n_ibd = 0u64;
    let mut n_hbd = 0u64;
    for (ib, hb, ni, nh) in &blocks {
        ibd.write_all(ib).unwrap();
        hbd.write_all(hb).unwrap();
        n_ibd += ni;
        n_hbd += nh;
    }
    ibd.flush().unwrap();
    hbd.flush().unwrap();
    (n_ibd, n_hbd)
}
