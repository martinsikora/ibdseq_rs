//! Sparse detection kernel. Score-cell tables are keyed by ALT-allele dose (the
//! minor-allele transform is folded into the tables). Only the small set of
//! *retained* (LD-kept) markers scores for every genotype; the large majority of
//! *exclusion* (LD-pruned) markers score only at an opposite homozygote (IBD) or
//! a heterozygote (HBD). For each outer sample s1 we therefore do the full
//! partner sweep only on retained markers and handle exclusion markers as sparse
//! events, while exactly reproducing the dense kernel's segment boundaries.

use crate::scorer::IbdScorer;
use crate::vcf::Markers;
use rayon::prelude::*;

pub struct Segment {
    pub id1: u32,
    pub id2: u32,
    pub start: u32,
    pub end: u32,
    pub score: f32,
    pub hbd: bool,
}

pub struct ScoreTables {
    pub ibd_cell: Vec<f64>,       // n*16, index (j<<4)+(altA<<2)+altB (alt doses, missing=3)
    pub hbd_cell: Vec<f64>,       // n*4,  index (j<<2)+alt
    pub is_retained: Vec<bool>,   // n: true == LD-kept (scores densely)
    pub excl_ibd: Vec<f64>,       // n: opposite-homozygote IBD score (exclusion markers)
    pub excl_hbd: Vec<f64>,       // n: heterozygote HBD score (exclusion markers)
    pub hom_minor: Vec<Vec<u32>>, // n: sorted homozygous-minor (dose 2) sample idx (exclusion markers)
    pub hom_major: Vec<Vec<u32>>, // n: sorted homozygous-major (dose 0) sample idx (exclusion markers)
    pub minor_is_alt: Vec<bool>,  // n: ALT->minor dose transform (for s1's dose at exclusion markers)
}

#[inline]
fn dose_pair_index(da: usize, db: usize) -> usize {
    if (da == 0 && db == 2) || (da == 2 && db == 0) {
        5
    } else {
        da + db
    }
}

/// minor-allele dose from ALT dose (missing 3 stays 3).
#[inline]
fn minor_dose(alt: usize, minor_is_alt: bool) -> usize {
    if alt == 3 {
        3
    } else if minor_is_alt {
        alt
    } else {
        2 - alt
    }
}

/// Builds ALT-keyed score tables plus the sparse-exclusion side tables.
/// `correlated[j]` markers are scored exclusion-only (opposite homozygotes for
/// IBD, het for HBD).
pub fn build_tables(m: &Markers, correlated: &[bool], scorer: &IbdScorer) -> ScoreTables {
    let n = m.n_markers();
    let mut ibd_cell = vec![0.0f64; n * 16];
    let mut hbd_cell = vec![0.0f64; n * 4];
    let mut is_retained = vec![false; n];
    let mut excl_ibd = vec![0.0f64; n];
    let mut excl_hbd = vec![0.0f64; n];
    let mut hom_minor: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut hom_major: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut minor_is_alt = vec![false; n];

    for j in 0..n {
        let fb = m.minor_freq[j] as f64; // float widened to double, matching Java
        let idx = [
            scorer.ibd_score(0, 0, fb),
            scorer.ibd_score(0, 1, fb),
            scorer.ibd_score(1, 1, fb),
            scorer.ibd_score(1, 2, fb),
            scorer.ibd_score(2, 2, fb),
            scorer.ibd_score(0, 2, fb),
        ];
        let h = [
            scorer.hbd_score(0, fb),
            scorer.hbd_score(1, fb),
            scorer.hbd_score(2, fb),
        ];
        let is_cor = correlated[j];
        let minor_alt = m.minor_is_alt[j];
        is_retained[j] = !is_cor;
        minor_is_alt[j] = minor_alt;
        excl_ibd[j] = idx[5]; // opposite-homozygote (minor doses 0,2)
        excl_hbd[j] = h[1]; // heterozygote

        let base = j << 4;
        for alt_a in 0..4usize {
            for alt_b in 0..4usize {
                let da = minor_dose(alt_a, minor_alt);
                let db = minor_dose(alt_b, minor_alt);
                let mut val = 0.0;
                if da < 3 && db < 3 {
                    let not_ibs = (da == 0 && db == 2) || (da == 2 && db == 0);
                    if not_ibs || !is_cor {
                        val = idx[dose_pair_index(da, db)];
                    }
                }
                ibd_cell[base + (alt_a << 2) + alt_b] = val;
            }
        }
        let hbase = j << 2;
        for alt in 0..4usize {
            let d = minor_dose(alt, minor_alt);
            let mut val = 0.0;
            if d < 3 && (!is_cor || d == 1) {
                val = h[d];
            }
            hbd_cell[hbase + alt] = val;
        }

        // homozygous-minor (dose 2) / homozygous-major (dose 0) sample lists for
        // exclusion markers only. Both are built ascending by sample index.
        if is_cor {
            let row = &m.alt_dose[j];
            let mut mn: Vec<u32> = Vec::new();
            let mut mj: Vec<u32> = Vec::new();
            for (s, &ad) in row.iter().enumerate() {
                match minor_dose(ad as usize, minor_alt) {
                    2 => mn.push(s as u32),
                    0 => mj.push(s as u32),
                    _ => {}
                }
            }
            hom_minor[j] = mn;
            hom_major[j] = mj;
        }
    }
    ScoreTables {
        ibd_cell,
        hbd_cell,
        is_retained,
        excl_ibd,
        excl_hbd,
        hom_minor,
        hom_major,
        minor_is_alt,
    }
}

pub fn detect(
    alt_dose: &[Vec<u8>],
    tables: &ScoreTables,
    n_markers: usize,
    n_samples: usize,
    focus: Option<&[bool]>,
    ibd_lod: f32,
    ibd_trim: f32,
) -> Vec<Segment> {
    (0..n_samples)
        .into_par_iter()
        .map(|s1| {
            detect_one_s1(
                s1, alt_dose, tables, n_markers, n_samples, focus, ibd_lod, ibd_trim,
            )
        })
        .reduce(Vec::new, |mut a, mut b| {
            if a.len() >= b.len() {
                a.append(&mut b);
                a
            } else {
                b.append(&mut a);
                b
            }
        })
}

#[allow(clippy::too_many_arguments)]
fn detect_one_s1(
    s1: usize,
    alt_dose: &[Vec<u8>],
    t: &ScoreTables,
    n_markers: usize,
    n_samples: usize,
    focus: Option<&[bool]>,
    ibd_lod: f32,
    ibd_trim: f32,
) -> Vec<Segment> {
    let focus_s1 = focus.map_or(true, |f| f[s1]);
    let has_focus = focus.is_some();
    let fref: &[bool] = focus.unwrap_or(&[]);

    let n_part = n_samples - s1; // partner p in 0..n_part, sample s2 = s1 + p (p==0 => self/HBD)
    let mut this_sum = vec![0f32; n_part];
    let mut max_sum = vec![0f32; n_part];
    let mut start = vec![0u32; n_part];
    let mut end = vec![0u32; n_part];
    let mut out: Vec<Segment> = Vec::new();

    let do_self = focus_s1; // self pair (HBD) accepted iff s1 is a focus sample
    let ibd_cell = &t.ibd_cell;
    let hbd_cell = &t.hbd_cell;

    for j in 0..n_markers {
        let djrow = &alt_dose[j];
        let alt1 = djrow[s1] as usize;

        if t.is_retained[j] {
            // ---- retained marker: full partner sweep (dense), same as stock ----
            if do_self {
                let sc = hbd_cell[(j << 2) + alt1];
                update(0, sc, j as u32, ibd_lod, &mut this_sum, &mut max_sum, &mut start, &mut end,
                       s1, s1, true, ibd_trim, alt_dose, ibd_cell, hbd_cell, &mut out);
            }
            let base = (j << 4) + (alt1 << 2);
            let row = [
                ibd_cell[base],
                ibd_cell[base + 1],
                ibd_cell[base + 2],
                ibd_cell[base + 3],
            ];
            if has_focus {
                for p in 1..n_part {
                    let s2 = s1 + p;
                    if focus_s1 || fref[s2] {
                        let sc = row[djrow[s2] as usize];
                        update(p, sc, j as u32, ibd_lod, &mut this_sum, &mut max_sum, &mut start,
                               &mut end, s1, s2, false, ibd_trim, alt_dose, ibd_cell, hbd_cell, &mut out);
                    }
                }
            } else {
                for p in 1..n_part {
                    let s2 = s1 + p;
                    let sc = row[djrow[s2] as usize];
                    // inlined update (hot path, no focus); lazy-start folds in skipped zeros
                    if this_sum[p] == 0.0 {
                        start[p] = j as u32;
                    }
                    let mut ts = ((this_sum[p] as f64) + sc) as f32;
                    if ts > max_sum[p] {
                        max_sum[p] = ts;
                        end[p] = j as u32;
                    } else if (ts as f64) <= 0.0 {
                        if max_sum[p] >= ibd_lod {
                            emit(&mut out, s1, s2, start[p], end[p], max_sum[p], false, ibd_trim,
                                 alt_dose, ibd_cell, hbd_cell);
                        }
                        start[p] = j as u32 + 1;
                        end[p] = start[p];
                        ts = 0.0;
                        max_sum[p] = 0.0;
                    }
                    this_sum[p] = ts;
                }
            }
        } else {
            // ---- exclusion marker: only opposite-homozygote / het events ----
            let minor_alt = t.minor_is_alt[j];
            let md1 = minor_dose(alt1, minor_alt);

            // self / HBD: nonzero only when s1 is heterozygous.
            if do_self && md1 == 1 {
                let sc = t.excl_hbd[j];
                update(0, sc, j as u32, ibd_lod, &mut this_sum, &mut max_sum, &mut start, &mut end,
                       s1, s1, true, ibd_trim, alt_dose, ibd_cell, hbd_cell, &mut out);
            }

            // IBD: opposite homozygote => one homMajor (dose 0), one homMinor (dose 2).
            // For a homozygous s1, the events are exactly the *opposite*-homozygote
            // list, partners with idx > s1. Partner updates at one marker are
            // independent, so list order is irrelevant to the result.
            let opp: Option<&Vec<u32>> = match md1 {
                0 => Some(&t.hom_minor[j]), // s1 homMajor -> homMinor partners
                2 => Some(&t.hom_major[j]), // s1 homMinor -> homMajor partners
                _ => None,                  // het / missing -> no IBD events
            };
            if let Some(list) = opp {
                let sc = t.excl_ibd[j];
                let lo = list.partition_point(|&x| (x as usize) <= s1);
                for &s2u in &list[lo..] {
                    let s2 = s2u as usize;
                    if !has_focus || focus_s1 || fref[s2] {
                        let p = s2 - s1;
                        update(p, sc, j as u32, ibd_lod, &mut this_sum, &mut max_sum, &mut start,
                               &mut end, s1, s2, false, ibd_trim, alt_dose, ibd_cell, hbd_cell, &mut out);
                    }
                }
            }
        }
    }

    // finalize trailing segments
    let finalize = |p: usize, s2: usize, hbd: bool, out: &mut Vec<Segment>| {
        if max_sum[p] >= ibd_lod {
            emit(out, s1, s2, start[p], end[p], max_sum[p], hbd, ibd_trim, alt_dose, ibd_cell, hbd_cell);
        }
    };
    if do_self {
        finalize(0, s1, true, &mut out);
    }
    for p in 1..n_part {
        let s2 = s1 + p;
        if !has_focus || focus_s1 || fref[s2] {
            finalize(p, s2, false, &mut out);
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn update(
    p: usize,
    sc: f64,
    j: u32,
    ibd_lod: f32,
    this_sum: &mut [f32],
    max_sum: &mut [f32],
    start: &mut [u32],
    end: &mut [u32],
    s1: usize,
    s2: usize,
    hbd: bool,
    ibd_trim: f32,
    alt_dose: &[Vec<u8>],
    ibd_cell: &[f64],
    hbd_cell: &[f64],
    out: &mut Vec<Segment>,
) {
    // lazy-start: an empty segment's start equals the first nonzero-contributing
    // marker, reproducing the dense kernel's per-zero-marker start advancement.
    if this_sum[p] == 0.0 {
        start[p] = j;
    }
    let mut ts = ((this_sum[p] as f64) + sc) as f32;
    if ts > max_sum[p] {
        max_sum[p] = ts;
        end[p] = j;
    } else if (ts as f64) <= 0.0 {
        if max_sum[p] >= ibd_lod {
            emit(out, s1, s2, start[p], end[p], max_sum[p], hbd, ibd_trim, alt_dose, ibd_cell, hbd_cell);
        }
        start[p] = j + 1;
        end[p] = start[p];
        ts = 0.0;
        max_sum[p] = 0.0;
    }
    this_sum[p] = ts;
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn emit(
    out: &mut Vec<Segment>,
    s1: usize,
    s2: usize,
    start: u32,
    end: u32,
    score: f32,
    hbd: bool,
    ibd_trim: f32,
    alt_dose: &[Vec<u8>],
    ibd_cell: &[f64],
    hbd_cell: &[f64],
) {
    let (st, en) = if ibd_trim > 0.0 {
        let st = trim_start(s1, s2, hbd, start, end, ibd_trim, alt_dose, ibd_cell, hbd_cell);
        let en = trim_end(s1, s2, hbd, st, end, ibd_trim, alt_dose, ibd_cell, hbd_cell);
        (st, en)
    } else {
        (start, end)
    };
    if en > st {
        out.push(Segment {
            id1: s1 as u32,
            id2: s2 as u32,
            start: st,
            end: en,
            score,
            hbd,
        });
    }
}

#[inline]
fn cell_at(
    s1: usize,
    s2: usize,
    hbd: bool,
    k: usize,
    alt_dose: &[Vec<u8>],
    ibd_cell: &[f64],
    hbd_cell: &[f64],
) -> f64 {
    let a1 = alt_dose[k][s1] as usize;
    if hbd {
        hbd_cell[(k << 2) + a1]
    } else {
        let a2 = alt_dose[k][s2] as usize;
        ibd_cell[(k << 4) + (a1 << 2) + a2]
    }
}

#[allow(clippy::too_many_arguments)]
fn trim_start(
    s1: usize,
    s2: usize,
    hbd: bool,
    start: u32,
    end: u32,
    ibd_trim: f32,
    alt_dose: &[Vec<u8>],
    ibd_cell: &[f64],
    hbd_cell: &[f64],
) -> u32 {
    let mut sum: f32 = 0.0;
    let mut index = start as i64;
    let end = end as i64;
    while index <= end && sum < ibd_trim {
        sum = ((sum as f64) + cell_at(s1, s2, hbd, index as usize, alt_dose, ibd_cell, hbd_cell)) as f32;
        index += 1;
    }
    (index - 1) as u32
}

#[allow(clippy::too_many_arguments)]
fn trim_end(
    s1: usize,
    s2: usize,
    hbd: bool,
    start: u32,
    end: u32,
    ibd_trim: f32,
    alt_dose: &[Vec<u8>],
    ibd_cell: &[f64],
    hbd_cell: &[f64],
) -> u32 {
    let mut sum: f32 = 0.0;
    let mut index = end as i64;
    let start = start as i64;
    while index >= start && sum < ibd_trim {
        sum = ((sum as f64) + cell_at(s1, s2, hbd, index as usize, alt_dose, ibd_cell, hbd_cell)) as f32;
        index -= 1;
    }
    (index + 1) as u32
}
