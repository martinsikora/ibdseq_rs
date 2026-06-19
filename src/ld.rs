//! Port of the LD-thinning in `VcfMarkerData.checkLastMarkerR2` + `MarkerData.r2`.
//! Returns a per-marker `correlated` flag (true == LD-pruned). Markers are kept
//! (exclusion-only scoring), matching the stock-equivalent full-run behavior.

use crate::vcf::Markers;

/// Squared dosage correlation between two markers using their major alleles.
fn r2(m: &Markers, a: usize, b: usize) -> f64 {
    let n_samples = m.n_samples;
    let mut cnt: i64 = 0;
    let mut cnt_a: i64 = 0;
    let mut cnt_b: i64 = 0;
    let mut cnt_ab: i64 = 0;
    let mut cnt_aa: i64 = 0;
    let mut cnt_bb: i64 = 0;
    let da = &m.alt_dose[a];
    let db = &m.alt_dose[b];
    let a_minor_alt = m.minor_is_alt[a];
    let b_minor_alt = m.minor_is_alt[b];
    for s in 0..n_samples {
        let dose_a = dose_from(da[s], a_minor_alt);
        let dose_b = dose_from(db[s], b_minor_alt);
        if dose_a >= 0 && dose_b >= 0 {
            cnt += 1;
            cnt_a += dose_a as i64;
            cnt_b += dose_b as i64;
            cnt_ab += (dose_a * dose_b) as i64;
            cnt_aa += (dose_a * dose_a) as i64;
            cnt_bb += (dose_b * dose_b) as i64;
        }
    }
    let n = cnt as f64;
    let mean_a = cnt_a as f64 / n;
    let mean_b = cnt_b as f64 / n;
    let cov = cnt_ab as f64 / n - mean_a * mean_b;
    let num = cov * cov;
    let den = (cnt_aa as f64 / n - mean_a * mean_a) * (cnt_bb as f64 / n - mean_b * mean_b);
    num / den
}

#[inline]
fn dose_from(ad: u8, minor_is_alt: bool) -> i32 {
    if ad == 3 {
        -1
    } else if minor_is_alt {
        (2 - ad) as i32 // major == REF
    } else {
        ad as i32 // major == ALT
    }
}

pub fn ld_prune(m: &Markers, r2_window: i32, r2_max: f32) -> Vec<bool> {
    let n = m.n_markers();
    let mut excluded = vec![false; n];
    if r2_window <= 0 {
        return excluded;
    }
    let radius = 1 + (r2_window / 2) as usize;
    let r2_max = r2_max as f64;
    for last in 0..n {
        // startIndex = max(0, (last+1) - radius)
        let start_index = (last + 1).saturating_sub(radius);
        let mut finished = false;
        let mut k = last;
        while k > start_index && !finished {
            k -= 1;
            if !excluded[k] {
                if r2(m, k, last) > r2_max {
                    // alleleFrequency(majorB) > alleleFrequency(majorA) ?
                    if m.major_freq[last] > m.major_freq[k] {
                        excluded[k] = true;
                    } else {
                        excluded[last] = true;
                        finished = true;
                    }
                }
            }
        }
    }
    excluded
}
