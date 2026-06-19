//! Port of `ibdseq.IbdScorer`: per-locus IBD/HBD LOD scores.
//! Operations mirror the Java arithmetic (Math.pow vs explicit multiply) to
//! maximize cross-language concordance.

pub struct IbdScorer {
    error_max: f64,
    error_prop: f64,
    max_error_array: [f64; 5],
}

impl IbdScorer {
    pub fn new(error_max: f64, error_prop: f64) -> Self {
        let max_error_array = error_array(error_max);
        IbdScorer { error_max, error_prop, max_error_array }
    }

    fn error_rate(&self, f_b: f64) -> f64 {
        let err = self.error_prop * f_b;
        if err <= self.error_max { err } else { self.error_max }
    }

    fn est_true_maf(f_b: f64, error_rate: f64) -> f64 {
        (f_b - error_rate) / (1.0 - 2.0 * error_rate)
    }

    /// log10 IBD likelihood ratio for a pair of allele doses (each 0,1,2).
    pub fn ibd_score(&self, dose1: i32, dose2: i32, f_b: f64) -> f64 {
        let e = self.error_rate(f_b);
        let p_b = Self::est_true_maf(f_b, e);
        let r = self.ibd_like(dose1, dose2, e, p_b) / Self::null_like(dose1, dose2, f_b);
        r.log10()
    }

    fn null_like(dose1: i32, dose2: i32, f_b: f64) -> f64 {
        let f_a = 1.0 - f_b;
        match dose1 + dose2 {
            0 => f_a.powf(4.0),
            1 => 4.0 * f_a.powf(3.0) * f_b,
            2 => {
                if dose1 == dose2 {
                    (2.0 * f_a * f_b).powf(2.0)
                } else {
                    2.0 * (f_a * f_b).powf(2.0)
                }
            }
            3 => 4.0 * f_a * f_b.powf(3.0),
            4 => f_b.powf(4.0),
            _ => f64::NAN,
        }
    }

    fn ibd_like(&self, dose1: i32, dose2: i32, err: f64, p_b: f64) -> f64 {
        let p_a = 1.0 - p_b;
        let e = if err == self.error_max {
            self.max_error_array
        } else {
            error_array(err)
        };
        match dose1 + dose2 {
            0 => e[0] * p_a.powf(3.0) + 2.0 * e[1] * p_a * p_a * p_b + e[2] * p_a * p_b,
            1 => {
                2.0 * (e[0] * p_a * p_a * p_b
                    + e[1] * (p_a * p_b + 2.0 * p_a.powf(3.0))
                    + 3.0 * e[2] * p_a * p_b)
            }
            2 => {
                if dose1 == dose2 {
                    (e[0] + 4.0 * e[1] + 2.0 * e[2]) * p_a * p_b
                        + 4.0 * e[2] * (p_a.powf(3.0) + p_b.powf(3.0))
                } else {
                    2.0 * ((e[1] + e[2] + e[3]) * p_a * p_b
                        + e[2] * (p_a.powf(3.0) + p_b.powf(3.0)))
                }
            }
            3 => {
                2.0 * (e[0] * p_a * p_b * p_b
                    + 2.0 * e[1] * p_b.powf(3.0)
                    + (e[1] + 3.0 * e[2] + e[3]) * p_a * p_b
                    + 2.0 * e[3] * p_a.powf(3.0))
            }
            4 => {
                e[0] * p_b.powf(3.0)
                    + 2.0 * e[1] * p_a * p_b * p_b
                    + e[2] * p_a * p_b
                    + 2.0 * e[3] * p_a * p_a * p_b
                    + e[4] * p_a.powf(3.0)
            }
            _ => f64::NAN,
        }
    }

    /// log10 HBD likelihood ratio for an allele dose (0,1,2).
    pub fn hbd_score(&self, dose: i32, f_b: f64) -> f64 {
        let e = self.error_rate(f_b);
        let f_a = 1.0 - f_b;
        let p_b = Self::est_true_maf(f_b, e);
        let p_a = 1.0 - p_b;
        match dose {
            0 => ((p_a + e * e * p_b) / (f_a * f_a)).log10(),
            1 => (e * (1.0 - e) / (f_a * f_b)).log10(),
            2 => ((e * e * p_a + p_b) / (f_b * f_b)).log10(),
            _ => 0.0,
        }
    }
}

fn error_array(err: f64) -> [f64; 5] {
    let mut a = [0.0f64; 5];
    a[0] = (1.0 - err).powf(4.0);
    for j in 1..5 {
        a[j] = a[j - 1] * err / (1.0 - err);
    }
    a
}
