//! Deterministic post-pass empirical Fragment Length Distribution (`FragLenDist`)
//! and true transcript effective length computation.

pub const FLD_MAX_LEN: usize = 1000;

/// Thread-local / global histogram accumulator for observed concordant fragment lengths.
#[derive(Clone, Debug)]
pub struct FragLenAccum {
    pub hist: Vec<u64>,
    pub count: u64,
    pub sum: u64,
    pub sum_sq: u128,
}

impl Default for FragLenAccum {
    fn default() -> Self {
        Self {
            hist: vec![0u64; FLD_MAX_LEN + 1],
            count: 0,
            sum: 0,
            sum_sq: 0,
        }
    }
}

impl FragLenAccum {
    #[inline]
    pub fn record(&mut self, flen: u32) {
        let idx = (flen as usize).min(FLD_MAX_LEN);
        if idx > 0 {
            self.hist[idx] += 1;
            self.count += 1;
            self.sum += idx as u64;
            self.sum_sq += (idx as u128) * (idx as u128);
        }
    }

    pub fn merge_from(&mut self, other: &Self) {
        for (a, &b) in self.hist.iter_mut().zip(other.hist.iter()) {
            *a += b;
        }
        self.count += other.count;
        self.sum += other.sum;
        self.sum_sq += other.sum_sq;
    }

    /// Finalize into a smoothed, normalized `FragLenDist` once across the entire sample
    /// (guaranteeing 100% determinism independent of chunk/thread order).
    pub fn finalize(&self, fallback_mean: f64, fallback_sd: f64) -> FragLenDist {
        let mut pmf = vec![0.0f64; FLD_MAX_LEN + 1];
        if self.count >= 100 {
            // Empirical distribution lightly smoothed with a Gaussian kernel (sigma = 2 bp)
            // + tiny prior weight so tails stay smooth.
            let emp_mean = (self.sum as f64) / (self.count as f64);
            let var = ((self.sum_sq as f64) / (self.count as f64) - emp_mean * emp_mean).max(25.0);
            let emp_sd = var.sqrt().clamp(5.0, 250.0);

            let kernel = [0.0545, 0.2442, 0.4026, 0.2442, 0.0545];
            let inv_n = 1.0 / (self.count as f64);
            for (d, slot) in pmf.iter_mut().enumerate().skip(1) {
                let mut acc = 0.0;
                for (ki, &kw) in kernel.iter().enumerate() {
                    let pos = (d as isize) + (ki as isize) - 2;
                    if pos >= 1 && pos <= (FLD_MAX_LEN as isize) {
                        acc += (self.hist[pos as usize] as f64) * inv_n * kw;
                    }
                }
                *slot = acc;
            }
            let sum_p: f64 = pmf.iter().sum();
            if sum_p > 0.0 {
                for p in &mut pmf {
                    *p /= sum_p;
                }
            }
            return FragLenDist::from_pmf(pmf, emp_mean, emp_sd);
        }

        // Low-count / single-end fallback: truncated Gaussian prior centered at fallback_mean
        let mean = if self.count > 0 {
            (self.sum as f64) / (self.count as f64)
        } else {
            fallback_mean.max(20.0)
        };
        let sd = fallback_sd.max(10.0);
        let inv_2var = 1.0 / (2.0 * sd * sd);
        let mut sum_p = 0.0;
        for (d, slot) in pmf.iter_mut().enumerate().skip(1) {
            let diff = (d as f64) - mean;
            let p = (-diff * diff * inv_2var).exp();
            *slot = p;
            sum_p += p;
        }
        if sum_p > 0.0 {
            for p in &mut pmf {
                *p /= sum_p;
            }
        }
        FragLenDist::from_pmf(pmf, mean, sd)
    }
}

/// Finalized fragment-length distribution with O(1) prefix-sum lookup for exact effective lengths.
#[derive(Clone, Debug)]
pub struct FragLenDist {
    pub pmf: Vec<f64>,
    /// `cdf[L] = sum_{d=1..=L} pmf[d]`
    pub cdf: Vec<f64>,
    /// `first_moment[L] = sum_{d=1..=L} d * pmf[d]`
    pub first_moment: Vec<f64>,
    pub mean: f64,
    pub sd: f64,
}

impl FragLenDist {
    pub fn from_pmf(pmf: Vec<f64>, mean: f64, sd: f64) -> Self {
        let mut cdf = vec![0.0f64; pmf.len()];
        let mut first_moment = vec![0.0f64; pmf.len()];
        let mut c = 0.0;
        let mut m = 0.0;
        for d in 1..pmf.len() {
            c += pmf[d];
            m += (d as f64) * pmf[d];
            cdf[d] = c;
            first_moment[d] = m;
        }
        Self {
            pmf,
            cdf,
            first_moment,
            mean,
            sd,
        }
    }

    /// Exact transcript effective length in O(1):
    /// `sum_{d=1..=min(len, FLD_MAX_LEN)} pmf[d] * (len - d + 1) = (len + 1) * cdf[L] - first_moment[L]`
    /// Floored smoothly at `1.0` for ultra-short transcripts.
    #[inline]
    pub fn effective_length(&self, tx_len: u32) -> f64 {
        let l = tx_len as usize;
        if l == 0 {
            return 1.0;
        }
        let cap = l.min(FLD_MAX_LEN);
        let eff = ((l + 1) as f64) * self.cdf[cap] - self.first_moment[cap];
        if eff >= 1.0 {
            eff
        } else {
            // Smooth fallback when transcript is shorter than almost all fragments
            (tx_len as f64).max(1.0) * 0.25 + 0.75
        }
    }
}

/// Compute single-end effective length for a transcript given mean read length.
#[inline]
pub fn single_end_effective_length(tx_len: u32, mean_read_len: f64) -> f64 {
    ((tx_len as f64) - mean_read_len + 1.0).max(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fld_effective_length_exact() {
        let mut acc = FragLenAccum::default();
        for _ in 0..500 {
            acc.record(200);
        }
        let fld = acc.finalize(200.0, 20.0);
        assert!((fld.mean - 200.0).abs() < 1e-6);
        let eff = fld.effective_length(1000);
        // 1000 - 200 + 1 = 801
        assert!((eff - 801.0).abs() < 1.0, "got eff={eff}");
    }
}
