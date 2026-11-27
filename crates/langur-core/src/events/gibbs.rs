//! Targeted Collapsed Gibbs Sampler for complex (`K_n >= 2`) local splice-graph events.
//!
//! - Triggered ONLY when `K_n >= 2 && !ambig_groups.is_empty()` and `gibbs_samples > 0`,
//!   keeping binary `K_1` events exact closed-form `O(1)`.
//! - Deterministic per-event `Xoshiro256++` RNG seeded via `SplitMix64(seed, gene_id, node_1based)`
//!   so outputs are 100% invariant across thread counts.
//! - Group-balanced priors (`alpha_p = 1/n_inc` for inclusion paths, `1/n_exc` for exclusion paths,
//!   unless `--complex-vb-prior` is explicitly set) preserve a symmetric `Beta(1, 1)` marginal prior
//!   on `Psi` regardless of complexity `K_n`.
//! - Exact Gamma rate parameterization (`theta_p = Gamma(alpha_p + z_p, scale=1) / eff_len[p]`).
//! - Fractional ambiguous counts `c_g` and `O(K)` depth-independent sequential Binomial-Normal/Poisson
//!   multinomial draws for `c_g > 50` so high-depth events never stall.

use crate::events::reduce::ReducedPath;
use std::f64::consts::TAU;

/// Default burn-in sweeps before collecting thinned posterior draws.
pub const DEFAULT_GIBBS_BURN_IN: usize = 32;
/// Default thinning interval between collected draws.
pub const DEFAULT_GIBBS_THIN: usize = 2;

/// Summary returned by the local-event collapsed Gibbs sampler.
#[derive(Copy, Clone, Debug)]
pub struct GibbsPosteriorSummary {
    /// Empirical posterior mean of `Psi = sum_{p in inc} pi_p`.
    pub mean_psi: f64,
    /// Empirical posterior variance of `Psi` (`Inf_Var`, column 15 of `.psi.gz`).
    pub inf_var: f64,
    /// Empirical 5th percentile of `Psi` (`CI_Lo`).
    pub ci_lo: f64,
    /// Empirical 95th percentile of `Psi` (`CI_Hi`).
    pub ci_hi: f64,
}

/// Fast deterministic `Xoshiro256++` PRNG seeded via `SplitMix64`.
struct EventRng {
    s: [u64; 4],
    spare_normal: Option<f64>,
}

impl EventRng {
    #[inline]
    fn splitmix64(x: &mut u64) -> u64 {
        *x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *x;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn from_event(global_seed: u64, gene_id: u32, node_1based: u32) -> Self {
        let mut sm = global_seed
            ^ ((gene_id as u64).wrapping_mul(0xD1B5_4A32_D192_ED03))
            ^ ((node_1based as u64).wrapping_mul(0x94D0_49BB_1331_11EB).rotate_left(17));
        let mut s = [0u64; 4];
        for slot in &mut s {
            *slot = Self::splitmix64(&mut sm);
        }
        if s == [0, 0, 0, 0] {
            s[0] = 0x9E37_79B9_7F4A_7C15;
        }
        Self {
            s,
            spare_normal: None,
        }
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let res = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        res
    }

    /// Uniform float strictly in `(0.0, 1.0)`.
    #[inline]
    fn uniform_open(&mut self) -> f64 {
        // 53 random bits mapped into (0, 1)
        let u = (self.next_u64() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64));
        u.clamp(f64::MIN_POSITIVE, 1.0 - f64::EPSILON)
    }

    /// Standard Normal `N(0, 1)` via Box-Muller with cached spare draw.
    #[inline]
    fn normal(&mut self) -> f64 {
        if let Some(z) = self.spare_normal.take() {
            return z;
        }
        let u1 = self.uniform_open();
        let u2 = self.uniform_open();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = TAU * u2;
        self.spare_normal = Some(r * theta.sin());
        r * theta.cos()
    }

    /// Standard `Gamma(shape, scale = 1.0)` via Marsaglia-Tsang (with Ahrens-Dieter boost for `shape < 1`).
    fn gamma_unit_scale(&mut self, shape: f64) -> f64 {
        if !shape.is_finite() || shape <= 0.0 {
            return 0.0;
        }
        if shape < 1.0 {
            // Boost: Gamma(a) = Gamma(a + 1) * U^(1/a)
            let g1 = self.gamma_unit_scale(shape + 1.0);
            let u = self.uniform_open();
            let log_scale = u.ln() / shape;
            if log_scale < -700.0 {
                return 0.0;
            }
            return g1 * log_scale.exp();
        }
        let d = shape - (1.0 / 3.0);
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v_root = 1.0 + c * x;
            if v_root <= 0.0 {
                continue;
            }
            let v = v_root * v_root * v_root;
            let u = self.uniform_open();
            let x2 = x * x;
            if u < 1.0 - 0.0331 * x2 * x2 {
                return d * v;
            }
            if u.ln() < 0.5 * x2 + d * (1.0 - v + v.ln()) {
                return d * v;
            }
        }
    }

    /// Exact Knuth Poisson sampler for small rate `lambda` (used only when `lambda < 10.0`).
    fn poisson_small(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        let l = (-lambda).exp();
        let mut k = 0u64;
        let mut p = 1.0f64;
        loop {
            k += 1;
            p *= self.uniform_open();
            if p <= l {
                return k - 1;
            }
        }
    }

    /// Draw a single categorical index in `0..probs.len()` given normalized `probs` summing to 1.
    #[inline]
    fn categorical(&mut self, probs: &[f64]) -> usize {
        let u = self.uniform_open();
        let mut acc = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            acc += p;
            if u <= acc {
                return i;
            }
        }
        probs.len() - 1
    }

    /// Draw from `Binomial(n, p)` (allowing non-negative fractional `n`) in `O(1)` time for `n > 50`.
    ///
    /// - For `n <= 50`: exact Bernoulli sum over `floor(n)` plus fractional residual `frac` via Bernoulli(`p`).
    /// - For `n > 50`: Normal approximation when `min(n*p, n*(1-p)) >= 10`, otherwise exact Knuth Poisson
    ///   on the rare side (`n*p` or `n*(1-p)`), clamped strictly to `[0.0, n]`.
    fn binomial_continuous(&mut self, n: f64, p: f64) -> f64 {
        if n <= 0.0 || p <= 0.0 {
            return 0.0;
        }
        if p >= 1.0 {
            return n;
        }
        if n <= 50.0 {
            let n_int = n.floor() as usize;
            let frac = n - (n_int as f64);
            let mut hits = 0.0;
            for _ in 0..n_int {
                if self.uniform_open() < p {
                    hits += 1.0;
                }
            }
            if frac > 1e-12 && self.uniform_open() < p {
                hits += frac;
            }
            return hits.clamp(0.0, n);
        }

        let mu = n * p;
        let q = 1.0 - p;
        let mu_comp = n * q;
        if mu >= 10.0 && mu_comp >= 10.0 {
            let sd = (mu * q).sqrt();
            (mu + sd * self.normal()).clamp(0.0, n)
        } else if mu < 10.0 {
            (self.poisson_small(mu) as f64).min(n)
        } else {
            let comp = (self.poisson_small(mu_comp) as f64).min(n);
            (n - comp).max(0.0)
        }
    }

    /// Allocate mass `c_g` (integer or fractional) across `paths` with normalized probabilities `probs`.
    ///
    /// - Preserves exact total mass `sum(alloc) == c_g`.
    /// - Runs in `O(K)` time independent of `c_g` depth when `c_g > 50` via sequential Binomial reduction.
    fn multinomial_allocate(&mut self, c_g: f64, probs: &[f64], paths: &[usize], counts: &mut [f64]) {
        if c_g <= 0.0 || paths.is_empty() {
            return;
        }
        let k = paths.len();
        if k == 1 {
            counts[paths[0]] += c_g;
            return;
        }

        if c_g <= 50.0 {
            let n_int = c_g.floor() as usize;
            let frac = c_g - (n_int as f64);
            for _ in 0..n_int {
                let idx = self.categorical(probs);
                counts[paths[idx]] += 1.0;
            }
            if frac > 1e-12 {
                let idx = self.categorical(probs);
                counts[paths[idx]] += frac;
            }
            return;
        }

        // Depth-independent O(K) sequential conditional Binomial reduction
        let mut rem_n = c_g;
        let mut rem_p = 1.0f64;
        for j in 0..(k - 1) {
            if rem_n <= 0.0 {
                break;
            }
            let pj = probs[j];
            if pj <= 0.0 {
                continue;
            }
            let p_cond = if rem_p > 0.0 {
                (pj / rem_p).clamp(0.0, 1.0)
            } else {
                1.0
            };
            let draw = self.binomial_continuous(rem_n, p_cond);
            counts[paths[j]] += draw;
            rem_n = (rem_n - draw).max(0.0);
            rem_p = (rem_p - pj).max(0.0);
        }
        if rem_n > 0.0 {
            counts[paths[k - 1]] += rem_n;
        }
    }
}

/// Type-7 quantile interpolation on a sorted non-empty slice (`q in [0, 1]`).
#[inline]
fn quantile_sorted(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return sorted[0];
    }
    let pos = q.clamp(0.0, 1.0) * ((n - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = pos - (lo as f64);
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

/// Run the targeted collapsed Gibbs sampler for a complex (`K_n >= 2`) local splice-graph event.
///
/// Returns `None` if `num_samples < 2` or if `inc`/`exc` is empty (caller should fall back to closed-form).
#[allow(clippy::too_many_arguments)]
pub fn sample_complex_event_gibbs(
    inc: &[ReducedPath],
    exc: &[ReducedPath],
    ambig: &[(Vec<usize>, f64)],
    init_psi: &[f64],
    vb_prior: Option<f64>,
    num_samples: usize,
    burn_in: usize,
    thin: usize,
    global_seed: u64,
    gene_id: u32,
    node_1based: u32,
) -> Option<GibbsPosteriorSummary> {
    let n_inc = inc.len();
    let n_exc = exc.len();
    let n_tot = n_inc + n_exc;
    if num_samples < 2 || n_inc == 0 || n_exc == 0 || init_psi.len() != n_tot {
        return None;
    }

    let mut uniq = vec![0.0f64; n_tot];
    let mut eff_len = vec![1.0f64; n_tot];
    let mut alpha = vec![0.0f64; n_tot];

    // Group-balanced prior preserves symmetric Beta(1, 1) marginal prior on Psi unless
    // `--complex-vb-prior` is explicitly supplied.
    let inc_alpha = vb_prior.map(|a| a.max(1e-6)).unwrap_or(1.0 / (n_inc as f64));
    let exc_alpha = vb_prior.map(|a| a.max(1e-6)).unwrap_or(1.0 / (n_exc as f64));

    for (i, p) in inc.iter().enumerate() {
        uniq[i] = p.count.max(0.0);
        eff_len[i] = p.eff_len.max(1e-6);
        alpha[i] = inc_alpha;
    }
    for (i, p) in exc.iter().enumerate() {
        uniq[n_inc + i] = p.count.max(0.0);
        eff_len[n_inc + i] = p.eff_len.max(1e-6);
        alpha[n_inc + i] = exc_alpha;
    }

    let mut rng = EventRng::from_event(global_seed, gene_id, node_1based);
    let mut pi = init_psi.to_vec();
    let pi_sum: f64 = pi.iter().sum();
    if pi_sum <= 0.0 || !pi_sum.is_finite() {
        let u = 1.0 / (n_tot as f64);
        pi.fill(u);
    } else {
        for v in &mut pi {
            *v = (*v).max(0.0) / pi_sum;
        }
    }

    let mut counts = vec![0.0f64; n_tot];
    let mut theta = vec![0.0f64; n_tot];
    let mut group_probs = Vec::with_capacity(n_tot);
    let mut psi_draws = Vec::with_capacity(num_samples);

    let thin_step = thin.max(1);
    let total_sweeps = burn_in + num_samples * thin_step;

    for sweep in 0..total_sweeps {
        // Step 1: Allocate each ambiguous edge group c_g across compatible paths ~ Multinomial(c_g, pi_G)
        counts.copy_from_slice(&uniq);
        for (paths, w) in ambig {
            if *w <= 0.0 || paths.is_empty() {
                continue;
            }
            group_probs.clear();
            let mut denom = 0.0;
            for &pid in paths {
                let v = pi[pid].max(0.0);
                group_probs.push(v);
                denom += v;
            }
            if denom > 0.0 {
                let inv = 1.0 / denom;
                for v in &mut group_probs {
                    *v *= inv;
                }
            } else {
                let unif = 1.0 / (paths.len() as f64);
                group_probs.fill(unif);
            }
            rng.multinomial_allocate(*w, &group_probs, paths, &mut counts);
        }

        // Step 2: Sample path expression rates with exact rate = eff_len[p]:
        //   theta_p = Gamma(alpha_p + z_p, scale = 1.0) / eff_len[p]
        let mut theta_sum = 0.0;
        for p in 0..n_tot {
            let shape = alpha[p] + counts[p].max(0.0);
            let t = rng.gamma_unit_scale(shape) / eff_len[p];
            theta[p] = t;
            theta_sum += t;
        }

        if theta_sum > 0.0 && theta_sum.is_finite() {
            let inv = 1.0 / theta_sum;
            for p in 0..n_tot {
                pi[p] = theta[p] * inv;
            }
        } else {
            let unif = 1.0 / (n_tot as f64);
            pi.fill(unif);
        }

        if sweep >= burn_in && (sweep - burn_in).is_multiple_of(thin_step) {
            let psi_b: f64 = pi[..n_inc].iter().sum::<f64>().clamp(0.0, 1.0);
            psi_draws.push(psi_b);
        }
    }

    let b = psi_draws.len();
    if b < 2 {
        return None;
    }

    let mean_psi = psi_draws.iter().sum::<f64>() / (b as f64);
    let ss: f64 = psi_draws
        .iter()
        .map(|&x| {
            let d = x - mean_psi;
            d * d
        })
        .sum();
    let inf_var = (ss / ((b - 1) as f64)).max(0.0);

    psi_draws.sort_unstable_by(|a, b| a.total_cmp(b));
    let ci_lo = quantile_sorted(&psi_draws, 0.05).clamp(0.0, 1.0);
    let ci_hi = quantile_sorted(&psi_draws, 0.95).clamp(0.0, 1.0);

    Some(GibbsPosteriorSummary {
        mean_psi,
        inf_var,
        ci_lo,
        ci_hi,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::events::psi::beta_posterior_var;

    #[test]
    fn test_gibbs_balanced_prior_matches_beta_when_unambiguous() {
        // Equal effective lengths & zero ambiguous reads across a K_3 event (4 inc paths, 4 exc paths):
        // With group-balanced priors (1/4 each), marginal posterior on Psi is EXACTLY Beta(1 + 30, 1 + 70)!
        let inc: Vec<ReducedPath> = (0..4)
            .map(|i| ReducedPath {
                nodes: vec![0, 1 + i, 9],
                count: 7.5, // sum = 30.0
                eff_len: 100.0,
                psi: 0.075,
            })
            .collect();
        let exc: Vec<ReducedPath> = (0..4)
            .map(|i| ReducedPath {
                nodes: vec![0, 5 + i, 9],
                count: 17.5, // sum = 70.0
                eff_len: 100.0,
                psi: 0.175,
            })
            .collect();
        let init_psi = vec![0.075, 0.075, 0.075, 0.075, 0.175, 0.175, 0.175, 0.175];

        let summary =
            sample_complex_event_gibbs(&inc, &exc, &[], &init_psi, None, 2000, 32, 1, 42, 7, 3)
                .unwrap();

        let expected_mean = 31.0 / 102.0; // Beta(31, 71)
        let expected_var = beta_posterior_var(0.30, 100.0);
        assert!(
            (summary.mean_psi - expected_mean).abs() < 0.01,
            "balanced-prior K_3 mean {} deviated from exact Beta(31,71) mean {}",
            summary.mean_psi,
            expected_mean
        );
        assert!(
            (summary.inf_var / expected_var - 1.0).abs() < 0.15,
            "balanced-prior K_3 variance {} deviated from exact Beta(31,71) var {}",
            summary.inf_var,
            expected_var
        );
    }

    #[test]
    fn test_gibbs_rate_parameterization_and_high_depth_o1_speed() {
        // Verify eff_len acts as Gamma RATE (theta = Gamma / eff_len):
        // 1 inc path with count 1000, eff_len 200 -> rate 5
        // 1 exc path with count 1000, eff_len 100 -> rate 10
        // => true Psi = 5 / (5 + 10) = 1/3, NOT 2/3!
        // Also test very large c_g = 100_000.5 to confirm O(K) depth-independent execution.
        let inc = vec![
            ReducedPath {
                nodes: vec![0, 1, 3],
                count: 500.0,
                eff_len: 200.0,
                psi: 0.1667,
            },
            ReducedPath {
                nodes: vec![0, 1, 2, 3],
                count: 500.0,
                eff_len: 200.0,
                psi: 0.1667,
            },
        ];
        let exc = vec![
            ReducedPath {
                nodes: vec![0, 3],
                count: 1000.0,
                eff_len: 100.0,
                psi: 0.6667,
            },
            ReducedPath {
                nodes: vec![0, 4, 3],
                count: 0.0,
                eff_len: 100.0,
                psi: 0.0,
            },
        ];
        let ambig = vec![(vec![0, 1, 2], 100_000.5)];
        let init_psi = vec![0.17, 0.17, 0.66, 0.0];
        let s = sample_complex_event_gibbs(&inc, &exc, &ambig, &init_psi, None, 100, 32, 2, 42, 1, 2)
            .unwrap();
        assert!(s.mean_psi.is_finite() && s.inf_var.is_finite());
        assert!(s.ci_lo <= s.mean_psi && s.mean_psi <= s.ci_hi);

        // Pure rate check (zero ambig)
        let s_rate =
            sample_complex_event_gibbs(&inc, &exc, &[], &init_psi, None, 500, 32, 1, 42, 1, 2)
                .unwrap();
        assert!(
            (s_rate.mean_psi - (1.0 / 3.0)).abs() < 0.02,
            "Gamma rate parameterization wrong: expected ~0.333, got {}",
            s_rate.mean_psi
        );
    }
}
