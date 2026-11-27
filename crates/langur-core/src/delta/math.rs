//! Numerical core for Beta-Binomial GLM: exact `trigamma`, continuous Beta-Binomial
//! log-likelihood, damped Fisher scoring (IRLS) with light ridge regularization,
//! dense `Matrix` with Cholesky & column-pivoted Householder QR, Cox-Reid APL,
//! 1D grid + golden-section dispersion maximizer, Empirical Bayes prior estimation
//! with Smyth finite-df log-chi2 bias correction, and Benjamini-Hochberg FDR.

use statrs::distribution::{ChiSquared, ContinuousCDF};
use statrs::function::gamma::{digamma, ln_gamma};

/// `theta = log(M)` bounds matching `Whippet2.jl/src/diff_glm.jl`:
/// `M = exp(-4) ≈ 0.0183` (`phi ≈ 0.982`) to `M = exp(11.5) ≈ 98,715` (`phi ≈ 1e-5`).
pub const THETA_MIN: f64 = -4.0;
pub const THETA_MAX: f64 = 11.5;
pub const MU_EPS: f64 = 1e-9;

/// Lightweight row-major dense matrix (`rows × cols`) for small GLM design / Fisher matrices.
#[derive(Clone, Debug)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl Matrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    pub fn from_row_major(rows: usize, cols: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), rows * cols);
        Self { rows, cols, data }
    }

    #[inline]
    pub fn get(&self, r: usize, c: usize) -> f64 {
        self.data[r * self.cols + c]
    }

    #[inline]
    pub fn set(&mut self, r: usize, c: usize, val: f64) {
        self.data[r * self.cols + c] = val;
    }

    #[inline]
    pub fn add_assign(&mut self, r: usize, c: usize, val: f64) {
        self.data[r * self.cols + c] += val;
    }

    pub fn fill(&mut self, val: f64) {
        self.data.fill(val);
    }

    /// Return column `c` as a freshly allocated `Vec<f64>`.
    pub fn col_vec(&self, c: usize) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.rows);
        for r in 0..self.rows {
            out.push(self.get(r, c));
        }
        out
    }

    /// Build a new matrix from a list of equal-length columns.
    pub fn from_cols(rows: usize, cols_data: &[Vec<f64>]) -> Self {
        let cols = cols_data.len();
        let mut m = Self::zeros(rows, cols);
        for (c, col) in cols_data.iter().enumerate() {
            debug_assert_eq!(col.len(), rows);
            for (r, &v) in col.iter().enumerate() {
                m.set(r, c, v);
            }
        }
        m
    }
}

#[inline]
pub fn logistic(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let ex = x.exp();
        ex / (1.0 + ex)
    }
}

#[inline]
pub fn logit(p: f64) -> f64 {
    let pc = p.clamp(MU_EPS, 1.0 - MU_EPS);
    (pc / (1.0 - pc)).ln()
}

/// Exact trigamma function `psi_1(x) = d^2/dx^2 ln Gamma(x)` for `x > 0`
/// to `< 1e-14` relative accuracy via recurrence to `x >= 8` + 7-term Bernoulli series.
#[inline]
pub fn trigamma(mut x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    let mut sum = 0.0;
    while x < 8.0 {
        sum += 1.0 / (x * x);
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    // Asymptotic Bernoulli expansion:
    // 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + 1/(42x^7) - 1/(30x^9) + 5/(66x^11) - 691/(2730x^13)
    let asymp = inv
        + 0.5 * inv2
        + inv
            * inv2
            * (1.0 / 6.0
                + inv2
                    * (-1.0 / 30.0
                        + inv2
                            * (1.0 / 42.0
                                + inv2
                                    * (-1.0 / 30.0
                                        + inv2 * (5.0 / 66.0 - (691.0 / 2730.0) * inv2)))));
    sum + asymp
}

/// Lower-triangular Cholesky decomposition `A = L L^T` for symmetric positive definite `p × p` matrix.
/// Returns `Some(L)` (row-major `p × p` lower triangular) or `None` if not positive definite.
pub fn cholesky_decompose(a: &Matrix) -> Option<Matrix> {
    let p = a.rows;
    debug_assert_eq!(a.cols, p);
    let mut l = Matrix::zeros(p, p);
    for i in 0..p {
        for j in 0..=i {
            let mut sum = a.get(i, j);
            for k in 0..j {
                sum -= l.get(i, k) * l.get(j, k);
            }
            if i == j {
                if !sum.is_finite() || sum <= 0.0 {
                    return None;
                }
                l.set(i, j, sum.sqrt());
            } else {
                let ljj = l.get(j, j);
                if ljj == 0.0 {
                    return None;
                }
                l.set(i, j, sum / ljj);
            }
        }
    }
    Some(l)
}

/// Solve `L L^T x = b` given lower-triangular Cholesky factor `L`.
pub fn cholesky_solve(l: &Matrix, b: &[f64]) -> Vec<f64> {
    let p = l.rows;
    let mut y = vec![0.0; p];
    for i in 0..p {
        let mut s = b[i];
        for k in 0..i {
            s -= l.get(i, k) * y[k];
        }
        y[i] = s / l.get(i, i);
    }
    let mut x = vec![0.0; p];
    for i in (0..p).rev() {
        let mut s = y[i];
        for k in (i + 1)..p {
            s -= l.get(k, i) * x[k];
        }
        x[i] = s / l.get(i, i);
    }
    x
}

/// Compute `log det(A) = 2 * sum_i ln(L[i, i])` from Cholesky factor `L`.
#[inline]
pub fn cholesky_logdet(l: &Matrix) -> f64 {
    let mut s = 0.0;
    for i in 0..l.rows {
        s += l.get(i, i).ln();
    }
    2.0 * s
}

/// Numerical rank of an `m × p` matrix via column-pivoted Householder QR decomposition
/// (matching Julia `qr(A, ColumnNorm())` rank tolerance `max(m, p) * eps * max(|diag(R)|, 1.0)`).
pub fn matrix_rank_pivoted_qr(a: &Matrix) -> usize {
    let m = a.rows;
    let p = a.cols;
    if m == 0 || p == 0 {
        return 0;
    }
    // Work in column-major layout for fast column swapping & Householder updates
    let mut cols: Vec<Vec<f64>> = (0..p).map(|c| a.col_vec(c)).collect();
    let kmax = m.min(p);
    let mut r_diag = Vec::with_capacity(kmax);

    for k in 0..kmax {
        // Pivot: choose column >= k with maximum residual 2-norm squared in rows k..m
        let mut best_c = k;
        let mut best_norm2 = -1.0f64;
        for c in k..p {
            let mut n2 = 0.0;
            for r in k..m {
                let v = cols[c][r];
                n2 += v * v;
            }
            if n2 > best_norm2 {
                best_norm2 = n2;
                best_c = c;
            }
        }
        if best_c != k {
            cols.swap(k, best_c);
        }

        let norm = best_norm2.sqrt();
        if norm == 0.0 || !norm.is_finite() {
            r_diag.push(0.0);
            continue;
        }

        let x0 = cols[k][k];
        let alpha = if x0 >= 0.0 { -norm } else { norm };
        r_diag.push(alpha.abs());

        // Householder vector v stored in cols[k][k..m]
        let mut v = vec![0.0; m - k];
        v[0] = x0 - alpha;
        for r in (k + 1)..m {
            v[r - k] = cols[k][r];
        }
        let v_norm2: f64 = v.iter().map(|&val| val * val).sum();
        if v_norm2 <= 1e-30 {
            continue;
        }
        let beta = 2.0 / v_norm2;

        // Apply reflection to remaining columns k..p
        for c in k..p {
            let mut dot = 0.0;
            for idx in 0..(m - k) {
                dot += v[idx] * cols[c][k + idx];
            }
            let scale = beta * dot;
            for idx in 0..(m - k) {
                cols[c][k + idx] -= scale * v[idx];
            }
        }
    }

    let rmax = r_diag.iter().copied().fold(0.0f64, f64::max);
    let tol = (m.max(p) as f64) * f64::EPSILON * rmax.max(1.0);
    r_diag.iter().filter(|&&d| d > tol).count()
}

/// Continuous Beta-Binomial log-likelihood (omitting constant log-binomial coefficient
/// that cancels identically across all `beta`/`M` derivatives and LRTs).
pub fn bb_loglik(k: &[f64], r: &[f64], n: &[f64], x: &Matrix, beta: &[f64], m_prec: f64) -> f64 {
    let m = k.len();
    let p = beta.len();
    let lg_m = ln_gamma(m_prec);
    let mut ll = 0.0;
    for j in 0..m {
        let mut eta = 0.0;
        let row_off = j * p;
        for c in 0..p {
            eta += x.data[row_off + c] * beta[c];
        }
        let mu = logistic(eta).clamp(MU_EPS, 1.0 - MU_EPS);
        let a = m_prec * mu;
        let b = m_prec * (1.0 - mu);
        ll += ln_gamma(k[j] + a) + ln_gamma(r[j] + b) - ln_gamma(n[j] + m_prec) + lg_m
            - ln_gamma(a)
            - ln_gamma(b);
    }
    ll
}

/// Result of fitting `beta` at fixed precision `M`.
#[derive(Clone, Debug)]
pub struct BetaFit {
    pub beta: Vec<f64>,
    pub loglik: f64,
    pub fisher_info: Matrix,
    pub converged: bool,
}

/// Fit `beta` for fixed precision `M` using damped Fisher scoring (IRLS) with light ridge
/// regularization for stability under complete separation (`0.0` vs `1.0` PSI).
pub fn fit_bb_beta(
    k: &[f64],
    r: &[f64],
    n: &[f64],
    x: &Matrix,
    m_prec: f64,
    beta_init: Option<&[f64]>,
    ridge: f64,
    maxiter: usize,
    tol: f64,
) -> BetaFit {
    let m = x.rows;
    let p = x.cols;
    let mut beta = match beta_init {
        Some(init) => init.to_vec(),
        None => {
            let mut b = vec![0.0; p];
            if p >= 1 {
                let sk: f64 = k.iter().sum();
                let sn: f64 = n.iter().sum();
                if sn > 0.0 {
                    b[0] = logit((sk / sn).clamp(1e-4, 1.0 - 1e-4));
                }
            }
            b
        }
    };

    let pen_val = |b: &[f64]| -> f64 {
        let mut pv = 0.0;
        for (c, &bc) in b.iter().enumerate() {
            let lam = if c == 0 { ridge * 0.01 } else { ridge };
            pv += 0.5 * lam * bc * bc;
        }
        pv
    };

    let mut g = vec![0.0; p];
    let mut h = Matrix::zeros(p, p);
    let mut beta_try = vec![0.0; p];

    let mut ll_cur = bb_loglik(k, r, n, x, &beta, m_prec);
    let mut obj_cur = ll_cur - pen_val(&beta);
    let mut converged = false;

    for _iter in 0..maxiter {
        g.fill(0.0);
        h.fill(0.0);

        for j in 0..m {
            let row_off = j * p;
            let mut eta = 0.0;
            for c in 0..p {
                eta += x.data[row_off + c] * beta[c];
            }
            let mu = logistic(eta).clamp(MU_EPS, 1.0 - MU_EPS);
            let mu1 = mu * (1.0 - mu);
            let a = m_prec * mu;
            let b = m_prec * (1.0 - mu);

            let score_j =
                m_prec * mu1 * (digamma(k[j] + a) - digamma(a) - digamma(r[j] + b) + digamma(b));
            let t_term = (trigamma(a) - trigamma(k[j] + a)) + (trigamma(b) - trigamma(r[j] + b));
            let w_j = (m_prec * m_prec * mu1 * mu1 * t_term).max(1e-12);

            for c in 0..p {
                let xjc = x.data[row_off + c];
                g[c] += xjc * score_j;
                let wx = w_j * xjc;
                for d in c..p {
                    h.add_assign(c, d, wx * x.data[row_off + d]);
                }
            }
        }

        for c in 0..p {
            let lam = if c == 0 { ridge * 0.01 } else { ridge };
            g[c] -= lam * beta[c];
            h.add_assign(c, c, lam);
            for d in (c + 1)..p {
                let val = h.get(c, d);
                h.set(d, c, val);
            }
        }

        // Cholesky solve with damped-diagonal fallback
        let step_opt = if let Some(l) = cholesky_decompose(&h) {
            Some(cholesky_solve(&l, &g))
        } else {
            let mut h_reg = h.clone();
            let mut bump = 1e-6;
            let mut res = None;
            for _ in 0..8 {
                for c in 0..p {
                    h_reg.set(c, c, h.get(c, c) + bump);
                }
                if let Some(l) = cholesky_decompose(&h_reg) {
                    res = Some(cholesky_solve(&l, &g));
                    break;
                }
                bump *= 10.0;
            }
            res
        };

        let Some(mut step) = step_opt else {
            break;
        };

        // Clamp extreme Newton step for stability
        let max_step = step.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        if max_step > 5.0 {
            let scale = 5.0 / max_step;
            for s in &mut step {
                *s *= scale;
            }
        }

        // Backtracking line search ensuring monotone ascent
        let mut alpha = 1.0;
        let mut accepted = false;
        let mut ll_new = ll_cur;
        let mut obj_new = obj_cur;
        for _ in 0..20 {
            for c in 0..p {
                beta_try[c] = beta[c] + alpha * step[c];
            }
            ll_new = bb_loglik(k, r, n, x, &beta_try, m_prec);
            obj_new = ll_new - pen_val(&beta_try);
            if obj_new.is_finite() && obj_new >= obj_cur - 1e-12 {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            converged = true;
            break;
        }

        let mut max_delta = 0.0f64;
        for c in 0..p {
            let d = (beta_try[c] - beta[c]).abs();
            if d > max_delta {
                max_delta = d;
            }
            beta[c] = beta_try[c];
        }
        let d_obj = (obj_new - obj_cur).abs();
        ll_cur = ll_new;
        obj_cur = obj_new;

        if max_delta < tol || d_obj < tol * (1.0 + obj_cur.abs()) {
            converged = true;
            break;
        }
    }

    // Recompute unpenalized Fisher information H at converged beta (for Cox-Reid logdet)
    h.fill(0.0);
    for j in 0..m {
        let row_off = j * p;
        let mut eta = 0.0;
        for c in 0..p {
            eta += x.data[row_off + c] * beta[c];
        }
        let mu = logistic(eta).clamp(MU_EPS, 1.0 - MU_EPS);
        let mu1 = mu * (1.0 - mu);
        let a = m_prec * mu;
        let b = m_prec * (1.0 - mu);
        let t_term = (trigamma(a) - trigamma(k[j] + a)) + (trigamma(b) - trigamma(r[j] + b));
        let w_j = (m_prec * m_prec * mu1 * mu1 * t_term).max(1e-12);
        for c in 0..p {
            let wx = w_j * x.data[row_off + c];
            for d in c..p {
                h.add_assign(c, d, wx * x.data[row_off + d]);
            }
        }
    }
    for c in 0..p {
        for d in (c + 1)..p {
            let val = h.get(c, d);
            h.set(d, c, val);
        }
    }

    BetaFit {
        beta,
        loglik: ll_cur,
        fisher_info: h,
        converged,
    }
}

/// Cox-Reid Adjusted Profile Log-Likelihood at `theta = log(M)`:
/// `APL(theta) = ll(beta_hat(theta), exp(theta)) - 0.5 * logdet(X' W X + 1e-8 I)`.
pub fn bb_cox_reid_apl(
    k: &[f64],
    r: &[f64],
    n: &[f64],
    x: &Matrix,
    theta: f64,
    beta_init: Option<&[f64]>,
    ridge: f64,
) -> (f64, Vec<f64>) {
    let m_prec = theta.exp();
    let fit = fit_bb_beta(k, r, n, x, m_prec, beta_init, ridge, 40, 1e-8);
    let p = x.cols;
    let mut h_cr = fit.fisher_info;
    for c in 0..p {
        h_cr.add_assign(c, c, 1e-8);
    }
    let cr_pen = match cholesky_decompose(&h_cr) {
        Some(l) => 0.5 * cholesky_logdet(&l),
        None => 0.0,
    };
    (fit.loglik - cr_pen, fit.beta)
}

/// 1D maximizer on `[lo, hi]` via coarse grid search + golden-section refinement.
pub fn maximize_theta<F: FnMut(f64) -> f64>(
    mut f: F,
    lo: f64,
    hi: f64,
    ngrid: usize,
    tol: f64,
) -> (f64, f64) {
    let ngrid = ngrid.max(3);
    let mut grid = Vec::with_capacity(ngrid);
    let step = (hi - lo) / ((ngrid - 1) as f64);
    for i in 0..ngrid {
        grid.push(lo + (i as f64) * step);
    }

    let mut best_t = lo;
    let mut best_v = f64::NEG_INFINITY;
    let mut best_idx = 0usize;
    for (i, &t) in grid.iter().enumerate() {
        let v = f(t);
        if v > best_v {
            best_v = v;
            best_t = t;
            best_idx = i;
        }
    }

    let mut a = if best_idx > 0 { grid[best_idx - 1] } else { lo };
    let mut b = if best_idx + 1 < ngrid {
        grid[best_idx + 1]
    } else {
        hi
    };

    const INVPHI: f64 = 0.618_033_988_749_894_9;
    const INVPHI2: f64 = 0.381_966_011_250_105_1;

    let mut h = b - a;
    let mut c = a + INVPHI2 * h;
    let mut d = a + INVPHI * h;
    let mut fc = f(c);
    let mut fd = f(d);

    for _ in 0..35 {
        if (b - a) < tol {
            break;
        }
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            h = b - a;
            c = a + INVPHI2 * h;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            h = b - a;
            d = a + INVPHI * h;
            fd = f(d);
        }
    }

    let t_opt = 0.5 * (a + b);
    let v_opt = f(t_opt);
    if best_v > v_opt {
        (best_t, best_v)
    } else {
        (t_opt, v_opt)
    }
}

/// Estimate raw Cox-Reid dispersion `theta_hat = log(M_hat)` and its local sampling variance
/// `1 / (-APL''(theta_hat))` via 3-point finite difference if interior.
pub fn estimate_raw_theta(
    k: &[f64],
    r: &[f64],
    n: &[f64],
    x: &Matrix,
    ridge: f64,
) -> (f64, bool, f64) {
    let mut beta_warm: Option<Vec<f64>> = None;
    let mut eval_apl = |th: f64| -> f64 {
        let (val, b) = bb_cox_reid_apl(k, r, n, x, th, beta_warm.as_deref(), ridge);
        beta_warm = Some(b);
        val
    };

    let (th_hat, apl_hat) = maximize_theta(&mut eval_apl, THETA_MIN, THETA_MAX, 9, 2e-3);
    let is_interior = th_hat > (THETA_MIN + 0.05) && th_hat < (THETA_MAX - 0.05);
    let mut samp_var = f64::NAN;
    if is_interior {
        let h = 0.1;
        let v_lo = eval_apl(th_hat - h);
        let v_hi = eval_apl(th_hat + h);
        let curv = -(v_hi - 2.0 * apl_hat + v_lo) / (h * h);
        if curv > 1e-4 {
            samp_var = 1.0 / curv;
        }
    }
    (th_hat, is_interior, samp_var)
}

/// Estimate Empirical Bayes MAP shrunk `theta = log(M)` under `Normal(prior_mu, prior_var)`.
pub fn estimate_eb_theta(
    k: &[f64],
    r: &[f64],
    n: &[f64],
    x: &Matrix,
    prior_mu: f64,
    prior_var: f64,
    ridge: f64,
) -> f64 {
    let mut beta_warm: Option<Vec<f64>> = None;
    let inv_2v = 0.5 / prior_var;
    let mut eval_post = |th: f64| -> f64 {
        let (val, b) = bb_cox_reid_apl(k, r, n, x, th, beta_warm.as_deref(), ridge);
        beta_warm = Some(b);
        let diff = th - prior_mu;
        val - diff * diff * inv_2v
    };
    let (th_map, _) = maximize_theta(&mut eval_post, THETA_MIN, THETA_MAX, 9, 2e-3);
    th_map
}

fn median_in_place(vals: &mut [f64]) -> f64 {
    if vals.is_empty() {
        return f64::NAN;
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    let n = vals.len();
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        0.5 * (vals[n / 2 - 1] + vals[n / 2])
    }
}

/// Fit genome-wide Empirical Bayes prior `Normal(mu_theta, var_theta)` across raw event estimates,
/// applying the exact Smyth finite-df `log(median(Chisq(d_res)) / d_res)` median bias correction.
pub fn fit_eb_dispersion_prior(
    raw_thetas: &[f64],
    is_interior: &[bool],
    samp_vars: &[f64],
    res_dfs: &[usize],
) -> (f64, f64, usize) {
    let mut interior_th = Vec::new();
    let mut interior_sv = Vec::new();
    let mut interior_df = Vec::new();
    let mut all_valid = Vec::new();
    let has_dfs = res_dfs.len() == raw_thetas.len();

    for i in 0..raw_thetas.len() {
        let t = raw_thetas[i];
        if !t.is_finite() {
            continue;
        }
        all_valid.push(t);
        if is_interior[i] {
            interior_th.push(t);
            if samp_vars[i].is_finite() {
                interior_sv.push(samp_vars[i]);
            }
            if has_dfs && res_dfs[i] >= 1 {
                interior_df.push(res_dfs[i] as f64);
            }
        }
    }

    let mut bias_corr = 0.0;
    if !interior_df.is_empty() {
        let med_df = median_in_place(&mut interior_df).max(1.0);
        if let Ok(chisq) = ChiSquared::new(med_df) {
            bias_corr = (chisq.inverse_cdf(0.5) / med_df).ln();
        }
    }

    if interior_th.len() >= 20 {
        let n_used = interior_th.len();
        let raw_med = median_in_place(&mut interior_th);
        let mu_prior = (raw_med + bias_corr).clamp(THETA_MIN, THETA_MAX);
        let mut abs_dev: Vec<f64> = interior_th.iter().map(|&x| (x - raw_med).abs()).collect();
        let mad = median_in_place(&mut abs_dev);
        let obs_sd = 1.4826 * mad;
        let med_sv = if interior_sv.is_empty() {
            0.0
        } else {
            median_in_place(&mut interior_sv)
        };
        let var_prior = (obs_sd * obs_sd - med_sv).clamp(0.08, 9.0);
        (mu_prior, var_prior, n_used)
    } else if all_valid.len() >= 10 {
        let n_used = all_valid.len();
        let raw_med = median_in_place(&mut all_valid);
        let mu_prior = (raw_med + bias_corr).clamp(THETA_MIN, THETA_MAX);
        let mut abs_dev: Vec<f64> = all_valid.iter().map(|&x| (x - raw_med).abs()).collect();
        let mad = median_in_place(&mut abs_dev);
        let obs_sd = 1.4826 * mad;
        let var_prior = (obs_sd * obs_sd).clamp(0.25, 9.0);
        (mu_prior, var_prior, n_used)
    } else {
        // Weakly informative fallback prior (M ≈ 20 => phi ≈ 0.047)
        (3.0, 4.0, all_valid.len())
    }
}

/// Benjamini-Hochberg FDR adjustment preserving `NaN` entries identically.
pub fn bh_fdr(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let mut fdr = vec![f64::NAN; n];
    let mut valid_idx: Vec<usize> = pvals
        .iter()
        .enumerate()
        .filter(|(_, &p)| p.is_finite())
        .map(|(i, _)| i)
        .collect();
    let m = valid_idx.len();
    if m == 0 {
        return fdr;
    }

    valid_idx.sort_by(|&a, &b| pvals[a].total_cmp(&pvals[b]));
    let mut running_min = 1.0f64;
    for rank_0 in (0..m).rev() {
        let idx = valid_idx[rank_0];
        let rank = (rank_0 + 1) as f64;
        let q = pvals[idx] * ((m as f64) / rank);
        if q < running_min {
            running_min = q;
        }
        fdr[idx] = running_min.clamp(0.0, 1.0);
    }
    fdr
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trigamma_high_accuracy() {
        // Known exact values: trigamma(1) = pi^2 / 6, trigamma(2) = pi^2/6 - 1, trigamma(0.5) = pi^2 / 2
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        assert!((trigamma(1.0) - pi2 / 6.0).abs() < 1e-13);
        assert!((trigamma(2.0) - (pi2 / 6.0 - 1.0)).abs() < 1e-13);
        assert!((trigamma(0.5) - pi2 / 2.0).abs() < 1e-13);
        // Check recurrence trigamma(x) - trigamma(x + 1) == 1 / x^2 at x = 10.5
        let x = 10.5;
        assert!(((trigamma(x) - trigamma(x + 1.0)) - 1.0 / (x * x)).abs() < 1e-14);
    }

    #[test]
    fn test_pivoted_qr_rank_and_cholesky() {
        // Rank-2 matrix with 3 columns (col 3 = col 1 + col 2)
        let a = Matrix::from_cols(
            4,
            &[
                vec![1.0, 1.0, 1.0, 1.0],
                vec![0.0, 1.0, 0.0, 1.0],
                vec![1.0, 2.0, 1.0, 2.0],
            ],
        );
        assert_eq!(matrix_rank_pivoted_qr(&a), 2);

        let spd = Matrix::from_row_major(2, 2, vec![4.0, 2.0, 2.0, 3.0]);
        let l = cholesky_decompose(&spd).unwrap();
        let x = cholesky_solve(&l, &[8.0, 7.0]);
        // 4*1.25 + 2*1.5 = 8, 2*1.25 + 3*1.5 = 7
        assert!((x[0] - 1.25).abs() < 1e-12);
        assert!((x[1] - 1.50).abs() < 1e-12);
        assert!((cholesky_logdet(&l) - (8.0f64).ln()).abs() < 1e-12);
    }
}
