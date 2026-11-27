//! 3-pass `rayon`-parallel Cox-Reid Empirical Bayes Beta-Binomial GLM engine,
//! finite-sample F(df, d_res + d_0) LRT, covariate-adjusted marginal Psi_A/Psi_B/DeltaPsi,
//! Benjamini-Hochberg FDR, and parallel RFC 1952 `.diff.gz` writer.

use crate::delta::cohort::{read_psi_cohort, CohortReadOptions, GlmEventData};
use crate::delta::design::{build_event_matrices, DesignSpec};
use crate::delta::math::{
    bb_loglik, bh_fdr, cholesky_decompose, cholesky_solve, estimate_eb_theta, estimate_raw_theta,
    fit_bb_beta, fit_eb_dispersion_prior, logistic, Matrix, THETA_MAX, THETA_MIN,
};
use crate::error::{Result, LangurError};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use statrs::distribution::{ChiSquared, ContinuousCDF, FisherSnedecor};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct GlmOptions {
    pub cohort: CohortReadOptions,
    /// Enable genome-wide Empirical Bayes dispersion shrinkage (default `true`).
    pub shrink_dispersion: bool,
    /// Stratify Empirical Bayes dispersion prior into Spliced (`CE/AA/AD/RI/AF/AL/BS`) vs.
    /// Tandem UTR (`TS/TE`) events, falling back to combined polymorphic prior if a stratum
    /// has `< 20` interior polymorphic events (default `true`).
    pub stratify_eb: bool,
    /// If `Some(phi)` with `0 < phi < 1`, use fixed Beta-Binomial overdispersion `phi` (ChiSquared LRT).
    pub fixed_phi: Option<f64>,
    /// Light ridge penalty on non-intercept coefficients for complete separation stability (default `1e-4`).
    pub ridge: f64,
}

impl Default for GlmOptions {
    fn default() -> Self {
        Self {
            cohort: CohortReadOptions::default(),
            shrink_dispersion: true,
            stratify_eb: true,
            fixed_phi: None,
            ridge: 1e-4,
        }
    }
}

/// Fitted Beta-Binomial GLM result for a single splicing event.
#[derive(Clone, Debug)]
pub struct GlmResult {
    pub tested: bool,
    pub psi_a: f64,
    pub psi_b: f64,
    pub delta_psi: f64,
    pub log_odds: f64,
    pub se: f64,
    pub dispersion: f64,
    pub lrt: f64,
    pub df: usize,
    pub pvalue: f64,
}

impl GlmResult {
    pub const fn empty() -> Self {
        Self {
            tested: false,
            psi_a: f64::NAN,
            psi_b: f64::NAN,
            delta_psi: f64::NAN,
            log_odds: f64::NAN,
            se: f64::NAN,
            dispersion: f64::NAN,
            lrt: f64::NAN,
            df: 0,
            pvalue: f64::NAN,
        }
    }
}

/// Summary of one Empirical Bayes dispersion prior stratum.
#[derive(Copy, Clone, Debug)]
pub struct EbStratumPrior {
    pub prior_mu: f64,
    pub prior_var: f64,
    pub n_prior_events: usize,
    pub df_prior: f64,
}

/// Summary of the fitted Empirical Bayes dispersion prior (primary/spliced + optional `TS`/`TE` UTR stratum).
#[derive(Copy, Clone, Debug)]
pub struct EbPriorSummary {
    pub prior_mu: f64,
    pub prior_var: f64,
    pub n_prior_events: usize,
    pub df_prior: f64,
    /// Separate Tandem UTR (`TS`/`TE`) EB prior when `stratify_eb` is enabled and `>= 20` interior UTR events exist.
    pub utr_prior: Option<EbStratumPrior>,
}

/// Returns `true` iff an event is polymorphic (non-constitutive) and carries statistical
/// information about Beta-Binomial overdispersion `phi`.
///
/// Constitutive/boundary events (`Psi = 0.0` or `1.0` across all samples) have a flat
/// Beta-Binomial log-likelihood in `M`, causing Cox-Reid `-0.5 * logdet(X'WX)` to create
/// an artificial interior peak near `phi ~ 0.64` that collapses the genome-wide EB prior.
#[inline]
pub fn is_polymorphic_event(ev: &GlmEventData) -> bool {
    let mut sum_k = 0.0f64;
    let mut sum_r = 0.0f64;
    let mut sum_n = 0.0f64;
    let mut min_psi = f64::INFINITY;
    let mut max_psi = f64::NEG_INFINITY;
    for (&p, &n) in ev.psi.iter().zip(ev.eff_reads.iter()) {
        if !p.is_finite() || !n.is_finite() || n <= 0.0 {
            continue;
        }
        let kj = p * n;
        sum_k += kj;
        sum_r += n - kj;
        sum_n += n;
        if p < min_psi {
            min_psi = p;
        }
        if p > max_psi {
            max_psi = p;
        }
    }
    if sum_n <= 0.0 || sum_k < 1.0 || sum_r < 1.0 {
        return false;
    }
    let pooled_psi = sum_k / sum_n;
    if pooled_psi <= 0.01 || pooled_psi >= 0.99 {
        return false;
    }
    (max_psi - min_psi) >= 1e-4
}

#[inline]
fn is_utr_event(ev: &GlmEventData) -> bool {
    matches!(ev.etype.as_str(), "TS" | "TE")
}

#[allow(clippy::too_many_arguments)]
fn test_single_event_with_matrices(
    ev: &GlmEventData,
    design: &DesignSpec,
    x_full: &Matrix,
    x_null: &Matrix,
    p0: usize,
    theta_use: f64,
    df_denom: f64,
    ridge: f64,
) -> GlmResult {
    let p_full = x_full.cols;
    let df = p_full.saturating_sub(p0);
    if df < 1 {
        return GlmResult::empty();
    }

    let m = ev.valid_idx.len();
    let n = &ev.eff_reads;
    let mut k = Vec::with_capacity(m);
    let mut r = Vec::with_capacity(m);
    for j in 0..m {
        let kj = ev.psi[j] * n[j];
        k.push(kj);
        r.push(n[j] - kj);
    }

    let m_prec = theta_use.exp();
    let phi = 1.0 / (1.0 + m_prec);

    // Full model fit at shrunk precision M
    let fit_full = fit_bb_beta(&k, &r, n, x_full, m_prec, None, ridge, 40, 1e-8);

    // Null model fit at the same M
    let ll_null = if p0 > 0 {
        let init_null = &fit_full.beta[..p0];
        let fit_null = fit_bb_beta(&k, &r, n, x_null, m_prec, Some(init_null), ridge, 40, 1e-8);
        fit_null.loglik
    } else {
        bb_loglik(&k, &r, n, &Matrix::zeros(m, 1), &[0.0], m_prec)
    };

    let lrt = (2.0 * (fit_full.loglik - ll_null)).max(0.0);
    let pvalue = if !lrt.is_finite() {
        f64::NAN
    } else if df_denom.is_finite() && df_denom > 0.0 {
        FisherSnedecor::new(df as f64, df_denom)
            .map(|dist| dist.sf(lrt / (df as f64)))
            .unwrap_or(f64::NAN)
    } else {
        ChiSquared::new(df as f64)
            .map(|dist| dist.sf(lrt))
            .unwrap_or(f64::NAN)
    };

    // Extract Wald log-odds ratio & SE for 1-df tests
    let mut log_odds = f64::NAN;
    let mut se = f64::NAN;
    if df == 1 {
        let t_pos = p_full - 1;
        log_odds = fit_full.beta[t_pos];
        let mut h_inv = fit_full.fisher_info.clone();
        for c in 0..p_full {
            h_inv.add_assign(c, c, if c == 0 { ridge * 0.01 } else { ridge });
        }
        if let Some(l) = cholesky_decompose(&h_inv) {
            let mut e_vec = vec![0.0; p_full];
            e_vec[t_pos] = 1.0;
            let v_col = cholesky_solve(&l, &e_vec);
            let v_diag = v_col[t_pos];
            if v_diag >= 0.0 && v_diag.is_finite() {
                se = v_diag.sqrt();
            }
        }
    }

    // Covariate-adjusted marginal Psi_A, Psi_B, and DeltaPsi (identical sign to log_odds!)
    let mut psi_a = f64::NAN;
    let mut psi_b = f64::NAN;
    let mut delta_psi = f64::NAN;

    if design.is_two_group && df == 1 {
        let t_pos = p_full - 1;
        let mut sum_mu_a = 0.0;
        let mut sum_mu_b = 0.0;
        let mut n_ab = 0usize;
        for (j, &s_orig) in ev.valid_idx.iter().enumerate() {
            if !design.group_a_mask[s_orig] && !design.group_b_mask[s_orig] {
                continue;
            }
            let mut eta_base = 0.0;
            for (c, &beta_c) in fit_full.beta.iter().enumerate().take(p0) {
                eta_base += x_full.get(j, c) * beta_c;
            }
            sum_mu_b += logistic(eta_base);
            sum_mu_a += logistic(eta_base + fit_full.beta[t_pos]);
            n_ab += 1;
        }
        if n_ab > 0 {
            psi_a = sum_mu_a / (n_ab as f64);
            psi_b = sum_mu_b / (n_ab as f64);
            delta_psi = psi_a - psi_b;
        }
    } else {
        let sk: f64 = k.iter().sum();
        let sn: f64 = n.iter().sum();
        psi_a = if sn > 0.0 { sk / sn } else { f64::NAN };
        if design.is_continuous_test && df == 1 && log_odds.is_finite() {
            // Average marginal effect d(Psi)/dx = mean_j(mu_j * (1 - mu_j) * beta)
            let mut sum_slope = 0.0;
            for j in 0..m {
                let mut eta = 0.0;
                for (c, &beta_c) in fit_full.beta.iter().enumerate().take(p_full) {
                    eta += x_full.get(j, c) * beta_c;
                }
                let mu = logistic(eta);
                sum_slope += mu * (1.0 - mu) * log_odds;
            }
            delta_psi = sum_slope / (m as f64);
        }
    }

    GlmResult {
        tested: true,
        psi_a,
        psi_b,
        delta_psi,
        log_odds,
        se,
        dispersion: phi,
        lrt,
        df,
        pvalue,
    }
}

fn test_single_event(
    ev: &GlmEventData,
    design: &DesignSpec,
    theta_use: f64,
    df_denom: f64,
    ridge: f64,
) -> GlmResult {
    let Some((x_full, x_null, p0, _keep_test)) =
        build_event_matrices(&design.x, &ev.valid_idx, &design.test_cols)
    else {
        return GlmResult::empty();
    };
    test_single_event_with_matrices(ev, design, &x_full, &x_null, p0, theta_use, df_denom, ridge)
}

/// Run the full 3-pass Cox-Reid Empirical Bayes Beta-Binomial GLM across all cohort events in parallel.
pub fn fit_glm_cohort(
    events: &[GlmEventData],
    design: &DesignSpec,
    opts: &GlmOptions,
) -> (Vec<GlmResult>, Vec<f64>, EbPriorSummary) {
    let nev = events.len();
    if nev == 0 {
        return (
            Vec::new(),
            Vec::new(),
            EbPriorSummary {
                prior_mu: f64::NAN,
                prior_var: f64::NAN,
                n_prior_events: 0,
                df_prior: 0.0,
                utr_prior: None,
            },
        );
    }

    // Fixed user-specified dispersion mode (ChiSquared LRT)
    if let Some(phi) = opts.fixed_phi {
        if phi > 0.0 && phi < 1.0 {
            let th_fixed = ((1.0 - phi) / phi).ln().clamp(THETA_MIN, THETA_MAX);
            let results: Vec<GlmResult> = events
                .par_iter()
                .map(|ev| test_single_event(ev, design, th_fixed, f64::INFINITY, opts.ridge))
                .collect();
            let pvals: Vec<f64> = results.iter().map(|r| r.pvalue).collect();
            let fdr = bh_fdr(&pvals);
            return (
                results,
                fdr,
                EbPriorSummary {
                    prior_mu: th_fixed,
                    prior_var: 0.0,
                    n_prior_events: 0,
                    df_prior: 0.0,
                    utr_prior: None,
                },
            );
        }
    }

    // Check whether any event has residual degrees of freedom (m > p)
    let has_replicates = events.iter().any(|ev| {
        build_event_matrices(&design.x, &ev.valid_idx, &design.test_cols)
            .is_some_and(|(xf, _, _, _)| ev.valid_idx.len() > xf.cols)
    });

    if !has_replicates {
        eprintln!(
            "warning: no residual degrees of freedom (unreplicated design); using default dispersion phi = 0.01"
        );
        let th_def = (99.0f64).ln();
        let results: Vec<GlmResult> = events
            .par_iter()
            .map(|ev| test_single_event(ev, design, th_def, f64::INFINITY, opts.ridge))
            .collect();
        let pvals: Vec<f64> = results.iter().map(|r| r.pvalue).collect();
        let fdr = bh_fdr(&pvals);
        return (
            results,
            fdr,
            EbPriorSummary {
                prior_mu: th_def,
                prior_var: 0.0,
                n_prior_events: 0,
                df_prior: 0.0,
                utr_prior: None,
            },
        );
    }

    // Pass 2: Parallel raw Cox-Reid dispersion estimation strictly on polymorphic events
    // (skipping constitutive/boundary events whose flat likelihood creates artificial phi ~ 0.64 peaks).
    let raw_est: Vec<(f64, bool, f64, usize, bool)> = events
        .par_iter()
        .map(|ev| {
            let poly = is_polymorphic_event(ev);
            if !poly {
                return (f64::NAN, false, f64::NAN, 0, false);
            }
            let Some((x_full, _, _, _)) =
                build_event_matrices(&design.x, &ev.valid_idx, &design.test_cols)
            else {
                return (f64::NAN, false, f64::NAN, 0, true);
            };
            let m = ev.valid_idx.len();
            if m <= x_full.cols {
                return (f64::NAN, false, f64::NAN, 0, true);
            }
            let d_res = m - x_full.cols;
            let n = &ev.eff_reads;
            let mut k = Vec::with_capacity(m);
            let mut r = Vec::with_capacity(m);
            for j in 0..m {
                let kj = ev.psi[j] * n[j];
                k.push(kj);
                r.push(n[j] - kj);
            }
            let (th, intr, sv) = estimate_raw_theta(&k, &r, n, &x_full, opts.ridge);
            (th, intr, sv, d_res, true)
        })
        .collect();

    let raw_thetas: Vec<f64> = raw_est.iter().map(|x| x.0).collect();
    let is_interior: Vec<bool> = raw_est.iter().map(|x| x.1).collect();
    let samp_vars: Vec<f64> = raw_est.iter().map(|x| x.2).collect();
    let res_dfs: Vec<usize> = raw_est.iter().map(|x| x.3).collect();
    let is_poly: Vec<bool> = raw_est.iter().map(|x| x.4).collect();

    let mk_stratum = |(mu, var, n_ev): (f64, f64, usize)| -> EbStratumPrior {
        let d0 = if opts.shrink_dispersion {
            (2.0 / var).clamp(1.0, 30.0)
        } else {
            0.0
        };
        EbStratumPrior {
            prior_mu: mu,
            prior_var: var,
            n_prior_events: n_ev,
            df_prior: d0,
        }
    };

    let combined = mk_stratum(fit_eb_dispersion_prior(
        &raw_thetas,
        &is_interior,
        &samp_vars,
        &res_dfs,
    ));

    let (spliced_prior, utr_prior_opt) = if opts.stratify_eb {
        let mut sp_th = Vec::new();
        let mut sp_int = Vec::new();
        let mut sp_sv = Vec::new();
        let mut sp_df = Vec::new();
        let mut n_sp_int = 0usize;

        let mut u_th = Vec::new();
        let mut u_int = Vec::new();
        let mut u_sv = Vec::new();
        let mut u_df = Vec::new();
        let mut n_u_int = 0usize;

        for (i, ev) in events.iter().enumerate() {
            if !is_poly[i] {
                continue;
            }
            if is_utr_event(ev) {
                u_th.push(raw_thetas[i]);
                u_int.push(is_interior[i]);
                u_sv.push(samp_vars[i]);
                u_df.push(res_dfs[i]);
                if is_interior[i] && raw_thetas[i].is_finite() {
                    n_u_int += 1;
                }
            } else {
                sp_th.push(raw_thetas[i]);
                sp_int.push(is_interior[i]);
                sp_sv.push(samp_vars[i]);
                sp_df.push(res_dfs[i]);
                if is_interior[i] && raw_thetas[i].is_finite() {
                    n_sp_int += 1;
                }
            }
        }

        let sp = if n_sp_int >= 20 {
            mk_stratum(fit_eb_dispersion_prior(&sp_th, &sp_int, &sp_sv, &sp_df))
        } else {
            combined
        };
        let utr = if n_u_int >= 20 {
            Some(mk_stratum(fit_eb_dispersion_prior(
                &u_th, &u_int, &u_sv, &u_df,
            )))
        } else {
            None
        };
        (sp, utr)
    } else {
        (combined, None)
    };

    // Pass 3: Parallel EB MAP dispersion shrinkage + finite-sample F(df, d_res + d_0) LRT
    let results: Vec<GlmResult> = events
        .par_iter()
        .enumerate()
        .map(|(i, ev)| {
            let Some((x_full, x_null, p0, _keep_test)) =
                build_event_matrices(&design.x, &ev.valid_idx, &design.test_cols)
            else {
                return GlmResult::empty();
            };
            let m = ev.valid_idx.len();
            let d_res = m.saturating_sub(x_full.cols);

            let stratum = if is_utr_event(ev) {
                utr_prior_opt.unwrap_or(spliced_prior)
            } else {
                spliced_prior
            };

            let mut th_use = stratum.prior_mu;
            if is_poly[i] && d_res > 0 {
                if opts.shrink_dispersion {
                    let n = &ev.eff_reads;
                    let mut k = Vec::with_capacity(m);
                    let mut r = Vec::with_capacity(m);
                    for j in 0..m {
                        let kj = ev.psi[j] * n[j];
                        k.push(kj);
                        r.push(n[j] - kj);
                    }
                    th_use = estimate_eb_theta(
                        &k,
                        &r,
                        n,
                        &x_full,
                        stratum.prior_mu,
                        stratum.prior_var,
                        opts.ridge,
                    );
                } else if raw_thetas[i].is_finite() {
                    th_use = raw_thetas[i];
                }
            }

            let df_denom = if d_res > 0 || stratum.df_prior > 0.0 {
                (d_res as f64) + stratum.df_prior
            } else {
                f64::INFINITY
            };
            test_single_event_with_matrices(
                ev, design, &x_full, &x_null, p0, th_use, df_denom, opts.ridge,
            )
        })
        .collect();

    let pvals: Vec<f64> = results.iter().map(|r| r.pvalue).collect();
    let fdr = bh_fdr(&pvals);

    (
        results,
        fdr,
        EbPriorSummary {
            prior_mu: spliced_prior.prior_mu,
            prior_var: spliced_prior.prior_var,
            n_prior_events: spliced_prior.n_prior_events,
            df_prior: spliced_prior.df_prior,
            utr_prior: utr_prior_opt,
        },
    )
}

#[inline]
fn fmt_num(x: f64, digits: usize) -> String {
    if x.is_finite() {
        format!("{x:.digits$}")
    } else {
        "NA".into()
    }
}

#[inline]
fn fmt_pval(x: f64) -> String {
    if !x.is_finite() {
        "NA".into()
    } else if x == 0.0 {
        "0.0".into()
    } else if x < 1e-4 {
        format!("{x:.4e}")
    } else {
        format!("{x:.6}")
    }
}

/// Write differential splicing results to `.diff.gz` using parallel RFC 1952 gzip members,
/// matching `Whippet2.jl/src/diff_glm.jl::output_diff_glm` column schema identically.
pub fn output_diff_glm(
    outfile: &Path,
    events: &[GlmEventData],
    results: &[GlmResult],
    fdr: &[f64],
    design: &DesignSpec,
) -> Result<()> {
    let header: &[u8] = if design.is_two_group {
        b"Gene\tNode\tCoord\tStrand\tType\tComplexity\tEntropy\tSampleCnt_A\tSampleCnt_B\tMeanReads_A\tMeanReads_B\tPsi_A\tPsi_B\tDeltaPsi\tLogOddsRatio\tSE\tDispersion\tLRT\tPvalue\tFDR\n"
    } else {
        b"Gene\tNode\tCoord\tStrand\tType\tComplexity\tEntropy\tSampleCnt\tMeanReads\tMeanPsi\tDeltaPsi\tLogOddsRatio\tSE\tDispersion\tLRT\tDF\tPvalue\tFDR\n"
    };

    let indices: Vec<usize> = (0..events.len()).collect();
    let chunk_size = 2048usize;
    let blocks: Vec<Vec<u8>> = indices
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut buf = String::with_capacity(chunk.len() * 160);
            for &i in chunk {
                let res = &results[i];
                if !res.tested {
                    continue;
                }
                let ev = &events[i];
                if design.is_two_group {
                    let _ = writeln!(
                        buf,
                        "{}\t{}\t{}\t{}\t{}\tK{}\t{:.3}\t{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        ev.gene,
                        ev.node,
                        ev.coord,
                        ev.strand,
                        ev.etype,
                        ev.complexity,
                        ev.entropy,
                        ev.samp_a,
                        ev.samp_b,
                        ev.mean_reads_a,
                        ev.mean_reads_b,
                        fmt_num(res.psi_a, 5),
                        fmt_num(res.psi_b, 5),
                        fmt_num(res.delta_psi, 5),
                        fmt_num(res.log_odds, 5),
                        fmt_num(res.se, 5),
                        fmt_pval(res.dispersion),
                        fmt_num(res.lrt, 4),
                        fmt_pval(res.pvalue),
                        fmt_pval(fdr[i]),
                    );
                } else {
                    let _ = writeln!(
                        buf,
                        "{}\t{}\t{}\t{}\t{}\tK{}\t{:.3}\t{}\t{:.1}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        ev.gene,
                        ev.node,
                        ev.coord,
                        ev.strand,
                        ev.etype,
                        ev.complexity,
                        ev.entropy,
                        ev.samp_a,
                        ev.mean_reads_a,
                        fmt_num(res.psi_a, 5),
                        fmt_num(res.delta_psi, 5),
                        fmt_num(res.log_odds, 5),
                        fmt_num(res.se, 5),
                        fmt_pval(res.dispersion),
                        fmt_num(res.lrt, 4),
                        res.df,
                        fmt_pval(res.pvalue),
                        fmt_pval(fdr[i]),
                    );
                }
            }
            buf.into_bytes()
        })
        .collect();

    let file = File::create(outfile).map_err(|source| LangurError::IoPath {
        path: outfile.to_path_buf(),
        source,
    })?;
    let mut out = BufWriter::with_capacity(2 * 1024 * 1024, file);

    let mut hdr_enc = GzEncoder::new(Vec::new(), Compression::fast());
    hdr_enc.write_all(header)?;
    out.write_all(&hdr_enc.finish()?)?;

    let compressed: Vec<std::io::Result<Vec<u8>>> = blocks
        .into_par_iter()
        .filter(|b| !b.is_empty())
        .map(|raw| {
            let mut enc = GzEncoder::new(Vec::with_capacity(raw.len() / 3), Compression::fast());
            enc.write_all(&raw)?;
            enc.finish()
        })
        .collect();
    for res in compressed {
        out.write_all(&res?)?;
    }
    out.flush()?;
    Ok(())
}

/// Return bundle for [`process_psi_files_glm`].
pub type CohortGlmOutput = (Vec<GlmEventData>, Vec<GlmResult>, Vec<f64>, EbPriorSummary);

/// High-level entry point called by `langur delta`.
pub fn process_psi_files_glm(
    outfile: &Path,
    design: &DesignSpec,
    opts: &GlmOptions,
) -> Result<CohortGlmOutput> {
    let events = read_psi_cohort(design, &opts.cohort)?;
    let (results, fdr, prior) = fit_glm_cohort(&events, design, opts);
    output_diff_glm(outfile, &events, &results, &fdr, design)?;
    Ok((events, results, fdr, prior))
}
