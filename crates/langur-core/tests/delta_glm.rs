//! Phase 3 validation suite for `langur delta` (Covariate-Aware Beta-Binomial GLM):
//! - Whippet2.jl parity (fractional Beta-Binomial MLE, complete separation 1.0 vs 0.0, BH FDR)
//! - Exact sign parity between DeltaPsi and LogOddsRatio in `-a`/`-b` and `--test var:A_vs_B`
//! - Simulated 400-event Beta-Binomial cohort Type-I error, EB dispersion recovery, and FDR calibration
//! - Batch-confounded design (`~ batch + condition`)
//! - Interactions (`~ batch * condition`), continuous covariates (`~ age`), CSV/TSV auto-delimiter,
//!   `--directory` resolution, `NA` sample dropping, and per-event pivoted-QR collinearity
//! - `Inf_Var` scientific notation & effective-depth scaling (with `--ignore-inf-var` escape hatch).

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use tempfile::tempdir;
use langur_core::delta::{
    bh_fdr, build_design_from_table, build_two_group_design, fit_bb_beta, fit_glm_cohort,
    logistic, parse_formula_terms, process_psi_files_glm, read_psi_cohort, CohortReadOptions,
    GlmEventData, GlmOptions, Matrix,
};
use langur_core::events::psi::beta_posterior_var;

/// Simple deterministic xorshift64 + Marsaglia-Tsang Gamma / Beta-Binomial simulator (no extra RNG deps).
struct SimpleRng(u64);

impl SimpleRng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn uniform(&mut self) -> f64 {
        let v = (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64);
        v.clamp(1e-12, 1.0 - 1e-12)
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn gamma(&mut self, alpha: f64) -> f64 {
        if alpha < 1.0 {
            let u = self.uniform();
            return self.gamma(alpha + 1.0) * u.powf(1.0 / alpha);
        }
        let d = alpha - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = 1.0 + c * x;
            if v <= 0.0 {
                continue;
            }
            let v3 = v * v * v;
            let u = self.uniform();
            if u < 1.0 - 0.0331 * x * x * x * x {
                return d * v3;
            }
            if u.ln() < 0.5 * x * x + d * (1.0 - v3 + v3.ln()) {
                return d * v3;
            }
        }
    }

    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        (x / (x + y)).clamp(0.0, 1.0)
    }

    fn binomial(&mut self, n: u32, p: f64) -> u32 {
        let mut k = 0u32;
        for _ in 0..n {
            if self.uniform() < p {
                k += 1;
            }
        }
        k
    }

    fn beta_binomial_psi(&mut self, n: u32, mu: f64, m_prec: f64) -> f64 {
        let p = self.beta(mu * m_prec, (1.0 - mu) * m_prec);
        (self.binomial(n, p) as f64) / (n as f64)
    }
}

fn write_psi_gz_fixture(path: &Path, rows: &[(&str, u32, &str, Option<f64>, f64, Option<f64>)]) {
    let f = std::fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::fast());
    writeln!(
        enc,
        "Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\tComplexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\tInf_Var"
    )
    .unwrap();
    for &(gene, node, etype, psi_opt, reads, iv_opt) in rows {
        match psi_opt {
            Some(psi) => {
                let iv = iv_opt.unwrap_or_else(|| beta_posterior_var(psi, reads));
                writeln!(
                    enc,
                    "{gene}\t{node}\tchr1:100-200\t+\t{etype}\t{psi:.5}\t0.1000\t0.4000,0.5000\t{reads:.1}\tK1\t1.000\t1-2-3:{psi:.4}\t1-3:{:.4}\t1-2:10.0\t{iv:.6e}",
                    1.0 - psi
                )
                .unwrap();
            }
            None => {
                writeln!(
                    enc,
                    "{gene}\t{node}\tchr1:100-200\t+\t{etype}\tNA\tNA\tNA\t0.0\tK0\t0.000\tNA\tNA\tNA\tNA"
                )
                .unwrap();
            }
        }
    }
    enc.finish().unwrap();
}

fn read_gz_to_string(path: &Path) -> String {
    let mut s = String::new();
    MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

#[test]
fn test_whippet2_glm_parity_mle_complete_separation_and_bh_fdr() {
    // 1. Exact test from Whippet2.jl/test/runtests.jl L778-788 (fractional counts)
    let k = vec![80.5, 76.2, 82.0, 20.1, 18.4, 21.5];
    let n = vec![100.0; 6];
    let r: Vec<f64> = n.iter().zip(k.iter()).map(|(ni, ki)| ni - ki).collect();
    let x = Matrix::from_row_major(
        6,
        2,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
    );
    let x0 = Matrix::from_row_major(6, 1, vec![1.0; 6]);

    let fit = fit_bb_beta(&k, &r, &n, &x, 50.0, None, 1e-4, 40, 1e-8);
    assert!(fit.converged);
    assert!((logistic(fit.beta[0]) - 0.20).abs() < 0.02);
    assert!((logistic(fit.beta[0] + fit.beta[1]) - 0.80).abs() < 0.02);
    let fit0 = fit_bb_beta(&k, &r, &n, &x0, 50.0, None, 1e-4, 40, 1e-8);
    let lrt = 2.0 * (fit.loglik - fit0.loglik);
    assert!(lrt > 30.0, "expected LRT > 30, got {lrt}");

    // 2. Complete separation (1.0 vs 0.0 PSI) converges cleanly without NaN/Inf
    let k_sep = vec![50.0, 40.0, 60.0, 0.0, 0.0, 0.0];
    let n_sep = vec![50.0, 40.0, 60.0, 50.0, 40.0, 60.0];
    let r_sep: Vec<f64> = n_sep.iter().zip(k_sep.iter()).map(|(a, b)| a - b).collect();
    let fit_sep = fit_bb_beta(&k_sep, &r_sep, &n_sep, &x, 20.0, None, 1e-4, 40, 1e-8);
    assert!(fit_sep.beta.iter().all(|v| v.is_finite()) && fit_sep.loglik.is_finite());
    assert!(logistic(fit_sep.beta[0] + fit_sep.beta[1]) > 0.98);
    assert!(logistic(fit_sep.beta[0]) < 0.02);

    // 3. Benjamini-Hochberg FDR preserves NaN cleanly
    let fdr = bh_fdr(&[0.001, 0.04, f64::NAN, 0.02, 0.5]);
    assert!((fdr[0] - 0.004).abs() < 1e-12);
    assert!(fdr[2].is_nan());
    assert!((fdr[3] - 0.04).abs() < 1e-12);
    assert!((fdr[1] - 0.053333333333).abs() < 1e-9);
    assert!((fdr[4] - 0.5).abs() < 1e-12);
}

#[test]
fn test_simulated_cohort_type1_error_and_fdr_calibration() {
    // 400 events across 6 Group-A vs 6 Group-B samples (true overdispersion phi = 0.04 -> M = 24):
    // - 360 true null events (same baseline PSI across A and B)
    // - 40 true differential events (DeltaPsi = ±0.35)
    let mut rng = SimpleRng(0xDEAD_BEEF_CAFE_1234);
    let m_true = 24.0; // phi = 1/25 = 0.04

    let files_a: Vec<_> = (0..6).map(|i| format!("a{i}.psi.gz").into()).collect();
    let files_b: Vec<_> = (0..6).map(|i| format!("b{i}.psi.gz").into()).collect();
    let design = build_two_group_design(files_a, files_b, "A", "B").unwrap();

    let mut events = Vec::with_capacity(400);
    for idx in 0..400usize {
        let is_alt = idx >= 360;
        let base_mu = 0.2 + 0.6 * rng.uniform();
        let (mu_a, mu_b) = if is_alt {
            if idx % 2 == 0 {
                (0.75, 0.35)
            } else {
                (0.25, 0.65)
            }
        } else {
            (base_mu, base_mu)
        };

        let mut psi = Vec::with_capacity(12);
        let mut eff_reads = Vec::with_capacity(12);
        for s in 0..12 {
            let n = 80u32;
            let mu = if s < 6 { mu_a } else { mu_b };
            psi.push(rng.beta_binomial_psi(n, mu, m_true));
            eff_reads.push(n as f64);
        }

        events.push(GlmEventData {
            gene: format!("G{idx}"),
            node: "2".into(),
            coord: "chr1:100-200".into(),
            strand: "+".into(),
            etype: "CE".into(),
            complexity: 1,
            entropy: 1.0,
            valid_idx: (0..12).collect(),
            psi,
            eff_reads,
            samp_a: 6,
            samp_b: 6,
            mean_reads_a: 80.0,
            mean_reads_b: 80.0,
        });
    }

    let (results, fdr, prior) = fit_glm_cohort(&events, &design, &GlmOptions::default());
    let est_phi = 1.0 / (1.0 + prior.prior_mu.exp());
    assert!(
        est_phi > 0.01 && est_phi < 0.12,
        "EB prior phi should recover ~0.04, got {est_phi:.4}"
    );

    // Null Type-I error at alpha = 0.05 across the 360 null events
    let null_fp = results[..360]
        .iter()
        .filter(|r| r.tested && r.pvalue < 0.05)
        .count();
    let fpr = (null_fp as f64) / 360.0;
    assert!(
        fpr <= 0.10,
        "null Type-I error must be controlled (got FPR={fpr:.3}, {null_fp}/360)"
    );

    // Power & realized empirical FDR at BH FDR <= 0.05
    let tp = (360..400)
        .filter(|&i| results[i].tested && fdr[i] <= 0.05)
        .count();
    let fp_fdr = (0..360)
        .filter(|&i| results[i].tested && fdr[i] <= 0.05)
        .count();
    let emp_fdr = (fp_fdr as f64) / ((tp + fp_fdr).max(1) as f64);
    assert!(tp >= 35, "expected high power (>=35/40 TP), got {tp}/40");
    assert!(
        emp_fdr <= 0.15,
        "realized FDR must be well-controlled (got {emp_fdr:.3}, fp={fp_fdr}, tp={tp})"
    );
}

#[test]
fn test_sign_parity_batch_confounding_interactions_continuous_na_and_collinearity() {
    let dir = tempdir().unwrap();

    // Create 8 samples across 2 batches x 2 conditions (balanced) + 1 extra NA-covariate sample
    // Event 1 ("G_COND_POS"): True positive condition shift (treated ~0.80, ctrl ~0.20), confounded with batch (+0.15 in B2)
    // Event 2 ("G_COND_NEG"): True negative condition shift (treated ~0.20, ctrl ~0.80)
    // Event 3 ("G_BATCH_ONLY"): Pure batch effect (B1 ~0.20, B2 ~0.80; zero condition effect!)
    // Event 4 ("G_COLLINEAR"): Missing in ALL batch B2 samples -> batch column becomes collinear in valid_idx, must still test condition cleanly!
    // Event 5 ("G_HIGH_DEPTH"): High depth (Total_Reads=8000, Psi=0.995 vs 0.995) testing Inf_Var scientific notation
    // Event 6 ("G_HIGH_INFVAR"): Elevated Inf_Var (5x nominal variance) testing effective-depth down-weighting
    let samples = [
        ("s1", "B1", "ctrl", 20.0, 0.18, 0.82, 0.20, Some(0.20)),
        ("s2", "B1", "ctrl", 25.0, 0.22, 0.78, 0.22, Some(0.22)),
        ("s3", "B2", "ctrl", 30.0, 0.33, 0.76, 0.80, None),
        ("s4", "B2", "ctrl", 35.0, 0.37, 0.80, 0.78, None),
        ("s5", "B1", "treated", 50.0, 0.78, 0.22, 0.19, Some(0.80)),
        ("s6", "B1", "treated", 55.0, 0.82, 0.18, 0.21, Some(0.82)),
        ("s7", "B2", "treated", 60.0, 0.88, 0.20, 0.82, None),
        ("s8", "B2", "treated", 65.0, 0.90, 0.24, 0.79, None),
        ("s9_na", "NA", "treated", 70.0, 0.50, 0.50, 0.50, Some(0.50)), // MUST be dropped due to batch=NA!
    ];

    let mut tsv = String::from("sample\tpsi_file\tbatch\tcondition\tage\n");
    for &(sid, batch, cond, age, p_pos, p_neg, p_batch, p_col) in &samples {
        let fname = format!("{sid}.psi.gz");
        let fpath = dir.path().join(&fname);
        let p_iv = if cond == "treated" { 0.65 } else { 0.45 };
        let v_nom_iv = beta_posterior_var(p_iv, 100.0);
        write_psi_gz_fixture(
            &fpath,
            &[
                ("G_COND_POS", 2, "CE", Some(p_pos), 100.0, None),
                ("G_COND_NEG", 2, "CE", Some(p_neg), 100.0, None),
                ("G_BATCH_ONLY", 2, "CE", Some(p_batch), 100.0, None),
                ("G_COLLINEAR", 2, "CE", p_col, 100.0, None),
                ("G_HIGH_DEPTH", 2, "CE", Some(0.995), 8000.0, None),
                ("G_HIGH_INFVAR", 2, "CE", Some(p_iv), 100.0, Some(v_nom_iv * 5.0)),
            ],
        );
        // Use bare basename to verify `--directory` resolution & TSV auto-delimiter!
        tsv.push_str(&format!("{sid}\t{fname}\t{batch}\t{cond}\t{age}\n"));
    }

    let design_path = dir.path().join("meta.tsv");
    std::fs::write(&design_path, &tsv).unwrap();

    // 1. Covariate GLM: `~ batch + condition` testing `condition:treated_vs_ctrl`
    let design = build_design_from_table(
        &design_path,
        "~ batch + condition",
        "condition:treated_vs_ctrl",
        dir.path(),
        "psi_file",
    )
    .unwrap();
    // s9_na must have been dropped -> exactly 8 samples, batch has 2 levels (NOT 3 with "NA"!)
    assert_eq!(design.files.len(), 8);
    assert_eq!(design.x.cols, 3); // (Intercept), batch:B2, condition:treated

    let out_diff = dir.path().join("cov.diff.gz");
    let opts = GlmOptions {
        fixed_phi: Some(0.02),
        ..Default::default()
    };
    let (events, results, _fdr, _) = process_psi_files_glm(&out_diff, &design, &opts).unwrap();
    assert_eq!(events.len(), 6);

    // Verify exact sign parity: DeltaPsi > 0 <=> LogOddsRatio > 0 (treated vs ctrl)
    let r_pos = &results[0];
    assert!(r_pos.tested);
    assert!(r_pos.delta_psi > 0.45 && r_pos.log_odds > 1.5);
    assert!(r_pos.pvalue < 1e-6);

    let r_neg = &results[1];
    assert!(r_neg.tested);
    assert!(r_neg.delta_psi < -0.45 && r_neg.log_odds < -1.5);
    assert!(r_neg.pvalue < 1e-6);

    // Pure batch effect MUST be suppressed (p > 0.2, |DeltaPsi| < 0.05)
    let r_batch = &results[2];
    assert!(r_batch.tested);
    assert!(
        r_batch.delta_psi.abs() < 0.05 && r_batch.pvalue > 0.2,
        "pure batch event leaked into condition test: dpsi={}, p={}",
        r_batch.delta_psi,
        r_batch.pvalue
    );

    // Per-event pivoted QR collinearity: G_COLLINEAR missing in all B2 samples still tests condition cleanly!
    let r_col = &results[3];
    assert!(
        r_col.tested && r_col.delta_psi > 0.50 && r_col.pvalue < 1e-4,
        "collinear nuisance batch should be dropped cleanly per event: {r_col:?}"
    );

    // High-depth event (Total_Reads=8000, Psi=0.995) must keep full eff_reads == 8000 (scientific Inf_Var never underflowed!)
    assert!(
        (events[4].eff_reads[0] - 8000.0).abs() < 1e-3,
        "high-depth scientific Inf_Var should not underflow n_eff (got {})",
        events[4].eff_reads[0]
    );

    // High-Inf_Var event (5x nominal variance) must scale eff_reads down to 100 * 0.2 = 20 and yield larger p-value than --ignore-inf-var!
    assert!(
        (events[5].eff_reads[0] - 20.0).abs() < 1e-3,
        "expected n_eff=20.0 for 5x Inf_Var, got {}",
        events[5].eff_reads[0]
    );
    let mut opts_ign = opts.clone();
    opts_ign.cohort.ignore_inf_var = true;
    let out_ign = dir.path().join("ign.diff.gz");
    let (_, res_ign, _, _) = process_psi_files_glm(&out_ign, &design, &opts_ign).unwrap();
    assert!(
        results[5].pvalue > res_ign[5].pvalue,
        "Inf_Var down-weighting must produce a more conservative p-value ({} vs {})",
        results[5].pvalue,
        res_ign[5].pvalue
    );

    // Verify `.diff.gz` written cleanly with 20 columns in 2-group mode
    let diff_txt = read_gz_to_string(&out_diff);
    let hdr_cols: Vec<&str> = diff_txt.lines().next().unwrap().split('\t').collect();
    assert_eq!(hdr_cols.len(), 20);
    assert_eq!(hdr_cols[13], "DeltaPsi");
    assert_eq!(hdr_cols[14], "LogOddsRatio");

    // 2. Classic `-a`/`-b` sign parity test (A = treated, B = ctrl -> DeltaPsi = Psi_A - Psi_B > 0 & LogOddsRatio > 0)
    let files_a: Vec<_> = ["s5", "s6", "s7", "s8"]
        .iter()
        .map(|s| dir.path().join(format!("{s}.psi.gz")))
        .collect();
    let files_b: Vec<_> = ["s1", "s2", "s3", "s4"]
        .iter()
        .map(|s| dir.path().join(format!("{s}.psi.gz")))
        .collect();
    let d_ab = build_two_group_design(files_a, files_b, "A", "B").unwrap();
    let out_ab = dir.path().join("ab.diff.gz");
    let (_, res_ab, _, _) = process_psi_files_glm(&out_ab, &d_ab, &opts).unwrap();
    assert!(res_ab[0].delta_psi > 0.0 && res_ab[0].log_odds > 0.0);
    assert!(res_ab[1].delta_psi < 0.0 && res_ab[1].log_odds < 0.0);

    // 3. Continuous numeric covariate (`~ age`) & interaction formula (`~ batch * condition`)
    assert_eq!(
        parse_formula_terms("~ batch * condition").unwrap(),
        vec![
            vec!["batch".to_string()],
            vec!["condition".to_string()],
            vec!["batch".to_string(), "condition".to_string()],
        ]
    );

    let d_age = build_design_from_table(&design_path, "~ age", "age", dir.path(), "psi_file").unwrap();
    assert!(d_age.is_continuous_test);
    assert!(!d_age.is_two_group);
    let out_age = dir.path().join("age.diff.gz");
    let (_, res_age, _, _) = process_psi_files_glm(&out_age, &d_age, &opts).unwrap();
    assert!(res_age[0].tested && res_age[0].log_odds > 0.0 && res_age[0].delta_psi > 0.0);
    let age_txt = read_gz_to_string(&out_age);
    assert_eq!(age_txt.lines().next().unwrap().split('\t').count(), 18); // general GLM schema with DF

    let d_int = build_design_from_table(
        &design_path,
        "~ batch * condition",
        "batch:condition",
        dir.path(),
        "psi_file",
    )
    .unwrap();
    assert_eq!(d_int.x.cols, 4); // (Intercept), batch:B2, condition:treated, batch:B2:condition:treated
    assert_eq!(d_int.test_cols, vec![3]);
}

fn write_gz(path: &Path, text: &str) {
    let f = std::fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::fast());
    enc.write_all(text.as_bytes()).unwrap();
    enc.finish().unwrap();
}

#[test]
fn test_sparse_circ_bs_events_and_truncated_psi_detection() {
    let dir = tempdir().unwrap();
    let hdr = "Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\tComplexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\tInf_Var\n";

    let row_ce = |gene: &str, node: u32, coord: &str, psi: f64, reads: u32| -> String {
        let iv = beta_posterior_var(psi, reads as f64);
        format!("{gene}\t{node}\t{coord}\t+\tCE\t{psi:.4}\t0.100\t0.0,1.0\t{reads}\tK1\t0.500\t1\t1\t1-2:1\t{iv:.6e}\n")
    };
    let row_bs = |gene: &str, node: &str, coord: &str, psi: f64, reads: u32| -> String {
        let iv = beta_posterior_var(psi, reads as f64);
        format!("{gene}\t{node}\t{coord}\t+\tBS\t{psi:.4}\t0.100\t0.0,1.0\t{reads}\tK1\t0.400\t1\t1\t2-1:1\t{iv:.6e}\n")
    };

    // Sample a1: G1 CE, BS_SHARED (interleaved), BS_ONLY_A (only in a1), G2 CE
    let a1_path = dir.path().join("a1.psi.gz");
    write_gz(
        &a1_path,
        &format!(
            "{hdr}{}{}{}{}",
            row_ce("G1", 2, "chr1:100-200", 0.50, 50),
            row_bs("G1", "3-2", "chr1:200-300", 0.80, 40),
            row_bs("G1", "4-2", "chr1:200-400", 0.60, 30),
            row_ce("G2", 2, "chr1:500-600", 0.50, 50),
        ),
    );
    // Sample a2: G1 CE, BS_SHARED, G2 CE, BS_TRAIL (trailing after last canonical row!)
    let a2_path = dir.path().join("a2.psi.gz");
    write_gz(
        &a2_path,
        &format!(
            "{hdr}{}{}{}{}",
            row_ce("G1", 2, "chr1:100-200", 0.52, 50),
            row_bs("G1", "3-2", "chr1:200-300", 0.78, 45),
            row_ce("G2", 2, "chr1:500-600", 0.48, 50),
            row_bs("G2", "3-2", "chr1:600-700", 0.70, 35),
        ),
    );
    // Sample b1: G1 CE, BS_SHARED (lower PSI), G2 CE
    let b1_path = dir.path().join("b1.psi.gz");
    write_gz(
        &b1_path,
        &format!(
            "{hdr}{}{}{}",
            row_ce("G1", 2, "chr1:100-200", 0.50, 50),
            row_bs("G1", "3-2", "chr1:200-300", 0.15, 40),
            row_ce("G2", 2, "chr1:500-600", 0.51, 50),
        ),
    );
    // Sample b2: G1 CE, G2 CE, BS_TRAIL (no BS_SHARED at all!)
    let b2_path = dir.path().join("b2.psi.gz");
    write_gz(
        &b2_path,
        &format!(
            "{hdr}{}{}{}",
            row_ce("G1", 2, "chr1:100-200", 0.49, 50),
            row_ce("G2", 2, "chr1:500-600", 0.50, 50),
            row_bs("G2", "3-2", "chr1:600-700", 0.20, 35),
        ),
    );

    let design = build_two_group_design(
        vec![a1_path.clone(), a2_path.clone()],
        vec![b1_path.clone(), b2_path.clone()],
        "A",
        "B",
    )
    .unwrap();

    // 1. Default (include_circ = false): sparse BS rows are skipped cleanly without mismatch error
    let opts_no_circ = CohortReadOptions {
        min_reads: 5,
        min_samples: 1,
        include_circ: false,
        ignore_inf_var: false,
    };
    let ev_no_circ = read_psi_cohort(&design, &opts_no_circ).unwrap();
    assert_eq!(ev_no_circ.len(), 2);
    assert!(ev_no_circ.iter().all(|e| e.etype == "CE"));

    // 2. With include_circ = true: canonical events stay lockstep, sparse BS events merge & sort cleanly
    let mut glm_opts = GlmOptions::default();
    glm_opts.cohort.include_circ = true;
    glm_opts.cohort.min_samples = 1;
    let out_diff = dir.path().join("circ.diff.gz");
    let (ev_circ, res_circ, _, _) = process_psi_files_glm(&out_diff, &design, &glm_opts).unwrap();
    assert_eq!(ev_circ.len(), 4);
    assert_eq!((ev_circ[0].gene.as_str(), ev_circ[0].etype.as_str()), ("G1", "CE"));
    assert_eq!((ev_circ[1].gene.as_str(), ev_circ[1].etype.as_str()), ("G2", "CE"));
    // Merged BS events appended in deterministic sorted order by (Gene, Node, Coord, Strand)
    assert_eq!((ev_circ[2].gene.as_str(), ev_circ[2].etype.as_str()), ("G1", "BS"));
    assert_eq!((ev_circ[2].samp_a, ev_circ[2].samp_b), (2, 1));
    assert_eq!(ev_circ[2].valid_idx, vec![0, 1, 2]);
    assert!(res_circ[2].tested && res_circ[2].delta_psi > 0.50);

    assert_eq!((ev_circ[3].gene.as_str(), ev_circ[3].etype.as_str()), ("G2", "BS"));
    assert_eq!((ev_circ[3].samp_a, ev_circ[3].samp_b), (1, 1));
    assert_eq!(ev_circ[3].valid_idx, vec![1, 3]);
    assert!(res_circ[3].tested && res_circ[3].delta_psi > 0.35);

    // 3. Truncated PSI detection: sample s > 0 missing trailing canonical event G2 CE
    let b2_trunc = dir.path().join("b2_trunc.psi.gz");
    write_gz(&b2_trunc, &format!("{hdr}{}", row_ce("G1", 2, "chr1:100-200", 0.49, 50)));
    let d_trunc_b = build_two_group_design(
        vec![a1_path.clone(), a2_path.clone()],
        vec![b1_path.clone(), b2_trunc],
        "A",
        "B",
    )
    .unwrap();
    let err_b = read_psi_cohort(&d_trunc_b, &glm_opts.cohort).unwrap_err();
    let msg_b = err_b.to_string();
    assert!(
        msg_b.contains("truncated") && msg_b.contains("b2_trunc"),
        "expected truncated file error naming b2_trunc, got: {msg_b}"
    );

    // 4. Truncated PSI detection: sample s == 0 missing trailing canonical event G2 CE
    let a1_trunc = dir.path().join("a1_trunc.psi.gz");
    write_gz(&a1_trunc, &format!("{hdr}{}", row_ce("G1", 2, "chr1:100-200", 0.50, 50)));
    let d_trunc_a = build_two_group_design(
        vec![a1_trunc, a2_path],
        vec![b1_path, b2_path],
        "A",
        "B",
    )
    .unwrap();
    let err_a = read_psi_cohort(&d_trunc_a, &glm_opts.cohort).unwrap_err();
    let msg_a = err_a.to_string();
    assert!(
        msg_a.contains("truncated") && msg_a.contains("a1_trunc"),
        "expected truncated file error naming a1_trunc, got: {msg_a}"
    );
}

#[test]
fn test_constitutive_majority_does_not_collapse_eb_prior_and_stratifies_utr() {
    // Regression test for human-scale E-MTAB-3766 finding:
    // Mix 80% constitutive (Psi = 1.0 / 0.0) exons with 20% polymorphic events
    // (10% spliced CE at low overdispersion phi = 0.015 + 10% tandem TE at high phi = 0.12)
    // and verify:
    //   1. Constitutive exons are excluded from EB prior estimation (never collapsing to phi ~ 0.57).
    //   2. Spliced vs. TS/TE stratified EB prior cleanly separates low spliced phi from high UTR phi.
    //   3. True DS CE events remain strongly significant at FDR <= 0.05.
    let mut rng = SimpleRng(20261126);
    let n_a = 3usize;
    let n_b = 3usize;
    let n_tot = n_a + n_b;

    let design = build_two_group_design(
        (0..n_a).map(|i| PathBuf::from(format!("a{i}.psi.gz"))).collect(),
        (0..n_b).map(|i| PathBuf::from(format!("b{i}.psi.gz"))).collect(),
        "A",
        "B",
    )
    .unwrap();

    let sample_bb_psi = |rng: &mut SimpleRng, mu: f64, phi: f64, n_reads: u32| -> f64 {
        let m_prec = (1.0 - phi) / phi;
        rng.beta_binomial_psi(n_reads, mu, m_prec)
    };

    let mut events: Vec<GlmEventData> = Vec::with_capacity(1000);
    // 800 constitutive events (80%): 400 all-1.0 CE + 400 all-0.0 TE
    for i in 0..800 {
        let is_one = i < 400;
        events.push(GlmEventData {
            gene: format!("CONST_{i}"),
            node: "1".into(),
            coord: format!("chr1:{}-{}", i * 100 + 1, i * 100 + 80),
            strand: "+".into(),
            etype: if is_one { "CE".into() } else { "TE".into() },
            complexity: 0,
            entropy: 0.0,
            valid_idx: (0..n_tot).collect(),
            psi: vec![if is_one { 1.0 } else { 0.0 }; n_tot],
            eff_reads: vec![150.0; n_tot],
            samp_a: n_a,
            samp_b: n_b,
            mean_reads_a: 150.0,
            mean_reads_b: 150.0,
        });
    }

    // 100 polymorphic spliced CE events at true phi = 0.015 (first 15 are true DS: 0.85 vs 0.25)
    for i in 0..100 {
        let is_ds = i < 15;
        let mu_a = if is_ds { 0.85 } else { 0.50 };
        let mu_b = if is_ds { 0.25 } else { 0.50 };
        let mut psi = Vec::with_capacity(n_tot);
        for j in 0..n_tot {
            let mu = if j < n_a { mu_a } else { mu_b };
            psi.push(sample_bb_psi(&mut rng, mu, 0.015, 160));
        }
        events.push(GlmEventData {
            gene: format!("POLY_CE_{i}"),
            node: "2".into(),
            coord: format!("chr2:{}-{}", i * 100 + 1, i * 100 + 80),
            strand: "+".into(),
            etype: "CE".into(),
            complexity: 1,
            entropy: 0.9,
            valid_idx: (0..n_tot).collect(),
            psi,
            eff_reads: vec![160.0; n_tot],
            samp_a: n_a,
            samp_b: n_b,
            mean_reads_a: 160.0,
            mean_reads_b: 160.0,
        });
    }

    // 100 polymorphic tandem UTR TE events at high true phi = 0.15 (all null: 0.50 vs 0.50)
    for i in 0..100 {
        let mut psi = Vec::with_capacity(n_tot);
        for _ in 0..n_tot {
            psi.push(sample_bb_psi(&mut rng, 0.50, 0.15, 160));
        }
        events.push(GlmEventData {
            gene: format!("POLY_TE_{i}"),
            node: "3".into(),
            coord: format!("chr3:{}-{}", i * 100 + 1, i * 100 + 80),
            strand: "+".into(),
            etype: "TE".into(),
            complexity: 2,
            entropy: 0.9,
            valid_idx: (0..n_tot).collect(),
            psi,
            eff_reads: vec![160.0; n_tot],
            samp_a: n_a,
            samp_b: n_b,
            mean_reads_a: 160.0,
            mean_reads_b: 160.0,
        });
    }

    let opts = GlmOptions::default();
    let (results, fdr, prior) = fit_glm_cohort(&events, &design, &opts);

    let sp_phi = 1.0 / (1.0 + prior.prior_mu.exp());
    let utr = prior
        .utr_prior
        .expect("stratified TS/TE EB prior should be fitted when >= 20 interior UTR events exist");
    let utr_phi = 1.0 / (1.0 + utr.prior_mu.exp());

    assert!(
        prior.n_prior_events <= 100 && utr.n_prior_events <= 100,
        "constitutive events must be excluded from EB prior (got sp={}, utr={})",
        prior.n_prior_events,
        utr.n_prior_events
    );
    assert!(
        sp_phi < 0.06,
        "spliced EB median phi collapsed instead of recovering ~0.015 (got {sp_phi:.4})"
    );
    assert!(
        utr_phi > sp_phi * 2.0,
        "stratified UTR phi ({utr_phi:.4}) should cleanly exceed spliced phi ({sp_phi:.4})"
    );

    // At least 13/15 true DS CE events must be detected at FDR <= 0.05!
    let tp_ds = (800..815)
        .filter(|&i| results[i].tested && fdr[i].is_finite() && fdr[i] <= 0.05)
        .count();
    assert!(
        tp_ds >= 13,
        "expected >= 13/15 true DS CE events at FDR <= 0.05 despite 80% constitutive background, got {tp_ds}"
    );
}
