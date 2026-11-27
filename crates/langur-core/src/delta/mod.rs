//! Covariate-aware Beta-Binomial GLM differential splicing engine (`langur delta`):
//! 3-pass Cox-Reid Adjusted Profile Likelihood + Empirical Bayes dispersion shrinkage,
//! Smyth finite-df log-chi2 median bias correction, finite-sample F-LRT, BH FDR,
//! DESeq2-style formula parser, per-event pivoted-QR collinearity handling, and
//! inferential variance (`Inf_Var`) propagation.

pub mod cohort;
pub mod design;
pub mod glm;
pub mod math;

pub use cohort::{read_psi_cohort, CohortReadOptions, GlmEventData};
pub use design::{
    build_design_from_table, build_two_group_design, parse_formula_terms, resolve_glob_or_list,
    DesignSpec,
};
pub use glm::{
    fit_glm_cohort, output_diff_glm, process_psi_files_glm, EbPriorSummary, GlmOptions, GlmResult,
};
pub use math::{
    bb_cox_reid_apl, bb_loglik, bh_fdr, fit_bb_beta, logistic, logit, trigamma, Matrix,
};
