use clap::Args;
use langur_core::delta::{
    build_design_from_table, build_two_group_design, process_psi_files_glm, resolve_glob_or_list,
    CohortReadOptions, GlmOptions,
};
use langur_core::{LangurError, Result};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args)]
pub struct DeltaArgs {
    /// Replicates for Group A — comma-delimited paths/basenames or pattern matching *<pat>*.psi* (2-group mode)
    #[arg(short = 'a', long = "a")]
    pub a: Option<String>,

    /// Replicates for Group B — same rules as -a (2-group mode; DeltaPsi = Psi_A - Psi_B)
    #[arg(short = 'b', long = "b")]
    pub b: Option<String>,

    /// Path to a CSV/TSV sample metadata table for covariate-aware GLM modeling
    #[arg(long = "design")]
    pub design: Option<PathBuf>,

    /// Column name in --design containing paths to .psi/.psi.gz files
    #[arg(long = "file-col", default_value = "file")]
    pub file_col: String,

    /// DESeq2-style model formula for --design mode (e.g. "~ batch + donor + condition" or "~ batch * condition")
    #[arg(short = 'f', long = "formula")]
    pub formula: Option<String>,

    /// Term or contrast to test (e.g. "condition", "condition:treated_vs_ctrl", or continuous column "age"; defaults to last term)
    #[arg(long = "test", default_value = "")]
    pub test: String,

    /// Output file prefix (writes <out>.diff.gz)
    #[arg(short = 'o', long = "out", default_value = "./output")]
    pub out: PathBuf,

    /// Directory used to resolve relative file paths or -a/-b patterns
    #[arg(short = 'd', long = "directory", default_value = ".")]
    pub directory: PathBuf,

    /// Minimum raw Total_Reads for an event in a sample to be included
    #[arg(short = 'r', long = "min-reads", default_value_t = 5)]
    pub min_reads: u32,

    /// Minimum valid samples (per group in 2-group/contrast mode, or total in general GLM mode)
    #[arg(short = 's', long = "min-samples", default_value_t = 1)]
    pub min_samples: usize,

    /// Disable Empirical Bayes dispersion shrinkage (use per-event Cox-Reid MLE)
    #[arg(long = "no-shrink")]
    pub no_shrink: bool,

    /// Disable Spliced vs. Tandem UTR (TS/TE) stratification of the Empirical Bayes dispersion prior
    #[arg(long = "no-stratify-eb")]
    pub no_stratify_eb: bool,

    /// Fixed Beta-Binomial overdispersion phi in (0, 1) (disables estimation; uses Chi-squared LRT)
    #[arg(long = "fixed-dispersion")]
    pub fixed_dispersion: Option<f64>,

    /// Ignore column 15 Inf_Var inferential variance scaling (use raw Total_Reads as in Whippet2)
    #[arg(long = "ignore-inf-var")]
    pub ignore_inf_var: bool,

    /// Include circular back-splicing (BS) events in differential analysis
    #[arg(long = "include-circ")]
    pub include_circ: bool,
}

pub fn run(args: DeltaArgs) -> Result<()> {
    let t0 = Instant::now();
    let has_design = args.design.is_some();
    let has_ab = args.a.is_some() || args.b.is_some();

    if has_design && has_ab {
        return Err(LangurError::InvalidParam(
            "specify either --design (with --formula) OR -a/-b, not both".into(),
        ));
    }
    if !has_design && !has_ab {
        return Err(LangurError::InvalidParam(
            "must specify either -a and -b (2-group mode) OR --design and --formula (covariate GLM)".into(),
        ));
    }

    let design = if let Some(ref design_path) = args.design {
        let formula = args.formula.as_deref().ok_or_else(|| {
            LangurError::InvalidParam(
                "--formula is required with --design (e.g. --formula '~ batch + condition')".into(),
            )
        })?;
        let d = build_design_from_table(
            design_path,
            formula,
            &args.test,
            &args.directory,
            &args.file_col,
        )?;
        let tcols: Vec<&str> = d.test_cols.iter().map(|&c| d.col_names[c].as_str()).collect();
        eprintln!(
            "langur delta: design '{}' ({} samples), formula '{}', testing '{}' (cols: {})",
            design_path.display(),
            d.files.len(),
            formula,
            d.test_label,
            tcols.join(", ")
        );
        d
    } else {
        let pa = args
            .a
            .as_deref()
            .ok_or_else(|| LangurError::InvalidParam("both -a and -b are required".into()))?;
        let pb = args
            .b
            .as_deref()
            .ok_or_else(|| LangurError::InvalidParam("both -a and -b are required".into()))?;
        let list_a = resolve_glob_or_list(pa, &args.directory)?;
        let list_b = resolve_glob_or_list(pb, &args.directory)?;
        if list_a.is_empty() || list_b.is_empty() {
            return Err(LangurError::InvalidParam(format!(
                "unable to match PSI files (matched {} for -a, {} for -b in {})",
                list_a.len(),
                list_b.len(),
                args.directory.display()
            )));
        }
        eprintln!(
            "langur delta: 2-group mode ({} samples in A vs {} samples in B)",
            list_a.len(),
            list_b.len()
        );
        build_two_group_design(list_a, list_b, "A", "B")?
    };

    let out_str = args.out.to_string_lossy();
    let outfile = if out_str.ends_with(".diff.gz") {
        args.out.clone()
    } else {
        PathBuf::from(format!("{out_str}.diff.gz"))
    };

    let opts = GlmOptions {
        cohort: CohortReadOptions {
            min_reads: args.min_reads,
            min_samples: args.min_samples,
            include_circ: args.include_circ,
            ignore_inf_var: args.ignore_inf_var,
        },
        shrink_dispersion: !args.no_shrink,
        stratify_eb: !args.no_stratify_eb,
        fixed_phi: args.fixed_dispersion,
        ridge: 1e-4,
    };

    let (events, results, fdr, prior) = process_psi_files_glm(&outfile, &design, &opts)?;
    let n_tested = results.iter().filter(|r| r.tested).count();
    let n_sig = results
        .iter()
        .zip(fdr.iter())
        .filter(|(r, &q)| r.tested && q.is_finite() && q <= 0.05)
        .count();

    if prior.prior_mu.is_finite() && prior.n_prior_events > 0 {
        let prior_phi = 1.0 / (1.0 + prior.prior_mu.exp());
        if let Some(u) = prior.utr_prior {
            let u_phi = 1.0 / (1.0 + u.prior_mu.exp());
            eprintln!(
                "  EB dispersion prior [spliced] ({} events): median phi = {:.4}, sd(log M) = {:.3}, d0 = {:.2}",
                prior.n_prior_events,
                prior_phi,
                prior.prior_var.sqrt(),
                prior.df_prior
            );
            eprintln!(
                "  EB dispersion prior [TS/TE]   ({} events): median phi = {:.4}, sd(log M) = {:.3}, d0 = {:.2}",
                u.n_prior_events,
                u_phi,
                u.prior_var.sqrt(),
                u.df_prior
            );
        } else {
            eprintln!(
                "  EB dispersion prior ({} events): median phi = {:.4}, sd(log M) = {:.3}, d0 = {:.2}",
                prior.n_prior_events,
                prior_phi,
                prior.prior_var.sqrt(),
                prior.df_prior
            );
        }
    }
    eprintln!(
        "langur delta: wrote {} in {:.2}s — {} loaded, {} tested, {} significant (FDR <= 0.05).",
        outfile.display(),
        t0.elapsed().as_secs_f64(),
        events.len(),
        n_tested,
        n_sig
    );
    Ok(())
}
