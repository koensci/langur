use clap::Args;
use langur_core::align::{AlignConfig, StrandMode};
use langur_core::events::EventConfig;
use langur_core::index::MmapIndex;
use langur_core::io::{quantify_batch_manifest, quantify_sample, QuantConfig};
use langur_core::quant::VbemConfig;
use langur_core::{LangurError, Result};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args)]
pub struct QuantArgs {
    /// Input FASTQ file (R1 for paired-end, or single-end FASTQ; supports .fq / .fq.gz)
    pub r1: Option<PathBuf>,

    /// Optional mate-2 FASTQ file (R2 for paired-end)
    pub r2: Option<PathBuf>,

    /// Path to .lidx splice-graph index (also accepts legacy .w3idx indexes)
    #[arg(short = 'x', long = "index")]
    pub index: PathBuf,

    /// Output file prefix (writes .psi.gz, .gene.tpm.gz, .isoform.tpm.gz, .jnc.gz, .map.gz)
    #[arg(short = 'o', long = "out", default_value = "./output")]
    pub out: PathBuf,

    /// Batch manifest CSV (columns: sample,r1,r2,out_prefix or sample,r1,out_prefix; loads index once)
    #[arg(long = "batch")]
    pub batch: Option<PathBuf>,

    /// Forward-stranded library (FR: R1 sense, R2 antisense)
    #[arg(long = "stranded", conflicts_with = "stranded_reverse")]
    pub stranded: bool,

    /// Reverse-stranded library (RF / Illumina TruSeq dUTP: R1 antisense, R2 sense)
    #[arg(long = "stranded-reverse", conflicts_with = "stranded")]
    pub stranded_reverse: bool,

    /// Paired reads are on the same strand (fwd/fwd or rev/rev instead of fwd/rev)
    #[arg(long = "pair-same-strand")]
    pub pair_same_strand: bool,

    /// Enable unified 5'-hexamer + 20-bin fragment GC bias correction (feeds into both TPM and PSI)
    #[arg(long = "biascorrect")]
    pub biascorrect: bool,

    /// Enable circular back-splicing (`BS`) detection and quantification
    #[arg(long = "circ")]
    pub circ: bool,

    /// Stream SAM alignments to stdout
    #[arg(long = "sam")]
    pub sam: bool,

    /// Force gzip decompression regardless of file extension
    #[arg(long = "force-gz")]
    pub force_gz: bool,

    /// Sparse Dirichlet prior alpha_0 for Stage 1 transcript VBEM (default 0.01)
    #[arg(long = "vb-prior", default_value_t = 0.01)]
    pub vb_prior: f64,

    /// Use standard MLE EM instead of Variational Bayes EM in Stage 1
    #[arg(long = "use-em")]
    pub use_em: bool,

    /// Optional sparse Dirichlet prior for complex (K_n >= 2) local PSI events
    #[arg(long = "complex-vb-prior")]
    pub complex_vb_prior: Option<f64>,

    /// Number of targeted collapsed Gibbs posterior draws for complex (K_n >= 2) events (0 = disable)
    #[arg(long = "gibbs-samples", default_value_t = 100)]
    pub gibbs_samples: usize,

    /// Deterministic RNG seed for per-event complex Gibbs sampling
    #[arg(long = "seed", default_value_t = 42)]
    pub seed: u64,

    /// Disable MAJIQ-style single-position PCR-stack soft-capping on novel/RI junctions
    #[arg(long = "no-stack-filter")]
    pub no_stack_filter: bool,

    /// Effective read cap applied when single-position PCR-stack soft-capping triggers
    #[arg(long = "stack-cap", default_value_t = 2.0)]
    pub stack_cap: f64,

    /// Maximum transcriptomic fragment length for concordant mate pairs
    #[arg(long = "max-frag-len", default_value_t = 2000)]
    pub max_frag_len: u32,

    /// Minimum alignment score fraction of read_len * match_score
    #[arg(short = 'S', long = "score-min", default_value_t = 0.65)]
    pub score_min: f32,
}

pub fn run(args: QuantArgs) -> Result<()> {
    let t_load = Instant::now();
    let idx_path = {
        let s = args.index.to_string_lossy();
        if s.ends_with(".lidx") || s.ends_with(".w3idx") || args.index.exists() {
            args.index.clone()
        } else {
            PathBuf::from(format!("{s}.lidx"))
        }
    };
    let idx = MmapIndex::open(&idx_path)?;
    eprintln!(
        "langur quant: loaded {} in {:.2} ms ({} genes, {} transcripts, {} decoys)",
        idx_path.display(),
        t_load.elapsed().as_secs_f64() * 1000.0,
        idx.header().num_genes,
        idx.header().num_transcripts,
        idx.header().num_decoys,
    );

    let strand_mode = if args.stranded {
        StrandMode::ForwardFR
    } else if args.stranded_reverse {
        StrandMode::ReverseRF
    } else {
        StrandMode::Unstranded
    };

    let mut align_cfg = AlignConfig::default();
    align_cfg.strand_mode = strand_mode;
    align_cfg.allow_circ = args.circ;
    align_cfg.pair_is_rc = !args.pair_same_strand;
    align_cfg.max_frag_len = args.max_frag_len;
    align_cfg.min_score_frac = args.score_min;

    let vbem_cfg = VbemConfig {
        vb_prior: args.vb_prior,
        use_em: args.use_em,
        ..Default::default()
    };

    let event_cfg = EventConfig {
        allow_circ: args.circ,
        stack_filter: !args.no_stack_filter,
        stack_cap_reads: args.stack_cap,
        complex_vb_prior: args.complex_vb_prior,
        gibbs_samples: args.gibbs_samples,
        gibbs_seed: args.seed,
        ..Default::default()
    };

    let base_cfg = QuantConfig {
        r1: args.r1.clone().unwrap_or_default(),
        r2: args.r2.clone(),
        out_prefix: args.out.clone(),
        index_label: idx_path.display().to_string(),
        force_gz: args.force_gz,
        bias_correct: args.biascorrect,
        emit_sam: args.sam,
        align: align_cfg,
        vbem: vbem_cfg,
        events: event_cfg,
        ..Default::default()
    };

    if let Some(ref manifest) = args.batch {
        let t0 = Instant::now();
        let results = quantify_batch_manifest(&idx, manifest, &base_cfg)?;
        eprintln!(
            "langur quant --batch: quantified {} samples in {:.2}s",
            results.len(),
            t0.elapsed().as_secs_f64()
        );
        for (sample, s) in results {
            eprintln!(
                "  {sample}: {}/{} mapped ({:.1}%), {} decoy-filtered, {} events [align={:.2}s, vbem={:.2}s, events={:.2}s, write={:.2}s]",
                s.mapped_reads,
                s.total_reads,
                (s.mapped_reads as f64 / s.total_reads.max(1) as f64) * 100.0,
                s.decoy_filtered_reads,
                s.num_events,
                s.align_secs,
                s.vbem_secs,
                s.events_secs,
                s.write_secs,
            );
        }
        return Ok(());
    }

    if args.r1.is_none() {
        return Err(LangurError::InvalidParam(
            "must provide input FASTQ (<R1> [R2]) or --batch <manifest.csv>".into(),
        ));
    }

    let t0 = Instant::now();
    let s = quantify_sample(&idx, &base_cfg)?;
    eprintln!(
        "langur quant: finished in {:.2}s [align={:.2}s, vbem={:.2}s, events={:.2}s, write={:.2}s] — {}/{} mapped ({:.1}%), {} multimap, {} decoy-filtered, {} events (VBEM: {} comps, pass1={} iters, pass2={} iters)",
        t0.elapsed().as_secs_f64(),
        s.align_secs,
        s.vbem_secs,
        s.events_secs,
        s.write_secs,
        s.mapped_reads,
        s.total_reads,
        (s.mapped_reads as f64 / s.total_reads.max(1) as f64) * 100.0,
        s.multimap_reads,
        s.decoy_filtered_reads,
        s.num_events,
        s.num_vbem_components,
        s.vbem_pass1_iters,
        s.vbem_pass2_iters,
    );
    Ok(())
}
