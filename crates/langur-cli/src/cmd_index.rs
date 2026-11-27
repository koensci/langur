use clap::Args;
use std::path::PathBuf;
use std::time::Instant;
use langur_core::annotation::{BamSupplementOptions, GtfOptions};
use langur_core::index::{
    build_index_from_files, write_exons_tsv_gz, write_paths_tsv_gz, DecoyMode, IndexBuildConfig,
    MmapIndex,
};
use langur_core::Result;

#[derive(Args)]
pub struct IndexArgs {
    /// Genome FASTA file (.fa, .fa.gz)
    #[arg(long = "fasta")]
    pub fasta: PathBuf,

    /// Gene annotation file in GTF2.2 format (.gtf, .gtf.gz)
    #[arg(long = "gtf")]
    pub gtf: PathBuf,

    /// Output index prefix (writes <prefix>.lidx, plus .exons.tsv.gz and .paths.tsv.gz)
    #[arg(short = 'x', long = "index", default_value = "langur_graph")]
    pub index: PathBuf,

    /// Odd minimizer k-mer size (13..=31; k <= 27 is mathematically collision-free at 54 bits)
    #[arg(short = 'k', long = "kmer", default_value_t = 17)]
    pub kmer: usize,

    /// Minimizer window size (exact-match seed guarantee = k + w - 1 bp)
    #[arg(short = 'w', long = "window", default_value_t = 8)]
    pub window: usize,

    /// Optional radix bucket bits override (10..=24; auto-selected by default)
    #[arg(long = "bucket-bits")]
    pub bucket_bits: Option<u16>,

    /// Maximum CSG postings retained per minimizer
    #[arg(long = "max-csg-occ", default_value_t = 256)]
    pub max_csg_occ: usize,

    /// Maximum decoy postings retained per minimizer
    #[arg(long = "max-decoy-occ", default_value_t = 32)]
    pub max_decoy_occ: usize,

    /// Disable annotated junction-bridge minimizer seeding
    #[arg(long = "no-jnc-seeds")]
    pub no_jnc_seeds: bool,

    /// Disable genomic decoy indexing
    #[arg(long = "no-decoys")]
    pub no_decoys: bool,

    /// Optional explicit external decoy FASTA file (instead of genome-derived selective decoys)
    #[arg(long = "decoys")]
    pub decoys: Option<PathBuf>,

    /// Flanking padding (bp) around CSG-colliding genomic decoy clusters
    #[arg(long = "decoy-pad", default_value_t = 120)]
    pub decoy_pad: usize,

    /// Restrict transcripts to GENCODE tag "basic" while using other transcripts for splice sites
    #[arg(long = "basic-or-ss-only")]
    pub basic_or_ss_only: bool,

    /// Ignore low-quality transcripts with TSL 3+
    #[arg(long = "suppress-low-tsl")]
    pub suppress_low_tsl: bool,

    /// Optional coordinate-sorted and indexed (.bai) BAM file for de novo splice-site supplementation
    #[arg(long = "bam")]
    pub bam: Option<PathBuf>,

    /// Minimum spliced read support required to incorporate a novel BAM splice site
    #[arg(long = "bam-min-reads", default_value_t = 1)]
    pub bam_min_reads: usize,

    /// Allow BAM splice junctions where neither donor nor acceptor is annotated
    #[arg(long = "bam-both-novel")]
    pub bam_both_novel: bool,

    /// Skip writing auxiliary .exons.tsv.gz and .paths.tsv.gz tables
    #[arg(long = "skip-aux-tsv")]
    pub skip_aux_tsv: bool,
}

pub fn run(args: IndexArgs) -> Result<()> {
    let t0 = Instant::now();

    let decoy_mode = if args.no_decoys {
        DecoyMode::None
    } else if let Some(ref dfa) = args.decoys {
        DecoyMode::ExternalFasta(dfa.clone())
    } else {
        DecoyMode::SelectiveGenome
    };

    let cfg = IndexBuildConfig {
        minimizer_k: args.kmer,
        minimizer_w: args.window,
        bucket_bits: args.bucket_bits,
        max_csg_occ: args.max_csg_occ,
        max_decoy_occ: args.max_decoy_occ,
        include_jnc_seeds: !args.no_jnc_seeds,
        decoy_mode,
        decoy_window_pad: args.decoy_pad,
        gtf: GtfOptions {
            basic_or_ss_only: args.basic_or_ss_only,
            suppress_low_tsl: args.suppress_low_tsl,
        },
        bam_path: args.bam,
        bam_opts: BamSupplementOptions {
            min_reads: args.bam_min_reads,
            allow_both_novel: args.bam_both_novel,
            unstranded_fallback: false,
        },
        ..Default::default()
    };

    let prefix_str = args.index.to_string_lossy();
    let idx_path = if prefix_str.ends_with(".lidx") {
        args.index.clone()
    } else {
        PathBuf::from(format!("{prefix_str}.lidx"))
    };
    let aux_prefix = prefix_str.strip_suffix(".lidx").unwrap_or(&prefix_str);

    eprintln!(
        "langur index: building {} (k={}, w={})...",
        idx_path.display(),
        cfg.minimizer_k,
        cfg.minimizer_w
    );

    let stats = build_index_from_files(&args.fasta, &args.gtf, &idx_path, &cfg)?;

    if !args.skip_aux_tsv {
        let mmap_idx = MmapIndex::open(&idx_path)?;
        let exons_path = PathBuf::from(format!("{aux_prefix}.exons.tsv.gz"));
        let paths_path = PathBuf::from(format!("{aux_prefix}.paths.tsv.gz"));
        write_exons_tsv_gz(&exons_path, &mmap_idx)?;
        write_paths_tsv_gz(&paths_path, &mmap_idx)?;
    }

    let elapsed = t0.elapsed().as_secs_f64();
    eprintln!(
        "Built {} in {:.2}s ({:.2} MB): {} chroms, {} genes, {} nodes, {} transcripts, {} anno junctions, {} decoys ({:.2} Mb), {} unique minimizers ({} postings, {} hash collisions).",
        idx_path.display(),
        elapsed,
        (stats.file_size_bytes as f64) / (1024.0 * 1024.0),
        stats.num_chroms,
        stats.num_genes,
        stats.num_nodes,
        stats.num_transcripts,
        stats.num_anno_junctions,
        stats.num_decoys,
        (stats.decoy_seq_bases as f64) / 1e6,
        stats.num_min_keys,
        stats.num_min_postings,
        stats.num_hash_collisions,
    );
    Ok(())
}
