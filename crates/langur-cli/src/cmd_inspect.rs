use clap::Args;
use langur_core::index::{
    MmapIndex, FLAG_HAS_BAM_NOVEL, FLAG_HAS_DECOYS, FLAG_HAS_JNC_SEEDS, FLAG_HASH_MAY_COLLIDE,
    LEGACY_LIDX_MAGIC,
};
use langur_core::Result;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args)]
pub struct InspectArgs {
    /// Path to .lidx index file (also accepts legacy .w3idx indexes)
    #[arg(short = 'x', long = "index")]
    pub index: PathBuf,
}

pub fn run(args: InspectArgs) -> Result<()> {
    let t0 = Instant::now();
    let idx = MmapIndex::open(&args.index)?;
    let load_us = t0.elapsed().as_micros();
    let h = idx.header();
    let magic_str = std::str::from_utf8(&h.magic).unwrap_or("????????");
    let magic_note = if h.magic == LEGACY_LIDX_MAGIC {
        " (legacy pre-rename .w3idx)"
    } else {
        ""
    };

    println!("Langur Index: {}", args.index.display());
    println!("  Load time (mmap + CRC32): {:.3} ms", (load_us as f64) / 1000.0);
    println!(
        "  File size:                {:.2} MB ({} bytes)",
        (idx.file_size_bytes() as f64) / (1024.0 * 1024.0),
        idx.file_size_bytes()
    );
    println!("  Magic / Format version:   {magic_str}{magic_note}, v{}", h.format_version);
    println!(
        "  Minimizer (k, w):         ({}, {}) [bucket_bits={}, hash_bits={}, bijective={}]",
        h.minimizer_k,
        h.minimizer_w,
        h.min_bucket_bits,
        idx.hash_bits(),
        (h.flags & FLAG_HASH_MAY_COLLIDE) == 0
    );
    println!(
        "  Flags:                    decoys={}, jnc_seeds={}, bam_novel={}, may_collide={}",
        (h.flags & FLAG_HAS_DECOYS) != 0,
        (h.flags & FLAG_HAS_JNC_SEEDS) != 0,
        (h.flags & FLAG_HAS_BAM_NOVEL) != 0,
        (h.flags & FLAG_HASH_MAY_COLLIDE) != 0,
    );
    println!("  Chromosomes:              {}", h.num_chroms);
    println!("  Genes:                    {}", h.num_genes);
    println!("  CSG Nodes:                {}", h.num_nodes);
    println!("  Annotated Transcripts:    {}", h.num_transcripts);
    println!("  Annotated Junctions:      {}", h.num_anno_junctions);
    println!("  CSG Sequence:             {} bp", h.csg_seq_len);
    println!(
        "  Selective Decoys:         {} contigs ({} bp, {} bytes 2-bit)",
        h.num_decoys,
        h.decoy_seq_bases,
        idx.decoy_seq_2bit().len()
    );
    println!(
        "  Minimizer Table:          {} unique keys, {} postings, {} hash collisions",
        h.num_min_keys, h.num_min_postings, h.num_hash_collisions
    );
    Ok(())
}
