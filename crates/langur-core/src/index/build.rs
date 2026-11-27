//! End-to-end `langur index` builder pipeline:
//! GTF2.2 (+ optional BAM supplementation) + Genome FASTA -> `.lidx` + `.exons.tsv.gz` + `.paths.tsv.gz`.

use crate::annotation::bam::{supplement_refset_from_bam, BamSupplementOptions};
use crate::annotation::{load_gtf_file, GtfOptions, RefGene, RefSet};
use crate::error::{Result, LangurError};
use crate::graph::{build_splice_graph, BuiltSpliceGraph, GC_BINS, GC_WINDOW_LEN};
use crate::index::decoy::{
    extract_chrom_decoys, extract_external_decoy_seq, DecoyExtractConfig,
};
use crate::index::format::{
    write_index_file, AnnoJunctionRecord, ChromRecord, DecoyRecord, GeneRecord, IndexBuffers,
    IndexHeader, MmapIndex, NodeRecord, StringPoolBuilder, TranscriptRecord, FLAG_HAS_BAM_NOVEL,
    FLAG_HAS_DECOYS, FLAG_HAS_JNC_SEEDS, FLAG_HASH_MAY_COLLIDE, LIDX_MAGIC, LIDX_VERSION,
};
use crate::index::minimizer::{
    build_minimizer_table, choose_bucket_bits, extract_csg_postings, RawPosting,
};
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecoyMode {
    /// Automatically extract selective non-CSG genomic decoys from the genome FASTA (default).
    SelectiveGenome,
    /// Use only explicitly supplied external decoy FASTA file(s).
    ExternalFasta(PathBuf),
    /// Disable decoy indexing entirely.
    None,
}

#[derive(Clone, Debug)]
pub struct IndexBuildConfig {
    /// Odd minimizer k-mer size (`13..=31`, default `17`; tests may use smaller odd `k >= 3`).
    pub minimizer_k: usize,
    /// Minimizer window size (`>= 1`, default `8`).
    pub minimizer_w: usize,
    /// Optional explicit bucket-bits override (`None` = auto-select preserving bijection for `k <= 27`).
    pub bucket_bits: Option<u16>,
    /// Maximum CSG postings stored per minimizer (default `256`).
    pub max_csg_occ: usize,
    /// Maximum decoy postings stored per minimizer (default `32`).
    pub max_decoy_occ: usize,
    /// Include `(k-1)+(k-1)` annotated junction-bridge minimizers (default `true`).
    pub include_jnc_seeds: bool,
    /// Decoy construction mode (default `SelectiveGenome`).
    pub decoy_mode: DecoyMode,
    /// Padding around CSG-colliding decoy hits in bp (default `250`).
    pub decoy_window_pad: usize,
    /// Exclude high-frequency CSG repeat minimizers (`occ > decoy_max_csg_occ`) from seeding
    /// selective genomic decoy windows (default `16`).
    pub decoy_max_csg_occ: usize,
    /// Minimum distinct CSG-colliding minimizers within `decoy_cluster_span_bp` to emit a
    /// selective genomic decoy window (default `2`).
    pub decoy_min_cluster_hits: usize,
    /// Maximum span (bp) between distinct CSG-colliding minimizers to form a decoy cluster (default `150`).
    pub decoy_cluster_span_bp: u32,
    /// GTF parsing/filtering options.
    pub gtf: GtfOptions,
    /// Optional indexed BAM file for de novo splice-site supplementation.
    pub bam_path: Option<PathBuf>,
    /// BAM supplementation options.
    pub bam_opts: BamSupplementOptions,
}

impl Default for IndexBuildConfig {
    fn default() -> Self {
        Self {
            minimizer_k: 17,
            minimizer_w: 8,
            bucket_bits: None,
            max_csg_occ: 256,
            max_decoy_occ: 32,
            include_jnc_seeds: true,
            decoy_mode: DecoyMode::SelectiveGenome,
            decoy_window_pad: 120,
            decoy_max_csg_occ: 16,
            decoy_min_cluster_hits: 4,
            decoy_cluster_span_bp: 80,
            gtf: GtfOptions::default(),
            bam_path: None,
            bam_opts: BamSupplementOptions::default(),
        }
    }
}

impl IndexBuildConfig {
    pub fn validate(&self) -> Result<()> {
        if self.minimizer_k < 3 || self.minimizer_k > 31 || self.minimizer_k % 2 == 0 {
            return Err(LangurError::InvalidParam(format!(
                "minimizer_k must be an odd integer in 3..=31 (got {})",
                self.minimizer_k
            )));
        }
        if self.minimizer_w == 0 || self.minimizer_w > 64 {
            return Err(LangurError::InvalidParam(format!(
                "minimizer_w must be in 1..=64 (got {})",
                self.minimizer_w
            )));
        }
        if self.max_csg_occ == 0 || self.max_csg_occ > (u16::MAX as usize) {
            return Err(LangurError::InvalidParam(format!(
                "max_csg_occ must be in 1..={} (got {})",
                u16::MAX,
                self.max_csg_occ
            )));
        }
        if self.max_decoy_occ > (u8::MAX as usize) {
            return Err(LangurError::InvalidParam(format!(
                "max_decoy_occ must be <= {} (got {})",
                u8::MAX,
                self.max_decoy_occ
            )));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct IndexBuildStats {
    pub num_chroms: usize,
    pub num_genes: usize,
    pub num_nodes: usize,
    pub num_transcripts: usize,
    pub num_anno_junctions: usize,
    pub num_decoys: usize,
    pub csg_seq_len: u64,
    pub decoy_seq_bases: u64,
    pub num_min_keys: u64,
    pub num_min_postings: u64,
    pub num_hash_collisions: u64,
    pub bam_novel_ss: usize,
    pub file_size_bytes: u64,
}

/// Build a `.lidx` index from genome FASTA + GTF2.2 (+ optional BAM) files on disk.
/// Uses a 2-pass FASTA stream so peak RAM never holds the entire genome simultaneously.
pub fn build_index_from_files(
    fasta_path: &Path,
    gtf_path: &Path,
    out_idx_path: &Path,
    cfg: &IndexBuildConfig,
) -> Result<IndexBuildStats> {
    cfg.validate()?;

    let mut refset = load_gtf_file(gtf_path, &cfg.gtf)?;
    let bam_novel_ss = if let Some(ref bam_path) = cfg.bam_path {
        supplement_refset_from_bam(&mut refset, bam_path, &cfg.bam_opts)?
    } else {
        0
    };

    // Group genes by chromosome name
    let mut chrom_genes: FxHashMap<String, Vec<&RefGene>> = FxHashMap::default();
    for g in &refset.genes {
        chrom_genes.entry(g.chrom.clone()).or_default().push(g);
    }

    // Pass 1 over genome FASTA: build CSGs per chromosome in parallel, drop chrom seq immediately.
    let mut chrom_list: Vec<(String, u32)> = Vec::new();
    let mut all_graphs: Vec<BuiltSpliceGraph> = Vec::new();
    let mut chrom_gene_ranges: Vec<(u32, u32)> = Vec::new();

    let mut fa_reader = parse_fastx_file(fasta_path)
        .map_err(|e| LangurError::Fasta(format!("failed to open FASTA {}: {e}", fasta_path.display())))?;

    while let Some(rec_res) = fa_reader.next() {
        let rec = rec_res.map_err(|e| LangurError::Fasta(format!("FASTA record error: {e}")))?;
        let chrom_name = extract_fasta_id(rec.id());
        let seq = rec.seq();
        let chrom_len = seq.len() as u32;

        let Some(genes_on_chr) = chrom_genes.get(&chrom_name) else {
            continue;
        };

        let built_res: Result<Vec<BuiltSpliceGraph>> = genes_on_chr
            .par_iter()
            .map(|gene| build_splice_graph(gene, seq.as_ref()))
            .collect();
        let mut built = built_res?;
        // Sort genes deterministically on each chromosome by (genomic_start, genomic_end, gene_id)
        built.sort_by(|a, b| {
            (a.genomic_start, a.genomic_end)
                .cmp(&(b.genomic_start, b.genomic_end))
                .then_with(|| a.gene_id.cmp(&b.gene_id))
        });

        let first_gene = all_graphs.len() as u32;
        let num_genes = built.len() as u32;
        chrom_list.push((chrom_name, chrom_len));
        chrom_gene_ranges.push((first_gene, num_genes));
        all_graphs.extend(built);
    }

    if all_graphs.is_empty() {
        return Err(LangurError::NoMatchingGenes);
    }

    // Estimate unique minimizers to choose optimal bucket_bits BEFORE extracting hashes
    let total_csg_bases: usize = all_graphs.iter().map(|g| g.seq.len()).sum();
    let est_unique = (total_csg_bases * 2) / (cfg.minimizer_w + 1);
    let bucket_bits = choose_bucket_bits(cfg.minimizer_k, est_unique, cfg.bucket_bits);
    let hash_bits = (bucket_bits as u32) + 32;

    // Extract CSG (+ annotated junction bridge) minimizer postings
    let csg_postings = extract_csg_postings(
        &all_graphs,
        cfg.minimizer_k,
        cfg.minimizer_w,
        hash_bits,
        cfg.include_jnc_seeds,
    );

    // Build selective decoys (Pass 2 over FASTA or external FASTA)
    let mut decoys: Vec<DecoyRecord> = Vec::new();
    let mut decoy_seq_2bit: Vec<u8> = Vec::new();
    let mut decoy_postings: Vec<RawPosting> = Vec::new();

    if cfg.decoy_mode != DecoyMode::None {
        let csg_hash_set = build_selective_csg_seed_set(&csg_postings, cfg.decoy_max_csg_occ);
        let decoy_cfg = DecoyExtractConfig {
            decoy_window_pad: cfg.decoy_window_pad,
            min_decoy_len: cfg.minimizer_k,
            min_cluster_hits: cfg.decoy_min_cluster_hits,
            cluster_span_bp: cfg.decoy_cluster_span_bp,
        };

        match &cfg.decoy_mode {
            DecoyMode::SelectiveGenome => {
                let chrom_to_idx: FxHashMap<&str, usize> = chrom_list
                    .iter()
                    .enumerate()
                    .map(|(i, (name, _))| (name.as_str(), i))
                    .collect();

                let mut fa_reader2 = parse_fastx_file(fasta_path).map_err(|e| {
                    LangurError::Fasta(format!(
                        "failed to reopen FASTA {} for decoy pass: {e}",
                        fasta_path.display()
                    ))
                })?;

                while let Some(rec_res) = fa_reader2.next() {
                    let rec = rec_res
                        .map_err(|e| LangurError::Fasta(format!("FASTA record error: {e}")))?;
                    let chrom_name = extract_fasta_id(rec.id());
                    let seq = rec.seq();

                    if let Some(&cidx) = chrom_to_idx.get(chrom_name.as_str()) {
                        let (fg, ng) = chrom_gene_ranges[cidx];
                        let graphs_ref: Vec<&BuiltSpliceGraph> = all_graphs
                            [fg as usize..(fg + ng) as usize]
                            .iter()
                            .collect();
                        extract_chrom_decoys(
                            cidx as u32,
                            seq.as_ref(),
                            &graphs_ref,
                            &csg_hash_set,
                            cfg.minimizer_k,
                            cfg.minimizer_w,
                            hash_bits,
                            &decoy_cfg,
                            &mut decoys,
                            &mut decoy_seq_2bit,
                            &mut decoy_postings,
                        );
                    } else {
                        // Unannotated chromosome/contig (e.g. decoy/alt/unplaced contig) -> scan fully!
                        extract_chrom_decoys(
                            u32::MAX,
                            seq.as_ref(),
                            &[],
                            &csg_hash_set,
                            cfg.minimizer_k,
                            cfg.minimizer_w,
                            hash_bits,
                            &decoy_cfg,
                            &mut decoys,
                            &mut decoy_seq_2bit,
                            &mut decoy_postings,
                        );
                    }
                }
            }
            DecoyMode::ExternalFasta(decoy_fa) => {
                let mut d_reader = parse_fastx_file(decoy_fa).map_err(|e| {
                    LangurError::Fasta(format!(
                        "failed to open decoy FASTA {}: {e}",
                        decoy_fa.display()
                    ))
                })?;
                while let Some(rec_res) = d_reader.next() {
                    let rec = rec_res.map_err(|e| {
                        LangurError::Fasta(format!("decoy FASTA record error: {e}"))
                    })?;
                    extract_external_decoy_seq(
                        rec.seq().as_ref(),
                        &csg_hash_set,
                        cfg.minimizer_k,
                        cfg.minimizer_w,
                        hash_bits,
                        &mut decoys,
                        &mut decoy_seq_2bit,
                        &mut decoy_postings,
                    );
                }
            }
            DecoyMode::None => {}
        }
    }

    finalize_and_write_index(
        out_idx_path,
        &chrom_list,
        &chrom_gene_ranges,
        &all_graphs,
        csg_postings,
        decoys,
        decoy_seq_2bit,
        decoy_postings,
        bucket_bits,
        bam_novel_ss,
        cfg,
    )
}

/// Collect CSG minimizer hashes occurring at most `max_occ` times in the transcriptome
/// so high-frequency repeats (Alu/LINE/poly-A) never seed false genomic decoy windows.
fn build_selective_csg_seed_set(csg_postings: &[RawPosting], max_occ: usize) -> FxHashSet<u64> {
    let mut freq: FxHashMap<u64, u32> =
        FxHashMap::with_capacity_and_hasher(csg_postings.len() / 2, Default::default());
    for p in csg_postings {
        let e = freq.entry(p.hash).or_insert(0);
        *e = e.saturating_add(1);
    }
    let cap = max_occ.min(u32::MAX as usize) as u32;
    freq.into_iter()
        .filter_map(|(h, c)| if c <= cap { Some(h) } else { None })
        .collect()
}

/// Build a `.lidx` index directly from in-memory chromosomes `(chrom_name, seq)` and a `RefSet`.
/// Ideal for unit/integration tests and programmatic API callers.
pub fn build_index_from_memory(
    chroms: &[(&str, &[u8])],
    refset: &RefSet,
    out_idx_path: &Path,
    cfg: &IndexBuildConfig,
) -> Result<IndexBuildStats> {
    cfg.validate()?;

    let mut chrom_genes: FxHashMap<&str, Vec<&RefGene>> = FxHashMap::default();
    for g in &refset.genes {
        chrom_genes.entry(g.chrom.as_str()).or_default().push(g);
    }

    let mut chrom_list: Vec<(String, u32)> = Vec::new();
    let mut chrom_gene_ranges: Vec<(u32, u32)> = Vec::new();
    let mut all_graphs: Vec<BuiltSpliceGraph> = Vec::new();

    for &(cname, cseq) in chroms {
        let Some(genes_on_chr) = chrom_genes.get(cname) else {
            continue;
        };
        let mut built: Vec<BuiltSpliceGraph> = genes_on_chr
            .iter()
            .map(|g| build_splice_graph(g, cseq))
            .collect::<Result<Vec<_>>>()?;
        built.sort_by(|a, b| {
            (a.genomic_start, a.genomic_end)
                .cmp(&(b.genomic_start, b.genomic_end))
                .then_with(|| a.gene_id.cmp(&b.gene_id))
        });

        let fg = all_graphs.len() as u32;
        let ng = built.len() as u32;
        chrom_list.push((cname.to_owned(), cseq.len() as u32));
        chrom_gene_ranges.push((fg, ng));
        all_graphs.extend(built);
    }

    if all_graphs.is_empty() {
        return Err(LangurError::NoMatchingGenes);
    }

    let total_csg_bases: usize = all_graphs.iter().map(|g| g.seq.len()).sum();
    let est_unique = (total_csg_bases * 2) / (cfg.minimizer_w + 1);
    let bucket_bits = choose_bucket_bits(cfg.minimizer_k, est_unique, cfg.bucket_bits);
    let hash_bits = (bucket_bits as u32) + 32;

    let csg_postings = extract_csg_postings(
        &all_graphs,
        cfg.minimizer_k,
        cfg.minimizer_w,
        hash_bits,
        cfg.include_jnc_seeds,
    );

    let mut decoys: Vec<DecoyRecord> = Vec::new();
    let mut decoy_seq_2bit: Vec<u8> = Vec::new();
    let mut decoy_postings: Vec<RawPosting> = Vec::new();

    if cfg.decoy_mode == DecoyMode::SelectiveGenome {
        let csg_hash_set = build_selective_csg_seed_set(&csg_postings, cfg.decoy_max_csg_occ);
        let decoy_cfg = DecoyExtractConfig {
            decoy_window_pad: cfg.decoy_window_pad,
            min_decoy_len: cfg.minimizer_k,
            min_cluster_hits: cfg.decoy_min_cluster_hits,
            cluster_span_bp: cfg.decoy_cluster_span_bp,
        };
        let chrom_to_idx: FxHashMap<&str, usize> = chrom_list
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.as_str(), i))
            .collect();

        for &(cname, cseq) in chroms {
            if let Some(&cidx) = chrom_to_idx.get(cname) {
                let (fg, ng) = chrom_gene_ranges[cidx];
                let graphs_ref: Vec<&BuiltSpliceGraph> =
                    all_graphs[fg as usize..(fg + ng) as usize].iter().collect();
                extract_chrom_decoys(
                    cidx as u32,
                    cseq,
                    &graphs_ref,
                    &csg_hash_set,
                    cfg.minimizer_k,
                    cfg.minimizer_w,
                    hash_bits,
                    &decoy_cfg,
                    &mut decoys,
                    &mut decoy_seq_2bit,
                    &mut decoy_postings,
                );
            } else {
                extract_chrom_decoys(
                    u32::MAX,
                    cseq,
                    &[],
                    &csg_hash_set,
                    cfg.minimizer_k,
                    cfg.minimizer_w,
                    hash_bits,
                    &decoy_cfg,
                    &mut decoys,
                    &mut decoy_seq_2bit,
                    &mut decoy_postings,
                );
            }
        }
    }

    finalize_and_write_index(
        out_idx_path,
        &chrom_list,
        &chrom_gene_ranges,
        &all_graphs,
        csg_postings,
        decoys,
        decoy_seq_2bit,
        decoy_postings,
        bucket_bits,
        0,
        cfg,
    )
}

#[allow(clippy::too_many_arguments)]
fn finalize_and_write_index(
    out_idx_path: &Path,
    chrom_list: &[(String, u32)],
    chrom_gene_ranges: &[(u32, u32)],
    all_graphs: &[BuiltSpliceGraph],
    csg_postings: Vec<RawPosting>,
    decoys: Vec<DecoyRecord>,
    decoy_seq_2bit: Vec<u8>,
    decoy_postings: Vec<RawPosting>,
    bucket_bits: u16,
    bam_novel_ss: usize,
    cfg: &IndexBuildConfig,
) -> Result<IndexBuildStats> {
    let min_table = build_minimizer_table(
        csg_postings,
        decoy_postings,
        cfg.minimizer_k,
        bucket_bits,
        cfg.max_csg_occ,
        cfg.max_decoy_occ,
    );

    let mut str_pool = StringPoolBuilder::default();
    let mut chrom_records = Vec::with_capacity(chrom_list.len());
    for ((cname, clen), &(fg, ng)) in chrom_list.iter().zip(chrom_gene_ranges.iter()) {
        chrom_records.push(ChromRecord {
            name: str_pool.push(cname),
            length: *clen,
            first_gene: fg,
            num_genes: ng,
            _pad: 0,
        });
    }

    let mut genes = Vec::with_capacity(all_graphs.len());
    let mut gene_prefix_max_end = vec![0u32; all_graphs.len()];
    let mut nodes: Vec<NodeRecord> = Vec::new();
    let mut anno_junctions: Vec<AnnoJunctionRecord> = Vec::new();
    let mut transcripts: Vec<TranscriptRecord> = Vec::new();
    let mut tx_nodes: Vec<u16> = Vec::new();
    let mut tx_node_prefix_len: Vec<u32> = Vec::new();
    let mut tx_expected_gc: Vec<f32> = Vec::new();
    let mut csg_seq: Vec<u8> = Vec::new();

    for (cidx, &(fg, ng)) in chrom_gene_ranges.iter().enumerate() {
        let mut running_max_end = 0u32;
        for gidx in (fg as usize)..((fg + ng) as usize) {
            let sg = &all_graphs[gidx];
            if sg.genomic_end > running_max_end {
                running_max_end = sg.genomic_end;
            }
            gene_prefix_max_end[gidx] = running_max_end;

            let seq_offset = csg_seq.len() as u64;
            csg_seq.extend_from_slice(&sg.seq);

            let node_offset = nodes.len() as u32;
            for nd in &sg.nodes {
                let edges_packed = nd.left_edge.as_u8() | (nd.right_edge.as_u8() << 4);
                let splice_meta = (nd.acc_dinuc as u8)
                    | ((nd.don_dinuc as u8) << 2)
                    | ((nd.is_novel_acc as u8) << 4)
                    | ((nd.is_novel_don as u8) << 5)
                    | ((nd.has_ambig as u8) << 6);
                nodes.push(NodeRecord {
                    gene_seq_offset: nd.gene_seq_offset,
                    genomic_coord: nd.genomic_coord,
                    len: nd.len,
                    acceptor_flank_2bit: nd.acceptor_flank_2bit,
                    donor_flank_2bit: nd.donor_flank_2bit,
                    edges_packed,
                    splice_meta,
                    acceptor_flank_len: nd.acceptor_flank_len,
                    donor_flank_len: nd.donor_flank_len,
                });
            }

            let anno_jnc_offset = anno_junctions.len() as u32;
            for j in &sg.anno_junctions {
                anno_junctions.push(AnnoJunctionRecord {
                    src_node: j.src_node,
                    dst_node: j.dst_node,
                    flags: j.flags,
                    tx_support: j.tx_support,
                });
            }

            let tx_offset = transcripts.len() as u32;
            for tx in &sg.transcripts {
                let tx_node_off = tx_nodes.len() as u32;
                tx_nodes.extend_from_slice(&tx.nodes);
                tx_node_prefix_len.extend_from_slice(&tx.node_prefix_len);
                tx_expected_gc.extend_from_slice(&tx.expected_gc);

                let tx_flags = (tx.is_basic as u32) | (((tx.tsl as u32) & 0x07) << 1);
                transcripts.push(TranscriptRecord {
                    tx_id: str_pool.push(&tx.tx_id),
                    gene_id: gidx as u32,
                    tx_node_offset: tx_node_off,
                    num_nodes: tx.nodes.len() as u32,
                    seq_len: tx.seq_len,
                    genomic_start: tx.tx_start,
                    genomic_end: tx.tx_end,
                    flags: tx_flags,
                });
            }

            let mut g_flags = 0u32;
            if sg.strand {
                g_flags |= GeneRecord::FLAG_POS_STRAND;
            }
            if sg.has_ambig {
                g_flags |= GeneRecord::FLAG_HAS_AMBIG;
            }
            if sg.has_bam_novel {
                g_flags |= GeneRecord::FLAG_HAS_BAM_NOVEL;
            }

            genes.push(GeneRecord {
                gene_id: str_pool.push(&sg.gene_id),
                gene_name: str_pool.push(&sg.gene_name),
                seq_offset,
                chrom_id: cidx as u32,
                genomic_start: sg.genomic_start,
                genomic_end: sg.genomic_end,
                seq_len: sg.seq.len() as u32,
                node_offset,
                num_nodes: sg.nodes.len() as u32,
                tx_offset,
                num_tx: sg.transcripts.len() as u32,
                anno_jnc_offset,
                num_anno_jnc: sg.anno_junctions.len() as u32,
                mean_tx_len: sg.mean_tx_len,
                flags: g_flags,
            });
        }
    }

    let decoy_seq_bases: u64 = decoys.iter().map(|d| d.seq_len as u64).sum();
    let mut header_flags = 0u32;
    if !decoys.is_empty() {
        header_flags |= FLAG_HAS_DECOYS;
    }
    if cfg.include_jnc_seeds {
        header_flags |= FLAG_HAS_JNC_SEEDS;
    }
    if bam_novel_ss > 0 {
        header_flags |= FLAG_HAS_BAM_NOVEL;
    }
    if min_table.hash_may_collide {
        header_flags |= FLAG_HASH_MAY_COLLIDE;
    }

    let header = IndexHeader {
        magic: LIDX_MAGIC,
        format_version: LIDX_VERSION,
        flags: header_flags,
        minimizer_k: cfg.minimizer_k as u16,
        minimizer_w: cfg.minimizer_w as u16,
        min_bucket_bits: min_table.bucket_bits,
        gc_bins: GC_BINS as u16,
        gc_window_len: GC_WINDOW_LEN as u16,
        max_csg_occ: cfg.max_csg_occ as u16,
        max_decoy_occ: cfg.max_decoy_occ as u16,
        section_count: 17,
        num_chroms: chrom_records.len() as u32,
        num_genes: genes.len() as u32,
        num_nodes: nodes.len() as u32,
        num_transcripts: transcripts.len() as u32,
        num_anno_junctions: anno_junctions.len() as u32,
        num_decoys: decoys.len() as u32,
        csg_seq_len: csg_seq.len() as u64,
        decoy_seq_bases,
        num_min_keys: min_table.key_lo.len() as u64,
        num_min_postings: min_table.postings.len() as u64,
        num_hash_collisions: min_table.num_hash_collisions,
        header_crc32: 0,
        _reserved: [0u32; 7],
    };

    let buffers = IndexBuffers {
        header,
        string_pool: str_pool.into_bytes(),
        chroms: chrom_records,
        genes,
        gene_prefix_max_end,
        nodes,
        anno_junctions,
        transcripts,
        tx_nodes,
        tx_node_prefix_len,
        tx_expected_gc,
        csg_seq,
        min_bucket_offsets: min_table.bucket_offsets,
        min_key_lo: min_table.key_lo,
        min_slots: min_table.slots,
        min_postings: min_table.postings,
        decoys,
        decoy_seq_2bit,
    };

    let file_size_bytes = write_index_file(out_idx_path, &buffers)?;

    Ok(IndexBuildStats {
        num_chroms: buffers.chroms.len(),
        num_genes: buffers.genes.len(),
        num_nodes: buffers.nodes.len(),
        num_transcripts: buffers.transcripts.len(),
        num_anno_junctions: buffers.anno_junctions.len(),
        num_decoys: buffers.decoys.len(),
        csg_seq_len: buffers.csg_seq.len() as u64,
        decoy_seq_bases,
        num_min_keys: buffers.min_key_lo.len() as u64,
        num_min_postings: buffers.min_postings.len() as u64,
        num_hash_collisions: buffers.header.num_hash_collisions,
        bam_novel_ss,
        file_size_bytes,
    })
}

fn extract_fasta_id(raw: &[u8]) -> String {
    let s = std::str::from_utf8(raw).unwrap_or("");
    s.split_whitespace().next().unwrap_or(s).to_owned()
}

/// Write Whippet2-compatible `<prefix>.exons.tsv.gz` directly from an open `MmapIndex`.
pub fn write_exons_tsv_gz(path: &Path, idx: &MmapIndex) -> Result<()> {
    let file = File::create(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut w = BufWriter::new(GzEncoder::new(file, Compression::fast()));
    writeln!(w, "Gene\tPotential_Exon\tWhippet_Nodes\tIs_Annotated")?;

    for (gidx, gene) in idx.genes().iter().enumerate() {
        let gname = idx.get_str(gene.gene_id);
        let chrom = idx.get_str(idx.chroms()[gene.chrom_id as usize].name);
        let strand_ch = if gene.is_pos_strand() { '+' } else { '-' };
        let nodes = idx.gene_nodes(gidx as u32);
        let txs = idx.gene_transcripts(gidx as u32);

        for i in 0..nodes.len() {
            if nodes[i].left_edge().is_exon_start() {
                for j in i..nodes.len() {
                    if nodes[j].right_edge().is_exon_stop() {
                        let is_anno = is_range_in_any_tx(idx, gene.tx_offset as usize, txs.len(), i as u16, j as u16);
                        let (coord_s, coord_e) = if gene.is_pos_strand() {
                            (nodes[i].genomic_coord, nodes[j].genomic_coord + nodes[j].len - 1)
                        } else {
                            (nodes[j].genomic_coord, nodes[i].genomic_coord + nodes[i].len - 1)
                        };
                        write!(w, "{gname}\t{chrom}:{coord_s}-{coord_e}:{strand_ch}\t")?;
                        for k in i..=j {
                            if k > i {
                                write!(w, ",")?;
                            }
                            write!(w, "{}", k + 1)?; // 1-based node numbering for Whippet2 parity
                        }
                        writeln!(w, "\t{}", if is_anno { "Y" } else { "N" })?;
                    }
                    if nodes[j].right_edge().is_hard_edge() {
                        break;
                    }
                }
            }
        }
    }
    w.flush()?;
    Ok(())
}

fn is_range_in_any_tx(idx: &MmapIndex, tx_off: usize, num_tx: usize, lo: u16, hi: u16) -> bool {
    for t in 0..num_tx {
        let tnodes = idx.transcript_nodes(tx_off + t);
        if let Ok(pos) = tnodes.binary_search(&lo) {
            let span = (hi - lo + 1) as usize;
            if pos + span <= tnodes.len()
                && (0..span).all(|d| tnodes[pos + d] == lo + (d as u16))
            {
                return true;
            }
        }
    }
    false
}

/// Write Whippet2-compatible `<prefix>.paths.tsv.gz` directly from an open `MmapIndex`.
pub fn write_paths_tsv_gz(path: &Path, idx: &MmapIndex) -> Result<()> {
    let file = File::create(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut w = BufWriter::new(GzEncoder::new(file, Compression::fast()));
    writeln!(w, "Gene\tIsoform\tPath\tEdges")?;

    for (gidx, gene) in idx.genes().iter().enumerate() {
        let gname = idx.get_str(gene.gene_id);
        let txs = idx.gene_transcripts(gidx as u32);
        for (rel, tx) in txs.iter().enumerate() {
            let tx_name = idx.get_str(tx.tx_id);
            let tnodes = idx.transcript_nodes(gene.tx_offset as usize + rel);
            write!(w, "{gname}\t{tx_name}\t")?;
            for (k, &n) in tnodes.iter().enumerate() {
                if k > 0 {
                    write!(w, "-")?;
                }
                write!(w, "{}", (n as u32) + 1)?;
            }
            write!(w, "\t")?;
            for (k, win) in tnodes.windows(2).enumerate() {
                if k > 0 {
                    write!(w, ",")?;
                }
                write!(w, "{}-{}", (win[0] as u32) + 1, (win[1] as u32) + 1)?;
            }
            writeln!(w)?;
        }
    }
    w.flush()?;
    Ok(())
}
