//! Short-read spliced alignment engine: two-tier minimizer seeding, hybrid ungapped +
//! zero-allocation stack-banded Gotoh affine-gap X-drop extension, SNP-tolerant splice-edge
//! jumping, strict-inequality selective decoy validation, and paired-end
//! concordant matching with exact transcriptomic fragment lengths.

pub mod extend;
pub mod paired;
pub mod seed;

use smallvec::SmallVec;

/// Library strandedness protocol for read orientation.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum StrandMode {
    /// Unstranded library (tests both sense and antisense orientations).
    Unstranded,
    /// Forward-stranded (`--stranded`, FR: SE sense, PE R1 sense / R2 antisense).
    ForwardFR,
    /// Reverse-stranded (`--stranded-reverse`, RF / TruSeq dUTP: SE antisense, PE R1 antisense / R2 sense).
    ReverseRF,
}

/// Alignment and scoring configuration shared identically across CSG and decoy scoring.
#[derive(Clone, Debug)]
pub struct AlignConfig {
    pub strand_mode: StrandMode,
    /// Allow circular / back-splicing edges (`dst_node <= src_node`).
    pub allow_circ: bool,
    /// Match reward (default `+2`).
    pub match_score: i32,
    /// Mismatch penalty (default `-4`).
    pub mismatch_penalty: i32,
    /// Ambiguous `'N'` score (default `-1`).
    pub ambig_score: i32,
    /// Affine gap open penalty (default `-6`; 1-bp indel costs `gap_open`).
    pub gap_open: i32,
    /// Affine gap extend penalty (default `-1`).
    pub gap_extend: i32,
    /// Penalty for non-annotated canonical (`GT-AG`/`GC-AG`/`AT-AC`) splice jump (default `0`).
    /// Annotated junctions ALWAYS have `0` penalty so processed pseudogene decoys cannot tie-break against real mRNA.
    pub novel_canonical_penalty: i32,
    /// Penalty for non-annotated non-canonical splice jump (default `-2`).
    pub novel_noncanonical_penalty: i32,
    /// Minimum overhang (bp) across an annotated splice junction (default `4`).
    pub min_anno_overhang: u32,
    /// Minimum overhang (bp) across a novel/combinatorial splice junction (default `8`).
    pub min_novel_overhang: u32,
    /// Minimum overhang (bp) across a back-splice (`BS`) circular junction (default `10`).
    pub min_circ_overhang: u32,
    /// Minimum alignment score as a fraction of `read_len * match_score` (default `0.65`).
    pub min_score_frac: f32,
    /// Score tolerance fraction relative to best score for multi-mapping retention (default `0.05`).
    pub score_range_frac: f32,
    /// Maximum candidate `(gene_id, strand)` clusters evaluated per read (default `32`).
    pub max_gene_candidates: usize,
    /// Maximum candidate decoy contigs evaluated per read (default `8`).
    pub max_decoy_candidates: usize,
    /// Maximum transcriptomic fragment length for paired-end concordance (default `2000`).
    pub max_frag_len: u32,
    /// Whether mate 2 is reverse-complemented relative to mate 1 (default `true`; false with `--pair-same-strand`).
    pub pair_is_rc: bool,
    /// Maximum multi-mapping alignments retained per read/pair (default `32`).
    pub max_multi_hits: usize,
    /// Enable zero-allocation stack-banded Gotoh affine-gap X-drop fallback on small-indel clusters.
    pub enable_gapped: bool,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            strand_mode: StrandMode::Unstranded,
            allow_circ: false,
            match_score: 2,
            mismatch_penalty: -4,
            ambig_score: -1,
            gap_open: -6,
            gap_extend: -1,
            novel_canonical_penalty: 0,
            novel_noncanonical_penalty: -2,
            min_anno_overhang: 4,
            min_novel_overhang: 8,
            min_circ_overhang: 10,
            min_score_frac: 0.65,
            score_range_frac: 0.05,
            max_gene_candidates: 32,
            max_decoy_candidates: 8,
            max_frag_len: 2000,
            pair_is_rc: true,
            max_multi_hits: 32,
            enable_gapped: true,
        }
    }
}

/// One aligned exonic node segment within a `ReadAlignment` (supports arbitrary `u32` read lengths).
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct AlignedNodeSeg {
    /// 0-based local node index within `gene_id`.
    pub node_id: u16,
    /// 0-based start offset inside this node's sequence.
    pub node_offset: u32,
    /// Reference bases consumed inside this node.
    pub ref_len: u32,
    /// Read bases consumed inside this node.
    pub read_len: u32,
    /// Number of base matches in this segment.
    pub matches: u32,
    /// Number of base mismatches in this segment.
    pub mismatches: u32,
    /// Number of indel bases (insertions + deletions) in this segment.
    pub indels: u32,
}

/// Complete single-read alignment onto a gene's Contiguous Splice Graph.
#[derive(Clone, Debug)]
pub struct ReadAlignment {
    pub gene_id: u32,
    /// `true` iff the original read aligned directly to the gene's 5'->3' CSG `+` strand
    /// (`false` if the read was reverse-complemented to align to the gene's CSG `+` strand).
    pub read_is_sense: bool,
    /// 0-based offset of the first aligned base in the oriented read (equals left soft-clip length).
    pub soft_left: u32,
    /// Right soft-clip length in the oriented read.
    pub soft_right: u32,
    /// Alignment score under `AlignConfig` (same scale as decoy alignment).
    pub score: i32,
    /// Ordered sequence of aligned CSG node segments in 5'->3' transcript order.
    pub segs: SmallVec<[AlignedNodeSeg; 4]>,
    /// True if any junction in `segs` is a circular back-splice (`dst <= src`).
    pub has_circ: bool,
}

impl ReadAlignment {
    #[inline]
    pub fn first_node(&self) -> u16 {
        self.segs[0].node_id
    }

    #[inline]
    pub fn last_node(&self) -> u16 {
        self.segs[self.segs.len() - 1].node_id
    }

    /// 0-based start offset inside `first_node`.
    #[inline]
    pub fn start_offset_in_first_node(&self) -> u32 {
        self.segs[0].node_offset
    }

    /// 0-based exclusive end offset inside `last_node`.
    #[inline]
    pub fn end_offset_in_last_node(&self) -> u32 {
        let last = &self.segs[self.segs.len() - 1];
        last.node_offset + last.ref_len
    }

    /// Node path key (sequence of `node_id`s) for equivalence-class & deduplication checks.
    #[inline]
    pub fn node_path(&self) -> SmallVec<[u16; 4]> {
        self.segs.iter().map(|s| s.node_id).collect()
    }

    #[inline]
    pub fn total_matches(&self) -> u32 {
        self.segs.iter().map(|s| s.matches).sum()
    }

    #[inline]
    pub fn total_mismatches(&self) -> u32 {
        self.segs.iter().map(|s| s.mismatches).sum()
    }
}

/// Outcome of aligning a single-end read (including decoy screening status).
#[derive(Clone, Debug)]
pub enum SingleAlignOutcome {
    Unmapped,
    DecoyFiltered,
    Mapped(SmallVec<[ReadAlignment; 2]>),
}

/// Concordant paired-end alignment on a single gene.
#[derive(Clone, Debug)]
pub struct PairedAlignment {
    pub gene_id: u32,
    /// Upstream (5') mate alignment on gene's 5'->3' CSG.
    pub mate_5p: ReadAlignment,
    /// Downstream (3') mate alignment on gene's 5'->3' CSG.
    pub mate_3p: ReadAlignment,
    /// `true` iff R1 was the 5' mate (FR), `false` iff R2 was the 5' mate (RF).
    pub r1_is_5p: bool,
    /// Combined alignment score (`mate_5p.score + mate_3p.score`).
    pub score: i32,
    /// Exact transcriptomic fragment length if compatible with ≥ 1 annotated transcript,
    /// or CSG-spliced estimate for novel-junction pairs.
    pub frag_len: Option<u32>,
    /// Sorted global `tx_idx` values of all annotated transcripts compatible with both mates
    /// (fused during `compute_pair_fragment_length` to avoid a second per-gene transcript scan).
    pub compat_txs: SmallVec<[u32; 4]>,
}

/// Outcome of aligning a paired-end read pair.
#[derive(Clone, Debug)]
pub enum PairedAlignOutcome {
    Unmapped,
    DecoyFiltered,
    Mapped(SmallVec<[PairedAlignment; 2]>),
}

pub use extend::{align_single_read, AlignScratch};
pub use paired::align_paired_reads;
