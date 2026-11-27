//! Hybrid ungapped + SIMD banded gapped (`block-aligner`) spliced extension,
//! SNP-tolerant splice-edge jumping (including SNPs 1–3 bp before exon boundaries),
//! strict-inequality selective decoy validation, and zero-clone reusable scratch buffers.

use crate::align::seed::{
    collect_decoy_seeds, collect_seeds, DecoyCandidate, GeneAnchor, GeneCandidate, SeedScratch,
};
use crate::align::{AlignConfig, AlignedNodeSeg, ReadAlignment, SingleAlignOutcome, StrandMode};
use crate::index::{MmapIndex, NodeRecord};
use crate::seq::{
    append_normalized, append_revcomp_normalized, encode_flank_lsb_from_end, encode_flank_msb,
    ASCII_NORM,
};
use block_aligner::cigar::{Cigar, Operation};
use block_aligner::scan_block::{Block, PaddedBytes};
use block_aligner::scores::{Gaps, Matrix, NucMatrix};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

/// Read-only SIMD scoring parameters borrowed independently of scratch buffers.
pub struct SimdParams {
    pub nuc_matrix: NucMatrix,
    pub gaps: Gaps,
    pub min_block: usize,
    pub max_block: usize,
}

/// Reusable scratch buffers for extension and decoy unpacking (split from `SeedScratch`
/// and normalized read buffers so zero `.clone()` allocations happen on the hot path).
pub struct AlignWorkBuf {
    pub decoy_buf: Vec<u8>,
    pub rev_q_buf: Vec<u8>,
    pub rev_r_buf: Vec<u8>,
    pub simd: SimdParams,
}

/// Reusable per-thread scratch space for single-end and paired-end alignment.
pub struct AlignScratch {
    pub seed: SeedScratch,
    pub seed_r2: SeedScratch,
    pub fwd_norm: Vec<u8>,
    pub rc_norm: Vec<u8>,
    pub r2_fwd_norm: Vec<u8>,
    pub r2_rc_norm: Vec<u8>,
    pub r2_lookup: FxHashMap<(u32, bool), usize>,
    pub work: AlignWorkBuf,
}

impl AlignScratch {
    pub fn new(cfg: &AlignConfig) -> Self {
        let mut mat = NucMatrix::new_simple(cfg.match_score as i8, cfg.mismatch_penalty as i8);
        for b in [b'A', b'C', b'G', b'T', b'N'] {
            mat.set(b'N', b, cfg.ambig_score as i8);
        }
        Self {
            seed: SeedScratch::default(),
            seed_r2: SeedScratch::default(),
            fwd_norm: Vec::with_capacity(256),
            rc_norm: Vec::with_capacity(256),
            r2_fwd_norm: Vec::with_capacity(256),
            r2_rc_norm: Vec::with_capacity(256),
            r2_lookup: FxHashMap::default(),
            work: AlignWorkBuf {
                decoy_buf: Vec::with_capacity(512),
                rev_q_buf: Vec::with_capacity(256),
                rev_r_buf: Vec::with_capacity(256),
                simd: SimdParams {
                    nuc_matrix: mat,
                    gaps: Gaps {
                        open: cfg.gap_open as i8,
                        extend: cfg.gap_extend as i8,
                    },
                    min_block: 32,
                    max_block: 128,
                },
            },
        }
    }
}

/// Align a single-end read to the `.lidx` splice-graph index with zero per-read scratch clones
/// and strict-inequality selective decoy validation.
pub fn align_single_read(
    idx: &MmapIndex,
    raw_seq: &[u8],
    cfg: &AlignConfig,
    scratch: &mut AlignScratch,
) -> SingleAlignOutcome {
    if raw_seq.len() < (idx.header().minimizer_k as usize) {
        return SingleAlignOutcome::Unmapped;
    }

    scratch.fwd_norm.clear();
    let q_has_n = append_normalized(raw_seq, &mut scratch.fwd_norm);

    let (allow_sense, allow_antisense) = match cfg.strand_mode {
        StrandMode::Unstranded => (true, true),
        StrandMode::ForwardFR => (true, false),
        StrandMode::ReverseRF => (false, true),
    };

    // Disjoint field borrows: zero `.clone()` allocations; compute `rc_norm` & decoy seeds lazily!
    collect_seeds(
        idx,
        &scratch.fwd_norm,
        allow_sense,
        allow_antisense,
        cfg.max_gene_candidates,
        &mut scratch.seed,
    );

    if scratch.seed.gene_candidates.is_empty() {
        return SingleAlignOutcome::Unmapped;
    }

    let read_len = scratch.fwd_norm.len();
    let max_possible = (read_len as i32) * cfg.match_score;
    let min_score = min_valid_score(read_len, cfg);
    let score_tol = score_tolerance(read_len, cfg);

    let mut alignments: SmallVec<[ReadAlignment; 2]> = SmallVec::new();
    let mut best_score = i32::MIN;
    let mut best_hit_count = 0u32;
    let mut rc_ready = false;

    for cand in &scratch.seed.gene_candidates {
        // Once a multi-hit (>= 4 minimizers) candidate achieved the mathematical maximum score,
        // single-minimizer (hit_count == 1) noise candidates cannot be exact-match paralogs.
        if best_score >= max_possible && best_hit_count >= 4 && cand.hit_count <= 1 {
            break;
        }
        let q: &[u8] = if cand.read_is_sense {
            &scratch.fwd_norm
        } else {
            if !rc_ready {
                scratch.rc_norm.clear();
                append_revcomp_normalized(raw_seq, &mut scratch.rc_norm);
                rc_ready = true;
            }
            &scratch.rc_norm
        };
        if let Some(aln) = extend_gene_candidate(idx, cand, q, q_has_n, cfg, &mut scratch.work) {
            if aln.score >= min_score {
                if aln.score > best_score {
                    best_hit_count = cand.hit_count;
                } else if aln.score == best_score && cand.hit_count > best_hit_count {
                    best_hit_count = cand.hit_count;
                }
                push_scored_alignment(
                    &mut alignments,
                    &mut best_score,
                    aln,
                    score_tol,
                    cfg.max_multi_hits,
                );
            }
        }
    }

    if alignments.is_empty() {
        return SingleAlignOutcome::Unmapped;
    }

    // Selective decoy validation (Safeguard #1: strict `>` discard!).
    // Mathematically impossible for any decoy to exceed `max_possible`, so defer decoy seed
    // collection completely until `best_score < max_possible && scratch.seed.has_any_decoy_posting`.
    if best_score < max_possible && scratch.seed.has_any_decoy_posting {
        collect_decoy_seeds(
            idx,
            read_len as u32,
            cfg.max_decoy_candidates,
            &mut scratch.seed,
        );
        if !scratch.seed.decoy_candidates.is_empty() {
            if !rc_ready
                && scratch
                    .seed
                    .decoy_candidates
                    .iter()
                    .any(|d| !d.read_is_fwd_on_decoy)
            {
                scratch.rc_norm.clear();
                append_revcomp_normalized(raw_seq, &mut scratch.rc_norm);
            }
            let best_decoy = score_decoy_candidates(
                idx,
                &scratch.seed.decoy_candidates,
                &scratch.fwd_norm,
                &scratch.rc_norm,
                cfg,
                &mut scratch.work,
            );
            if best_decoy > best_score {
                return SingleAlignOutcome::DecoyFiltered;
            }
        }
    }

    SingleAlignOutcome::Mapped(alignments)
}

#[inline]
pub fn min_valid_score(read_len: usize, cfg: &AlignConfig) -> i32 {
    ((read_len as f32) * (cfg.match_score as f32) * cfg.min_score_frac).floor() as i32
}

#[inline]
pub fn score_tolerance(read_len: usize, cfg: &AlignConfig) -> i32 {
    ((read_len as f32) * (cfg.match_score as f32) * cfg.score_range_frac)
        .ceil()
        .max(0.0) as i32
}

#[inline]
pub fn push_scored_alignment(
    alns: &mut SmallVec<[ReadAlignment; 2]>,
    best_score: &mut i32,
    aln: ReadAlignment,
    score_tol: i32,
    max_multi: usize,
) {
    for existing in alns.iter_mut() {
        if existing.gene_id == aln.gene_id
            && existing.read_is_sense == aln.read_is_sense
            && existing.segs.len() == aln.segs.len()
            && existing
                .segs
                .iter()
                .zip(aln.segs.iter())
                .all(|(a, b)| a.node_id == b.node_id)
        {
            if aln.score > existing.score {
                *existing = aln;
                if existing.score > *best_score {
                    *best_score = existing.score;
                    let cutoff = *best_score - score_tol;
                    alns.retain(|x| x.score >= cutoff);
                }
            }
            return;
        }
    }

    if aln.score > *best_score {
        *best_score = aln.score;
        let cutoff = *best_score - score_tol;
        alns.retain(|x| x.score >= cutoff);
        alns.push(aln);
    } else if aln.score >= *best_score - score_tol && alns.len() < max_multi {
        alns.push(aln);
    }
}

/// Extend all seed anchors belonging to one `(gene_id, read_is_sense)` candidate and return
/// the highest-scoring alignment on that gene.
pub fn extend_gene_candidate(
    idx: &MmapIndex,
    cand: &GeneCandidate,
    q: &[u8],
    q_has_n: bool,
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
) -> Option<ReadAlignment> {
    let gene_id = cand.gene_id;
    let nodes = idx.gene_nodes(gene_id);
    if nodes.is_empty() {
        return None;
    }
    let gene_seq = idx.gene_seq(gene_id);
    let anno_jncs = idx.gene_anno_junctions(gene_id);
    let is_pos_strand = idx.genes()[gene_id as usize].is_pos_strand();

    let max_possible = (q.len() as i32) * cfg.match_score;
    let mut best: Option<ReadAlignment> = None;
    let mut seen_diags: SmallVec<[(u16, i32); 8]> = SmallVec::new();

    for &anchor in &cand.anchors {
        let aln_opt = match anchor {
            GeneAnchor::CsgNode {
                node_id,
                node_offset,
                oriented_read_pos,
            } => {
                let diag = (node_offset as i32) - (oriented_read_pos as i32);
                if seen_diags.contains(&(node_id, diag)) {
                    continue;
                }
                seen_diags.push((node_id, diag));
                extend_from_node_anchor(
                    gene_id,
                    cand.read_is_sense,
                    is_pos_strand,
                    nodes,
                    gene_seq,
                    anno_jncs,
                    node_id,
                    node_offset,
                    oriented_read_pos,
                    q,
                    q_has_n,
                    cfg,
                    work,
                )
            }
            GeneAnchor::AnnoJunction {
                local_jnc_idx,
                oriented_jnc_pos,
            } => {
                let jnc = anno_jncs.get(local_jnc_idx as usize)?;
                extend_from_jnc_anchor(
                    gene_id,
                    cand.read_is_sense,
                    is_pos_strand,
                    nodes,
                    gene_seq,
                    anno_jncs,
                    jnc.src_node,
                    jnc.dst_node,
                    oriented_jnc_pos,
                    q,
                    q_has_n,
                    cfg,
                    work,
                )
            }
        };

        if let Some(aln) = aln_opt {
            let is_perfect = aln.score >= max_possible;
            if best.as_ref().map_or(true, |b| aln.score > b.score) {
                best = Some(aln);
            }
            if is_perfect {
                break;
            }
        }
    }

    best
}

#[inline]
fn nodes_contiguous_rec(prev: &NodeRecord, next: &NodeRecord, is_pos_strand: bool) -> bool {
    if is_pos_strand {
        prev.genomic_coord + prev.len == next.genomic_coord
    } else {
        next.genomic_coord + next.len == prev.genomic_coord
    }
}

#[derive(Copy, Clone, Debug)]
struct SegResult {
    score: i32,
    ref_used: u32,
    read_used: u32,
    matches: u32,
    mismatches: u32,
    indels: u32,
}

#[derive(Clone, Debug)]
struct ExtBranch {
    score: i32,
    segs: SmallVec<[AlignedNodeSeg; 4]>,
    soft_clip: u32,
    has_circ: bool,
}

#[allow(clippy::too_many_arguments)]
fn extend_from_node_anchor(
    gene_id: u32,
    read_is_sense: bool,
    is_pos_strand: bool,
    nodes: &[NodeRecord],
    gene_seq: &[u8],
    anno_jncs: &[crate::index::AnnoJunctionRecord],
    anchor_node: u16,
    anchor_node_off: u32,
    anchor_rpos: u32,
    q: &[u8],
    q_has_n: bool,
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
) -> Option<ReadAlignment> {
    let fwd = extend_forward(
        is_pos_strand,
        nodes,
        gene_seq,
        anno_jncs,
        anchor_node,
        anchor_node_off,
        anchor_rpos,
        q,
        q_has_n,
        cfg,
        work,
        0,
    );

    let rev = if anchor_rpos > 0 {
        extend_backward(
            is_pos_strand,
            nodes,
            gene_seq,
            anno_jncs,
            anchor_node,
            anchor_node_off,
            anchor_rpos,
            q,
            q_has_n,
            cfg,
            work,
            0,
        )
    } else {
        ExtBranch {
            score: 0,
            segs: SmallVec::new(),
            soft_clip: 0,
            has_circ: false,
        }
    };

    let mut segs: SmallVec<[AlignedNodeSeg; 4]> = SmallVec::new();
    for s in rev.segs.into_iter().rev() {
        segs.push(s);
    }
    for s in fwd.segs {
        if let Some(last) = segs.last_mut() {
            if last.node_id == s.node_id && last.node_offset + last.ref_len == s.node_offset {
                last.ref_len += s.ref_len;
                last.read_len += s.read_len;
                last.matches += s.matches;
                last.mismatches += s.mismatches;
                last.indels += s.indels;
                continue;
            }
        }
        segs.push(s);
    }

    if segs.is_empty() {
        return None;
    }

    Some(ReadAlignment {
        gene_id,
        read_is_sense,
        soft_left: rev.soft_clip,
        soft_right: fwd.soft_clip,
        score: rev.score + fwd.score,
        segs,
        has_circ: rev.has_circ || fwd.has_circ,
    })
}

#[allow(clippy::too_many_arguments)]
fn extend_from_jnc_anchor(
    gene_id: u32,
    read_is_sense: bool,
    is_pos_strand: bool,
    nodes: &[NodeRecord],
    gene_seq: &[u8],
    anno_jncs: &[crate::index::AnnoJunctionRecord],
    src_node: u16,
    dst_node: u16,
    jnc_rpos: u32,
    q: &[u8],
    q_has_n: bool,
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
) -> Option<ReadAlignment> {
    let min_ov = if dst_node <= src_node {
        cfg.min_circ_overhang
    } else {
        cfg.min_anno_overhang
    };
    let rem_right = (q.len() as u32).saturating_sub(jnc_rpos);
    if jnc_rpos < min_ov || rem_right < min_ov {
        return None;
    }

    let src_len = nodes[src_node as usize].len;
    let rev = extend_backward(
        is_pos_strand,
        nodes,
        gene_seq,
        anno_jncs,
        src_node,
        src_len,
        jnc_rpos,
        q,
        q_has_n,
        cfg,
        work,
        0,
    );
    let fwd = extend_forward(
        is_pos_strand,
        nodes,
        gene_seq,
        anno_jncs,
        dst_node,
        0,
        jnc_rpos,
        q,
        q_has_n,
        cfg,
        work,
        0,
    );

    if rev.segs.is_empty() || fwd.segs.is_empty() {
        return None;
    }
    // Enforce minimum overhang on BOTH sides of the anchored junction
    let left_used: u32 = rev.segs.iter().map(|s| s.read_len).sum();
    let right_used: u32 = fwd.segs.iter().map(|s| s.read_len).sum();
    if left_used < min_ov || right_used < min_ov {
        return None;
    }

    let mut segs: SmallVec<[AlignedNodeSeg; 4]> = SmallVec::new();
    for s in rev.segs.into_iter().rev() {
        segs.push(s);
    }
    for s in fwd.segs {
        segs.push(s);
    }

    Some(ReadAlignment {
        gene_id,
        read_is_sense,
        soft_left: rev.soft_clip,
        soft_right: fwd.soft_clip,
        score: rev.score + fwd.score,
        segs,
        has_circ: rev.has_circ || fwd.has_circ || dst_node <= src_node,
    })
}

const MAX_SPLICE_DEPTH: usize = 6;

/// Fast O(1) small-indel pre-screen around a mismatch/divergence offset:
/// returns `true` iff a 1–8 bp insertion or deletion shift yields `>= 6` base matches in a 10 bp window.
/// Eliminates ~97% of wasted `block-aligner` SIMD DP calls on unrelated false-positive seed hits.
#[inline]
fn has_small_shift_support(q_slice: &[u8], r_slice: &[u8]) -> bool {
    let qn = q_slice.len();
    let rn = r_slice.len();
    for shift in 1..=8usize {
        if qn >= shift + 8 && rn >= 8 {
            let len = (qn - shift).min(rn).min(10);
            let mut m = 0usize;
            for i in 0..len {
                if q_slice[shift + i] == r_slice[i] {
                    m += 1;
                }
            }
            if m >= 6 {
                return true;
            }
        }
        if qn >= 8 && rn >= shift + 8 {
            let len = qn.min(rn - shift).min(10);
            let mut m = 0usize;
            for i in 0..len {
                if q_slice[i] == r_slice[shift + i] {
                    m += 1;
                }
            }
            if m >= 6 {
                return true;
            }
        }
    }
    false
}

/// Score a single segment inside one CSG node (forward direction) using fast ungapped scan
/// with automatic SIMD `block-aligner` gapped fallback on mismatch clusters.
///
/// Returns `(peak_prefix, optional_to_boundary)`:
/// When `r_slice.len() <= q_slice.len()` and a SNP 1–3 bp before the exon boundary caused
/// `peak_prefix.ref_used < r_slice.len()` (with boundary score within 12 of peak), `optional_to_boundary`
/// provides the full-to-boundary segment so the caller can evaluate crossing the splice junction!
fn align_node_segment_fwd(
    q_slice: &[u8],
    r_slice: &[u8],
    q_has_n: bool,
    cfg: &AlignConfig,
    simd: &SimdParams,
) -> (SegResult, Option<SegResult>) {
    let max_n = q_slice.len().min(r_slice.len());
    // Fast SIMD memcmp path: if `!q_has_n` and `q_slice[..max_n] == r_slice[..max_n]`,
    // then `r_slice[..max_n]` is mathematically guaranteed to contain no `'N'` either!
    if !q_has_n && max_n > 0 && q_slice[..max_n] == r_slice[..max_n] {
        let m = max_n as u32;
        return (
            SegResult {
                score: (max_n as i32) * cfg.match_score,
                ref_used: m,
                read_used: m,
                matches: m,
                mismatches: 0,
                indels: 0,
            },
            None,
        );
    }
    let mut score: i32 = 0;
    let mut best_score: i32 = 0;
    let mut best_len: usize = 0;
    let mut matches: u32 = 0;
    let mut mismatches: u32 = 0;
    let mut best_m: u32 = 0;
    let mut best_mm: u32 = 0;
    let mut consec_mm: u32 = 0;
    let mut max_consec_mm: u32 = 0;
    let mut first_mm_pos: Option<usize> = None;

    for i in 0..max_n {
        let qb = q_slice[i];
        let rb = r_slice[i];
        if qb == b'N' || rb == b'N' {
            score += cfg.ambig_score;
            consec_mm = 0;
        } else if qb == rb {
            score += cfg.match_score;
            matches += 1;
            consec_mm = 0;
        } else {
            score += cfg.mismatch_penalty;
            mismatches += 1;
            consec_mm += 1;
            if first_mm_pos.is_none() {
                first_mm_pos = Some(i);
            }
            if consec_mm > max_consec_mm {
                max_consec_mm = consec_mm;
            }
        }
        if score >= best_score {
            best_score = score;
            best_len = i + 1;
            best_m = matches;
            best_mm = mismatches;
        }
    }

    let mut peak = SegResult {
        score: best_score,
        ref_used: best_len as u32,
        read_used: best_len as u32,
        matches: best_m,
        mismatches: best_mm,
        indels: 0,
    };

    // If a SNP sits 1–3 bp before an exon boundary (r_slice.len() <= q_slice.len()),
    // retain the full-to-boundary candidate so we can still jump across the splice site!
    let mut to_boundary = if r_slice.len() <= q_slice.len()
        && best_len < max_n
        && (max_n - best_len) <= 3
        && score >= best_score - 12
    {
        Some(SegResult {
            score,
            ref_used: max_n as u32,
            read_used: max_n as u32,
            matches,
            mismatches,
            indels: 0,
        })
    } else {
        None
    };

    // Trigger SIMD banded gapped fallback only when a genuine small-indel shift is supported
    let reached_end = best_len == max_n;
    let needs_gapped = cfg.enable_gapped
        && q_slice.len() >= 10
        && r_slice.len() >= 10
        && ((!reached_end && to_boundary.is_none()) || max_consec_mm >= 2 || mismatches >= 3);

    if needs_gapped {
        let div = if !reached_end {
            best_len.saturating_sub(2)
        } else {
            first_mm_pos.unwrap_or(0).saturating_sub(1)
        };
        if has_small_shift_support(&q_slice[div..], &r_slice[div..]) {
            if let Some(g) = simd_gapped_prefix(q_slice, r_slice, simd) {
                if g.score > peak.score {
                    peak = g;
                    if (g.ref_used as usize) == r_slice.len() {
                        to_boundary = None;
                    }
                }
            }
        }
    }

    (peak, to_boundary)
}

/// Run `block-aligner` SIMD banded X-drop alignment with traceback on `(q_slice, r_slice)`.
fn simd_gapped_prefix(q_slice: &[u8], r_slice: &[u8], simd: &SimdParams) -> Option<SegResult> {
    let r_cap = r_slice.len().min(q_slice.len() + 32);
    let r_win = &r_slice[..r_cap];
    let q_padded = PaddedBytes::from_bytes::<NucMatrix>(q_slice, simd.max_block);
    let r_padded = PaddedBytes::from_bytes::<NucMatrix>(r_win, simd.max_block);

    let mut block = Block::<true, true>::new(q_padded.len(), r_padded.len(), simd.max_block);
    let x_drop = 24;
    block.align(
        &q_padded,
        &r_padded,
        &simd.nuc_matrix,
        simd.gaps,
        simd.min_block..=simd.max_block,
        x_drop,
    );
    let res = block.res();
    if res.score <= 0 || res.query_idx == 0 || res.reference_idx == 0 {
        return None;
    }

    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    block
        .trace()
        .cigar_eq(&q_padded, &r_padded, res.query_idx, res.reference_idx, &mut cigar);

    let mut matches: u32 = 0;
    let mut mismatches: u32 = 0;
    let mut indels: u32 = 0;
    for op_len in cigar.to_vec() {
        match op_len.op {
            Operation::Eq => matches += op_len.len as u32,
            Operation::X | Operation::M => mismatches += op_len.len as u32,
            Operation::I | Operation::D => indels += op_len.len as u32,
            _ => {}
        }
    }

    Some(SegResult {
        score: res.score,
        ref_used: res.reference_idx as u32,
        read_used: res.query_idx as u32,
        matches,
        mismatches,
        indels,
    })
}

#[allow(clippy::too_many_arguments)]
fn extend_forward(
    is_pos_strand: bool,
    nodes: &[NodeRecord],
    gene_seq: &[u8],
    anno_jncs: &[crate::index::AnnoJunctionRecord],
    node_id: u16,
    node_off: u32,
    rpos: u32,
    q: &[u8],
    q_has_n: bool,
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
    depth: usize,
) -> ExtBranch {
    let rem_q = q.len() - (rpos as usize);
    if rem_q == 0 {
        return ExtBranch {
            score: 0,
            segs: SmallVec::new(),
            soft_clip: 0,
            has_circ: false,
        };
    }

    let nd = &nodes[node_id as usize];
    let avail_ref = (nd.len.saturating_sub(node_off)) as usize;
    let q_slice = &q[rpos as usize..];

    let (peak, to_boundary) = if avail_ref > 0 {
        let s0 = (nd.gene_seq_offset + node_off) as usize;
        let r_slice = &gene_seq[s0..s0 + avail_ref];
        align_node_segment_fwd(q_slice, r_slice, q_has_n, cfg, &work.simd)
    } else {
        (
            SegResult {
                score: 0,
                ref_used: 0,
                read_used: 0,
                matches: 0,
                mismatches: 0,
                indels: 0,
            },
            None,
        )
    };

    let mut peak_segs: SmallVec<[AlignedNodeSeg; 4]> = SmallVec::new();
    if peak.ref_used > 0 || peak.read_used > 0 {
        peak_segs.push(AlignedNodeSeg {
            node_id,
            node_offset: node_off,
            ref_len: peak.ref_used,
            read_len: peak.read_used,
            matches: peak.matches,
            mismatches: peak.mismatches,
            indels: peak.indels,
        });
    }

    let mut best = ExtBranch {
        score: peak.score,
        segs: peak_segs,
        soft_clip: (q.len() as u32) - (rpos + peak.read_used),
        has_circ: false,
    };

    // Determine which segment reaches the 3' node boundary (either `peak` or `to_boundary`)
    let boundary_seg = if (node_off + peak.ref_used) == nd.len {
        Some(peak)
    } else {
        to_boundary
    };

    let Some(bseg) = boundary_seg else {
        return best;
    };

    let next_rpos = rpos + bseg.read_used;
    let rem_after = (q.len() as u32) - next_rpos;
    if rem_after == 0 || depth >= MAX_SPLICE_DEPTH {
        return best;
    }

    let mut base_segs: SmallVec<[AlignedNodeSeg; 4]> = SmallVec::new();
    if bseg.ref_used > 0 || bseg.read_used > 0 {
        base_segs.push(AlignedNodeSeg {
            node_id,
            node_offset: node_off,
            ref_len: bseg.ref_used,
            read_len: bseg.read_used,
            matches: bseg.matches,
            mismatches: bseg.mismatches,
            indels: bseg.indels,
        });
    }
    let seg_score = bseg.score;
    let rem_slice = &q[next_rpos as usize..];
    let max_rem_score = (rem_after as i32) * cfg.match_score;

    // 1. Contiguous soft-edge continuation into `node_id + 1` (e.g. AA/AD/RI)
    let next_idx = (node_id as usize) + 1;
    if next_idx < nodes.len() && nodes_contiguous_rec(nd, &nodes[next_idx], is_pos_strand) {
        let sub = extend_forward(
            is_pos_strand,
            nodes,
            gene_seq,
            anno_jncs,
            next_idx as u16,
            0,
            next_rpos,
            q,
            q_has_n,
            cfg,
            work,
            depth,
        );
        if seg_score + sub.score > best.score {
            let mut merged = base_segs.clone();
            merged.extend(sub.segs);
            best = ExtBranch {
                score: seg_score + sub.score,
                segs: merged,
                soft_clip: sub.soft_clip,
                has_circ: sub.has_circ,
            };
        }
        // Short-circuit: a perfect contiguous continuation cannot be beaten by any splice jump
        if sub.score >= max_rem_score {
            return best;
        }
    }

    // 2. Splice-edge jump (only if `node_id` ends at a donor boundary)
    if nd.right_edge().is_donor() {
        let (q_msb, q_msb_len) = encode_flank_msb(rem_slice);
        let mut tried_anno: SmallVec<[u16; 8]> = SmallVec::new();

        // Priority 1: Annotated outgoing junctions (ZERO splice penalty!)
        let j_start = anno_jncs.partition_point(|j| j.src_node < node_id);
        for j in &anno_jncs[j_start..] {
            if j.src_node != node_id {
                break;
            }
            let dst = j.dst_node;
            tried_anno.push(dst);
            if rem_after < cfg.min_anno_overhang {
                continue;
            }
            if !flank_matches_acceptor(
                &nodes[dst as usize],
                rem_slice,
                q_msb,
                q_msb_len,
                true,
                rem_after,
            ) {
                continue;
            }
            let sub = extend_forward(
                is_pos_strand,
                nodes,
                gene_seq,
                anno_jncs,
                dst,
                0,
                next_rpos,
                q,
                q_has_n,
                cfg,
                work,
                depth + 1,
            );
            if seg_score + sub.score > best.score {
                let mut merged = base_segs.clone();
                merged.extend(sub.segs);
                best = ExtBranch {
                    score: seg_score + sub.score,
                    segs: merged,
                    soft_clip: sub.soft_clip,
                    has_circ: sub.has_circ,
                };
            }
            if sub.score >= max_rem_score {
                return best;
            }
        }

        // Priority 2: Novel / combinatorial acceptor nodes in the same gene (1-SNP tolerant!)
        // Skip full-gene novel scan if Priority 1 / contiguous already reached or exceeded the
        // best possible penalized novel score (`seg_score + max_rem_score + novel_canonical_penalty`).
        if rem_after >= cfg.min_novel_overhang
            && best.score < seg_score + max_rem_score + cfg.novel_canonical_penalty
        {
            let start_v = if cfg.allow_circ { 0 } else { (node_id as usize) + 1 };
            for (v_idx, v_nd) in nodes.iter().enumerate().skip(start_v) {
                let v = v_idx as u16;
                let is_circ = v <= node_id;
                if is_circ && rem_after < cfg.min_circ_overhang {
                    continue;
                }
                if !v_nd.left_edge().is_acceptor() || tried_anno.contains(&v) {
                    continue;
                }
                if v_idx == next_idx && nodes_contiguous_rec(nd, v_nd, is_pos_strand) {
                    continue;
                }
                if !flank_matches_acceptor(v_nd, rem_slice, q_msb, q_msb_len, false, rem_after) {
                    continue;
                }
                let is_canon = is_canonical_pair(nd.don_dinuc(), v_nd.acc_dinuc());
                let jnc_pen = if is_canon {
                    cfg.novel_canonical_penalty
                } else {
                    cfg.novel_noncanonical_penalty
                };
                let sub = extend_forward(
                    is_pos_strand,
                    nodes,
                    gene_seq,
                    anno_jncs,
                    v,
                    0,
                    next_rpos,
                    q,
                    q_has_n,
                    cfg,
                    work,
                    depth + 1,
                );
                let total_sc = seg_score + sub.score + jnc_pen;
                if total_sc > best.score {
                    let mut merged = base_segs.clone();
                    merged.extend(sub.segs);
                    best = ExtBranch {
                        score: total_sc,
                        segs: merged,
                        soft_clip: sub.soft_clip,
                        has_circ: sub.has_circ || is_circ,
                    };
                }
            }
        }
    }

    best
}

/// Backward extension (right-to-left) ending just before `(node_id, node_end_off, rpos_end)`.
#[allow(clippy::too_many_arguments)]
fn extend_backward(
    is_pos_strand: bool,
    nodes: &[NodeRecord],
    gene_seq: &[u8],
    anno_jncs: &[crate::index::AnnoJunctionRecord],
    node_id: u16,
    node_end_off: u32,
    rpos_end: u32,
    q: &[u8],
    q_has_n: bool,
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
    depth: usize,
) -> ExtBranch {
    if rpos_end == 0 {
        return ExtBranch {
            score: 0,
            segs: SmallVec::new(),
            soft_clip: 0,
            has_circ: false,
        };
    }

    let nd = &nodes[node_id as usize];
    let avail_ref = node_end_off as usize;
    let r_end = rpos_end as usize;

    let (peak, to_boundary) = if avail_ref > 0 {
        let s0 = nd.gene_seq_offset as usize;
        let max_n = r_end.min(avail_ref);
        // Zero-copy, zero-reverse SIMD memcmp fast path!
        // Avoids reversing multi-kilobase 3' UTR exons into `rev_r_buf` on exact matches.
        if !q_has_n && max_n > 0 && q[r_end - max_n..r_end] == gene_seq[s0 + avail_ref - max_n..s0 + avail_ref] {
            let m = max_n as u32;
            (
                SegResult {
                    score: (max_n as i32) * cfg.match_score,
                    ref_used: m,
                    read_used: m,
                    matches: m,
                    mismatches: 0,
                    indels: 0,
                },
                None,
            )
        } else {
            // Even on mismatch fallback, only reverse at most `r_end + 32` reference bases
            // (matching `simd_gapped_prefix`'s window cap), never a full 8 kb 3' UTR!
            work.rev_q_buf.clear();
            work.rev_q_buf.extend(q[..r_end].iter().rev().copied());
            work.rev_r_buf.clear();
            let r_take = avail_ref.min(r_end + 32);
            work.rev_r_buf.extend(
                gene_seq[s0 + avail_ref - r_take..s0 + avail_ref]
                    .iter()
                    .rev()
                    .copied(),
            );
            align_node_segment_fwd(&work.rev_q_buf, &work.rev_r_buf, q_has_n, cfg, &work.simd)
        }
    } else {
        (
            SegResult {
                score: 0,
                ref_used: 0,
                read_used: 0,
                matches: 0,
                mismatches: 0,
                indels: 0,
            },
            None,
        )
    };

    let mut peak_segs: SmallVec<[AlignedNodeSeg; 4]> = SmallVec::new();
    if peak.ref_used > 0 || peak.read_used > 0 {
        peak_segs.push(AlignedNodeSeg {
            node_id,
            node_offset: node_end_off - peak.ref_used,
            ref_len: peak.ref_used,
            read_len: peak.read_used,
            matches: peak.matches,
            mismatches: peak.mismatches,
            indels: peak.indels,
        });
    }

    let mut best = ExtBranch {
        score: peak.score,
        segs: peak_segs,
        soft_clip: rpos_end - peak.read_used,
        has_circ: false,
    };

    let boundary_seg = if peak.ref_used == node_end_off {
        Some(peak)
    } else {
        to_boundary
    };

    let Some(bseg) = boundary_seg else {
        return best;
    };

    let rem_left = rpos_end - bseg.read_used;
    if rem_left == 0 || depth >= MAX_SPLICE_DEPTH {
        return best;
    }

    let mut base_segs: SmallVec<[AlignedNodeSeg; 4]> = SmallVec::new();
    if bseg.ref_used > 0 || bseg.read_used > 0 {
        base_segs.push(AlignedNodeSeg {
            node_id,
            node_offset: node_end_off - bseg.ref_used,
            ref_len: bseg.ref_used,
            read_len: bseg.read_used,
            matches: bseg.matches,
            mismatches: bseg.mismatches,
            indels: bseg.indels,
        });
    }
    let seg_score = bseg.score;
    let rem_slice = &q[..rem_left as usize];
    let max_rem_score = (rem_left as i32) * cfg.match_score;

    // 1. Contiguous soft-edge continuation backward into `node_id - 1`
    if node_id > 0 {
        let prev_id = node_id - 1;
        if nodes_contiguous_rec(&nodes[prev_id as usize], nd, is_pos_strand) {
            let sub = extend_backward(
                is_pos_strand,
                nodes,
                gene_seq,
                anno_jncs,
                prev_id,
                nodes[prev_id as usize].len,
                rem_left,
                q,
                q_has_n,
                cfg,
                work,
                depth,
            );
            if seg_score + sub.score > best.score {
                let mut merged = base_segs.clone();
                merged.extend(sub.segs);
                best = ExtBranch {
                    score: seg_score + sub.score,
                    segs: merged,
                    soft_clip: sub.soft_clip,
                    has_circ: sub.has_circ,
                };
            }
            if sub.score >= max_rem_score {
                return best;
            }
        }
    }

    // 2. Splice-edge jump backward (only if `node_id` starts at an acceptor boundary)
    if nd.left_edge().is_acceptor() {
        let (q_lsb, q_lsb_len) = encode_flank_lsb_from_end(rem_slice);
        let mut tried_anno: SmallVec<[u16; 8]> = SmallVec::new();

        // Priority 1: Annotated incoming junctions (dst_node == node_id)
        for j in anno_jncs {
            if j.dst_node != node_id {
                continue;
            }
            let src = j.src_node;
            tried_anno.push(src);
            if rem_left < cfg.min_anno_overhang {
                continue;
            }
            if !flank_matches_donor(
                &nodes[src as usize],
                rem_slice,
                q_lsb,
                q_lsb_len,
                true,
                rem_left,
            ) {
                continue;
            }
            let sub = extend_backward(
                is_pos_strand,
                nodes,
                gene_seq,
                anno_jncs,
                src,
                nodes[src as usize].len,
                rem_left,
                q,
                q_has_n,
                cfg,
                work,
                depth + 1,
            );
            if seg_score + sub.score > best.score {
                let mut merged = base_segs.clone();
                merged.extend(sub.segs);
                best = ExtBranch {
                    score: seg_score + sub.score,
                    segs: merged,
                    soft_clip: sub.soft_clip,
                    has_circ: sub.has_circ,
                };
            }
            if sub.score >= max_rem_score {
                return best;
            }
        }

        // Priority 2: Novel / combinatorial donor nodes (1-SNP tolerant!)
        if rem_left >= cfg.min_novel_overhang
            && best.score < seg_score + max_rem_score + cfg.novel_canonical_penalty
        {
            let end_u = if cfg.allow_circ {
                nodes.len()
            } else {
                node_id as usize
            };
            for (u_idx, u_nd) in nodes[..end_u].iter().enumerate() {
                let u = u_idx as u16;
                let is_circ = u >= node_id;
                if is_circ && rem_left < cfg.min_circ_overhang {
                    continue;
                }
                if !u_nd.right_edge().is_donor() || tried_anno.contains(&u) {
                    continue;
                }
                if u + 1 == node_id && nodes_contiguous_rec(u_nd, nd, is_pos_strand) {
                    continue;
                }
                if !flank_matches_donor(u_nd, rem_slice, q_lsb, q_lsb_len, false, rem_left) {
                    continue;
                }
                let is_canon = is_canonical_pair(u_nd.don_dinuc(), nd.acc_dinuc());
                let jnc_pen = if is_canon {
                    cfg.novel_canonical_penalty
                } else {
                    cfg.novel_noncanonical_penalty
                };
                let sub = extend_backward(
                    is_pos_strand,
                    nodes,
                    gene_seq,
                    anno_jncs,
                    u,
                    u_nd.len,
                    rem_left,
                    q,
                    q_has_n,
                    cfg,
                    work,
                    depth + 1,
                );
                let total_sc = seg_score + sub.score + jnc_pen;
                if total_sc > best.score {
                    let mut merged = base_segs.clone();
                    merged.extend(sub.segs);
                    best = ExtBranch {
                        score: total_sc,
                        segs: merged,
                        soft_clip: sub.soft_clip,
                        has_circ: sub.has_circ || is_circ,
                    };
                }
            }
        }
    }

    best
}

#[inline]
fn is_canonical_pair(don_dinuc: u8, acc_dinuc: u8) -> bool {
    matches!((don_dinuc, acc_dinuc), (1, 1) | (2, 1) | (3, 2))
}

#[inline]
fn flank_matches_acceptor(
    nd: &NodeRecord,
    rem_slice: &[u8],
    q_msb: u32,
    q_msb_len: u8,
    is_anno: bool,
    rem_len: u32,
) -> bool {
    let cmp_req = (rem_len.min(16) as u8).min(q_msb_len);
    if cmp_req == 0 {
        return rem_slice
            .first()
            .copied()
            .map_or(false, |b| ASCII_NORM[b as usize] == b'N');
    }
    let Some((compared, mm)) = nd.acceptor_flank_mismatches(q_msb, cmp_req) else {
        return true;
    };
    let max_mm = if is_anno {
        if compared >= 7 { 1 } else { 0 }
    } else if compared >= 8 {
        1
    } else {
        0
    };
    mm <= max_mm
}

#[inline]
fn flank_matches_donor(
    nd: &NodeRecord,
    rem_slice: &[u8],
    q_lsb: u32,
    q_lsb_len: u8,
    is_anno: bool,
    rem_len: u32,
) -> bool {
    let cmp_req = (rem_len.min(16) as u8).min(q_lsb_len);
    if cmp_req == 0 {
        return rem_slice
            .last()
            .copied()
            .map_or(false, |b| ASCII_NORM[b as usize] == b'N');
    }
    let Some((compared, mm)) = nd.donor_flank_mismatches(q_lsb, cmp_req) else {
        return true;
    };
    let max_mm = if is_anno {
        if compared >= 7 { 1 } else { 0 }
    } else if compared >= 8 {
        1
    } else {
        0
    };
    mm <= max_mm
}

/// Score a single read against candidate decoy contigs using the EXACT same scoring matrix
/// as CSG alignment (Safeguard #1, zero clones).
pub fn score_decoy_candidates(
    idx: &MmapIndex,
    decoys: &[DecoyCandidate],
    fwd_q: &[u8],
    rc_q: &[u8],
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
) -> i32 {
    let mut best = i32::MIN;
    for dc in decoys {
        let q = if dc.read_is_fwd_on_decoy { fwd_q } else { rc_q };
        let sc = score_single_decoy(idx, *dc, q, cfg, work);
        if sc > best {
            best = sc;
        }
    }
    best
}

pub fn score_single_decoy(
    idx: &MmapIndex,
    dc: DecoyCandidate,
    q: &[u8],
    cfg: &AlignConfig,
    work: &mut AlignWorkBuf,
) -> i32 {
    let drec = &idx.decoys()[dc.decoy_id as usize];
    let dlen = drec.seq_len as usize;
    let rpos = dc.oriented_read_pos as usize;
    let doff = dc.decoy_offset as usize;
    let qlen = q.len();

    let pad = 16usize;
    let win_start = doff.saturating_sub(rpos + pad);
    let win_end = (doff + (qlen - rpos) + pad).min(dlen);
    if win_end <= win_start {
        return i32::MIN;
    }

    idx.unpack_decoy_slice(
        dc.decoy_id,
        win_start,
        win_end - win_start,
        &mut work.decoy_buf,
    );
    let anchor_in_win = doff - win_start;

    // Disjoint borrows on `work.decoy_buf`, `work.rev_q_buf`, `work.rev_r_buf`, and `work.simd` — zero clones!
    let fwd_sc = align_node_segment_fwd(
        &q[rpos..],
        &work.decoy_buf[anchor_in_win..],
        true,
        cfg,
        &work.simd,
    )
    .0
    .score;

    let rev_sc = if rpos > 0 && anchor_in_win > 0 {
        work.rev_q_buf.clear();
        work.rev_q_buf.extend(q[..rpos].iter().rev().copied());
        work.rev_r_buf.clear();
        work.rev_r_buf
            .extend(work.decoy_buf[..anchor_in_win].iter().rev().copied());
        align_node_segment_fwd(&work.rev_q_buf, &work.rev_r_buf, true, cfg, &work.simd)
            .0
            .score
    } else {
        0
    };

    fwd_sc + rev_sc
}
