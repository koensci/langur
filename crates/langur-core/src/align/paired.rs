//! Paired-end concordant alignment across Unstranded, FR (`--stranded`), and
//! RF (`--stranded-reverse` / TruSeq dUTP) protocols, with O(log deg) exact
//! transcriptomic fragment-length calculation via `tx_node_prefix_len[]`,
//! zero per-pair scratch clones, and strict-inequality paired decoy validation.

use crate::align::extend::{
    extend_gene_candidate, min_valid_score, score_single_decoy, score_tolerance, AlignScratch,
};
use crate::align::seed::{collect_decoy_seeds, collect_seeds};
use crate::align::{
    AlignConfig, PairedAlignOutcome, PairedAlignment, ReadAlignment, StrandMode,
};
use crate::index::MmapIndex;
use crate::seq::{append_normalized, append_revcomp_normalized};
use smallvec::SmallVec;

/// Align a paired-end read pair `(r1_raw, r2_raw)` concordantly to the `.lidx` index
/// with zero per-pair scratch clones or hash-map allocations.
pub fn align_paired_reads(
    idx: &MmapIndex,
    r1_raw: &[u8],
    r2_raw: &[u8],
    cfg: &AlignConfig,
    scratch: &mut AlignScratch,
) -> PairedAlignOutcome {
    let k = idx.header().minimizer_k as usize;
    if r1_raw.len() < k || r2_raw.len() < k {
        return PairedAlignOutcome::Unmapped;
    }

    scratch.fwd_norm.clear();
    let r1_has_n = append_normalized(r1_raw, &mut scratch.fwd_norm);
    scratch.r2_fwd_norm.clear();
    let r2_has_n = append_normalized(r2_raw, &mut scratch.r2_fwd_norm);

    let (r1_allow_sense, r1_allow_anti, r2_allow_sense, r2_allow_anti) = match cfg.strand_mode {
        StrandMode::ForwardFR => {
            if cfg.pair_is_rc {
                (true, false, false, true)
            } else {
                (true, false, true, false)
            }
        }
        StrandMode::ReverseRF => {
            if cfg.pair_is_rc {
                (false, true, true, false)
            } else {
                (false, true, false, true)
            }
        }
        StrandMode::Unstranded => (true, true, true, true),
    };

    collect_seeds(
        idx,
        &scratch.fwd_norm,
        r1_allow_sense,
        r1_allow_anti,
        cfg.max_gene_candidates,
        &mut scratch.seed,
    );
    if scratch.seed.gene_candidates.is_empty() {
        return PairedAlignOutcome::Unmapped;
    }

    collect_seeds(
        idx,
        &scratch.r2_fwd_norm,
        r2_allow_sense,
        r2_allow_anti,
        cfg.max_gene_candidates,
        &mut scratch.seed_r2,
    );
    if scratch.seed_r2.gene_candidates.is_empty() {
        return PairedAlignOutcome::Unmapped;
    }

    // For the common case (<= 12 R2 candidates), linear-scan in L1 cache instead of building FxHashMap
    let use_r2_map = scratch.seed_r2.gene_candidates.len() > 12;
    if use_r2_map {
        scratch.r2_lookup.clear();
        for (i, c) in scratch.seed_r2.gene_candidates.iter().enumerate() {
            scratch.r2_lookup.insert((c.gene_id, c.read_is_sense), i);
        }
    }

    let r1_len = scratch.fwd_norm.len();
    let r2_len = scratch.r2_fwd_norm.len();
    let max_r2 = (r2_len as i32) * cfg.match_score;
    let max_possible_pair = ((r1_len + r2_len) as i32) * cfg.match_score;
    let min_r1 = min_valid_score(r1_len, cfg);
    let min_r2 = min_valid_score(r2_len, cfg);
    let pair_tol = score_tolerance(r1_len + r2_len, cfg);

    let mut pairs: SmallVec<[PairedAlignment; 2]> = SmallVec::new();
    let mut best_pair_score = i32::MIN;
    let mut best_c1_hits = 0u32;
    let mut r1_rc_ready = false;
    let mut r2_rc_ready = false;

    for c1 in &scratch.seed.gene_candidates {
        // Once a multi-hit (>= 4 minimizers on R1) candidate achieved the maximum possible pair score,
        // single-minimizer noise R1 candidates (hit_count == 1) cannot be exact-match paralogs.
        if best_pair_score >= max_possible_pair && best_c1_hits >= 4 && c1.hit_count <= 1 {
            break;
        }
        let expected_r2_sense = if cfg.pair_is_rc {
            !c1.read_is_sense
        } else {
            c1.read_is_sense
        };
        let i2 = if use_r2_map {
            let Some(&idx2) = scratch.r2_lookup.get(&(c1.gene_id, expected_r2_sense)) else {
                continue;
            };
            idx2
        } else {
            let Some(idx2) = scratch
                .seed_r2
                .gene_candidates
                .iter()
                .position(|c| c.gene_id == c1.gene_id && c.read_is_sense == expected_r2_sense)
            else {
                continue;
            };
            idx2
        };
        let c2 = &scratch.seed_r2.gene_candidates[i2];
        // If we already have a near-exact pair (<= 2 mismatches across 200 bp) supported by >= 5 R1
        // minimizers, skip trailing noise candidates where BOTH mates have only 1 minimizer hit.
        if best_pair_score >= max_possible_pair - 12
            && best_c1_hits >= 5
            && c1.hit_count <= 1
            && c2.hit_count <= 1
        {
            continue;
        }

        // Exact mathematical guard: once `best_pair_score` is high enough that even a single
        // 1-bp gap (`gap_open = -6`) cannot reach `best_pair_score - pair_tol`, gapped DP on
        // trailing candidates is mathematically guaranteed not to survive retention!
        let mut cand_cfg = cfg.clone();
        if best_pair_score != i32::MIN
            && max_possible_pair + cfg.gap_open < best_pair_score - pair_tol
        {
            cand_cfg.enable_gapped = false;
        }

        let q1: &[u8] = if c1.read_is_sense {
            &scratch.fwd_norm
        } else {
            if !r1_rc_ready {
                scratch.rc_norm.clear();
                append_revcomp_normalized(r1_raw, &mut scratch.rc_norm);
                r1_rc_ready = true;
            }
            &scratch.rc_norm
        };
        let Some(aln1) = extend_gene_candidate(idx, c1, q1, r1_has_n, &cand_cfg, &mut scratch.work) else {
            continue;
        };
        if aln1.score < min_r1 {
            continue;
        }
        // Exact mathematical upper-bound pruning: if even a perfect R2 alignment cannot reach
        // `best_pair_score - pair_tol`, skip extending R2 completely (zero output change).
        if best_pair_score != i32::MIN && aln1.score + max_r2 < best_pair_score - pair_tol {
            continue;
        }

        let q2: &[u8] = if c2.read_is_sense {
            &scratch.r2_fwd_norm
        } else {
            if !r2_rc_ready {
                scratch.r2_rc_norm.clear();
                append_revcomp_normalized(r2_raw, &mut scratch.r2_rc_norm);
                r2_rc_ready = true;
            }
            &scratch.r2_rc_norm
        };
        let Some(aln2) = extend_gene_candidate(idx, c2, q2, r2_has_n, &cand_cfg, &mut scratch.work) else {
            continue;
        };
        if aln2.score < min_r2 {
            continue;
        }

        let (mate_5p, mate_3p, r1_is_5p) = if c1.read_is_sense {
            (aln1, aln2, true)
        } else {
            (aln2, aln1, false)
        };

        let Some((frag_len, compat_txs)) =
            compute_pair_fragment_length(idx, c1.gene_id, &mate_5p, &mate_3p, cfg.max_frag_len)
        else {
            continue;
        };

        let pair_score = mate_5p.score + mate_3p.score;
        if pair_score > best_pair_score
            || (pair_score == best_pair_score && c1.hit_count > best_c1_hits)
        {
            best_c1_hits = c1.hit_count;
        }
        let pa = PairedAlignment {
            gene_id: c1.gene_id,
            mate_5p,
            mate_3p,
            r1_is_5p,
            score: pair_score,
            frag_len: Some(frag_len),
            compat_txs,
        };
        push_scored_pair(
            &mut pairs,
            &mut best_pair_score,
            pa,
            pair_tol,
            cfg.max_multi_hits,
        );
    }

    if pairs.is_empty() {
        return PairedAlignOutcome::Unmapped;
    }

    // Paired selective decoy validation (strict `best_decoy_pair_score > best_pair_score`).
    // Mathematically impossible for any decoy pair to exceed `max_possible_pair`, so defer decoy
    // seed collection completely until `best_pair_score < max_possible_pair` and both mates saw
    // at least one decoy posting!
    if best_pair_score < max_possible_pair
        && scratch.seed.has_any_decoy_posting
        && scratch.seed_r2.has_any_decoy_posting
    {
        collect_decoy_seeds(
            idx,
            r1_len as u32,
            cfg.max_decoy_candidates,
            &mut scratch.seed,
        );
        if !scratch.seed.decoy_candidates.is_empty() {
            collect_decoy_seeds(
                idx,
                r2_len as u32,
                cfg.max_decoy_candidates,
                &mut scratch.seed_r2,
            );
            if !scratch.seed_r2.decoy_candidates.is_empty() {
                let mut best_decoy_pair = i32::MIN;
                for d1 in &scratch.seed.decoy_candidates {
                    let want_d2_fwd = if cfg.pair_is_rc {
                        !d1.read_is_fwd_on_decoy
                    } else {
                        d1.read_is_fwd_on_decoy
                    };
                    for d2 in &scratch.seed_r2.decoy_candidates {
                        if d2.decoy_id != d1.decoy_id || d2.read_is_fwd_on_decoy != want_d2_fwd {
                            continue;
                        }
                        let dist = d1.decoy_offset.abs_diff(d2.decoy_offset);
                        if dist > cfg.max_frag_len {
                            continue;
                        }
                        let q1: &[u8] = if d1.read_is_fwd_on_decoy {
                            &scratch.fwd_norm
                        } else {
                            if !r1_rc_ready {
                                scratch.rc_norm.clear();
                                append_revcomp_normalized(r1_raw, &mut scratch.rc_norm);
                                r1_rc_ready = true;
                            }
                            &scratch.rc_norm
                        };
                        let q2: &[u8] = if d2.read_is_fwd_on_decoy {
                            &scratch.r2_fwd_norm
                        } else {
                            if !r2_rc_ready {
                                scratch.r2_rc_norm.clear();
                                append_revcomp_normalized(r2_raw, &mut scratch.r2_rc_norm);
                                r2_rc_ready = true;
                            }
                            &scratch.r2_rc_norm
                        };
                        // Exact mathematical guard: a decoy alignment containing any gap has score
                        // <= max_possible_pair + gap_open; if that cannot strictly exceed best_pair_score,
                        // gapped DP on decoys is mathematically superfluous (zero output change).
                        let decoy_gapped = max_possible_pair + cfg.gap_open > best_pair_score;
                        let sc1 = score_single_decoy(
                            idx,
                            *d1,
                            q1,
                            r1_has_n,
                            decoy_gapped,
                            cfg,
                            &mut scratch.work,
                        );
                        // Exact mathematical upper-bound pruning: if even a perfect R2 decoy hit
                        // cannot strictly beat `best_pair_score`, skip scoring d2 (zero output change).
                        if sc1.saturating_add(max_r2) <= best_pair_score {
                            continue;
                        }
                        let sc2 = score_single_decoy(
                            idx,
                            *d2,
                            q2,
                            r2_has_n,
                            decoy_gapped,
                            cfg,
                            &mut scratch.work,
                        );
                        let total = sc1.saturating_add(sc2);
                        if total > best_decoy_pair {
                            best_decoy_pair = total;
                            if best_decoy_pair > best_pair_score {
                                return PairedAlignOutcome::DecoyFiltered;
                            }
                        }
                    }
                }
            }
        }
    }

    PairedAlignOutcome::Mapped(pairs)
}

#[inline]
fn push_scored_pair(
    pairs: &mut SmallVec<[PairedAlignment; 2]>,
    best_score: &mut i32,
    pa: PairedAlignment,
    score_tol: i32,
    max_multi: usize,
) {
    for existing in pairs.iter_mut() {
        if existing.gene_id == pa.gene_id
            && existing.r1_is_5p == pa.r1_is_5p
            && existing.mate_5p.segs.len() == pa.mate_5p.segs.len()
            && existing.mate_3p.segs.len() == pa.mate_3p.segs.len()
            && existing
                .mate_5p
                .segs
                .iter()
                .zip(pa.mate_5p.segs.iter())
                .all(|(a, b)| a.node_id == b.node_id)
            && existing
                .mate_3p
                .segs
                .iter()
                .zip(pa.mate_3p.segs.iter())
                .all(|(a, b)| a.node_id == b.node_id)
        {
            if pa.score > existing.score {
                *existing = pa;
                if existing.score > *best_score {
                    *best_score = existing.score;
                    let cutoff = *best_score - score_tol;
                    pairs.retain(|x| x.score >= cutoff);
                }
            }
            return;
        }
    }

    if pa.score > *best_score {
        *best_score = pa.score;
        let cutoff = *best_score - score_tol;
        pairs.retain(|x| x.score >= cutoff);
        pairs.push(pa);
    } else if pa.score >= *best_score - score_tol && pairs.len() < max_multi {
        pairs.push(pa);
    }
}

/// Check whether a single mate alignment path is strictly compatible with a sorted `tx_nodes` slice.
#[inline]
pub fn is_mate_compatible_with_tx(aln: &ReadAlignment, tx_nodes: &[u16]) -> bool {
    if aln.segs.is_empty() {
        return false;
    }
    let Ok(mut prev_pos) = tx_nodes.binary_search(&aln.segs[0].node_id) else {
        return false;
    };
    for seg in &aln.segs[1..] {
        let next_pos = prev_pos + 1;
        if next_pos >= tx_nodes.len() || tx_nodes[next_pos] != seg.node_id {
            return false;
        }
        prev_pos = next_pos;
    }
    true
}

/// Compute exact transcriptomic fragment length of `(mate_5p, mate_3p)` on transcript `tx_idx`
/// in O(log deg) using precomputed `tx_node_prefix_len[]`.
#[inline]
pub fn exact_tx_fragment_length(
    idx: &MmapIndex,
    tx_idx: usize,
    mate_5p: &ReadAlignment,
    mate_3p: &ReadAlignment,
    max_frag_len: u32,
) -> Option<u32> {
    let tnodes = idx.transcript_nodes(tx_idx);
    if !is_mate_compatible_with_tx(mate_5p, tnodes) || !is_mate_compatible_with_tx(mate_3p, tnodes) {
        return None;
    }
    let prefix = idx.transcript_prefix_lens(tx_idx);
    let idx_5p = tnodes.binary_search(&mate_5p.first_node()).ok()?;
    let idx_3p = tnodes.binary_search(&mate_3p.last_node()).ok()?;
    if idx_3p < idx_5p {
        return None;
    }
    let pos_5p = prefix[idx_5p] + mate_5p.start_offset_in_first_node();
    let pos_3p = prefix[idx_3p] + mate_3p.end_offset_in_last_node();
    if pos_3p <= pos_5p {
        return None;
    }
    let flen = pos_3p - pos_5p;
    if flen <= max_frag_len {
        Some(flen)
    } else {
        None
    }
}

/// Return `(first_idx, last_idx)` in `tx_nodes` iff `aln` is strictly compatible with `tx_nodes`
/// (single binary search per mate per transcript).
#[inline]
fn mate_tx_node_range(aln: &ReadAlignment, tx_nodes: &[u16]) -> Option<(usize, usize)> {
    if aln.segs.is_empty() {
        return None;
    }
    let first_pos = tx_nodes.binary_search(&aln.segs[0].node_id).ok()?;
    let mut prev_pos = first_pos;
    for seg in &aln.segs[1..] {
        let next_pos = prev_pos + 1;
        if next_pos >= tx_nodes.len() || tx_nodes[next_pos] != seg.node_id {
            return None;
        }
        prev_pos = next_pos;
    }
    Some((first_pos, prev_pos))
}

/// Fused fragment-length calculation AND compatible-transcript discovery in a single pass
/// over `gene_id`'s transcripts (eliminates the second per-gene transcript scan in `QuantShard`).
pub fn compute_pair_fragment_length(
    idx: &MmapIndex,
    gene_id: u32,
    mate_5p: &ReadAlignment,
    mate_3p: &ReadAlignment,
    max_frag_len: u32,
) -> Option<(u32, SmallVec<[u32; 4]>)> {
    let g = &idx.genes()[gene_id as usize];
    let tx_off = g.tx_offset as usize;
    let num_tx = g.num_tx as usize;

    let mut compat_txs: SmallVec<[u32; 4]> = SmallVec::new();
    let mut sum: u64 = 0;
    let mut cnt: u64 = 0;
    let start_off_5p = mate_5p.start_offset_in_first_node();
    let end_off_3p = mate_3p.end_offset_in_last_node();

    for rel in 0..num_tx {
        let tid = tx_off + rel;
        let tnodes = idx.transcript_nodes(tid);
        let Some((idx_5p, _)) = mate_tx_node_range(mate_5p, tnodes) else {
            continue;
        };
        let Some((_, idx_3p)) = mate_tx_node_range(mate_3p, tnodes) else {
            continue;
        };
        compat_txs.push(tid as u32);
        if idx_3p >= idx_5p {
            let prefix = idx.transcript_prefix_lens(tid);
            let pos_5p = prefix[idx_5p] + start_off_5p;
            let pos_3p = prefix[idx_3p] + end_off_3p;
            if pos_3p > pos_5p {
                let flen = pos_3p - pos_5p;
                if flen <= max_frag_len {
                    sum += flen as u64;
                    cnt += 1;
                }
            }
        }
    }
    if let Some(avg) = sum.checked_div(cnt) {
        return Some((avg as u32, compat_txs));
    }

    let nodes = idx.gene_nodes(gene_id);
    let n5 = mate_5p.first_node() as usize;
    let n3 = mate_3p.last_node() as usize;
    if n3 < n5 {
        return None;
    }
    let off_5p = nodes[n5].gene_seq_offset + start_off_5p;
    let off_3p = nodes[n3].gene_seq_offset + end_off_3p;
    if off_3p <= off_5p {
        return None;
    }
    let span = off_3p - off_5p;
    if span <= max_frag_len {
        Some((span, compat_txs))
    } else {
        None
    }
}
