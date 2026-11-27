//! Selective genomic decoy extraction (`salmon --validateMappings` style) with
//! `± (k - 1)` CSG boundary masking so exon-boundary overhangs and annotated
//! retained-intron (`RI`) nodes are never captured as decoys, ambiguous (`N`)
//! run splitting so `decoy_seq_2bit` never encodes ambiguous bases, and `rayon`
//! parallelism across unmasked genomic gaps.

use crate::graph::BuiltSpliceGraph;
use crate::index::format::{DecoyRecord, PackedSeedHit};
use crate::index::minimizer::RawPosting;
use crate::seq::{extract_minimizers, is_ambiguous, pack_2bit_append, Minimizer};
use rayon::prelude::*;
use rustc_hash::FxHashSet;

/// Configuration for selective genomic decoy extraction.
#[derive(Clone, Debug)]
pub struct DecoyExtractConfig {
    /// Padding around each CSG-colliding decoy hit cluster when forming decoy contigs (default 250 bp).
    pub decoy_window_pad: usize,
    /// Minimum decoy contig length in bp (defaults to `k`).
    pub min_decoy_len: usize,
    /// Minimum number of distinct CSG-colliding minimizers required within `cluster_span_bp`
    /// to emit a selective genomic decoy window (default `2`, preventing isolated random 17-mer
    /// matches from pulling in 88% of the human genome!).
    pub min_cluster_hits: usize,
    /// Maximum genomic span (bp) between two distinct CSG-colliding minimizers to form a decoy
    /// cluster (default `150` bp).
    pub cluster_span_bp: u32,
}

impl Default for DecoyExtractConfig {
    fn default() -> Self {
        Self {
            decoy_window_pad: 120,
            min_decoy_len: 17,
            min_cluster_hits: 4,
            cluster_span_bp: 80,
        }
    }
}

/// Build sorted, merged, 1-based inclusive masked intervals covering all CSG nodes on a chromosome,
/// each padded by `± (k - 1)` bp (clamped to `[1, chrom_len]`).
pub fn build_masked_csg_intervals(
    chrom_graphs: &[&BuiltSpliceGraph],
    chrom_len: u32,
    k: usize,
) -> Vec<(u32, u32)> {
    if chrom_len == 0 {
        return Vec::new();
    }
    let pad = (k.saturating_sub(1)) as u32;
    let mut intervals: Vec<(u32, u32)> = Vec::new();
    for sg in chrom_graphs {
        for nd in &sg.nodes {
            if nd.len == 0 {
                continue;
            }
            let s = nd.genomic_coord.saturating_sub(pad).max(1);
            let e = nd
                .genomic_coord
                .saturating_add(nd.len - 1)
                .saturating_add(pad)
                .min(chrom_len);
            if s <= e {
                intervals.push((s, e));
            }
        }
    }
    merge_intervals(intervals)
}

/// Build sorted, merged 1-based inclusive gene-span intervals on a chromosome (used to tag
/// whether a decoy contig is intronic vs intergenic).
pub fn build_gene_span_intervals(chrom_graphs: &[&BuiltSpliceGraph]) -> Vec<(u32, u32)> {
    let intervals: Vec<(u32, u32)> = chrom_graphs
        .iter()
        .filter(|sg| sg.genomic_start > 0 && sg.genomic_end >= sg.genomic_start)
        .map(|sg| (sg.genomic_start, sg.genomic_end))
        .collect();
    merge_intervals(intervals)
}

fn merge_intervals(mut intervals: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    if intervals.is_empty() {
        return intervals;
    }
    intervals.sort_unstable();
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(intervals.len());
    let mut cur = intervals[0];
    for &(s, e) in &intervals[1..] {
        if s <= cur.1.saturating_add(1) {
            if e > cur.1 {
                cur.1 = e;
            }
        } else {
            merged.push(cur);
            cur = (s, e);
        }
    }
    merged.push(cur);
    merged
}

/// Compute the complement of sorted merged 1-based inclusive `masked` intervals within `[1, chrom_len]`.
pub fn unmasked_gaps(masked: &[(u32, u32)], chrom_len: u32) -> Vec<(u32, u32)> {
    if chrom_len == 0 {
        return Vec::new();
    }
    let mut gaps = Vec::with_capacity(masked.len() + 1);
    let mut cursor: u32 = 1;
    for &(ms, me) in masked {
        if cursor < ms {
            gaps.push((cursor, ms - 1));
        }
        if me >= chrom_len {
            return gaps;
        }
        cursor = me + 1;
    }
    if cursor <= chrom_len {
        gaps.push((cursor, chrom_len));
    }
    gaps
}

fn overlaps_any(sorted_intervals: &[(u32, u32)], start: u32, end: u32) -> bool {
    let idx = sorted_intervals.partition_point(|&(_, e)| e < start);
    idx < sorted_intervals.len() && sorted_intervals[idx].0 <= end
}

/// Split a slice at ambiguous (`N`/IUPAC) runs, yielding contiguous unambiguous sub-ranges
/// `[start..end)` of length `>= min_len`.
fn split_unambiguous_runs(seq: &[u8], min_len: usize) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut i = 0usize;
    let n = seq.len();
    while i < n {
        while i < n && is_ambiguous(seq[i]) {
            i += 1;
        }
        let start = i;
        while i < n && !is_ambiguous(seq[i]) {
            i += 1;
        }
        if i - start >= min_len {
            runs.push((start, i));
        }
    }
    runs
}

/// Per-gap decoy candidate discovered in parallel before deterministic assembly.
struct PendingDecoy {
    genomic_start: u32,
    genomic_end: u32,
    flags: u32,
    /// Minimizers on this clean (ambiguous-free) contig that collide with `csg_hash_set`.
    hits: Vec<Minimizer>,
}

/// Extract selective genomic decoys for a single chromosome outside its `± (k - 1)` padded
/// CSG node intervals, retaining only windows containing at least one CSG-colliding minimizer.
///
/// - Splits decoy contigs cleanly at `N`/ambiguous runs so `decoy_seq_2bit` never encodes `'N'`.
/// - Parallelizes gap scanning across cores with `rayon` while preserving deterministic output order.
#[allow(clippy::too_many_arguments)]
pub fn extract_chrom_decoys(
    chrom_id: u32,
    chrom_seq: &[u8],
    chrom_graphs: &[&BuiltSpliceGraph],
    csg_hash_set: &FxHashSet<u64>,
    k: usize,
    w: usize,
    hash_bits: u32,
    cfg: &DecoyExtractConfig,
    decoys_out: &mut Vec<DecoyRecord>,
    decoy_seq_2bit: &mut Vec<u8>,
    decoy_postings: &mut Vec<RawPosting>,
) {
    let chrom_len = chrom_seq.len() as u32;
    let masked = build_masked_csg_intervals(chrom_graphs, chrom_len, k);
    let gene_spans = build_gene_span_intervals(chrom_graphs);
    let gaps = unmasked_gaps(&masked, chrom_len);
    let min_len = cfg.min_decoy_len.max(k);

    // Parallel scan across unmasked genomic gaps (splitting large intergenic/centromeric gaps
    // would also work automatically via chunks if needed, and per-gap parallel preserves exact order).
    let per_gap_pending: Vec<Vec<PendingDecoy>> = gaps
        .par_iter()
        .map(|&(gap_start, gap_end)| {
            let gap_len = (gap_end - gap_start + 1) as usize;
            if gap_len < min_len {
                return Vec::new();
            }
            let g0 = (gap_start - 1) as usize;
            let g1 = gap_end as usize;
            let gap_slice = &chrom_seq[g0..g1];

            let mut min_buf: Vec<Minimizer> = Vec::new();
            extract_minimizers(gap_slice, k, w, hash_bits, &mut min_buf);

            // Collect all non-repetitive CSG-colliding minimizers inside this gap
            let mut raw_hits: Vec<(u32, u64)> = Vec::new();
            for m in &min_buf {
                if csg_hash_set.contains(&m.hash) {
                    raw_hits.push((gap_start + m.pos, m.hash));
                }
            }
            if raw_hits.len() < cfg.min_cluster_hits.max(1) {
                return Vec::new();
            }

            // Require at least `min_cluster_hits` (default 2) DISTINCT CSG-colliding minimizers
            // within `cluster_span_bp` (default 150 bp) before emitting a selective decoy window!
            let mut hit_windows: Vec<(u32, u32)> = Vec::new();
            let req = cfg.min_cluster_hits.max(1);
            let span = cfg.cluster_span_bp;
            if req == 1 {
                for &(abs_start, _) in &raw_hits {
                    let abs_end = abs_start + (k as u32) - 1;
                    let ws = abs_start
                        .saturating_sub(cfg.decoy_window_pad as u32)
                        .max(gap_start);
                    let we = abs_end
                        .saturating_add(cfg.decoy_window_pad as u32)
                        .min(gap_end);
                    hit_windows.push((ws, we));
                }
            } else {
                let mut left = 0usize;
                for right in 0..raw_hits.len() {
                    let (r_pos, r_hash) = raw_hits[right];
                    while left < right && r_pos - raw_hits[left].0 > span {
                        left += 1;
                    }
                    // Check if [left..=right] contains >= `req` distinct hashes
                    let mut distinct = 1usize;
                    if right > left {
                        for i in left..right {
                            let hi = raw_hits[i].1;
                            if hi != r_hash && !raw_hits[left..i].iter().any(|&(_, h)| h == hi) {
                                distinct += 1;
                                if distinct >= req {
                                    break;
                                }
                            }
                        }
                    }
                    if distinct >= req {
                        let abs_start = raw_hits[left].0;
                        let abs_end = r_pos + (k as u32) - 1;
                        let ws = abs_start
                            .saturating_sub(cfg.decoy_window_pad as u32)
                            .max(gap_start);
                        let we = abs_end
                            .saturating_add(cfg.decoy_window_pad as u32)
                            .min(gap_end);
                        hit_windows.push((ws, we));
                    }
                }
            }
            if hit_windows.is_empty() {
                return Vec::new();
            }

            let merged_windows = merge_intervals(hit_windows);
            let mut local_decoys = Vec::new();

            for (ws, we) in merged_windows {
                let w0 = (ws - 1) as usize;
                let w1 = we as usize;
                let win_slice = &chrom_seq[w0..w1];

                // Split window at ambiguous ('N'/IUPAC) runs so 2-bit decoy never encodes 'N'
                for (rs, re) in split_unambiguous_runs(win_slice, min_len) {
                    let clean_slice = &win_slice[rs..re];
                    extract_minimizers(clean_slice, k, w, hash_bits, &mut min_buf);
                    let hits: Vec<Minimizer> = min_buf
                        .iter()
                        .copied()
                        .filter(|m| csg_hash_set.contains(&m.hash))
                        .collect();
                    if hits.len() < req {
                        continue;
                    }
                    if req >= 2 {
                        let h0 = hits[0].hash;
                        if !hits[1..].iter().any(|m| m.hash != h0) {
                            continue;
                        }
                    }
                    let abs_s = ws + (rs as u32);
                    let abs_e = ws + (re as u32) - 1;
                    let is_intronic = overlaps_any(&gene_spans, abs_s, abs_e);
                    let flags = if is_intronic {
                        DecoyRecord::FLAG_INTRONIC
                    } else {
                        DecoyRecord::FLAG_INTERGENIC
                    };
                    local_decoys.push(PendingDecoy {
                        genomic_start: abs_s,
                        genomic_end: abs_e,
                        flags,
                        hits,
                    });
                }
            }

            local_decoys
        })
        .collect();

    // Deterministic sequential assembly in chromosome coordinate order
    for group in per_gap_pending {
        for pd in group {
            let decoy_id = decoys_out.len() as u32;
            if decoy_id >= PackedSeedHit::MAX_DECOY_ID {
                return;
            }
            let s0 = (pd.genomic_start - 1) as usize;
            let e0 = pd.genomic_end as usize;
            let clean_slice = &chrom_seq[s0..e0];
            let seq_byte_offset = decoy_seq_2bit.len() as u64;
            pack_2bit_append(clean_slice, decoy_seq_2bit);

            decoys_out.push(DecoyRecord {
                seq_byte_offset,
                chrom_id,
                genomic_start: pd.genomic_start,
                seq_len: clean_slice.len() as u32,
                flags: pd.flags,
            });
            for m in pd.hits {
                decoy_postings.push(RawPosting {
                    hash: m.hash,
                    canon_kmer: m.canon_kmer,
                    hit: PackedSeedHit::decoy(decoy_id, m.pos, m.is_rc),
                });
            }
        }
    }
}

/// Extract decoys from an explicit external decoy FASTA sequence (retaining only CSG-colliding hits),
/// splitting cleanly at `N`/ambiguous runs so `decoy_seq_2bit` never encodes ambiguous bases.
#[allow(clippy::too_many_arguments)]
pub fn extract_external_decoy_seq(
    seq: &[u8],
    csg_hash_set: &FxHashSet<u64>,
    k: usize,
    w: usize,
    hash_bits: u32,
    decoys_out: &mut Vec<DecoyRecord>,
    decoy_seq_2bit: &mut Vec<u8>,
    decoy_postings: &mut Vec<RawPosting>,
) {
    if seq.len() < k {
        return;
    }
    let mut min_buf: Vec<Minimizer> = Vec::new();
    for (rs, re) in split_unambiguous_runs(seq, k) {
        let clean_slice = &seq[rs..re];
        extract_minimizers(clean_slice, k, w, hash_bits, &mut min_buf);
        let matched: Vec<Minimizer> = min_buf
            .iter()
            .copied()
            .filter(|m| csg_hash_set.contains(&m.hash))
            .collect();
        if matched.is_empty() {
            continue;
        }
        let decoy_id = decoys_out.len() as u32;
        if decoy_id >= PackedSeedHit::MAX_DECOY_ID {
            return;
        }
        let seq_byte_offset = decoy_seq_2bit.len() as u64;
        pack_2bit_append(clean_slice, decoy_seq_2bit);
        decoys_out.push(DecoyRecord {
            seq_byte_offset,
            chrom_id: u32::MAX,
            genomic_start: (rs as u32) + 1,
            seq_len: clean_slice.len() as u32,
            flags: DecoyRecord::FLAG_EXTERNAL_FASTA,
        });
        for m in matched {
            decoy_postings.push(RawPosting {
                hash: m.hash,
                canon_kmer: m.canon_kmer,
                hit: PackedSeedHit::decoy(decoy_id, m.pos, m.is_rc),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unmasked_gaps_and_padding() {
        let masked = vec![(10, 25), (30, 40)];
        assert_eq!(unmasked_gaps(&masked, 50), vec![(1, 9), (26, 29), (41, 50)]);
        assert_eq!(unmasked_gaps(&[(1, 50)], 50), Vec::<(u32, u32)>::new());
    }

    #[test]
    fn test_split_unambiguous_runs() {
        let seq = b"ACGTACGTNNNNNNTTGCATGCNACGT";
        assert_eq!(split_unambiguous_runs(seq, 5), vec![(0, 8), (14, 22)]);
    }
}
