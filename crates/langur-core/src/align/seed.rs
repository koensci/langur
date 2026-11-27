//! Minimizer seeding and candidate `(gene_id, sense)` / decoy `(decoy_id, fwd, diag_bucket)` clustering.

use crate::index::{MinimizerLookup, MmapIndex, PackedSeedHit, SeedHitKind};
use crate::seq::{extract_minimizers, Minimizer};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

/// A single seed anchor on a gene in the oriented (gene 5'->3' sense) coordinate system.
#[derive(Copy, Clone, Debug)]
pub enum GeneAnchor {
    CsgNode {
        node_id: u16,
        node_offset: u32,
        /// 0-based start of k-mer in the **oriented** read (already revcomp'd if `!read_is_sense`).
        oriented_read_pos: u32,
    },
    AnnoJunction {
        local_jnc_idx: u32,
        /// Exact 0-based position of the splice junction inside the **oriented** read.
        oriented_jnc_pos: u32,
    },
}

/// Candidate `(gene_id, read_is_sense)` cluster with all supporting seed anchors.
#[derive(Clone, Debug)]
pub struct GeneCandidate {
    pub gene_id: u32,
    pub read_is_sense: bool,
    pub hit_count: u32,
    pub anchors: SmallVec<[GeneAnchor; 8]>,
}

/// Candidate decoy alignment anchor clustered by `(decoy_id, read_is_fwd_on_decoy, diag_bucket)`.
#[derive(Copy, Clone, Debug)]
pub struct DecoyCandidate {
    pub decoy_id: u32,
    pub read_is_fwd_on_decoy: bool,
    pub decoy_offset: u32,
    /// 0-based start of k-mer in the read oriented to the decoy's `+` strand.
    pub oriented_read_pos: u32,
    pub hit_count: u32,
}

/// Reusable thread-local scratch buffers for seeding (zero per-read heap allocations).
#[derive(Default)]
pub struct SeedScratch {
    pub mins: Vec<Minimizer>,
    pub gene_map: FxHashMap<(u32, bool), usize>,
    /// Keyed by `(decoy_id, read_is_fwd, diag_bucket)` so spurious hits on a different diagonal
    /// of a long decoy contig never shadow a real match downstream.
    pub decoy_map: FxHashMap<(u32, bool, i32), usize>,
    pub gene_candidates: Vec<GeneCandidate>,
    pub decoy_candidates: Vec<DecoyCandidate>,
    /// True iff at least one minimizer from the read had a non-empty decoy posting list.
    pub has_any_decoy_posting: bool,
}

impl SeedScratch {
    pub fn clear(&mut self) {
        self.mins.clear();
        if !self.gene_map.is_empty() {
            self.gene_map.clear();
        }
        self.gene_candidates.clear();
        self.has_any_decoy_posting = false;
    }
}

/// Extract canonical minimizers from `read_seq` and populate `scratch.gene_candidates`
/// (filtered by allowed orientations `allow_sense`, `allow_antisense`).
///
/// Decoy seed collection is intentionally deferred to `collect_decoy_seeds` so reads that achieve
/// `max_possible_score` (~80% of human reads) never touch decoy postings, hash maps, or sorts.
pub fn collect_seeds(
    idx: &MmapIndex,
    read_seq: &[u8],
    allow_sense: bool,
    allow_antisense: bool,
    max_gene_candidates: usize,
    scratch: &mut SeedScratch,
) {
    scratch.clear();
    let h = idx.header();
    let k = h.minimizer_k as usize;
    let w = h.minimizer_w as usize;
    if read_seq.len() < k {
        return;
    }
    let read_len = read_seq.len() as u32;

    extract_minimizers(read_seq, k, w, idx.hash_bits(), &mut scratch.mins);
    if scratch.mins.is_empty() {
        return;
    }

    // Two-tier repetitive minimizer screening:
    // On 60k-gene human annotations (GENCODE v38), poly-A/Alu/LINE 17-mers hit up to 256 CSG
    // genes per minimizer (~19% of minimizers on real RNA-seq). Look up all ~11 minimizers once
    // into a stack SmallVec; if >= 2 low/moderate-copy minimizers (<= 64 CSG postings) exist,
    // skip > 64-copy repetitive postings (eliminating ~95% of Alu/poly-A posting bloat while
    // preserving full fallback sensitivity for genuine high-copy paralog families).
    const LOW_OCC_CSG_CUTOFF: usize = 64;
    let mut lookups: SmallVec<[(Minimizer, MinimizerLookup<'_>); 24]> = SmallVec::new();
    let mut n_low_occ_csg = 0usize;
    for i in 0..scratch.mins.len() {
        let m = scratch.mins[i];
        let Some(lookup) = idx.lookup_minimizer(m.hash) else {
            continue;
        };
        match lookup {
            MinimizerLookup::Inline(_) => {
                n_low_occ_csg += 1;
            }
            MinimizerLookup::Postings { csg, decoys, .. } => {
                if !decoys.is_empty() {
                    scratch.has_any_decoy_posting = true;
                }
                if !csg.is_empty() && csg.len() <= LOW_OCC_CSG_CUTOFF {
                    n_low_occ_csg += 1;
                }
            }
        }
        lookups.push((m, lookup));
    }
    let skip_high_occ = n_low_occ_csg >= 2;

    for (m, lookup) in lookups {
        match lookup {
            MinimizerLookup::Inline(hit) => {
                record_csg_hit(
                    hit,
                    m,
                    k as u32,
                    read_len,
                    allow_sense,
                    allow_antisense,
                    &mut scratch.gene_map,
                    &mut scratch.gene_candidates,
                );
            }
            MinimizerLookup::Postings { csg, .. } => {
                if skip_high_occ && csg.len() > LOW_OCC_CSG_CUTOFF {
                    continue;
                }
                for &hit in csg {
                    record_csg_hit(
                        hit,
                        m,
                        k as u32,
                        read_len,
                        allow_sense,
                        allow_antisense,
                        &mut scratch.gene_map,
                        &mut scratch.gene_candidates,
                    );
                }
            }
        }
    }

    // Sort deterministically by descending hit_count, then ascending (gene_id, !read_is_sense)
    scratch.gene_candidates.sort_unstable_by(|a, b| {
        b.hit_count
            .cmp(&a.hit_count)
            .then_with(|| a.gene_id.cmp(&b.gene_id))
            .then_with(|| b.read_is_sense.cmp(&a.read_is_sense))
    });
    if scratch.gene_candidates.len() > max_gene_candidates {
        scratch.gene_candidates.truncate(max_gene_candidates);
    }
}

/// Lazily collect and sort decoy candidates from the already-extracted `scratch.mins`.
/// Called ONLY when a read/pair has an imperfect CSG best score (`best_score < max_possible_score`)
/// and `scratch.has_any_decoy_posting` is true.
pub fn collect_decoy_seeds(
    idx: &MmapIndex,
    read_len: u32,
    max_decoy_candidates: usize,
    scratch: &mut SeedScratch,
) {
    scratch.decoy_map.clear();
    scratch.decoy_candidates.clear();
    if !scratch.has_any_decoy_posting || scratch.mins.is_empty() {
        return;
    }
    let k = idx.header().minimizer_k as u32;

    for i in 0..scratch.mins.len() {
        let m = scratch.mins[i];
        if let Some(MinimizerLookup::Postings { decoys, .. }) = idx.lookup_minimizer(m.hash) {
            for &dhit in decoys {
                record_decoy_hit(
                    dhit,
                    m,
                    k,
                    read_len,
                    &mut scratch.decoy_map,
                    &mut scratch.decoy_candidates,
                );
            }
        }
    }

    scratch.decoy_candidates.sort_unstable_by(|a, b| {
        b.hit_count
            .cmp(&a.hit_count)
            .then_with(|| a.decoy_id.cmp(&b.decoy_id))
            .then_with(|| b.read_is_fwd_on_decoy.cmp(&a.read_is_fwd_on_decoy))
            .then_with(|| a.decoy_offset.cmp(&b.decoy_offset))
    });
    if scratch.decoy_candidates.len() > max_decoy_candidates {
        scratch.decoy_candidates.truncate(max_decoy_candidates);
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn record_csg_hit(
    hit: PackedSeedHit,
    m: Minimizer,
    k: u32,
    read_len: u32,
    allow_sense: bool,
    allow_antisense: bool,
    gene_map: &mut FxHashMap<(u32, bool), usize>,
    candidates: &mut Vec<GeneCandidate>,
) {
    match hit.kind() {
        SeedHitKind::CsgNode => {
            let (gene_id, node_id, node_offset, hit_rc) = hit.unpack_csg_node();
            let read_is_sense = m.is_rc == hit_rc;
            if (read_is_sense && !allow_sense) || (!read_is_sense && !allow_antisense) {
                return;
            }
            let oriented_read_pos = if read_is_sense {
                m.pos
            } else {
                read_len - (m.pos + k)
            };
            push_gene_anchor(
                gene_id,
                read_is_sense,
                GeneAnchor::CsgNode {
                    node_id,
                    node_offset,
                    oriented_read_pos,
                },
                gene_map,
                candidates,
            );
        }
        SeedHitKind::AnnoJunction => {
            let (gene_id, local_jnc_idx, left_overhang, hit_rc) = hit.unpack_anno_junction();
            let read_is_sense = m.is_rc == hit_rc;
            if (read_is_sense && !allow_sense) || (!read_is_sense && !allow_antisense) {
                return;
            }
            let oriented_pos = if read_is_sense {
                m.pos
            } else {
                read_len - (m.pos + k)
            };
            let oriented_jnc_pos = oriented_pos + (left_overhang as u32);
            if oriented_jnc_pos == 0 || oriented_jnc_pos >= read_len {
                return;
            }
            push_gene_anchor(
                gene_id,
                read_is_sense,
                GeneAnchor::AnnoJunction {
                    local_jnc_idx,
                    oriented_jnc_pos,
                },
                gene_map,
                candidates,
            );
        }
        SeedHitKind::Decoy => {}
    }
}

#[inline]
fn push_gene_anchor(
    gene_id: u32,
    read_is_sense: bool,
    anchor: GeneAnchor,
    gene_map: &mut FxHashMap<(u32, bool), usize>,
    candidates: &mut Vec<GeneCandidate>,
) {
    // Fast L1-cache linear scan for the common case (<= 8 distinct candidate genes per read),
    // eliminating ~75M FxHashMap insert/clear cycles per sample!
    if gene_map.is_empty() {
        for cand in candidates.iter_mut() {
            if cand.gene_id == gene_id && cand.read_is_sense == read_is_sense {
                cand.hit_count += 1;
                if cand.anchors.len() < 16 {
                    cand.anchors.push(anchor);
                }
                return;
            }
        }
        let mut anchors = SmallVec::new();
        anchors.push(anchor);
        candidates.push(GeneCandidate {
            gene_id,
            read_is_sense,
            hit_count: 1,
            anchors,
        });
        if candidates.len() >= 8 {
            for (i, c) in candidates.iter().enumerate() {
                gene_map.insert((c.gene_id, c.read_is_sense), i);
            }
        }
        return;
    }

    if let Some(&cidx) = gene_map.get(&(gene_id, read_is_sense)) {
        let cand = &mut candidates[cidx];
        cand.hit_count += 1;
        if cand.anchors.len() < 16 {
            cand.anchors.push(anchor);
        }
    } else {
        let cidx = candidates.len();
        gene_map.insert((gene_id, read_is_sense), cidx);
        let mut anchors = SmallVec::new();
        anchors.push(anchor);
        candidates.push(GeneCandidate {
            gene_id,
            read_is_sense,
            hit_count: 1,
            anchors,
        });
    }
}

#[inline]
fn record_decoy_hit(
    hit: PackedSeedHit,
    m: Minimizer,
    k: u32,
    read_len: u32,
    decoy_map: &mut FxHashMap<(u32, bool, i32), usize>,
    candidates: &mut Vec<DecoyCandidate>,
) {
    let (decoy_id, decoy_offset, hit_rc) = hit.unpack_decoy();
    let read_is_fwd = m.is_rc == hit_rc;
    let oriented_read_pos = if read_is_fwd {
        m.pos
    } else {
        read_len - (m.pos + k)
    };
    let diag_bucket = ((decoy_offset as i32) - (oriented_read_pos as i32)) >> 4;
    if let Some(&cidx) = decoy_map.get(&(decoy_id, read_is_fwd, diag_bucket)) {
        candidates[cidx].hit_count += 1;
    } else {
        let cidx = candidates.len();
        decoy_map.insert((decoy_id, read_is_fwd, diag_bucket), cidx);
        candidates.push(DecoyCandidate {
            decoy_id,
            read_is_fwd_on_decoy: read_is_fwd,
            decoy_offset,
            oriented_read_pos,
            hit_count: 1,
        });
    }
}
