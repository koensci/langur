//! Memory-bounded sparse per-worker `QuantShard` and deterministic wave-merged
//! `SampleCounts` (zero dense 650k-node allocations per thread/chunk).

use crate::align::paired::is_mate_compatible_with_tx;
use crate::align::{PairedAlignment, ReadAlignment};
use crate::index::MmapIndex;
use crate::quant::bias::{BiasAccum, BiasModel, FeatureBias};
use crate::quant::fld::FragLenAccum;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

/// Per-read bias observation mode:
/// - `None`: `--biascorrect` disabled
/// - `Wave0`: exact `(hexamer, weight)` list accumulated during Wave 0 and collapsed immediately after Wave 0
/// - `Fast`: O(1) locked primer weight scalar + GC bin used on Waves 1..N
#[derive(Copy, Clone)]
pub enum ReadBiasObs<'a> {
    None,
    Wave0 {
        hexamers: &'a [u16],
        gc_bin: Option<usize>,
    },
    Fast {
        primer_w: Option<f64>,
        gc_bin: Option<usize>,
    },
}

/// Observed count + inline compact bias summary for a single node, edge, or equivalence class
/// (zero `Box`, zero `SmallVec` spills on Waves 1..N).
#[derive(Clone, Debug, Default)]
pub struct FeatureCount {
    pub raw: f64,
    pub bias: FeatureBias,
}

impl FeatureCount {
    #[inline]
    pub fn add(&mut self, weight: f64, bias_obs: ReadBiasObs<'_>) {
        self.raw += weight;
        match bias_obs {
            ReadBiasObs::None => {}
            ReadBiasObs::Wave0 { hexamers, gc_bin } => {
                self.bias.record_wave0(hexamers, gc_bin, weight as f32);
            }
            ReadBiasObs::Fast { primer_w, gc_bin } => {
                self.bias.record_fast(primer_w, gc_bin, weight);
            }
        }
    }

    #[inline]
    pub fn merge_from(&mut self, other: &Self) {
        self.raw += other.raw;
        self.bias.merge_from(&other.bias);
    }

    #[inline]
    pub fn merge_owned(&mut self, other: Self) {
        self.raw += other.raw;
        self.bias.merge_owned(other.bias);
    }

    #[inline]
    pub fn primer_adjusted(&self, model: &BiasModel, bias_correct: bool) -> f64 {
        if !bias_correct {
            return self.raw;
        }
        self.raw * self.bias.primer_factor(model)
    }

    #[inline]
    pub fn joint_adjusted(&self, model: &BiasModel, bias_correct: bool) -> f64 {
        if !bias_correct {
            return self.raw;
        }
        self.raw * self.bias.joint_factor(model)
    }
}

/// Per-junction count + MAJIQ-style position-diversity / PCR-stack tracker.
#[derive(Clone, Debug, Default)]
pub struct JunctionStats {
    pub count: FeatureCount,
    /// 64-bit bitmask of distinct left-overhang offsets (`0..=63`) observed on this junction.
    pub overhang_mask: u64,
    /// Exact per-overhang read counts (`(overhang_clamped_u8, count)`), kept inline for `<= 4` positions.
    pub pos_counts: SmallVec<[(u8, u32); 4]>,
    /// Count of raw reads (before VBEM/bias scaling) supporting this junction.
    pub raw_reads: u32,
}

impl JunctionStats {
    #[inline]
    pub fn record_overhang(&mut self, left_overhang: u32) {
        let pos = left_overhang.min(63) as u8;
        self.overhang_mask |= 1u64 << pos;
        self.raw_reads = self.raw_reads.saturating_add(1);
        if let Some(entry) = self.pos_counts.iter_mut().find(|e| e.0 == pos) {
            entry.1 = entry.1.saturating_add(1);
        } else {
            self.pos_counts.push((pos, 1));
        }
    }

    pub fn merge_from(&mut self, other: &Self) {
        self.count.merge_from(&other.count);
        self.overhang_mask |= other.overhang_mask;
        self.raw_reads = self.raw_reads.saturating_add(other.raw_reads);
        for &(pos, cnt) in &other.pos_counts {
            if let Some(entry) = self.pos_counts.iter_mut().find(|e| e.0 == pos) {
                entry.1 = entry.1.saturating_add(cnt);
            } else {
                self.pos_counts.push((pos, cnt));
            }
        }
    }

    pub fn merge_owned(&mut self, other: Self) {
        self.count.merge_owned(other.count);
        self.overhang_mask |= other.overhang_mask;
        self.raw_reads = self.raw_reads.saturating_add(other.raw_reads);
        for &(pos, cnt) in &other.pos_counts {
            if let Some(entry) = self.pos_counts.iter_mut().find(|e| e.0 == pos) {
                entry.1 = entry.1.saturating_add(cnt);
            } else {
                self.pos_counts.push((pos, cnt));
            }
        }
    }

    #[inline]
    pub fn distinct_positions(&self) -> u32 {
        self.overhang_mask.count_ones()
    }

    #[inline]
    pub fn max_stack(&self) -> u32 {
        self.pos_counts.iter().map(|&(_, c)| c).max().unwrap_or(0)
    }
}

/// One candidate alignment inside a multi-mapping read/pair group.
#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct MultiAlignCandidate {
    pub gene_id: u32,
    pub path_5p: SmallVec<[u16; 4]>,
    pub path_3p: SmallVec<[u16; 4]>, // empty for single-end reads
}

/// Deterministic key representing a set of equivalent multi-mapping alignment loci.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct MultiGroupKey {
    pub candidates: SmallVec<[MultiAlignCandidate; 3]>,
}

/// Summary mapping counters emitted to `.map.gz`.
#[derive(Clone, Debug, Default)]
pub struct SampleMapStats {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub multimap_reads: u64,
    pub decoy_filtered_reads: u64,
    pub read_len_sum: u64,
    pub read_len_cnt: u64,
}

impl SampleMapStats {
    pub fn merge_from(&mut self, other: &Self) {
        self.total_reads += other.total_reads;
        self.mapped_reads += other.mapped_reads;
        self.multimap_reads += other.multimap_reads;
        self.decoy_filtered_reads += other.decoy_filtered_reads;
        self.read_len_sum += other.read_len_sum;
        self.read_len_cnt += other.read_len_cnt;
    }

    #[inline]
    pub fn mean_read_len(&self) -> f64 {
        if self.read_len_cnt == 0 {
            50.0
        } else {
            (self.read_len_sum as f64) / (self.read_len_cnt as f64)
        }
    }
}

/// Sparse per-worker shard accumulated during one parallel batch of reads.
/// Never allocates dense arrays over all 650k genomic nodes!
#[derive(Clone, Debug, Default)]
pub struct QuantShard {
    pub stats: SampleMapStats,
    pub fld: FragLenAccum,
    pub bias: BiasAccum,
    /// Uniquely-mapped single-node body counts: `(gene_id, node_id)` -> count.
    pub node_counts: FxHashMap<(u32, u16), FeatureCount>,
    /// Uniquely-mapped edge / splice-junction counts: `(gene_id, src_node, dst_node)` -> stats.
    pub edge_counts: FxHashMap<(u32, u16, u16), JunctionStats>,
    /// Transcript-level equivalence classes (sorted global `tx_idx`s) -> count.
    pub ec_counts: FxHashMap<SmallVec<[u32; 4]>, FeatureCount>,
    /// Multi-mapping read groups (resolved after Stage 1 VBEM and routed into Stage 2).
    pub multi_groups: FxHashMap<MultiGroupKey, FeatureCount>,
}

impl QuantShard {
    /// Lightweight shard for Waves 1..N (skips allocating 64 KB hexamer arrays in `BiasAccum`).
    pub fn new_gc_only() -> Self {
        Self {
            bias: BiasAccum::new_gc_only(),
            ..Self::default()
        }
    }

    /// Pre-sized shard for a chunk of `chunk_reads` reads/pairs (eliminates incremental FxHashMap rehashes).
    pub fn with_chunk_hint(gc_only: bool, chunk_reads: usize) -> Self {
        let est = (chunk_reads / 2).clamp(64, 16384);
        Self {
            stats: SampleMapStats::default(),
            fld: FragLenAccum::default(),
            bias: if gc_only {
                BiasAccum::new_gc_only()
            } else {
                BiasAccum::default()
            },
            node_counts: FxHashMap::with_capacity_and_hasher(est, Default::default()),
            edge_counts: FxHashMap::with_capacity_and_hasher(est, Default::default()),
            ec_counts: FxHashMap::with_capacity_and_hasher(est, Default::default()),
            multi_groups: FxHashMap::with_capacity_and_hasher(est / 8, Default::default()),
        }
    }

    pub fn clear(&mut self) {
        self.stats = SampleMapStats::default();
        self.fld = FragLenAccum::default();
        self.bias = BiasAccum::default();
        self.node_counts.clear();
        self.edge_counts.clear();
        self.ec_counts.clear();
        self.multi_groups.clear();
    }

    /// Collapse Wave 0 exact hexamer vectors into scalar `(primer_weighted, primer_den)`
    /// using the deterministically finalized Wave 0 `primer_weight`, freeing all hexamer heap memory.
    pub fn collapse_wave0_hexamers(&mut self, primer_weight: &[f64]) {
        for fc in self.node_counts.values_mut() {
            fc.bias.collapse_wave0(primer_weight);
        }
        for jst in self.edge_counts.values_mut() {
            jst.count.bias.collapse_wave0(primer_weight);
        }
        for fc in self.ec_counts.values_mut() {
            fc.bias.collapse_wave0(primer_weight);
        }
        for fc in self.multi_groups.values_mut() {
            fc.bias.collapse_wave0(primer_weight);
        }
    }

    /// Record a uniquely mapped single-end read alignment into this sparse shard.
    pub fn record_single_unique(
        &mut self,
        idx: &MmapIndex,
        aln: &ReadAlignment,
        bias_obs: ReadBiasObs<'_>,
    ) {
        record_single_path_into_nodes_edges(
            aln,
            1.0,
            bias_obs,
            &mut self.node_counts,
            &mut self.edge_counts,
        );

        let compat = compatible_txs_for_single(idx, aln);
        if !compat.is_empty() {
            self.ec_counts.entry(compat).or_default().add(1.0, bias_obs);
        }
    }

    /// Record a uniquely mapped concordant paired-end alignment into this sparse shard
    /// (reusing `pa.compat_txs` already computed during fragment-length validation).
    pub fn record_paired_unique(
        &mut self,
        pa: &PairedAlignment,
        bias_obs: ReadBiasObs<'_>,
    ) {
        if let Some(fl) = pa.frag_len {
            self.fld.record(fl);
        }

        // Record edges/nodes from both mates (avoiding double-counting an identical overlapping segment)
        record_paired_paths_into_nodes_edges(
            &pa.mate_5p,
            &pa.mate_3p,
            1.0,
            bias_obs,
            &mut self.node_counts,
            &mut self.edge_counts,
        );

        if !pa.compat_txs.is_empty() {
            if let Some(entry) = self.ec_counts.get_mut(&pa.compat_txs) {
                entry.add(1.0, bias_obs);
            } else {
                let mut fc = FeatureCount::default();
                fc.add(1.0, bias_obs);
                self.ec_counts.insert(pa.compat_txs.clone(), fc);
            }
        }
    }

    /// Record a multi-mapping single-end read (`alns.len() >= 2`).
    pub fn record_single_multi(
        &mut self,
        idx: &MmapIndex,
        alns: &[ReadAlignment],
        bias_obs: ReadBiasObs<'_>,
    ) {
        let mut cands: SmallVec<[MultiAlignCandidate; 3]> = SmallVec::with_capacity(alns.len());
        let mut union_txs: SmallVec<[u32; 4]> = SmallVec::new();
        let mut all_have_anno = true;

        for a in alns {
            cands.push(MultiAlignCandidate {
                gene_id: a.gene_id,
                path_5p: a.node_path(),
                path_3p: SmallVec::new(),
            });
            let txs = compatible_txs_for_single(idx, a);
            if txs.is_empty() {
                all_have_anno = false;
            } else {
                union_txs.extend_from_slice(&txs);
            }
        }
        cands.sort_unstable();
        cands.dedup();

        self.multi_groups
            .entry(MultiGroupKey { candidates: cands })
            .or_default()
            .add(1.0, bias_obs);

        if all_have_anno && !union_txs.is_empty() {
            union_txs.sort_unstable();
            union_txs.dedup();
            self.ec_counts
                .entry(union_txs)
                .or_default()
                .add(1.0, bias_obs);
        }
    }

    /// Record a multi-mapping paired-end read pair (`pairs.len() >= 2`),
    /// reusing precomputed `p.compat_txs` with zero redundant transcript scans.
    pub fn record_paired_multi(
        &mut self,
        pairs: &[PairedAlignment],
        bias_obs: ReadBiasObs<'_>,
    ) {
        let mut cands: SmallVec<[MultiAlignCandidate; 3]> = SmallVec::with_capacity(pairs.len());
        let mut union_txs: SmallVec<[u32; 4]> = SmallVec::new();
        let mut all_have_anno = true;

        for p in pairs {
            cands.push(MultiAlignCandidate {
                gene_id: p.gene_id,
                path_5p: p.mate_5p.node_path(),
                path_3p: p.mate_3p.node_path(),
            });
            if p.compat_txs.is_empty() {
                all_have_anno = false;
            } else {
                union_txs.extend_from_slice(&p.compat_txs);
            }
        }
        cands.sort_unstable();
        cands.dedup();

        self.multi_groups
            .entry(MultiGroupKey { candidates: cands })
            .or_default()
            .add(1.0, bias_obs);

        if all_have_anno && !union_txs.is_empty() {
            union_txs.sort_unstable();
            union_txs.dedup();
            self.ec_counts
                .entry(union_txs)
                .or_default()
                .add(1.0, bias_obs);
        }
    }

    /// Deterministically merge `other` into `self` by value (zero key cloning, used for parallel
    /// binary-tree wave reduction).
    pub fn merge_from(&mut self, other: QuantShard) {
        self.stats.merge_from(&other.stats);
        self.fld.merge_from(&other.fld);
        self.bias.merge_from(&other.bias);
        for (k, v) in other.node_counts {
            use std::collections::hash_map::Entry;
            match self.node_counts.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        for (k, v) in other.edge_counts {
            use std::collections::hash_map::Entry;
            match self.edge_counts.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        for (k, v) in other.ec_counts {
            use std::collections::hash_map::Entry;
            match self.ec_counts.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        for (k, v) in other.multi_groups {
            use std::collections::hash_map::Entry;
            match self.multi_groups.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
    }
}

/// Global merged sparse counts for one sample (assembled deterministically across waves).
#[derive(Clone, Debug, Default)]
pub struct SampleCounts {
    pub stats: SampleMapStats,
    pub fld: FragLenAccum,
    pub bias: BiasAccum,
    pub node_counts: FxHashMap<(u32, u16), FeatureCount>,
    pub edge_counts: FxHashMap<(u32, u16, u16), JunctionStats>,
    pub ec_counts: FxHashMap<SmallVec<[u32; 4]>, FeatureCount>,
    pub multi_groups: FxHashMap<MultiGroupKey, FeatureCount>,
}

impl SampleCounts {
    /// Pre-size global sample tables from index dimensions so the background merger thread
    /// avoids expensive multi-hundred-thousand-entry rehashes across waves.
    pub fn with_index_capacity(idx: &MmapIndex) -> Self {
        let n_nodes = (idx.nodes().len() / 2).clamp(256, 300_000);
        let n_edges = idx.anno_junctions().len().clamp(256, 300_000);
        let n_ecs = (idx.transcripts().len() * 2).clamp(256, 350_000);
        Self {
            stats: SampleMapStats::default(),
            fld: FragLenAccum::default(),
            bias: BiasAccum::default(),
            node_counts: FxHashMap::with_capacity_and_hasher(n_nodes, Default::default()),
            edge_counts: FxHashMap::with_capacity_and_hasher(n_edges, Default::default()),
            ec_counts: FxHashMap::with_capacity_and_hasher(n_ecs, Default::default()),
            multi_groups: FxHashMap::with_capacity_and_hasher(
                (n_nodes / 8).clamp(64, 65_536),
                Default::default(),
            ),
        }
    }

    pub fn merge_shard(&mut self, shard: QuantShard) {
        self.stats.merge_from(&shard.stats);
        self.fld.merge_from(&shard.fld);
        self.bias.merge_from(&shard.bias);
        for (k, v) in shard.node_counts {
            use std::collections::hash_map::Entry;
            match self.node_counts.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        for (k, v) in shard.edge_counts {
            use std::collections::hash_map::Entry;
            match self.edge_counts.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        for (k, v) in shard.ec_counts {
            use std::collections::hash_map::Entry;
            match self.ec_counts.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
        for (k, v) in shard.multi_groups {
            use std::collections::hash_map::Entry;
            match self.multi_groups.entry(k) {
                Entry::Occupied(mut e) => e.get_mut().merge_owned(v),
                Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
    }
}

/// Record a single mate's node/edge counts + exact junction left-overhang diversity.
fn record_single_path_into_nodes_edges(
    aln: &ReadAlignment,
    weight: f64,
    bias_obs: ReadBiasObs<'_>,
    node_counts: &mut FxHashMap<(u32, u16), FeatureCount>,
    edge_counts: &mut FxHashMap<(u32, u16, u16), JunctionStats>,
) {
    if aln.segs.is_empty() {
        return;
    }
    if aln.segs.len() == 1 {
        node_counts
            .entry((aln.gene_id, aln.segs[0].node_id))
            .or_default()
            .add(weight, bias_obs);
    } else {
        let mut cum_read_pos: u32 = 0;
        for win in aln.segs.windows(2) {
            cum_read_pos += win[0].read_len;
            let u = win[0].node_id;
            let v = win[1].node_id;
            let jst = edge_counts.entry((aln.gene_id, u, v)).or_default();
            jst.count.add(weight, bias_obs);
            jst.record_overhang(cum_read_pos);
        }
    }
}

/// Record paired mates into node/edge counts cleanly (matching Whippet 2's subset check so
/// overlapping mates don't double-count the same node or junction).
fn record_paired_paths_into_nodes_edges(
    m5: &ReadAlignment,
    m3: &ReadAlignment,
    weight: f64,
    bias_obs: ReadBiasObs<'_>,
    node_counts: &mut FxHashMap<(u32, u16), FeatureCount>,
    edge_counts: &mut FxHashMap<(u32, u16, u16), JunctionStats>,
) {
    // If both mates map to the exact same single node, count the pair once (like Whippet2)
    if m5.segs.len() == 1 && m3.segs.len() == 1 && m5.segs[0].node_id == m3.segs[0].node_id {
        node_counts
            .entry((m5.gene_id, m5.segs[0].node_id))
            .or_default()
            .add(weight, bias_obs);
        return;
    }

    // Track edges already counted by mate_5p so an overlapping mate_3p spanning the exact same
    // splice junction doesn't double-count that junction for one fragment!
    let mut seen_edges: SmallVec<[(u16, u16); 4]> = SmallVec::new();
    let mut seen_nodes: SmallVec<[u16; 2]> = SmallVec::new();

    if m5.segs.len() == 1 {
        let n = m5.segs[0].node_id;
        seen_nodes.push(n);
        node_counts
            .entry((m5.gene_id, n))
            .or_default()
            .add(weight, bias_obs);
    } else {
        let mut cum = 0u32;
        for win in m5.segs.windows(2) {
            cum += win[0].read_len;
            let key = (win[0].node_id, win[1].node_id);
            seen_edges.push(key);
            let jst = edge_counts.entry((m5.gene_id, key.0, key.1)).or_default();
            jst.count.add(weight, bias_obs);
            jst.record_overhang(cum);
        }
    }

    if m3.segs.len() == 1 {
        let n = m3.segs[0].node_id;
        // Count mate_3p's node if not already covered by mate_5p's segments
        if !seen_nodes.contains(&n) && !m5.segs.iter().any(|s| s.node_id == n) {
            node_counts
                .entry((m3.gene_id, n))
                .or_default()
                .add(weight, bias_obs);
        }
    } else {
        let mut cum = 0u32;
        for win in m3.segs.windows(2) {
            cum += win[0].read_len;
            let key = (win[0].node_id, win[1].node_id);
            if !seen_edges.contains(&key) {
                let jst = edge_counts.entry((m3.gene_id, key.0, key.1)).or_default();
                jst.count.add(weight, bias_obs);
                jst.record_overhang(cum);
            }
        }
    }
}

/// Find sorted global `tx_idx` values of all annotated transcripts compatible with a single alignment.
pub fn compatible_txs_for_single(idx: &MmapIndex, aln: &ReadAlignment) -> SmallVec<[u32; 4]> {
    let g = &idx.genes()[aln.gene_id as usize];
    let tx_off = g.tx_offset as usize;
    let num_tx = g.num_tx as usize;
    let mut out = SmallVec::new();
    for rel in 0..num_tx {
        let tid = tx_off + rel;
        if is_mate_compatible_with_tx(aln, idx.transcript_nodes(tid)) {
            out.push(tid as u32);
        }
    }
    out
}

/// Find sorted global `tx_idx` values of all annotated transcripts compatible with a paired alignment.
pub fn compatible_txs_for_pair(idx: &MmapIndex, pa: &PairedAlignment) -> SmallVec<[u32; 4]> {
    let g = &idx.genes()[pa.gene_id as usize];
    let tx_off = g.tx_offset as usize;
    let num_tx = g.num_tx as usize;
    let mut out = SmallVec::new();
    for rel in 0..num_tx {
        let tid = tx_off + rel;
        let tnodes = idx.transcript_nodes(tid);
        if is_mate_compatible_with_tx(&pa.mate_5p, tnodes)
            && is_mate_compatible_with_tx(&pa.mate_3p, tnodes)
        {
            out.push(tid as u32);
        }
    }
    out
}

/// Find global `tx_idx` values compatible with a `MultiAlignCandidate` (used during Stage 1 -> Stage 2 routing).
pub fn compatible_txs_for_candidate(idx: &MmapIndex, cand: &MultiAlignCandidate) -> SmallVec<[u32; 4]> {
    let g = &idx.genes()[cand.gene_id as usize];
    let tx_off = g.tx_offset as usize;
    let num_tx = g.num_tx as usize;
    let mut out = SmallVec::new();
    for rel in 0..num_tx {
        let tid = tx_off + rel;
        let tnodes = idx.transcript_nodes(tid);
        let ok5 = is_node_path_compatible(&cand.path_5p, tnodes);
        let ok3 = cand.path_3p.is_empty() || is_node_path_compatible(&cand.path_3p, tnodes);
        if ok5 && ok3 {
            out.push(tid as u32);
        }
    }
    out
}

#[inline]
pub fn is_node_path_compatible(path: &[u16], tx_nodes: &[u16]) -> bool {
    if path.is_empty() {
        return false;
    }
    let Ok(mut prev) = tx_nodes.binary_search(&path[0]) else {
        return false;
    };
    for &nid in &path[1..] {
        let next = prev + 1;
        if next >= tx_nodes.len() || tx_nodes[next] != nid {
            return false;
        }
        prev = next;
    }
    true
}
