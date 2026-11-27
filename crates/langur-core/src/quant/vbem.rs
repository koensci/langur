//! Stage 1: Connected-Component Transcript VBEM + SQUAREM acceleration
//! (with `--use-em` MLE fallback, warm-starting for 2-pass GC bias feedback,
//! and deterministic Stage 1 -> Stage 2 multi-mapper routing).

use crate::index::MmapIndex;
use crate::quant::bias::BiasModel;
use crate::quant::shard::{compatible_txs_for_candidate, MultiGroupKey, SampleCounts};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use smallvec::SmallVec;
use statrs::function::gamma::digamma;

#[derive(Clone, Debug)]
pub struct VbemConfig {
    /// Sparse Dirichlet prior `alpha_0` (default `0.01`, matching Salmon VBEM).
    pub vb_prior: f64,
    /// If `true`, run standard MLE EM (`theta_t = c_t / eff_len_t`) instead of VBEM.
    pub use_em: bool,
    /// Maximum SQUAREM iterations per connected component (default `1000`).
    pub max_iters: usize,
    /// Relative convergence tolerance (default `1e-5`).
    pub rel_tol: f64,
}

impl Default for VbemConfig {
    fn default() -> Self {
        Self {
            vb_prior: 0.01,
            use_em: false,
            max_iters: 1000,
            rel_tol: 1e-5,
        }
    }
}

/// Converged Stage 1 transcript & gene abundance estimates.
#[derive(Clone, Debug)]
pub struct VbemResult {
    /// Converged expected read counts per transcript (`len == num_transcripts`).
    pub tx_counts: Vec<f64>,
    /// Converged TPM per transcript (`len == num_transcripts`, sums to 1e6 when non-empty).
    pub tx_tpm: Vec<f64>,
    /// Converged expected read counts per gene (`len == num_genes`).
    pub gene_counts: Vec<f64>,
    /// Converged TPM per gene (`len == num_genes`).
    pub gene_tpm: Vec<f64>,
    /// Effective lengths used in this pass (`len == num_transcripts`).
    pub tx_eff_len: Vec<f64>,
    /// Number of independent connected components solved.
    pub num_components: usize,
    /// Total SQUAREM iterations summed across components.
    pub total_iters: usize,
}

/// Disjoint-set Union-Find for partitioning active transcripts into independent connected components.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}

struct ComponentInput {
    tx_ids: Vec<u32>,
    eff_lens: Vec<f64>,
    uniq_counts: Vec<f64>,
    init_counts: Option<Vec<f64>>,
    ecs: Vec<(SmallVec<[u32; 4]>, f64)>, // local tx indices + EC count
}

/// Solve Stage 1 transcript quantification across independent connected components in parallel.
///
/// - `ecs`: deterministic slice of `(sorted_tx_ids, adjusted_count)`.
/// - `tx_eff_len`: effective length for every transcript (`len == num_transcripts`).
/// - `warm_start`: optional Pass 1 `tx_counts` to accelerate Pass 2 GC-adjusted convergence.
pub fn solve_cc_vbem(
    idx: &MmapIndex,
    ecs: &[(SmallVec<[u32; 4]>, f64)],
    tx_eff_len: Vec<f64>,
    cfg: &VbemConfig,
    warm_start: Option<&[f64]>,
) -> VbemResult {
    let n_tx = idx.transcripts().len();
    let n_genes = idx.genes().len();
    let mut tx_counts = vec![0.0f64; n_tx];

    if n_tx == 0 || ecs.is_empty() {
        return VbemResult {
            tx_counts,
            tx_tpm: vec![0.0; n_tx],
            gene_counts: vec![0.0; n_genes],
            gene_tpm: vec![0.0; n_genes],
            tx_eff_len,
            num_components: 0,
            total_iters: 0,
        };
    }

    // 1. Collect all active transcripts (appearing in at least one positive EC)
    let mut active_tx_map: FxHashMap<u32, usize> = FxHashMap::default();
    let mut active_txs: Vec<u32> = Vec::new();
    for (txs, cnt) in ecs {
        if *cnt <= 0.0 {
            continue;
        }
        for &tid in txs {
            if let std::collections::hash_map::Entry::Vacant(e) = active_tx_map.entry(tid) {
                let lid = active_txs.len();
                e.insert(lid);
                active_txs.push(tid);
            }
        }
    }
    // Sort active transcripts deterministically so Union-Find / local IDs are 100% reproducible
    active_txs.sort_unstable();
    active_tx_map.clear();
    for (lid, &tid) in active_txs.iter().enumerate() {
        active_tx_map.insert(tid, lid);
    }

    let m_active = active_txs.len();
    let mut uf = UnionFind::new(m_active);
    let mut uniq_active = vec![0.0f64; m_active];
    let mut multi_ecs: Vec<(SmallVec<[usize; 4]>, f64)> = Vec::new();

    for (txs, cnt) in ecs {
        if *cnt <= 0.0 || txs.is_empty() {
            continue;
        }
        if txs.len() == 1 {
            let lid = active_tx_map[&txs[0]];
            uniq_active[lid] += *cnt;
        } else {
            let mut local_ids: SmallVec<[usize; 4]> = SmallVec::with_capacity(txs.len());
            for &tid in txs {
                local_ids.push(active_tx_map[&tid]);
            }
            let first = local_ids[0];
            for &other in &local_ids[1..] {
                uf.union(first, other);
            }
            multi_ecs.push((local_ids, *cnt));
        }
    }

    // 2. Partition into connected components
    let mut root_to_comp: FxHashMap<usize, usize> = FxHashMap::default();
    let mut comp_members: Vec<Vec<usize>> = Vec::new();
    for lid in 0..m_active {
        let r = uf.find(lid);
        let cidx = if let Some(&c) = root_to_comp.get(&r) {
            c
        } else {
            let c = comp_members.len();
            root_to_comp.insert(r, c);
            comp_members.push(Vec::new());
            c
        };
        comp_members[cidx].push(lid);
    }

    // Map each active transcript to its (comp_idx, pos_in_comp)
    let mut lid_to_comp_pos = vec![(0usize, 0u32); m_active];
    let mut comp_inputs: Vec<ComponentInput> = Vec::with_capacity(comp_members.len());
    for (cidx, members) in comp_members.into_iter().enumerate() {
        let mut tx_ids = Vec::with_capacity(members.len());
        let mut eff_lens = Vec::with_capacity(members.len());
        let mut uniq_counts = Vec::with_capacity(members.len());
        let mut init_c = warm_start.map(|_| Vec::with_capacity(members.len()));

        for (pos, lid) in members.into_iter().enumerate() {
            lid_to_comp_pos[lid] = (cidx, pos as u32);
            let tid = active_txs[lid];
            tx_ids.push(tid);
            eff_lens.push(tx_eff_len[tid as usize].max(1.0));
            uniq_counts.push(uniq_active[lid]);
            if let (Some(ref mut ic), Some(ws)) = (&mut init_c, warm_start) {
                ic.push(ws[tid as usize]);
            }
        }
        comp_inputs.push(ComponentInput {
            tx_ids,
            eff_lens,
            uniq_counts,
            init_counts: init_c,
            ecs: Vec::new(),
        });
    }

    for (local_ids, cnt) in multi_ecs {
        let cidx = lid_to_comp_pos[local_ids[0]].0;
        let comp_local: SmallVec<[u32; 4]> =
            local_ids.iter().map(|&lid| lid_to_comp_pos[lid].1).collect();
        comp_inputs[cidx].ecs.push((comp_local, cnt));
    }

    // 3. Solve every connected component independently in parallel via `rayon`
    let num_components = comp_inputs.len();
    let solved: Vec<(Vec<u32>, Vec<f64>, usize)> = comp_inputs
        .into_par_iter()
        .map(|comp| {
            let (counts, iters) = solve_single_component(&comp, cfg);
            (comp.tx_ids, counts, iters)
        })
        .collect();

    let mut total_iters = 0usize;
    for (tx_ids, counts, iters) in solved {
        total_iters += iters;
        for (tid, c) in tx_ids.into_iter().zip(counts) {
            tx_counts[tid as usize] = c;
        }
    }

    // 4. Compute isoform TPMs & aggregate to gene counts / gene TPMs
    let mut tx_tpm = vec![0.0f64; n_tx];
    let mut rpk_sum = 0.0f64;
    for t in 0..n_tx {
        if tx_counts[t] > 0.0 {
            let rpk = tx_counts[t] / tx_eff_len[t].max(1.0);
            tx_tpm[t] = rpk;
            rpk_sum += rpk;
        }
    }
    if rpk_sum > 0.0 {
        let scale = 1_000_000.0 / rpk_sum;
        for v in &mut tx_tpm {
            *v *= scale;
        }
    }

    let mut gene_counts = vec![0.0f64; n_genes];
    let mut gene_tpm = vec![0.0f64; n_genes];
    for (t, trec) in idx.transcripts().iter().enumerate() {
        let g = trec.gene_id as usize;
        gene_counts[g] += tx_counts[t];
        gene_tpm[g] += tx_tpm[t];
    }

    VbemResult {
        tx_counts,
        tx_tpm,
        gene_counts,
        gene_tpm,
        tx_eff_len,
        num_components,
        total_iters,
    }
}

/// Solve one independent connected component via VBEM (or MLE EM) + SQUAREM acceleration.
fn solve_single_component(comp: &ComponentInput, cfg: &VbemConfig) -> (Vec<f64>, usize) {
    let m = comp.tx_ids.len();
    if comp.ecs.is_empty() {
        return (comp.uniq_counts.clone(), 0);
    }

    // Initialize counts (warm-start if provided, else uniform split of multi-ECs + unique counts)
    let mut x = if let Some(ref init) = comp.init_counts {
        init.clone()
    } else {
        let mut init = comp.uniq_counts.clone();
        for (members, cnt) in &comp.ecs {
            let share = *cnt / (members.len() as f64);
            for &t in members {
                init[t as usize] += share;
            }
        }
        init
    };

    let mut x1 = vec![0.0f64; m];
    let mut x2 = vec![0.0f64; m];
    let mut x_cand = vec![0.0f64; m];
    let mut x_next = vec![0.0f64; m];
    let mut weights = vec![0.0f64; m];

    let em_step = |src: &[f64], dst: &mut [f64], wbuf: &mut [f64]| {
        // Compute E-step weights: exp(digamma(alpha_0 + c_t)) / eff_len_t (VBEM) or c_t / eff_len_t (MLE EM)
        for i in 0..m {
            let inv_l = 1.0 / comp.eff_lens[i];
            wbuf[i] = if cfg.use_em {
                src[i].max(0.0) * inv_l
            } else {
                let alpha = cfg.vb_prior + src[i].max(0.0);
                if alpha > 1e-12 {
                    digamma(alpha).exp() * inv_l
                } else {
                    0.0
                }
            };
        }

        dst.copy_from_slice(&comp.uniq_counts);
        for (members, cnt) in &comp.ecs {
            let mut denom = 0.0f64;
            for &t in members {
                denom += wbuf[t as usize];
            }
            if denom > 0.0 {
                let scale = *cnt / denom;
                for &t in members {
                    dst[t as usize] += wbuf[t as usize] * scale;
                }
            } else {
                let share = *cnt / (members.len() as f64);
                for &t in members {
                    dst[t as usize] += share;
                }
            }
        }
    };

    let mut iters = 0usize;
    for it in 1..=cfg.max_iters {
        iters = it;
        em_step(&x, &mut x1, &mut weights);
        em_step(&x1, &mut x2, &mut weights);

        // SQUAREM step length: alpha = -sqrt(||r||^2 / ||v||^2), clamped to [-100, -1]
        let mut r_sq = 0.0f64;
        let mut v_sq = 0.0f64;
        for i in 0..m {
            let r = x1[i] - x[i];
            let v = (x2[i] - x1[i]) - r;
            r_sq += r * r;
            v_sq += v * v;
        }

        let converged_inner = max_rel_diff(&x1, &x2) < cfg.rel_tol;
        if converged_inner || v_sq <= 1e-30 {
            x.copy_from_slice(&x2);
            if converged_inner {
                break;
            }
            continue;
        }

        let alpha = -(r_sq / v_sq).sqrt().clamp(1.0, 100.0);
        let mut valid_cand = true;
        for i in 0..m {
            let r = x1[i] - x[i];
            let v = (x2[i] - x1[i]) - r;
            let val = x[i] - 2.0 * alpha * r + alpha * alpha * v;
            if !val.is_finite() {
                valid_cand = false;
                break;
            }
            x_cand[i] = val.max(0.0);
        }

        if valid_cand {
            em_step(&x_cand, &mut x_next, &mut weights);
            // Safeguard: accept SQUAREM step iff it didn't diverge wildly from x2 mass
            let sum_cand: f64 = x_next.iter().sum();
            let sum_x2: f64 = x2.iter().sum();
            if sum_cand.is_finite() && (sum_cand - sum_x2).abs() <= 1e-3 * sum_x2.max(1.0) {
                let diff = max_rel_diff(&x, &x_next);
                x.copy_from_slice(&x_next);
                if diff < cfg.rel_tol {
                    break;
                }
                continue;
            }
        }

        let diff = max_rel_diff(&x, &x2);
        x.copy_from_slice(&x2);
        if diff < cfg.rel_tol {
            break;
        }
    }

    (x, iters)
}

#[inline]
fn max_rel_diff(a: &[f64], b: &[f64]) -> f64 {
    let mut md = 0.0f64;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let d = (x - y).abs() / x.max(y).max(1.0);
        if d > md {
            md = d;
        }
    }
    md
}

/// Stage 2 routed `(node_counts, edge_counts)` maps returned by [`route_multimappers_to_stage2`].
pub type Stage2RoutedCounts = (
    FxHashMap<(u32, u16), f64>,
    FxHashMap<(u32, u16, u16), (f64, u32, u32)>,
);

/// Deterministic Stage 1 -> Stage 2 routing:
/// Returns finalized per-gene `(node_counts, edge_counts)` maps containing BOTH uniquely-mapped
/// counts (bias-adjusted if enabled) AND VBEM-resolved multi-mapper assignments!
pub fn route_multimappers_to_stage2(
    idx: &MmapIndex,
    counts: &SampleCounts,
    vbem: &VbemResult,
    bias_model: &BiasModel,
    bias_correct: bool,
) -> Stage2RoutedCounts {
    // 1. Copy bias-adjusted unique node and edge counts
    let mut final_nodes: FxHashMap<(u32, u16), f64> =
        FxHashMap::with_capacity_and_hasher(counts.node_counts.len(), Default::default());
    for (&k, fc) in &counts.node_counts {
        let v = fc.joint_adjusted(bias_model, bias_correct);
        if v > 0.0 {
            *final_nodes.entry(k).or_insert(0.0) += v;
        }
    }

    let mut final_edges: FxHashMap<(u32, u16, u16), (f64, u32, u32)> =
        FxHashMap::with_capacity_and_hasher(counts.edge_counts.len(), Default::default());
    for (&k, jst) in &counts.edge_counts {
        let v = jst.count.joint_adjusted(bias_model, bias_correct);
        if v > 0.0 {
            final_edges.insert(k, (v, jst.distinct_positions(), jst.max_stack()));
        }
    }

    // 2. Route every multi-mapping group deterministically using converged transcript posterior
    //    weights (or converged gene TPM weights for unannotated novel-junction multi-mappers).
    let mut sorted_groups: Vec<(&MultiGroupKey, f64)> = counts
        .multi_groups
        .iter()
        .map(|(k, fc)| (k, fc.joint_adjusted(bias_model, bias_correct)))
        .filter(|&(_, v)| v > 0.0)
        .collect();
    sorted_groups.sort_by(|a, b| a.0.candidates.cmp(&b.0.candidates));

    let mut cand_weights: SmallVec<[f64; 8]> = SmallVec::new();
    for (group, total_w) in sorted_groups {
        let n_cands = group.candidates.len();
        if n_cands == 0 {
            continue;
        }
        cand_weights.clear();

        // Try transcript-level converged weights first (sum of tx_tpm / eff_len_rate for compatible txs)
        let mut all_tx_matched = true;
        let mut tx_w_sum = 0.0f64;
        for cand in &group.candidates {
            let txs = compatible_txs_for_candidate(idx, cand);
            if txs.is_empty() {
                all_tx_matched = false;
                break;
            }
            let w: f64 = txs.iter().map(|&t| vbem.tx_tpm[t as usize]).sum();
            cand_weights.push(w);
            tx_w_sum += w;
        }

        if !all_tx_matched || tx_w_sum <= 0.0 {
            // Fallback for novel-junction multi-mappers: use converged gene-level TPM weights
            cand_weights.clear();
            let mut gene_w_sum = 0.0f64;
            for cand in &group.candidates {
                let gw = vbem.gene_tpm[cand.gene_id as usize];
                cand_weights.push(gw);
                gene_w_sum += gw;
            }
            if gene_w_sum > 0.0 {
                for w in &mut cand_weights {
                    *w /= gene_w_sum;
                }
            } else {
                let unif = 1.0 / (n_cands as f64);
                cand_weights.fill(unif);
            }
        } else {
            for w in &mut cand_weights {
                *w /= tx_w_sum;
            }
        }

        // Distribute weighted counts into Stage 2 node/edge maps
        for (cand, &prop) in group.candidates.iter().zip(cand_weights.iter()) {
            let share = total_w * prop;
            if share <= 1e-9 {
                continue;
            }
            let mut seen_edges: SmallVec<[(u16, u16); 4]> = SmallVec::new();
            if cand.path_5p.len() == 1 {
                *final_nodes.entry((cand.gene_id, cand.path_5p[0])).or_insert(0.0) += share;
            } else {
                for win in cand.path_5p.windows(2) {
                    let ek = (win[0], win[1]);
                    seen_edges.push(ek);
                    let entry = final_edges
                        .entry((cand.gene_id, ek.0, ek.1))
                        .or_insert((0.0, 0, 0));
                    entry.0 += share;
                }
            }
            if !cand.path_3p.is_empty() {
                if cand.path_3p.len() == 1 {
                    let n = cand.path_3p[0];
                    if !cand.path_5p.contains(&n) {
                        *final_nodes.entry((cand.gene_id, n)).or_insert(0.0) += share;
                    }
                } else {
                    for win in cand.path_3p.windows(2) {
                        let ek = (win[0], win[1]);
                        if !seen_edges.contains(&ek) {
                            let entry = final_edges
                                .entry((cand.gene_id, ek.0, ek.1))
                                .or_insert((0.0, 0, 0));
                            entry.0 += share;
                        }
                    }
                }
            }
        }
    }

    (final_nodes, final_edges)
}
