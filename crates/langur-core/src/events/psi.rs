//! Stage 2: Local splice-graph event PSI quantification across all Whippet motifs
//! (`CE`, `AA`, `AD`, `RI`, `TS`, `TE`, `AF`, `AL`, and `BS` with `--circ`),
//! true junction-overhang effective lengths, MAJIQ-style PCR-stack soft-capping,
//! neutral `K_1` EM / optional complex VBEM, and Beta/Dirichlet `Inf_Var` + 90% CIs.

use crate::events::gibbs::{
    sample_complex_event_gibbs, DEFAULT_GIBBS_BURN_IN, DEFAULT_GIBBS_THIN,
};
use crate::events::reduce::{complexity, reduce_graph, shannon_entropy, PathEdge, ReducedPath};
use crate::graph::EdgeMotif;
use crate::index::{MmapIndex, NodeRecord};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use statrs::distribution::{Beta, ContinuousCDF};
use statrs::function::gamma::digamma;

#[derive(Clone, Debug)]
pub struct EventConfig {
    /// Emit circular back-splicing (`BS`) events.
    pub allow_circ: bool,
    /// Minimum relative edge weight vs max edge in gene (matches Whippet2 `0.0001`).
    pub min_edge_weight: f64,
    /// Enable MAJIQ-style PCR-stack soft-capping on novel/`RI` junctions (default `true`).
    pub stack_filter: bool,
    /// Minimum count threshold to trigger single-position PCR-stack soft-capping (default `10.0`).
    pub stack_min_count: f64,
    /// Effective read cap applied when single-position PCR-stack filter triggers (default `2.0`).
    pub stack_cap_reads: f64,
    /// Minimum overhang for annotated junctions (default `4`).
    pub min_anno_overhang: u32,
    /// Minimum overhang for novel junctions & stack geometry checks (default `8`).
    pub min_overhang: u32,
    /// Mean read length observed in the sample.
    pub mean_read_len: f64,
    /// Optional sparse VBEM prior for complex (`K_n >= 2`) multi-path events (default `None` = neutral EM).
    pub complex_vb_prior: Option<f64>,
    /// Number of targeted collapsed Gibbs posterior draws for complex (`K_n >= 2`) events with
    /// ambiguous shared edges (default `100`; `0` disables Gibbs and uses closed-form variance).
    pub gibbs_samples: usize,
    /// Deterministic RNG seed for per-event complex Gibbs sampling (default `42`).
    pub gibbs_seed: u64,
}

impl Default for EventConfig {
    fn default() -> Self {
        Self {
            allow_circ: false,
            min_edge_weight: 0.0001,
            stack_filter: true,
            stack_min_count: 10.0,
            stack_cap_reads: 2.0,
            min_anno_overhang: 4,
            min_overhang: 8,
            mean_read_len: 100.0,
            complex_vb_prior: None,
            gibbs_samples: 100,
            gibbs_seed: 42,
        }
    }
}

/// One quantified alternative-splicing event row corresponding to `.psi.gz`
/// (exact 14 Whippet2 columns + column 15 `Inf_Var`).
#[derive(Clone, Debug)]
pub struct EventRecord {
    pub gene_id: u32,
    /// 1-based CSG node index (matching Whippet2 `.psi.gz`).
    pub node_1based: u32,
    pub chrom_id: u32,
    pub coord_start: u32,
    pub coord_end: u32,
    pub is_pos_strand: bool,
    pub event_type: &'static str, // "CE", "AA", "AD", "RI", "TS", "TE", "AF", "AL", "BS", "NA"
    pub psi: Option<f64>,
    pub ci_lo: Option<f64>,
    pub ci_hi: Option<f64>,
    pub total_reads: f64,
    pub complexity: usize,
    pub entropy: f64,
    pub inc_paths: Vec<(Vec<u16>, f64)>, // 1-based node lists + path psi
    pub exc_paths: Vec<(Vec<u16>, f64)>,
    pub edges: Vec<(u16, u16, f64)>,     // 1-based (u, v, count)
    /// Posterior inferential variance on `Psi` (column 15 `Inf_Var`).
    pub inf_var: Option<f64>,
}

/// Quantify all splice-graph events across every gene in deterministic gene order (parallelized with `rayon`).
pub fn quantify_all_events(
    idx: &MmapIndex,
    final_nodes: &FxHashMap<(u32, u16), f64>,
    final_edges: &FxHashMap<(u32, u16, u16), (f64, u32, u32)>,
    cfg: &EventConfig,
) -> Vec<EventRecord> {
    let n_genes = idx.genes().len();

    // Group sparse node & edge counts per gene once in O(nnz)
    let mut gene_nodes_map: FxHashMap<u32, Vec<(u16, f64)>> = FxHashMap::default();
    for (&(gid, nid), &cnt) in final_nodes {
        if cnt > 0.0 {
            gene_nodes_map.entry(gid).or_default().push((nid, cnt));
        }
    }

    type GeneEdgeMap = FxHashMap<u32, Vec<(u16, u16, f64, u32, u32)>>;
    let mut gene_edges_map: GeneEdgeMap = FxHashMap::default();
    for (&(gid, u, v), &(cnt, dpos, mstack)) in final_edges {
        if cnt > 0.0 {
            gene_edges_map
                .entry(gid)
                .or_default()
                .push((u, v, cnt, dpos, mstack));
        }
    }

    let per_gene: Vec<Vec<EventRecord>> = (0..n_genes)
        .into_par_iter()
        .map(|gidx| {
            let gid = gidx as u32;
            let n_cnts = gene_nodes_map.get(&gid).map(Vec::as_slice).unwrap_or(&[]);
            let e_cnts = gene_edges_map.get(&gid).map(Vec::as_slice).unwrap_or(&[]);
            quantify_gene_events(idx, gid, n_cnts, e_cnts, cfg)
        })
        .collect();

    let total: usize = per_gene.iter().map(Vec::len).sum();
    let mut out = Vec::with_capacity(total);
    for mut v in per_gene {
        out.append(&mut v);
    }
    out
}

/// Apply MAJIQ-style single-position PCR-stack soft-capping ONLY to novel or RI-touching
/// junctions when read geometry allows multiple overhang positions (`mean_read_len >= 2*min_ov + 4`).
/// Annotated non-RI junctions are NEVER capped.
#[inline]
pub fn apply_pcr_stack_cap(
    count: f64,
    distinct_pos: u32,
    is_annotated: bool,
    touches_ri: bool,
    cfg: &EventConfig,
) -> f64 {
    if !cfg.stack_filter {
        return count;
    }
    if is_annotated && !touches_ri {
        return count;
    }
    let geometry_allows_multi = cfg.mean_read_len >= (2 * cfg.min_overhang + 4) as f64;
    if geometry_allows_multi && distinct_pos == 1 && count >= cfg.stack_min_count {
        count.min(cfg.stack_cap_reads)
    } else {
        count
    }
}

#[inline]
fn nodes_contiguous_psi(prev: &NodeRecord, next: &NodeRecord, is_pos_strand: bool) -> bool {
    if is_pos_strand {
        prev.genomic_coord + prev.len == next.genomic_coord
    } else {
        next.genomic_coord + next.len == prev.genomic_coord
    }
}

/// Compute true junction-overhang effective length for an edge `(u, v)`.
///
/// Walks contiguous upstream exonic nodes from `u` (stopping at `max_side` or an upstream splice
/// acceptor) and contiguous downstream exonic nodes from `v` (stopping at `max_side` or a
/// downstream splice donor). This ensures:
/// - Internal microexons (e.g. 15 bp `CE`) retain full mappable overhang space.
/// - Tandem `TS`/`TE` nodes inside multi-node exons (e.g. a 2 bp `TS` sub-node contiguous with
///   adjacent exonic nodes) never collapse `eff_len` to `~0.01`, which would otherwise amplify
///   1-read skipping paths by ~90x!
/// - True short terminal exons are floored at `0.15` relative effective length to guard against
///   extreme variance inflation from noisy terminal reads.
#[inline]
pub fn edge_effective_length(
    nodes: &[NodeRecord],
    is_pos_strand: bool,
    u: u16,
    v: u16,
    is_annotated: bool,
    cfg: &EventConfig,
) -> f64 {
    let rlen = cfg.mean_read_len;
    let min_ov = if is_annotated {
        cfg.min_anno_overhang as f64
    } else {
        cfg.min_overhang as f64
    };
    let ref_span = (rlen - 2.0 * (cfg.min_anno_overhang as f64) + 1.0).max(1.0);
    let max_span = (rlen - 2.0 * min_ov + 1.0).max(1.0);
    let max_side = (rlen - min_ov).max(1.0);

    // Walk contiguous upstream exonic extent ending at u's 3' boundary
    let mut lu = 0.0f64;
    let mut cur = u as usize;
    loop {
        lu += nodes[cur].len as f64;
        if lu >= max_side || nodes[cur].left_edge().is_acceptor() {
            lu = max_side;
            break;
        }
        if cur > 0 && nodes_contiguous_psi(&nodes[cur - 1], &nodes[cur], is_pos_strand) {
            cur -= 1;
        } else {
            break;
        }
    }

    // Walk contiguous downstream exonic extent starting at v's 5' boundary
    let mut lv = 0.0f64;
    let mut cur = v as usize;
    loop {
        lv += nodes[cur].len as f64;
        if lv >= max_side || nodes[cur].right_edge().is_donor() {
            lv = max_side;
            break;
        }
        if cur + 1 < nodes.len() && nodes_contiguous_psi(&nodes[cur], &nodes[cur + 1], is_pos_strand)
        {
            cur += 1;
        } else {
            break;
        }
    }

    let avail = (lu + lv - rlen + 1.0).clamp(1.0, max_span);
    (avail / ref_span).max(0.15)
}

fn quantify_gene_events(
    idx: &MmapIndex,
    gene_id: u32,
    node_counts_raw: &[(u16, f64)],
    edge_counts_raw: &[(u16, u16, f64, u32, u32)],
    cfg: &EventConfig,
) -> Vec<EventRecord> {
    let gene = &idx.genes()[gene_id as usize];
    let nodes = idx.gene_nodes(gene_id);
    let n_nodes = nodes.len();
    if n_nodes <= 1 {
        return Vec::new();
    }
    let anno_jncs = idx.gene_anno_junctions(gene_id);

    // Build local node count slice for O(1) lookup
    let mut local_node_cnt = vec![0.0f64; n_nodes];
    for &(nid, c) in node_counts_raw {
        if (nid as usize) < n_nodes {
            local_node_cnt[nid as usize] += c;
        }
    }

    // Separate linear forward edges (u < v) from circular back-splice edges (u >= v)
    let mut linear_edges: Vec<PathEdge> = Vec::new();
    let mut circ_edges: Vec<(u16, u16, f64)> = Vec::new();

    for &(u, v, raw_cnt, dpos, _mstack) in edge_counts_raw {
        if (u as usize) >= n_nodes || (v as usize) >= n_nodes {
            continue;
        }
        if u < v {
            let is_anno = anno_jncs
                .binary_search_by_key(&(u, v), |j| (j.src_node, j.dst_node))
                .is_ok();
            let touches_ri =
                nodes[u as usize].motif() == EdgeMotif::RI || nodes[v as usize].motif() == EdgeMotif::RI;
            let adj_cnt = apply_pcr_stack_cap(raw_cnt, dpos, is_anno, touches_ri, cfg);
            let eff_len = edge_effective_length(nodes, gene.is_pos_strand(), u, v, is_anno, cfg);
            linear_edges.push(PathEdge {
                src: u,
                dst: v,
                count: adj_cnt,
                eff_len,
            });
        } else if cfg.allow_circ {
            let adj_cnt = apply_pcr_stack_cap(raw_cnt, dpos, false, false, cfg);
            circ_edges.push((u, v, adj_cnt));
        }
    }
    linear_edges.sort_unstable_by_key(|e| (e.src, e.dst));
    circ_edges.sort_unstable_by_key(|e| (e.0, e.1));

    let max_edge_val = linear_edges
        .iter()
        .map(|e| e.count)
        .fold(0.0f64, f64::max);

    let mut out = Vec::new();
    let mut i = 0usize;

    while i < n_nodes {
        let motif = nodes[i].motif();
        let next_is_te = i + 1 < n_nodes && nodes[i + 1].motif() == EdgeMotif::TE;

        match motif {
            EdgeMotif::NA => {
                if !next_is_te {
                    out.push(empty_event(gene_id, gene, nodes, i, motif.as_str()));
                }
                i += 1;
            }
            EdgeMotif::TS | EdgeMotif::TE => {
                let (evts, next_i) = quantify_tandem_utr(
                    gene_id,
                    gene,
                    nodes,
                    &local_node_cnt,
                    &linear_edges,
                    i,
                    motif,
                    cfg,
                );
                out.extend(evts);
                i = next_i.max(i + 1);
            }
            _ => {
                // Spliced event: CE, AA, AD, RI, AF, AL
                let evt = quantify_spliced_node(
                    gene_id,
                    gene,
                    nodes,
                    &linear_edges,
                    max_edge_val,
                    i as u16,
                    motif,
                    cfg,
                );
                out.push(evt);
                i += 1;
            }
        }
    }

    // Emit back-splicing (BS) events when `--circ` is enabled
    if cfg.allow_circ && !circ_edges.is_empty() {
        for (u, v, bs_cnt) in circ_edges {
            // Linear reads flanking the circularized exon boundaries (v -> downstream or upstream -> u)
            let lin_cnt: f64 = linear_edges
                .iter()
                .filter(|e| (e.src == u && e.dst > u) || (e.dst == v && e.src < v))
                .map(|e| e.count)
                .sum();
            let total = bs_cnt + lin_cnt;
            let (coord_s, coord_e) = if gene.is_pos_strand() {
                (
                    nodes[v as usize].genomic_coord,
                    nodes[u as usize].genomic_coord + nodes[u as usize].len.saturating_sub(1),
                )
            } else {
                (
                    nodes[u as usize].genomic_coord,
                    nodes[v as usize].genomic_coord + nodes[v as usize].len.saturating_sub(1),
                )
            };
            let psi = if total > 0.0 { Some(bs_cnt / total) } else { None };
            let (ci_lo, ci_hi, inf_var) = match psi {
                Some(p) => {
                    let (lo, hi) = beta_posterior_ci_90(p, total);
                    (Some(lo), Some(hi), Some(beta_posterior_var(p, total)))
                }
                None => (None, None, None),
            };
            out.push(EventRecord {
                gene_id,
                node_1based: (v as u32) + 1,
                chrom_id: gene.chrom_id,
                coord_start: coord_s,
                coord_end: coord_e,
                is_pos_strand: gene.is_pos_strand(),
                event_type: "BS",
                psi,
                ci_lo,
                ci_hi,
                total_reads: total,
                complexity: if lin_cnt > 0.0 { 1 } else { 0 },
                entropy: psi.map_or(0.0, |p| shannon_entropy(&[p, 1.0 - p])),
                inc_paths: vec![(vec![u + 1, v + 1], psi.unwrap_or(0.0))],
                exc_paths: Vec::new(),
                edges: vec![(u + 1, v + 1, bs_cnt)],
                inf_var,
            });
        }
    }

    out
}

#[allow(clippy::too_many_arguments)]
fn quantify_spliced_node(
    gene_id: u32,
    gene: &crate::index::GeneRecord,
    nodes: &[NodeRecord],
    all_edges: &[PathEdge],
    max_edge_val: f64,
    node: u16,
    motif: EdgeMotif,
    cfg: &EventConfig,
) -> EventRecord {
    if all_edges.is_empty() || max_edge_val <= 0.0 {
        return empty_event(gene_id, gene, nodes, node as usize, motif.as_str());
    }

    let mut inc_edges: Vec<PathEdge> = Vec::new();
    let mut exc_edges: Vec<PathEdge> = Vec::new();
    let mut local_edges: Vec<PathEdge> = Vec::new();
    let mut connecting_val = 0.0f64;
    let mut spanning_val = 0.0f64;

    for &e in all_edges {
        if e.count / max_edge_val < cfg.min_edge_weight {
            continue;
        }
        if e.src == node || e.dst == node {
            inc_edges.push(e);
            local_edges.push(e);
            connecting_val += e.count;
        } else if e.src < node && node < e.dst {
            exc_edges.push(e);
            local_edges.push(e);
            spanning_val += e.count;
        }
    }

    let total_reads = connecting_val + spanning_val;
    let is_alt_ss = matches!(motif, EdgeMotif::AA | EdgeMotif::AD);

    // Pure inclusion (no skipping edges)
    if exc_edges.is_empty() {
        let valid = !inc_edges.is_empty()
            && (connecting_val >= 2.0 || (connecting_val >= 1.0 && is_alt_ss));
        if !valid {
            return empty_event(gene_id, gene, nodes, node as usize, motif.as_str());
        }
        let mut inc_paths = reduce_graph(&inc_edges, 64, 1024);
        normalize_single_side_path_psi(&mut inc_paths);
        let (ci_lo, ci_hi) = beta_posterior_ci_90(1.0, total_reads);
        let inf_var = beta_posterior_var(1.0, total_reads);
        return build_event_record(
            gene_id,
            gene,
            nodes,
            node as usize,
            motif.as_str(),
            1.0,
            ci_lo,
            ci_hi,
            inf_var,
            total_reads,
            &inc_paths,
            &[],
            &local_edges,
        );
    }

    // Pure exclusion (no inclusion edges)
    if inc_edges.is_empty() {
        if spanning_val < 1.0 {
            return empty_event(gene_id, gene, nodes, node as usize, motif.as_str());
        }
        let mut exc_paths = reduce_graph(&exc_edges, 64, 1024);
        normalize_single_side_path_psi(&mut exc_paths);
        let (ci_lo, ci_hi) = beta_posterior_ci_90(0.0, total_reads);
        let inf_var = beta_posterior_var(0.0, total_reads);
        return build_event_record(
            gene_id,
            gene,
            nodes,
            node as usize,
            motif.as_str(),
            0.0,
            ci_lo,
            ci_hi,
            inf_var,
            total_reads,
            &[],
            &exc_paths,
            &local_edges,
        );
    }

    // Both inclusion and exclusion edges present: perform two outward sweeps strictly
    // bounded by the initial [min_v, max_v] extent (matching Whippet2 `extend_edges!`).
    let min_v = inc_edges
        .iter()
        .chain(exc_edges.iter())
        .map(|e| e.src)
        .min()
        .unwrap_or(node);
    let max_v = inc_edges
        .iter()
        .chain(exc_edges.iter())
        .map(|e| e.dst)
        .max()
        .unwrap_or(node);

    let span = (max_v - min_v) as usize + 1;
    let mut inc_nodes = vec![false; span];
    let mut exc_nodes = vec![false; span];
    for e in &inc_edges {
        inc_nodes[(e.src - min_v) as usize] = true;
        inc_nodes[(e.dst - min_v) as usize] = true;
    }
    for e in &exc_edges {
        exc_nodes[(e.src - min_v) as usize] = true;
        exc_nodes[(e.dst - min_v) as usize] = true;
    }

    let mut left_cands: Vec<PathEdge> = Vec::new();
    let mut right_cands: Vec<PathEdge> = Vec::new();
    for &e in all_edges {
        if e.count / max_edge_val < cfg.min_edge_weight {
            continue;
        }
        if e.src >= min_v && e.dst > min_v && e.dst < node {
            left_cands.push(e);
        } else if e.src > node && e.src < max_v && e.dst <= max_v {
            right_cands.push(e);
        }
    }
    left_cands.sort_unstable_by(|a, b| b.dst.cmp(&a.dst).then_with(|| a.src.cmp(&b.src)));
    right_cands.sort_unstable_by(|a, b| a.src.cmp(&b.src).then_with(|| a.dst.cmp(&b.dst)));

    let mut bridge_edges: Vec<PathEdge> = Vec::new();
    // Left outward sweep: propagate from `node - 1` down to `min_v + 1`
    for e in left_cands {
        let dst_off = (e.dst - min_v) as usize;
        let touches_inc = inc_nodes[dst_off];
        let touches_exc = exc_nodes[dst_off];
        if touches_inc {
            inc_edges.push(e);
            inc_nodes[(e.src - min_v) as usize] = true;
        }
        if touches_exc {
            exc_edges.push(e);
            exc_nodes[(e.src - min_v) as usize] = true;
        }
        if touches_inc || touches_exc {
            bridge_edges.push(e);
        }
    }
    // Right outward sweep: propagate from `node + 1` up to `max_v - 1`
    for e in right_cands {
        let src_off = (e.src - min_v) as usize;
        let touches_inc = inc_nodes[src_off];
        let touches_exc = exc_nodes[src_off];
        if touches_inc {
            inc_edges.push(e);
            inc_nodes[(e.dst - min_v) as usize] = true;
        }
        if touches_exc {
            exc_edges.push(e);
            exc_nodes[(e.dst - min_v) as usize] = true;
        }
        if touches_inc || touches_exc {
            bridge_edges.push(e);
        }
    }

    local_edges.extend_from_slice(&bridge_edges);
    local_edges.sort_unstable_by_key(|e| (e.src, e.dst));
    local_edges.dedup_by_key(|e| (e.src, e.dst));

    let mut inc_paths = reduce_graph(&inc_edges, 64, 1024);
    let mut exc_paths = reduce_graph(&exc_edges, 64, 1024);
    // Keep only valid inclusion paths containing `node` and valid exclusion paths that
    // actually skip `node` (preventing dead-end upstream bridge fragments from polluting exc_paths).
    inc_paths.retain(|p| p.nodes.contains(&node));
    exc_paths.retain(|p| p.nodes.windows(2).any(|w| w[0] < node && node < w[1]));
    if inc_paths.is_empty() || exc_paths.is_empty() {
        return empty_event(gene_id, gene, nodes, node as usize, motif.as_str());
    }
    // Reset counts/lengths to accumulate exact per-edge assignments
    for p in inc_paths.iter_mut().chain(exc_paths.iter_mut()) {
        p.count = 0.0;
        p.eff_len = 0.0;
    }

    let n_inc = inc_paths.len();
    let n_exc = exc_paths.len();
    let n_tot = n_inc + n_exc;

    let mut ambig_groups: Vec<(Vec<usize>, f64)> = Vec::new();
    let mut ambig_read_sum = 0.0f64;

    for e in &local_edges {
        let mut compat: Vec<usize> = Vec::new();
        for (pi, p) in inc_paths.iter_mut().enumerate() {
            if p.contains_edge(e.src, e.dst) {
                compat.push(pi);
                p.eff_len += e.eff_len;
            }
        }
        for (pi, p) in exc_paths.iter_mut().enumerate() {
            if p.contains_edge(e.src, e.dst) {
                compat.push(n_inc + pi);
                p.eff_len += e.eff_len;
            }
        }
        if compat.len() == 1 {
            let pid = compat[0];
            if pid < n_inc {
                inc_paths[pid].count += e.count;
            } else {
                exc_paths[pid - n_inc].count += e.count;
            }
        } else if compat.len() > 1 && e.count > 0.0 {
            ambig_read_sum += e.count;
            if let Some(ag) = ambig_groups.iter_mut().find(|g| g.0 == compat) {
                ag.1 += e.count;
            } else {
                ambig_groups.push((compat, e.count));
            }
        }
    }

    // Ensure every path has positive effective length
    for p in inc_paths.iter_mut().chain(exc_paths.iter_mut()) {
        if p.eff_len <= 0.0 {
            p.eff_len = 1.0;
        }
    }

    // Run local EM (neutral EM for binary K_1 events; optional VBEM for complex K_n >= 2 events)
    let k_n = complexity(n_tot);
    let vb_prior = if k_n >= 2 { cfg.complex_vb_prior } else { None };
    let (all_psi, converged_counts) =
        solve_local_path_em(&inc_paths, &exc_paths, &ambig_groups, vb_prior, 750);

    for i in 0..n_inc {
        inc_paths[i].psi = all_psi[i];
    }
    for i in 0..n_exc {
        exc_paths[i].psi = all_psi[n_inc + i];
    }

    let psi: f64 = inc_paths.iter().map(|p| p.psi).sum::<f64>().clamp(0.0, 1.0);
    if !psi.is_finite() || total_reads <= 0.0 {
        return empty_event(gene_id, gene, nodes, node as usize, motif.as_str());
    }

    let (ci_lo, ci_hi, inf_var) = if k_n >= 2
        && !ambig_groups.is_empty()
        && cfg.gibbs_samples >= 2
    {
        if let Some(gs) = sample_complex_event_gibbs(
            &inc_paths,
            &exc_paths,
            &ambig_groups,
            &all_psi,
            vb_prior,
            cfg.gibbs_samples,
            DEFAULT_GIBBS_BURN_IN,
            DEFAULT_GIBBS_THIN,
            cfg.gibbs_seed,
            gene_id,
            (node as u32) + 1,
        ) {
            (gs.ci_lo, gs.ci_hi, gs.inf_var)
        } else {
            let (lo, hi) = beta_posterior_ci_90(psi, total_reads);
            let iv = compute_inferential_var(
                psi,
                total_reads,
                &converged_counts[..n_inc],
                &converged_counts[n_inc..],
                ambig_read_sum,
                k_n,
            );
            (lo, hi, iv)
        }
    } else {
        // Binary K_1 events (or unambiguous events, or --gibbs-samples 0): exact closed-form O(1)
        let (lo, hi) = beta_posterior_ci_90(psi, total_reads);
        let iv = compute_inferential_var(
            psi,
            total_reads,
            &converged_counts[..n_inc],
            &converged_counts[n_inc..],
            ambig_read_sum,
            k_n,
        );
        (lo, hi, iv)
    };

    build_event_record(
        gene_id,
        gene,
        nodes,
        node as usize,
        motif.as_str(),
        psi,
        ci_lo,
        ci_hi,
        inf_var,
        total_reads,
        &inc_paths,
        &exc_paths,
        &local_edges,
    )
}

/// Populate normalized path `psi` proportions (summing to 1.0) on pure-inclusion or
/// pure-exclusion reduced paths so `.psi.gz` never prints `:0.0000` for a 100% path.
fn normalize_single_side_path_psi(paths: &mut [ReducedPath]) {
    if paths.is_empty() {
        return;
    }
    let mut sum = 0.0f64;
    for p in paths.iter_mut() {
        let rate = p.count.max(0.0) / p.eff_len.max(1e-6);
        p.psi = rate;
        sum += rate;
    }
    if sum > 0.0 {
        for p in paths.iter_mut() {
            p.psi /= sum;
        }
    } else {
        let unif = 1.0 / (paths.len() as f64);
        for p in paths.iter_mut() {
            p.psi = unif;
        }
    }
}

/// Solve local splice-graph path proportions via neutral EM (or optional VBEM for `K_n >= 2`).
fn solve_local_path_em(
    inc: &[ReducedPath],
    exc: &[ReducedPath],
    ambig: &[(Vec<usize>, f64)],
    vb_prior: Option<f64>,
    max_iters: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n_inc = inc.len();
    let n_tot = n_inc + exc.len();
    let mut uniq = vec![0.0f64; n_tot];
    let mut eff_len = vec![1.0f64; n_tot];
    for (i, p) in inc.iter().enumerate() {
        uniq[i] = p.count;
        eff_len[i] = p.eff_len.max(1e-6);
    }
    for (i, p) in exc.iter().enumerate() {
        uniq[n_inc + i] = p.count;
        eff_len[n_inc + i] = p.eff_len.max(1e-6);
    }

    // Initialize counts uniformly across ambiguous groups
    let mut counts = uniq.clone();
    for (paths, w) in ambig {
        let share = *w / (paths.len() as f64);
        for &pid in paths {
            counts[pid] += share;
        }
    }

    let mut psi = vec![0.0f64; n_tot];
    let mut prev_psi = vec![0.0f64; n_tot];

    let update_psi = |c: &[f64], out: &mut [f64]| {
        let mut s = 0.0;
        for i in 0..n_tot {
            let v = if let Some(a0) = vb_prior {
                let alpha = a0 + c[i].max(0.0);
                if alpha > 1e-12 {
                    digamma(alpha).exp() / eff_len[i]
                } else {
                    0.0
                }
            } else {
                c[i].max(0.0) / eff_len[i]
            };
            out[i] = v;
            s += v;
        }
        if s > 0.0 {
            for v in out.iter_mut() {
                *v /= s;
            }
        }
    };

    update_psi(&counts, &mut psi);

    if ambig.is_empty() {
        return (psi, counts);
    }

    for _ in 0..max_iters {
        prev_psi.copy_from_slice(&psi);
        counts.copy_from_slice(&uniq);

        for (paths, w) in ambig {
            let denom: f64 = paths.iter().map(|&pid| prev_psi[pid]).sum();
            if denom > 0.0 {
                let scale = *w / denom;
                for &pid in paths {
                    counts[pid] += prev_psi[pid] * scale;
                }
            } else {
                let share = *w / (paths.len() as f64);
                for &pid in paths {
                    counts[pid] += share;
                }
            }
        }

        update_psi(&counts, &mut psi);
        let max_d = psi
            .iter()
            .zip(prev_psi.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        if max_d < 1e-6 {
            break;
        }
    }

    (psi, counts)
}

#[allow(clippy::too_many_arguments)]
fn quantify_tandem_utr(
    gene_id: u32,
    gene: &crate::index::GeneRecord,
    nodes: &[NodeRecord],
    local_node_cnt: &[f64],
    linear_edges: &[PathEdge],
    start_idx: usize,
    motif: EdgeMotif,
    cfg: &EventConfig,
) -> (Vec<EventRecord>, usize) {
    let n_nodes = nodes.len();
    let start_n = if motif == EdgeMotif::TE && start_idx > 0 {
        start_idx - 1
    } else {
        start_idx
    };

    let mut cur = start_idx;
    while cur < n_nodes && nodes[cur].motif() == motif {
        cur += 1;
    }
    let end_n = if motif == EdgeMotif::TS && cur < n_nodes {
        cur
    } else {
        cur.saturating_sub(1)
    };
    // Prevent duplicate node rows when `TS` incorporates the downstream anchor `nodes[cur]`
    // (or `TE` incorporates the upstream anchor `nodes[start_n]`):
    // - For `TS`: if `nodes[cur]` is `NA`, emit it inside the `TS` batch and advance `next_i = cur + 1`.
    //   If `nodes[cur]` has its own active splicing motif (e.g. `CE`), omit emitting it here and
    //   return `next_i = cur` so the outer loop quantifies it under its own motif.
    // - For `TE`: `nodes[start_n]` was skipped in the outer loop iff its motif was `NA` (`next_is_te`);
    //   if it had an active motif, it was already emitted at `i = start_idx - 1`, so omit duplicating it!
    let emit_anchor_ts = motif == EdgeMotif::TS && cur < n_nodes && nodes[cur].motif() == EdgeMotif::NA;
    let emit_anchor_te = motif == EdgeMotif::TE && start_idx > 0 && nodes[start_n].motif() == EdgeMotif::NA;
    let next_i = if emit_anchor_ts { cur + 1 } else { cur };

    let should_emit_node = |n_idx: usize| -> bool {
        if motif == EdgeMotif::TS && cur < n_nodes && n_idx == cur {
            emit_anchor_ts
        } else if motif == EdgeMotif::TE && start_idx > 0 && n_idx == start_n {
            emit_anchor_te
        } else {
            true
        }
    };

    let seg_len = end_n - start_n + 1;
    if seg_len < 2 {
        return (
            vec![empty_event(gene_id, gene, nodes, start_idx, motif.as_str())],
            next_i,
        );
    }

    // Build nested UTR paths (TS grows backward from end_n..start_n; TE grows forward from start_n..end_n)
    let mut paths: Vec<Vec<u16>> = Vec::with_capacity(seg_len);
    let mut cur_set: Vec<u16> = Vec::new();
    for step in 0..seg_len {
        let nid = if motif == EdgeMotif::TS {
            (end_n - step) as u16
        } else {
            (start_n + step) as u16
        };
        cur_set.push(nid);
        cur_set.sort_unstable();
        paths.push(cur_set.clone());
    }

    let mut eff_lens = vec![0.0f64; seg_len];
    let mut uniq_cnt = vec![0.0f64; seg_len];
    let mut ambig: Vec<(Vec<usize>, f64)> = Vec::new();
    let mut total_cnt = 0.0f64;

    for n in start_n..=end_n {
        let nid = n as u16;
        let c = local_node_cnt[n];
        total_cnt += c;
        let mut compat = Vec::new();
        for (pi, p) in paths.iter().enumerate() {
            if p.binary_search(&nid).is_ok() {
                compat.push(pi);
                eff_lens[pi] += nodes[n].len as f64;
            }
        }
        if compat.len() == 1 {
            uniq_cnt[compat[0]] += c;
        } else if compat.len() > 1 && c > 0.0 {
            ambig.push((compat, c));
        }
    }

    for e in linear_edges {
        if (e.src as usize) >= start_n && (e.dst as usize) <= end_n {
            total_cnt += e.count;
            let mut compat = Vec::new();
            for (pi, p) in paths.iter().enumerate() {
                if let Ok(pos) = p.binary_search(&e.src) {
                    if pos + 1 < p.len() && p[pos + 1] == e.dst {
                        compat.push(pi);
                    }
                }
            }
            if compat.len() == 1 {
                uniq_cnt[compat[0]] += e.count;
            } else if compat.len() > 1 && e.count > 0.0 {
                ambig.push((compat, e.count));
            }
        }
    }

    if total_cnt <= 2.0 {
        let evts = (start_n..=end_n)
            .filter(|&n| should_emit_node(n))
            .map(|n| empty_event(gene_id, gene, nodes, n, motif.as_str()))
            .collect();
        return (evts, next_i);
    }

    let reduced: Vec<ReducedPath> = paths
        .iter()
        .enumerate()
        .map(|(i, p)| ReducedPath {
            nodes: p.clone(),
            count: uniq_cnt[i],
            eff_len: (eff_lens[i] - cfg.mean_read_len + 1.0).max(1.0),
            psi: 0.0,
        })
        .collect();

    let (psi_vec, _counts) = solve_local_path_em(&reduced, &[], &ambig, None, 1000);
    let comp_k = complexity(seg_len);
    let ent = shannon_entropy(&psi_vec);

    let mut evts = Vec::with_capacity(seg_len);
    for step in 0..seg_len {
        let n_idx = if motif == EdgeMotif::TS {
            end_n - step
        } else {
            start_n + step
        };
        if !should_emit_node(n_idx) {
            continue;
        }
        let p = psi_vec[step].clamp(0.0, 1.0);
        let (ci_lo, ci_hi) = beta_posterior_ci_90(p, total_cnt);
        let inf_var = beta_posterior_var(p, total_cnt);
        let nd = &nodes[n_idx];
        let inc_1b: Vec<(Vec<u16>, f64)> = paths
            .iter()
            .zip(psi_vec.iter())
            .map(|(pt, &pv)| (pt.iter().map(|&x| x + 1).collect(), pv))
            .collect();
        evts.push(EventRecord {
            gene_id,
            node_1based: (n_idx as u32) + 1,
            chrom_id: gene.chrom_id,
            coord_start: nd.genomic_coord,
            coord_end: nd.genomic_coord + nd.len.saturating_sub(1),
            is_pos_strand: gene.is_pos_strand(),
            event_type: motif.as_str(),
            psi: Some(p),
            ci_lo: Some(ci_lo),
            ci_hi: Some(ci_hi),
            total_reads: total_cnt,
            complexity: comp_k,
            entropy: ent,
            inc_paths: inc_1b,
            exc_paths: Vec::new(),
            edges: Vec::new(),
            inf_var: Some(inf_var),
        });
    }
    // Sort UTR rows in ascending node order for deterministic `.psi.gz` output
    evts.sort_by_key(|e| e.node_1based);
    (evts, next_i)
}

#[allow(clippy::too_many_arguments)]
fn build_event_record(
    gene_id: u32,
    gene: &crate::index::GeneRecord,
    nodes: &[NodeRecord],
    node_idx: usize,
    event_type: &'static str,
    psi: f64,
    ci_lo: f64,
    ci_hi: f64,
    inf_var: f64,
    total_reads: f64,
    inc_paths: &[ReducedPath],
    exc_paths: &[ReducedPath],
    local_edges: &[PathEdge],
) -> EventRecord {
    let nd = &nodes[node_idx];
    let n_paths = inc_paths.len() + exc_paths.len();
    let all_probs: Vec<f64> = inc_paths
        .iter()
        .chain(exc_paths.iter())
        .map(|p| p.psi)
        .collect();
    let ent = if all_probs.iter().sum::<f64>() > 0.0 {
        shannon_entropy(&all_probs)
    } else {
        shannon_entropy(&[psi, 1.0 - psi])
    };

    EventRecord {
        gene_id,
        node_1based: (node_idx as u32) + 1,
        chrom_id: gene.chrom_id,
        coord_start: nd.genomic_coord,
        coord_end: nd.genomic_coord + nd.len.saturating_sub(1),
        is_pos_strand: gene.is_pos_strand(),
        event_type,
        psi: Some(psi),
        ci_lo: Some(ci_lo),
        ci_hi: Some(ci_hi),
        total_reads,
        complexity: complexity(n_paths),
        entropy: ent,
        inc_paths: inc_paths
            .iter()
            .map(|p| (p.nodes.iter().map(|&x| x + 1).collect(), p.psi))
            .collect(),
        exc_paths: exc_paths
            .iter()
            .map(|p| (p.nodes.iter().map(|&x| x + 1).collect(), p.psi))
            .collect(),
        edges: local_edges
            .iter()
            .map(|e| (e.src + 1, e.dst + 1, e.count))
            .collect(),
        inf_var: Some(inf_var),
    }
}

fn empty_event(
    gene_id: u32,
    gene: &crate::index::GeneRecord,
    nodes: &[NodeRecord],
    node_idx: usize,
    event_type: &'static str,
) -> EventRecord {
    let nd = &nodes[node_idx];
    EventRecord {
        gene_id,
        node_1based: (node_idx as u32) + 1,
        chrom_id: gene.chrom_id,
        coord_start: nd.genomic_coord,
        coord_end: nd.genomic_coord + nd.len.saturating_sub(1),
        is_pos_strand: gene.is_pos_strand(),
        event_type,
        psi: None,
        ci_lo: None,
        ci_hi: None,
        total_reads: 0.0,
        complexity: 0,
        entropy: 0.0,
        inc_paths: Vec::new(),
        exc_paths: Vec::new(),
        edges: Vec::new(),
        inf_var: None,
    }
}

/// 90% Beta posterior credible interval `[0.05, 0.95]` matching `Whippet2.jl::beta_posterior_ci`.
#[inline]
pub fn beta_posterior_ci_90(p: f64, n: f64) -> (f64, f64) {
    let a = (p * n + 1.0).max(1e-6);
    let b = ((1.0 - p) * n + 1.0).max(1e-6);
    if let Ok(dist) = Beta::new(a, b) {
        (dist.inverse_cdf(0.05), dist.inverse_cdf(0.95))
    } else {
        (0.0, 1.0)
    }
}

/// Beta posterior variance `alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))`.
#[inline]
pub fn beta_posterior_var(p: f64, n: f64) -> f64 {
    let a = (p * n + 1.0).max(1e-6);
    let b = ((1.0 - p) * n + 1.0).max(1e-6);
    let s = a + b;
    (a * b) / (s * s * (s + 1.0))
}

/// Variational Dirichlet/Beta inferential variance on `Psi` (propagated into `.psi.gz` column 15).
fn compute_inferential_var(
    psi: f64,
    total_reads: f64,
    inc_counts: &[f64],
    exc_counts: &[f64],
    ambig_reads: f64,
    k_n: usize,
) -> f64 {
    if k_n <= 1 || ambig_reads <= 0.0 {
        return beta_posterior_var(psi, total_reads);
    }
    let a_inc = inc_counts.iter().sum::<f64>() + 1.0;
    let a_exc = exc_counts.iter().sum::<f64>() + 1.0;
    let s = a_inc + a_exc;
    let base_var = (a_inc * a_exc) / (s * s * (s + 1.0));
    // Inflate variance modestly by fraction of path-ambiguous reads (up to 2x) to reflect
    // inferential assignment uncertainty across complex K_n >= 2 splice graphs.
    let ambig_frac = (ambig_reads / total_reads.max(1.0)).clamp(0.0, 1.0);
    base_var * (1.0 + ambig_frac)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::GeneRecord;
    use bytemuck::Zeroable;

    fn make_dummy_gene_and_nodes(n_nodes: usize) -> (GeneRecord, Vec<NodeRecord>) {
        let mut gene = GeneRecord::zeroed();
        gene.flags = GeneRecord::FLAG_POS_STRAND;
        let nodes: Vec<NodeRecord> = (0..n_nodes)
            .map(|i| {
                let mut nd = NodeRecord::zeroed();
                nd.genomic_coord = (i as u32) * 1000 + 100;
                nd.len = 100;
                nd
            })
            .collect();
        (gene, nodes)
    }

    #[test]
    fn test_ltbp2_multi_exon_skip_outward_bridge_sweeps() {
        // Model the LTBP2 (ENSG00000119681) multi-exon skip topology around 1-based node 18
        // (0-based index 17):
        //   - High-depth constitutive chain: 14->15->16->17->18->19->...->34 (~650 reads each)
        //   - 1-read multi-exon skipping junction: 15->34 (0-based 14->33)
        // Verify:
        //   1. Upstream edge (14, 15) outside [15, 34] is NOT attached to exc_paths,
        //   2. inc_paths chains 15-16-17-18-...-34 cleanly back to 15,
        //   3. PSI stays > 0.99 (instead of dropping to ~0.65 under the old single forward pass),
        //   4. Even when an alt-TSS edge (1, 18) expands min_v to 1 so (14, 15) lies inside [1, 34],
        //      BOTH inc_paths and exc_paths share (14, 15) symmetrically and PSI still stays > 0.99.
        let (gene, nodes) = make_dummy_gene_and_nodes(36);
        let cfg = EventConfig::default();

        // 0-based indices for 1-based nodes 14..=34 are 13..=33; target node 18 is 0-based 17.
        let mut edges: Vec<PathEdge> = Vec::new();
        for u in 13u16..33u16 {
            edges.push(PathEdge {
                src: u,
                dst: u + 1,
                count: 650.0,
                eff_len: 1.0,
            });
        }
        // 1-read multi-exon skip 15 -> 34 (0-based 14 -> 33)
        edges.push(PathEdge {
            src: 14,
            dst: 33,
            count: 1.0,
            eff_len: 1.0,
        });
        edges.sort_unstable_by_key(|e| (e.src, e.dst));

        let ev = quantify_spliced_node(0, &gene, &nodes, &edges, 650.0, 17, EdgeMotif::CE, &cfg);
        let psi = ev.psi.expect("valid PSI");
        assert!(
            psi > 0.99,
            "expected PSI > 0.99 for 650 vs 1-read skip, got {psi:.4} (inc={:?}, exc={:?})",
            ev.inc_paths,
            ev.exc_paths
        );
        // (14, 15) is outside [15, 34] and must NOT appear in exc_paths or inc_paths
        assert!(
            ev.exc_paths.iter().all(|(p, _)| !p.contains(&14)),
            "upstream edge (14,15) outside [15,34] leaked into exc_paths: {:?}",
            ev.exc_paths
        );
        // inc_paths must chain cleanly all the way back to 1-based node 15 and forward to 34
        assert_eq!(ev.inc_paths.len(), 1);
        assert_eq!(ev.inc_paths[0].0.first().copied(), Some(15));
        assert_eq!(ev.inc_paths[0].0.last().copied(), Some(34));
        assert_eq!(ev.exc_paths.len(), 1);
        assert_eq!(ev.exc_paths[0].0, vec![15, 34]);

        // Case 2: add alt-TSS edge 1 -> 18 (0-based 0 -> 17, 40 reads) so min_v = 1 and
        // (14, 15) lies inside [1, 34]. Both inc and exc must share (14, 15) symmetrically!
        let mut edges_tss = edges.clone();
        edges_tss.push(PathEdge {
            src: 0,
            dst: 17,
            count: 40.0,
            eff_len: 1.0,
        });
        edges_tss.sort_unstable_by_key(|e| (e.src, e.dst));

        let ev2 = quantify_spliced_node(0, &gene, &nodes, &edges_tss, 650.0, 17, EdgeMotif::CE, &cfg);
        let psi2 = ev2.psi.expect("valid PSI with alt-TSS");
        assert!(
            psi2 > 0.99,
            "expected PSI > 0.99 with alt-TSS + shared (14,15), got {psi2:.4} (inc={:?}, exc={:?})",
            ev2.inc_paths,
            ev2.exc_paths
        );
        assert!(
            ev2.inc_paths.iter().any(|(p, _)| p.first().copied() == Some(14) && p.last().copied() == Some(34)),
            "expected inclusion path starting at 14 sharing (14,15), got {:?}",
            ev2.inc_paths
        );
    }
}
