//! Contiguous Splice Graph (CSG) builder implementing exact Whippet2 topology,
//! canonical splice-site dinucleotide annotation, annotated junction extraction,
//! and contiguous-exonic 16 bp 2-bit splice-flank encoding.

use crate::annotation::{RefGene, RefTx};
use crate::error::{Result, LangurError};
use crate::graph::edge_type::EdgeType;
use crate::graph::gc::{expected_gc, GC_BINS, GC_WINDOW_LEN};
use crate::seq::{
    append_normalized, append_revcomp_normalized, encode_flank_lsb_from_end, encode_flank_msb,
    ASCII_COMPLEMENT, ASCII_NORM,
};
use rustc_hash::{FxHashMap, FxHashSet};

/// Maximum CSG nodes per gene (`u16` local node index, reserving `0xFFFF`).
pub const MAX_NODES_PER_GENE: usize = (u16::MAX as usize) - 1;

/// Canonical dinucleotide class for splice-site scoring (2-bit code).
#[repr(u8)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum AcceptorDinuc {
    Other = 0,
    Ag = 1,
    Ac = 2,
}

#[repr(u8)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum DonorDinuc {
    Other = 0,
    Gt = 1,
    Gc = 2,
    At = 3,
}

/// Intermediate built CSG node before packing into flat `.lidx` arrays.
#[derive(Clone, Debug)]
pub struct BuiltNode {
    pub gene_seq_offset: u32,
    pub genomic_coord: u32,
    pub len: u32,
    pub left_edge: EdgeType,
    pub right_edge: EdgeType,
    pub acceptor_flank_2bit: u32,
    pub acceptor_flank_len: u8,
    pub donor_flank_2bit: u32,
    pub donor_flank_len: u8,
    pub acc_dinuc: AcceptorDinuc,
    pub don_dinuc: DonorDinuc,
    pub is_novel_acc: bool,
    pub is_novel_don: bool,
    pub has_ambig: bool,
}

/// Intermediate built annotated transcript path.
#[derive(Clone, Debug)]
pub struct BuiltTx {
    pub tx_id: String,
    pub tx_start: u32,
    pub tx_end: u32,
    pub is_basic: bool,
    pub tsl: u8,
    /// Sorted 0-based local node IDs traversed by this transcript.
    pub nodes: Vec<u16>,
    /// Cumulative spliced length *before* each node in `nodes`.
    pub node_prefix_len: Vec<u32>,
    /// Exact spliced sequence length in bp.
    pub seq_len: u32,
    /// Normalized 20-bin expected GC distribution.
    pub expected_gc: [f32; GC_BINS],
}

/// Intermediate built annotated splice junction (local node pair `src_node < dst_node`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BuiltAnnoJunction {
    pub src_node: u16,
    pub dst_node: u16,
    pub tx_support: u16,
    pub flags: u16,
}

/// Complete built splice graph for a single gene.
#[derive(Clone, Debug)]
pub struct BuiltSpliceGraph {
    pub gene_id: String,
    pub gene_name: String,
    pub chrom: String,
    pub strand: bool,
    pub genomic_start: u32,
    pub genomic_end: u32,
    pub mean_tx_len: f32,
    pub has_ambig: bool,
    pub has_bam_novel: bool,
    /// Uppercase normalized (`A/C/G/T/N`) concatenated CSG sequence in transcript 5'->3' order.
    pub seq: Vec<u8>,
    pub nodes: Vec<BuiltNode>,
    /// Length `nodes.len() + 1` (matches Whippet2 `sg.edgetype`).
    pub edgetypes: Vec<EdgeType>,
    pub transcripts: Vec<BuiltTx>,
    pub anno_junctions: Vec<BuiltAnnoJunction>,
}

/// Build the Contiguous Splice Graph for a single `RefGene` given its chromosome sequence.
pub fn build_splice_graph(gene: &RefGene, chrom_seq: &[u8]) -> Result<BuiltSpliceGraph> {
    let strand = gene.strand;
    let chrom_len = chrom_seq.len() as u32;

    let novel_acc_set: FxHashSet<u32> = gene.novel_acc.iter().copied().collect();
    let novel_don_set: FxHashSet<u32> = gene.novel_don.iter().copied().collect();

    // 4-way merge indices: 0=txen, 1=acc, 2=don, 3=txst (exact Whippet2 tie-breaking order!)
    let mut idx = [0usize; 4];
    let get_min = |idx: &[usize; 4]| -> (usize, u32) {
        let vals = [
            gene.txen.get(idx[0]).copied().unwrap_or(u32::MAX),
            gene.acc.get(idx[1]).copied().unwrap_or(u32::MAX),
            gene.don.get(idx[2]).copied().unwrap_or(u32::MAX),
            gene.txst.get(idx[3]).copied().unwrap_or(u32::MAX),
        ];
        let mut best_i = 0;
        let mut best_v = vals[0];
        for (i, &v) in vals.iter().enumerate().skip(1) {
            if v < best_v {
                best_v = v;
                best_i = i;
            }
        }
        (best_i, best_v)
    };

    // Build nodes in genomic left-to-right order first, then reverse for `-` strand.
    // This produces identical order to Whippet2's `stranded_push!` (push_front on `-` strand) in O(N).
    struct RawNode {
        genomic_coord: u32,
        len: u32,
        edge: EdgeType, // edge on the genomic-left side (becomes left_edge on `+`, right_edge before shift on `-`)
        seq_slice: Vec<u8>,
        has_ambig: bool,
        is_novel_acc: bool,
        is_novel_don: bool,
    }

    let mut raw_nodes: Vec<RawNode> = Vec::new();
    let mut term_edge: Option<EdgeType> = None;

    while idx[0] < gene.txen.len()
        || idx[1] < gene.acc.len()
        || idx[2] < gene.don.len()
        || idx[3] < gene.txst.len()
    {
        let (min_idx, min_val) = get_min(&idx);
        idx[min_idx] += 1;
        let (sec_idx, sec_val) = get_min(&idx);

        if sec_val == u32::MAX {
            term_edge = Some(if strand { EdgeType::RS } else { EdgeType::SL });
            break;
        }

        let (push_val, node_size, edge, is_nov_acc, is_nov_don) =
            if is_subinterval(&gene.exons, min_val, sec_val)
                || (min_idx == 1 && novel_acc_set.contains(&min_val))
                || (sec_idx == 2 && novel_don_set.contains(&sec_val))
            {
                let left_adj: u32 = if (min_idx == 0 || min_idx == 2) && min_val != sec_val {
                    1
                } else {
                    0
                };
                let right_adj: u32 = if (sec_idx == 1 || sec_idx == 3) && min_val != sec_val {
                    1
                } else {
                    0
                };
                let g_start = min_val + left_adj;
                let size_i64 =
                    (sec_val as i64 - min_val as i64) - (left_adj as i64) - (right_adj as i64) + 1;
                let size = size_i64.max(0) as u32;
                let e = EdgeType::from_coord_indices(min_idx, sec_idx, true, strand);
                (
                    g_start,
                    size,
                    e,
                    min_idx == 1 && novel_acc_set.contains(&min_val),
                    sec_idx == 2 && novel_don_set.contains(&sec_val),
                )
            } else {
                // Sequence gap (intron): advance past sec_idx
                idx[sec_idx] += 1;
                let (thr_idx, thr_val) = get_min(&idx);
                if thr_val == u32::MAX {
                    term_edge = Some(if strand { EdgeType::RS } else { EdgeType::SL });
                    break;
                }
                let right_adj: u32 = if (thr_idx == 1 || thr_idx == 3) && sec_val != thr_val {
                    1
                } else {
                    0
                };
                let g_start = sec_val;
                let size_i64 = (thr_val as i64 - sec_val as i64) - (right_adj as i64) + 1;
                let size = size_i64.max(0) as u32;
                let e = EdgeType::from_coord_indices(min_idx, sec_idx, false, strand);
                (
                    g_start,
                    size,
                    e,
                    sec_idx == 1 && novel_acc_set.contains(&sec_val),
                    thr_idx == 2 && novel_don_set.contains(&thr_val),
                )
            };

        let mut node_seq = Vec::with_capacity(node_size as usize);
        let has_ambig = extract_genomic_slice(chrom_seq, push_val, node_size, strand, &mut node_seq)?;

        raw_nodes.push(RawNode {
            genomic_coord: push_val,
            len: node_size,
            edge,
            seq_slice: node_seq,
            has_ambig,
            is_novel_acc: is_nov_acc,
            is_novel_don: is_nov_don,
        });
    }

    if raw_nodes.len() > MAX_NODES_PER_GENE {
        return Err(LangurError::TooManyNodes {
            gene_id: gene.gene_id.clone(),
            count: raw_nodes.len(),
            max: MAX_NODES_PER_GENE,
        });
    }

    // Reconstruct exact Whippet2 `edgetype` vector (length = num_nodes + 1) and node order.
    // On `+` strand: push! appends each node's edge then term_edge at the end.
    // On `-` strand: pushfirst! reverses the node list and puts term_edge (SL) at index 0,
    // followed by the reversed node edges!
    let n_nodes = raw_nodes.len();
    let mut edgetypes: Vec<EdgeType> = Vec::with_capacity(n_nodes + 1);
    if strand {
        for rn in &raw_nodes {
            edgetypes.push(rn.edge);
        }
        edgetypes.push(term_edge.unwrap_or(EdgeType::RS));
    } else {
        edgetypes.push(term_edge.unwrap_or(EdgeType::SL));
        for rn in raw_nodes.iter().rev() {
            edgetypes.push(rn.edge);
        }
        raw_nodes.reverse();
    }

    // Concatenate CSG sequence and compute node offsets
    let total_seq_len: usize = raw_nodes.iter().map(|r| r.seq_slice.len()).sum();
    let mut seq = Vec::with_capacity(total_seq_len);
    let mut gene_has_ambig = false;
    let mut nodes = Vec::with_capacity(n_nodes);

    for (i, rn) in raw_nodes.into_iter().enumerate() {
        let seq_off = seq.len() as u32;
        seq.extend_from_slice(&rn.seq_slice);
        if rn.has_ambig {
            gene_has_ambig = true;
        }
        let left_edge = edgetypes[i];
        let right_edge = edgetypes[i + 1];

        let acc_dinuc = if left_edge.is_acceptor() {
            inspect_acceptor_dinuc(chrom_seq, chrom_len, rn.genomic_coord, rn.len, strand)
        } else {
            AcceptorDinuc::Other
        };
        let don_dinuc = if right_edge.is_donor() {
            inspect_donor_dinuc(chrom_seq, chrom_len, rn.genomic_coord, rn.len, strand)
        } else {
            DonorDinuc::Other
        };

        // Note: on `-` strand, genomic-left novel_acc corresponds to transcript-right donor,
        // and genomic-right novel_don corresponds to transcript-left acceptor.
        let (is_nov_acc, is_nov_don) = if strand {
            (rn.is_novel_acc, rn.is_novel_don)
        } else {
            (rn.is_novel_don, rn.is_novel_acc)
        };

        nodes.push(BuiltNode {
            gene_seq_offset: seq_off,
            genomic_coord: rn.genomic_coord,
            len: rn.len,
            left_edge,
            right_edge,
            acceptor_flank_2bit: 0,
            acceptor_flank_len: 0,
            donor_flank_2bit: 0,
            donor_flank_len: 0,
            acc_dinuc,
            don_dinuc,
            is_novel_acc: is_nov_acc,
            is_novel_don: is_nov_don,
            has_ambig: rn.has_ambig,
        });
    }

    // Compute contiguous-exonic 16 bp 2-bit splice flanks (Refinement 1):
    // Extend across genomically contiguous adjacent CSG nodes (e.g. 3 bp AA/AD nodes)
    // without bleeding across hard intron gaps, and store exact valid flank lengths.
    for i in 0..nodes.len() {
        if nodes[i].left_edge.is_acceptor() {
            let start = nodes[i].gene_seq_offset as usize;
            let mut end = start + (nodes[i].len as usize);
            let mut cur = i;
            while end - start < 16 && cur + 1 < nodes.len() && are_nodes_contiguous(&nodes[cur], &nodes[cur + 1], strand) {
                cur += 1;
                end += nodes[cur].len as usize;
            }
            let flank_slice = &seq[start..end.min(start + 16)];
            let (packed, flen) = encode_flank_msb(flank_slice);
            nodes[i].acceptor_flank_2bit = packed;
            nodes[i].acceptor_flank_len = flen;
        }
        if nodes[i].right_edge.is_donor() {
            let end = (nodes[i].gene_seq_offset + nodes[i].len) as usize;
            let mut start = nodes[i].gene_seq_offset as usize;
            let mut cur = i;
            while end - start < 16 && cur > 0 && are_nodes_contiguous(&nodes[cur - 1], &nodes[cur], strand) {
                cur -= 1;
                start = nodes[cur].gene_seq_offset as usize;
            }
            let flank_slice = &seq[start.max(end.saturating_sub(16))..end];
            let (packed, flen) = encode_flank_lsb_from_end(flank_slice);
            nodes[i].donor_flank_2bit = packed;
            nodes[i].donor_flank_len = flen;
        }
    }

    // Build annotated transcript paths, prefix lengths, ExpectedGC, and annotated splice junctions
    let mut transcripts = Vec::with_capacity(gene.reftx.len());
    let mut jnc_counts: FxHashMap<(u16, u16), u16> = FxHashMap::default();
    let mut tx_seq_buf: Vec<u8> = Vec::new();

    for tx in &gene.reftx {
        let tx_nodes = build_annotated_path(&nodes, tx, strand);
        let mut prefix_lens = Vec::with_capacity(tx_nodes.len());
        tx_seq_buf.clear();
        let mut running_len: u32 = 0;

        for &nid in &tx_nodes {
            prefix_lens.push(running_len);
            let nd = &nodes[nid as usize];
            let s = nd.gene_seq_offset as usize;
            let e = s + (nd.len as usize);
            tx_seq_buf.extend_from_slice(&seq[s..e]);
            running_len += nd.len;
        }

        // Extract annotated splice junctions between consecutive transcript nodes
        // that cross an intron (non-contiguous nodes, or hard donor->acceptor edge).
        for win in tx_nodes.windows(2) {
            let u = win[0];
            let v = win[1];
            if v > u
                && nodes[u as usize].right_edge.is_donor()
                && nodes[v as usize].left_edge.is_acceptor()
                && !are_nodes_contiguous(&nodes[u as usize], &nodes[v as usize], strand)
            {
                let entry = jnc_counts.entry((u, v)).or_insert(0);
                *entry = entry.saturating_add(1);
            }
        }

        let gc = expected_gc(&tx_seq_buf, GC_WINDOW_LEN);
        transcripts.push(BuiltTx {
            tx_id: tx.tx_id.clone(),
            tx_start: tx.tx_start,
            tx_end: tx.tx_end,
            is_basic: tx.is_basic,
            tsl: tx.tsl,
            nodes: tx_nodes,
            node_prefix_len: prefix_lens,
            seq_len: running_len,
            expected_gc: gc,
        });
    }

    let mut anno_junctions: Vec<BuiltAnnoJunction> = jnc_counts
        .into_iter()
        .map(|((src, dst), support)| {
            let don_d = nodes[src as usize].don_dinuc as u16;
            let acc_d = nodes[dst as usize].acc_dinuc as u16;
            let dinuc_class: u16 = match (don_d, acc_d) {
                (1, 1) => 1, // GT-AG
                (2, 1) => 2, // GC-AG
                (3, 2) => 3, // AT-AC
                _ => 0,
            };
            let flags = 0x0001u16 | (dinuc_class << 2);
            BuiltAnnoJunction {
                src_node: src,
                dst_node: dst,
                tx_support: support,
                flags,
            }
        })
        .collect();
    anno_junctions.sort_unstable_by_key(|j| (j.src_node, j.dst_node));

    let (genomic_start, genomic_end) = gene.genomic_span();
    let has_bam_novel = !gene.novel_acc.is_empty() || !gene.novel_don.is_empty();

    Ok(BuiltSpliceGraph {
        gene_id: gene.gene_id.clone(),
        gene_name: gene.gene_name.clone(),
        chrom: gene.chrom.clone(),
        strand,
        genomic_start,
        genomic_end,
        mean_tx_len: gene.mean_tx_len,
        has_ambig: gene_has_ambig,
        has_bam_novel,
        seq,
        nodes,
        edgetypes,
        transcripts,
        anno_junctions,
    })
}

/// Returns `true` iff two consecutive CSG nodes `(prev, next)` in 5'->3' transcript order
/// are immediately adjacent with zero intervening intronic bases on the genome.
#[inline]
pub fn are_nodes_contiguous(prev: &BuiltNode, next: &BuiltNode, strand: bool) -> bool {
    if strand {
        prev.genomic_coord + prev.len == next.genomic_coord
    } else {
        next.genomic_coord + next.len == prev.genomic_coord
    }
}

/// Returns `true` iff some exon interval `(s, e)` fully contains `[sub_start, sub_end]`.
fn is_subinterval(exons: &[(u32, u32)], sub_start: u32, sub_end: u32) -> bool {
    // `exons` is sorted by `(start, end)`. Only exons with `s <= sub_start` can contain it.
    let upper = exons.partition_point(|&(s, _)| s <= sub_start);
    for &(s, e) in &exons[..upper] {
        if s <= sub_start && e >= sub_end {
            return true;
        }
    }
    false
}

fn extract_genomic_slice(
    chrom_seq: &[u8],
    coord_1based: u32,
    len: u32,
    strand: bool,
    dst: &mut Vec<u8>,
) -> Result<bool> {
    if len == 0 {
        return Ok(false);
    }
    let start0 = (coord_1based as usize).saturating_sub(1);
    let end0 = start0 + (len as usize);
    if end0 > chrom_seq.len() {
        return Err(LangurError::Fasta(format!(
            "coordinate range {coord_1based}..={} exceeds chromosome length {}",
            coord_1based + len - 1,
            chrom_seq.len()
        )));
    }
    let raw = &chrom_seq[start0..end0];
    if strand {
        Ok(append_normalized(raw, dst))
    } else {
        Ok(append_revcomp_normalized(raw, dst))
    }
}

/// Build sorted 0-based local CSG node IDs for an annotated transcript,
/// matching `Whippet2.jl/src/graph.jl::build_annotated_path`.
fn build_annotated_path(nodes: &[BuiltNode], tx: &RefTx, strand: bool) -> Vec<u16> {
    let mut path: Vec<u16> = Vec::new();
    if nodes.is_empty() {
        return path;
    }

    for (&ex_acc, &ex_don) in tx.acc.iter().zip(tx.don.iter()) {
        if strand {
            // `nodes` is sorted ascending by genomic_coord
            let start_idx = nodes.partition_point(|n| n.genomic_coord < ex_acc);
            for (idx, nd) in nodes[start_idx..].iter().enumerate() {
                if nd.genomic_coord > ex_don {
                    break;
                }
                path.push((start_idx + idx) as u16);
            }
        } else {
            // `nodes` is sorted descending by genomic_coord
            let start_idx = nodes.partition_point(|n| n.genomic_coord > ex_don);
            for (idx, nd) in nodes[start_idx..].iter().enumerate() {
                if nd.genomic_coord < ex_acc {
                    break;
                }
                path.push((start_idx + idx) as u16);
            }
        }
    }

    path.sort_unstable();
    path.dedup();
    path
}

fn inspect_acceptor_dinuc(
    chrom_seq: &[u8],
    chrom_len: u32,
    coord_1based: u32,
    len: u32,
    strand: bool,
) -> AcceptorDinuc {
    if len == 0 {
        return AcceptorDinuc::Other;
    }
    let (b0, b1) = if strand {
        // 2 bp immediately upstream of coord_1based: 0-based [coord-3, coord-2]
        if coord_1based < 3 {
            return AcceptorDinuc::Other;
        }
        let s = (coord_1based - 3) as usize;
        (ASCII_NORM[chrom_seq[s] as usize], ASCII_NORM[chrom_seq[s + 1] as usize])
    } else {
        // On `-` strand, 5' end of node is genomic `coord + len - 1`; upstream intron is `[coord+len-1, coord+len]` (0-based)
        let e = coord_1based + len - 1;
        if e + 2 > chrom_len {
            return AcceptorDinuc::Other;
        }
        let s = e as usize; // 0-based index of (e + 1)
        let c0 = ASCII_NORM[ASCII_COMPLEMENT[chrom_seq[s + 1] as usize] as usize];
        let c1 = ASCII_NORM[ASCII_COMPLEMENT[chrom_seq[s] as usize] as usize];
        (c0, c1)
    };
    match (b0, b1) {
        (b'A', b'G') => AcceptorDinuc::Ag,
        (b'A', b'C') => AcceptorDinuc::Ac,
        _ => AcceptorDinuc::Other,
    }
}

fn inspect_donor_dinuc(
    chrom_seq: &[u8],
    chrom_len: u32,
    coord_1based: u32,
    len: u32,
    strand: bool,
) -> DonorDinuc {
    if len == 0 {
        return DonorDinuc::Other;
    }
    let (b0, b1) = if strand {
        let e = coord_1based + len - 1;
        if e + 2 > chrom_len {
            return DonorDinuc::Other;
        }
        let s = e as usize;
        (ASCII_NORM[chrom_seq[s] as usize], ASCII_NORM[chrom_seq[s + 1] as usize])
    } else {
        if coord_1based < 3 {
            return DonorDinuc::Other;
        }
        let s = (coord_1based - 3) as usize;
        let c0 = ASCII_NORM[ASCII_COMPLEMENT[chrom_seq[s + 1] as usize] as usize];
        let c1 = ASCII_NORM[ASCII_COMPLEMENT[chrom_seq[s] as usize] as usize];
        (c0, c1)
    };
    match (b0, b1) {
        (b'G', b'T') => DonorDinuc::Gt,
        (b'G', b'C') => DonorDinuc::Gc,
        (b'A', b'T') => DonorDinuc::At,
        _ => DonorDinuc::Other,
    }
}
