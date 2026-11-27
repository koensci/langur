//! Parallel minimizer extraction across CSG nodes and annotated splice-junction
//! bridges, plus 2-level radix-bucketed table construction with singleton inlining
//! and collision detection for `k > 27`.

use crate::graph::BuiltSpliceGraph;
use crate::index::format::{pack_multi_slot, PackedSeedHit};
use crate::seq::{extract_minimizers, is_bijective_hash, Minimizer};
use rayon::prelude::*;

/// Raw unsorted minimizer posting emitted during index construction.
#[derive(Copy, Clone, Debug)]
pub struct RawPosting {
    pub hash: u64,
    pub canon_kmer: u64,
    pub hit: PackedSeedHit,
}

/// Extract CSG node minimizers and (optionally) annotated junction-bridge minimizers
/// across all splice graphs in parallel.
pub fn extract_csg_postings(
    graphs: &[BuiltSpliceGraph],
    k: usize,
    w: usize,
    hash_bits: u32,
    include_jnc_seeds: bool,
) -> Vec<RawPosting> {
    let per_gene: Vec<Vec<RawPosting>> = graphs
        .par_iter()
        .enumerate()
        .map(|(gene_idx, sg)| {
            let gene_id = gene_idx as u32;
            let mut postings = Vec::new();
            let mut min_buf: Vec<Minimizer> = Vec::new();
            let mut bridge_buf: Vec<u8> = Vec::new();

            // 1. Extract minimizers per CSG node (and across contiguous exonic runs so
            //    short nodes like 3 bp AA/AD nodes don't miss boundary-spanning seeds!)
            let mut run_start = 0usize;
            while run_start < sg.nodes.len() {
                let mut run_end = run_start;
                while run_end + 1 < sg.nodes.len()
                    && crate::graph::builder::are_nodes_contiguous(
                        &sg.nodes[run_end],
                        &sg.nodes[run_end + 1],
                        sg.strand,
                    )
                {
                    run_end += 1;
                }

                let seq_start = sg.nodes[run_start].gene_seq_offset as usize;
                let last_nd = &sg.nodes[run_end];
                let seq_end = (last_nd.gene_seq_offset + last_nd.len) as usize;
                let run_seq = &sg.seq[seq_start..seq_end];

                extract_minimizers(run_seq, k, w, hash_bits, &mut min_buf);
                let mut cur_node = run_start;
                for m in &min_buf {
                    let gene_off = (seq_start as u32) + m.pos;
                    while cur_node < run_end
                        && gene_off
                            >= sg.nodes[cur_node].gene_seq_offset + sg.nodes[cur_node].len
                    {
                        cur_node += 1;
                    }
                    let node_off = gene_off - sg.nodes[cur_node].gene_seq_offset;
                    postings.push(RawPosting {
                        hash: m.hash,
                        canon_kmer: m.canon_kmer,
                        hit: PackedSeedHit::csg_node(gene_id, cur_node as u16, node_off, m.is_rc),
                    });
                }

                run_start = run_end + 1;
            }

            // 2. Optionally extract annotated junction-bridge minimizers:
            //    (k-1) bp ending at src_node 3' end + (k-1) bp starting at dst_node 5' start.
            if include_jnc_seeds && k > 1 {
                let flank = k - 1;
                for (jnc_idx, jnc) in sg.anno_junctions.iter().enumerate() {
                    let src = jnc.src_node as usize;
                    let dst = jnc.dst_node as usize;

                    // Walk backward from src 3' end across contiguous exonic nodes up to `flank` bp
                    let src_end = (sg.nodes[src].gene_seq_offset + sg.nodes[src].len) as usize;
                    let mut left_start = sg.nodes[src].gene_seq_offset as usize;
                    let mut cur = src;
                    while src_end - left_start < flank
                        && cur > 0
                        && crate::graph::builder::are_nodes_contiguous(
                            &sg.nodes[cur - 1],
                            &sg.nodes[cur],
                            sg.strand,
                        )
                    {
                        cur -= 1;
                        left_start = sg.nodes[cur].gene_seq_offset as usize;
                    }
                    let left_slice = &sg.seq[left_start.max(src_end.saturating_sub(flank))..src_end];

                    // Walk forward from dst 5' start across contiguous exonic nodes up to `flank` bp
                    let dst_start = sg.nodes[dst].gene_seq_offset as usize;
                    let mut right_end = dst_start + (sg.nodes[dst].len as usize);
                    let mut cur = dst;
                    while right_end - dst_start < flank
                        && cur + 1 < sg.nodes.len()
                        && crate::graph::builder::are_nodes_contiguous(
                            &sg.nodes[cur],
                            &sg.nodes[cur + 1],
                            sg.strand,
                        )
                    {
                        cur += 1;
                        right_end += sg.nodes[cur].len as usize;
                    }
                    let right_slice = &sg.seq[dst_start..right_end.min(dst_start + flank)];

                    if left_slice.len() + right_slice.len() < k {
                        continue;
                    }
                    bridge_buf.clear();
                    bridge_buf.extend_from_slice(left_slice);
                    bridge_buf.extend_from_slice(right_slice);

                    extract_minimizers(&bridge_buf, k, w, hash_bits, &mut min_buf);
                    let left_len = left_slice.len() as u32;
                    for m in &min_buf {
                        // Only keep minimizers that actually straddle the junction boundary;
                        // store `left_overhang = left_len - m.pos` (strictly in 1..k) so Phase 2
                        // recovers exact `jnc_pos = seed_pos + left_overhang` even when `left_len < k - 1`.
                        if m.pos < left_len && m.pos + (k as u32) > left_len {
                            let left_overhang = (left_len - m.pos) as u8;
                            postings.push(RawPosting {
                                hash: m.hash,
                                canon_kmer: m.canon_kmer,
                                hit: PackedSeedHit::anno_junction(
                                    gene_id,
                                    jnc_idx as u32,
                                    left_overhang,
                                    m.is_rc,
                                ),
                            });
                        }
                    }
                }
            }

            postings
        })
        .collect();

    let total: usize = per_gene.iter().map(Vec::len).sum();
    let mut flat = Vec::with_capacity(total);
    for mut v in per_gene {
        flat.append(&mut v);
    }
    flat
}

/// Built 2-level radix-bucketed minimizer table ready for `.lidx` serialization.
pub struct BuiltMinimizerTable {
    pub bucket_bits: u16,
    pub bucket_offsets: Vec<u32>,
    pub key_lo: Vec<u32>,
    pub slots: Vec<u64>,
    pub postings: Vec<PackedSeedHit>,
    pub num_hash_collisions: u64,
    pub hash_may_collide: bool,
}

/// Build the final 2-level radix-bucketed minimizer table from CSG and decoy postings.
///
/// - For `k <= 27` (when `2*k <= bucket_bits + 32`), hashes are strictly bijective (`num_hash_collisions == 0`).
/// - For `k > 27`, explicitly detects any distinct `canon_kmer` colliding on `(bucket, key_lo)`,
///   merges their postings cleanly (zero false negatives), counts collisions, and sets `hash_may_collide`.
pub fn build_minimizer_table(
    mut csg_postings: Vec<RawPosting>,
    mut decoy_postings: Vec<RawPosting>,
    k: usize,
    bucket_bits: u16,
    max_csg_occ: usize,
    max_decoy_occ: usize,
) -> BuiltMinimizerTable {
    let hash_bits = (bucket_bits as u32) + 32;
    let bijective = is_bijective_hash(k, hash_bits);
    let num_buckets = 1usize << (bucket_bits as usize);

    // Parallel sort CSG and decoy postings by `(hash, hit)`
    csg_postings.par_sort_unstable_by_key(|p| (p.hash, p.hit.0));
    decoy_postings.par_sort_unstable_by_key(|p| (p.hash, p.hit.0));

    let mut bucket_offsets = vec![0u32; num_buckets + 1];
    let mut key_lo: Vec<u32> = Vec::new();
    let mut slots: Vec<u64> = Vec::new();
    let mut postings: Vec<PackedSeedHit> = Vec::new();
    let mut num_hash_collisions: u64 = 0;

    let mut ci = 0usize;
    let mut di = 0usize;

    while ci < csg_postings.len() {
        let cur_hash = csg_postings[ci].hash;
        let c_start = ci;
        let first_canon = csg_postings[ci].canon_kmer;
        let mut saw_collision = false;

        while ci < csg_postings.len() && csg_postings[ci].hash == cur_hash {
            if !bijective && csg_postings[ci].canon_kmer != first_canon {
                saw_collision = true;
            }
            ci += 1;
        }

        // Advance decoy pointer to cur_hash (skip any decoy hash < cur_hash)
        while di < decoy_postings.len() && decoy_postings[di].hash < cur_hash {
            di += 1;
        }
        let d_start = di;
        while di < decoy_postings.len() && decoy_postings[di].hash == cur_hash {
            if !bijective && decoy_postings[di].canon_kmer != first_canon {
                saw_collision = true;
            }
            di += 1;
        }
        if saw_collision {
            num_hash_collisions += 1;
        }

        // Deduplicate exact identical PackedSeedHits in [c_start..ci] and [d_start..di]
        let c_raw = &csg_postings[c_start..ci];
        let d_raw = &decoy_postings[d_start..di];

        let mut c_uniq: Vec<PackedSeedHit> = Vec::with_capacity(c_raw.len().min(max_csg_occ));
        let mut c_trunc = false;
        let mut last_hit = u64::MAX;
        for p in c_raw {
            if p.hit.0 != last_hit {
                last_hit = p.hit.0;
                if c_uniq.len() < max_csg_occ {
                    c_uniq.push(p.hit);
                } else {
                    c_trunc = true;
                    break;
                }
            }
        }

        let mut d_uniq: Vec<PackedSeedHit> = Vec::with_capacity(d_raw.len().min(max_decoy_occ));
        let mut d_trunc = false;
        last_hit = u64::MAX;
        for p in d_raw {
            if p.hit.0 != last_hit {
                last_hit = p.hit.0;
                if d_uniq.len() < max_decoy_occ {
                    d_uniq.push(p.hit);
                } else {
                    d_trunc = true;
                    break;
                }
            }
        }

        let bucket = (cur_hash >> 32) as usize;
        let klo = cur_hash as u32;
        bucket_offsets[bucket + 1] += 1;
        key_lo.push(klo);

        // Singleton inline optimization: bit 63 == 0 when exactly 1 CSG hit and 0 decoy hits
        if c_uniq.len() == 1 && d_uniq.is_empty() && !c_trunc && !d_trunc {
            debug_assert!((c_uniq[0].0 >> 63) == 0);
            slots.push(c_uniq[0].0);
        } else {
            let off = postings.len() as u64;
            let c_cnt = c_uniq.len().min(u16::MAX as usize) as u16;
            let d_cnt = d_uniq.len().min(u8::MAX as usize) as u8;
            postings.extend_from_slice(&c_uniq[..c_cnt as usize]);
            postings.extend_from_slice(&d_uniq[..d_cnt as usize]);
            slots.push(pack_multi_slot(off, c_cnt, d_cnt, c_trunc, d_trunc));
        }
    }

    // Prefix sum bucket_offsets
    for i in 1..bucket_offsets.len() {
        bucket_offsets[i] += bucket_offsets[i - 1];
    }

    BuiltMinimizerTable {
        bucket_bits,
        bucket_offsets,
        key_lo,
        slots,
        postings,
        num_hash_collisions,
        hash_may_collide: !bijective,
    }
}

/// Choose an optimal `bucket_bits` (`10..=22`) that guarantees mathematical bijectivity
/// (`bucket_bits + 32 >= 2 * k`) for all `k <= 27`, while staying tiny for small test indexes.
pub fn choose_bucket_bits(k: usize, estimated_unique_keys: usize, override_bits: Option<u16>) -> u16 {
    if let Some(b) = override_bits {
        return b.clamp(10, 24);
    }
    let min_for_bijection: u16 = if k <= 27 {
        (2 * k).saturating_sub(32) as u16
    } else {
        22
    };
    let density_bits: u16 = if estimated_unique_keys < 4096 {
        10
    } else {
        // Target ~4 keys per bucket -> log2(keys / 4)
        (estimated_unique_keys.ilog2().saturating_sub(2)).clamp(10, 22) as u16
    };
    density_bits.max(min_for_bijection).min(22)
}
