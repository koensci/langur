//! Parallel RFC 1952 multi-member gzip writers (`.psi.gz`, `.gene.tpm.gz`,
//! `.isoform.tpm.gz`, `.jnc.gz`, `.map.gz`) and optional SAM alignment output.

use crate::align::ReadAlignment;
use crate::error::{Result, LangurError};
use crate::events::EventRecord;
use crate::index::MmapIndex;
use crate::quant::{FragLenDist, SampleMapStats, VbemResult};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Compress a sequence of raw byte blocks into valid concatenated RFC 1952 gzip members
/// in parallel via `rayon`, writing them sequentially in deterministic order.
fn write_parallel_gz(path: &Path, header: &[u8], blocks: Vec<Vec<u8>>) -> Result<()> {
    let file = File::create(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut out = BufWriter::with_capacity(2 * 1024 * 1024, file);

    // Header member
    let mut hdr_enc = GzEncoder::new(Vec::new(), Compression::fast());
    hdr_enc.write_all(header)?;
    let hdr_gz = hdr_enc.finish()?;
    out.write_all(&hdr_gz)?;

    // Parallel RFC 1952 body members
    let compressed: Vec<std::io::Result<Vec<u8>>> = blocks
        .into_par_iter()
        .filter(|b| !b.is_empty())
        .map(|raw| {
            let mut enc = GzEncoder::new(Vec::with_capacity(raw.len() / 3), Compression::fast());
            enc.write_all(&raw)?;
            enc.finish()
        })
        .collect();

    for res in compressed {
        let gz_bytes = res?;
        out.write_all(&gz_bytes)?;
    }
    out.flush()?;
    Ok(())
}

/// Write `.psi.gz` with exact 14 Whippet2 columns + column 15 `Inf_Var`.
pub fn write_psi_gz(path: &Path, idx: &MmapIndex, events: &[EventRecord]) -> Result<()> {
    let header = b"Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\tComplexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\tInf_Var\n";
    let chunk_size = 2048usize;

    let blocks: Vec<Vec<u8>> = events
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut buf = String::with_capacity(chunk.len() * 128);
            for ev in chunk {
                let gname = idx.get_str(idx.genes()[ev.gene_id as usize].gene_id);
                let chrom = idx.get_str(idx.chroms()[ev.chrom_id as usize].name);
                let strand = if ev.is_pos_strand { '+' } else { '-' };

                let _ = write!(
                    buf,
                    "{gname}\t{}\t{chrom}:{}-{}\t{strand}\t{}\t",
                    ev.node_1based, ev.coord_start, ev.coord_end, ev.event_type
                );

                match (ev.psi, ev.ci_lo, ev.ci_hi, ev.inf_var) {
                    (Some(psi), Some(lo), Some(hi), Some(iv)) => {
                        let width = (hi - lo).max(0.0);
                        let _ = write!(
                            buf,
                            "{psi:.4}\t{width:.4}\t{lo:.4},{hi:.4}\t{:.1}\tK{}\t{:.3}\t",
                            ev.total_reads, ev.complexity, ev.entropy
                        );
                        write_path_list(&mut buf, &ev.inc_paths);
                        buf.push('\t');
                        write_path_list(&mut buf, &ev.exc_paths);
                        buf.push('\t');
                        write_edge_list(&mut buf, &ev.edges);
                        let _ = writeln!(buf, "\t{iv:.6e}");
                    }
                    _ => {
                        buf.push_str("NA\tNA\tNA\t0.0\tK0\t0.000\tNA\tNA\tNA\tNA\n");
                    }
                }
            }
            buf.into_bytes()
        })
        .collect();

    write_parallel_gz(path, header, blocks)
}

fn write_path_list(buf: &mut String, paths: &[(Vec<u16>, f64)]) {
    if paths.is_empty() {
        buf.push_str("NA");
        return;
    }
    for (pi, (nodes, psi)) in paths.iter().enumerate() {
        if pi > 0 {
            buf.push(',');
        }
        for (ni, &n) in nodes.iter().enumerate() {
            if ni > 0 {
                buf.push('-');
            }
            let _ = write!(buf, "{n}");
        }
        let _ = write!(buf, ":{psi:.4}");
    }
}

fn write_edge_list(buf: &mut String, edges: &[(u16, u16, f64)]) {
    if edges.is_empty() {
        buf.push_str("NA");
        return;
    }
    for (ei, &(u, v, c)) in edges.iter().enumerate() {
        if ei > 0 {
            buf.push(',');
        }
        let _ = write!(buf, "{u}-{v}:{c:.1}");
    }
}

/// Write `.gene.tpm.gz` (`Gene\tTpM\tRead_Counts`).
pub fn write_gene_tpm_gz(path: &Path, idx: &MmapIndex, vbem: &VbemResult) -> Result<()> {
    let header = b"Gene\tTpM\tRead_Counts\n";
    let genes = idx.genes();
    let chunk_size = 4096usize;
    let blocks: Vec<Vec<u8>> = (0..genes.len())
        .collect::<Vec<_>>()
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut buf = String::with_capacity(chunk.len() * 48);
            for &g in chunk {
                let gname = idx.get_str(genes[g].gene_id);
                let _ = writeln!(buf, "{gname}\t{:.4}\t{:.2}", vbem.gene_tpm[g], vbem.gene_counts[g]);
            }
            buf.into_bytes()
        })
        .collect();
    write_parallel_gz(path, header, blocks)
}

/// Write `.isoform.tpm.gz` (`Isoform\tTpM\tRead_Counts`).
pub fn write_isoform_tpm_gz(path: &Path, idx: &MmapIndex, vbem: &VbemResult) -> Result<()> {
    let header = b"Isoform\tTpM\tRead_Counts\n";
    let txs = idx.transcripts();
    let chunk_size = 8192usize;
    let blocks: Vec<Vec<u8>> = (0..txs.len())
        .collect::<Vec<_>>()
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut buf = String::with_capacity(chunk.len() * 48);
            for &t in chunk {
                let tname = idx.get_str(txs[t].tx_id);
                let _ = writeln!(buf, "{tname}\t{:.4}\t{:.2}", vbem.tx_tpm[t], vbem.tx_counts[t]);
            }
            buf.into_bytes()
        })
        .collect();
    write_parallel_gz(path, header, blocks)
}

/// Write `.jnc.gz` with Whippet2 core columns + `Distinct_Pos` and `Max_Stack` PCR-stack diversity columns.
pub fn write_jnc_gz(
    path: &Path,
    idx: &MmapIndex,
    final_edges: &FxHashMap<(u32, u16, u16), (f64, u32, u32)>,
) -> Result<()> {
    let header = b"CHROM\tDONOR_COORD\tACCEPTOR_COORD\tGENEID:NODES:TYPE\tREADCOUNT\tSTRAND\tDistinct_Pos\tMax_Stack\n";

    let mut rows: Vec<(u32, u16, u16, f64, u32, u32)> = final_edges
        .iter()
        .filter(|(_, &(cnt, _, _))| cnt > 0.0)
        .map(|(&(gid, u, v), &(cnt, dpos, mstack))| (gid, u, v, cnt, dpos, mstack))
        .collect();
    rows.sort_unstable_by_key(|&(gid, u, v, _, _, _)| (gid, u, v));

    let chunk_size = 4096usize;
    let blocks: Vec<Vec<u8>> = rows
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut buf = String::with_capacity(chunk.len() * 80);
            for &(gid, u, v, cnt, dpos, mstack) in chunk {
                let gene = &idx.genes()[gid as usize];
                let nodes = idx.gene_nodes(gid);
                if (u as usize) >= nodes.len() || (v as usize) >= nodes.len() {
                    continue;
                }
                let u_nd = &nodes[u as usize];
                let v_nd = &nodes[v as usize];
                if !u_nd.right_edge().is_donor() || !v_nd.left_edge().is_acceptor() {
                    continue;
                }
                // Skip contiguous non-intron node transitions (e.g. RI exon-intron boundary)
                let contiguous = if gene.is_pos_strand() {
                    u_nd.genomic_coord + u_nd.len == v_nd.genomic_coord
                } else {
                    v_nd.genomic_coord + v_nd.len == u_nd.genomic_coord
                };
                if contiguous && u < v {
                    continue;
                }

                let (donor, acceptor) = if gene.is_pos_strand() {
                    (u_nd.genomic_coord + u_nd.len - 1, v_nd.genomic_coord)
                } else {
                    (v_nd.genomic_coord + v_nd.len - 1, u_nd.genomic_coord)
                };

                let anno_jncs = idx.gene_anno_junctions(gid);
                let jnc_type = if u >= v {
                    "CIRC_BS"
                } else if let Ok(jpos) = anno_jncs.binary_search_by_key(&(u, v), |j| (j.src_node, j.dst_node)) {
                    if (anno_jncs[jpos].tx_support as u32) >= gene.num_tx && gene.num_tx > 0 {
                        "CON_ANNO"
                    } else {
                        "ALT_ANNO"
                    }
                } else {
                    "ALT_UNIQ"
                };

                let chrom = idx.get_str(idx.chroms()[gene.chrom_id as usize].name);
                let gname = idx.get_str(gene.gene_id);
                let strand = if gene.is_pos_strand() { '+' } else { '-' };
                let _ = writeln!(
                    buf,
                    "{chrom}\t{donor}\t{acceptor}\t{gname}:{}-{}:{jnc_type}\t{cnt:.2}\t{strand}\t{dpos}\t{mstack}",
                    (u as u32) + 1,
                    (v as u32) + 1
                );
            }
            buf.into_bytes()
        })
        .collect();

    write_parallel_gz(path, header, blocks)
}

/// Write `.map.gz` summary mapping statistics (Whippet2 compatible + Langur metrics).
#[allow(clippy::too_many_arguments)]
pub fn write_map_gz(
    path: &Path,
    index_name: &str,
    stats: &SampleMapStats,
    fld: Option<&FragLenDist>,
    vbem: &VbemResult,
    final_nodes: &FxHashMap<(u32, u16), f64>,
    final_edges: &FxHashMap<(u32, u16, u16), (f64, u32, u32)>,
    idx: &MmapIndex,
) -> Result<()> {
    let total = stats.total_reads.max(1);
    let mapped = stats.mapped_reads;
    let mapped_pct = (mapped as f64 / total as f64) * 100.0;
    let multi_pct = if mapped > 0 {
        (stats.multimap_reads as f64 / mapped as f64) * 100.0
    } else {
        0.0
    };
    let decoy_pct = (stats.decoy_filtered_reads as f64 / total as f64) * 100.0;

    let body_reads: f64 = final_nodes.values().sum();
    let mut junc_reads = 0.0f64;
    let mut anno_junc_reads = 0.0f64;
    let mut junc_num = 0usize;
    let mut anno_junc_num = 0usize;

    for (&(gid, u, v), &(cnt, _, _)) in final_edges {
        if cnt <= 0.0 {
            continue;
        }
        junc_reads += cnt;
        junc_num += 1;
        let anno_jncs = idx.gene_anno_junctions(gid);
        if anno_jncs
            .binary_search_by_key(&(u, v), |j| (j.src_node, j.dst_node))
            .is_ok()
        {
            anno_junc_reads += cnt;
            anno_junc_num += 1;
        }
    }

    let mut s = String::with_capacity(1024);
    let _ = writeln!(s, "# -------- Langur {} ---------", env!("CARGO_PKG_VERSION"));
    let _ = writeln!(s, "Mapped_Index\t{index_name}");
    let _ = writeln!(s, "Read_Length\t{:.1}", stats.mean_read_len());
    let _ = writeln!(s, "Total_Reads\t{}", stats.total_reads);
    if let Some(f) = fld {
        let _ = writeln!(s, "Mean_Frag_Len\t{:.2}", f.mean);
        let _ = writeln!(s, "SD_Frag_Len\t{:.2}", f.sd);
    }
    let _ = writeln!(s, "# -------- Mapping Stats ---------");
    let _ = writeln!(s, "Mapped_Reads\t{mapped}");
    let _ = writeln!(s, "Mapped_Percent\t{mapped_pct:.2}%");
    let _ = writeln!(s, "Decoy_Filtered_Reads\t{}", stats.decoy_filtered_reads);
    let _ = writeln!(s, "Decoy_Filtered_Percent\t{decoy_pct:.2}%");
    let _ = writeln!(s, "# -------- Multi-loci Mapping ---------");
    let _ = writeln!(s, "Multimap_Reads\t{}", stats.multimap_reads);
    let _ = writeln!(s, "Multimap_Percent\t{multi_pct:.2}%");
    let _ = writeln!(s, "VBEM_Components\t{}", vbem.num_components);
    let _ = writeln!(s, "# -------- Unspliced / Spliced Reads ---------");
    let tot_bj = (body_reads + junc_reads).max(1.0);
    let _ = writeln!(s, "Exon_Body_Reads\t{:.0}", body_reads);
    let _ = writeln!(s, "Exon_Body_Percent\t{:.2}%", (body_reads / tot_bj) * 100.0);
    let _ = writeln!(s, "Junction_Reads\t{:.0}", junc_reads);
    let _ = writeln!(s, "Junction_Percent\t{:.2}%", (junc_reads / tot_bj) * 100.0);
    let _ = writeln!(s, "Junction_Number\t{junc_num}");
    let _ = writeln!(s, "Annotated_Junc_Reads\t{:.0}", anno_junc_reads);
    let _ = writeln!(s, "Annotated_Junc_Number\t{anno_junc_num}");
    let _ = writeln!(s, "Novel_Junc_Number\t{}", junc_num.saturating_sub(anno_junc_num));

    let file = File::create(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut enc = GzEncoder::new(BufWriter::new(file), Compression::fast());
    enc.write_all(s.as_bytes())?;
    enc.finish()?;
    Ok(())
}

/// Write SAM `@HD` and `@SQ` header lines to `w` (used when `--sam` is enabled).
pub fn write_sam_header<W: Write>(w: &mut W, idx: &MmapIndex) -> Result<()> {
    writeln!(w, "@HD\tVN:1.6\tSO:unsorted")?;
    for chr in idx.chroms() {
        let name = idx.get_str(chr.name);
        writeln!(w, "@SQ\tSN:{name}\tLN:{}", chr.length)?;
    }
    writeln!(w, "@PG\tID:langur\tPN:langur\tVN:{}", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

/// Format a single `ReadAlignment` as a valid SAMv1 record line.
/// Sets FLAG `0x10` based on true genomic strand (`pos_strand != aln.read_is_sense`) and
/// emits reverse-complemented `SEQ` + reversed `QUAL` whenever aligned to the genomic `-` strand.
pub fn format_sam_line(
    idx: &MmapIndex,
    qname: &[u8],
    seq: &[u8],
    qual: &[u8],
    aln: &ReadAlignment,
    is_paired: bool,
    is_first_mate: bool,
    is_secondary: bool,
) -> String {
    let gene = &idx.genes()[aln.gene_id as usize];
    let nodes = idx.gene_nodes(aln.gene_id);
    let chrom = idx.get_str(idx.chroms()[gene.chrom_id as usize].name);
    let pos_strand = gene.is_pos_strand();
    // Genomic `-` strand iff gene strand and read-vs-gene sense disagree
    let genomic_minus = pos_strand != aln.read_is_sense;

    let mut flag: u16 = 0;
    if is_paired {
        flag |= 0x1 | 0x2;
        flag |= if is_first_mate { 0x40 } else { 0x80 };
    }
    if genomic_minus {
        flag |= 0x10;
    }
    if is_secondary {
        flag |= 0x100;
    }

    let (pos, cigar) = build_sam_pos_and_cigar(aln, nodes, pos_strand);
    let qname_s = std::str::from_utf8(qname)
        .unwrap_or("read")
        .split_whitespace()
        .next()
        .unwrap_or("read");

    let rc_seq_buf;
    let seq_s = if genomic_minus {
        rc_seq_buf = crate::seq::revcomp(seq);
        std::str::from_utf8(&rc_seq_buf).unwrap_or("*")
    } else {
        std::str::from_utf8(seq).unwrap_or("*")
    };

    let rev_qual_buf;
    let qual_s = if qual.len() == seq.len() && !qual.is_empty() {
        if genomic_minus {
            rev_qual_buf = qual.iter().rev().copied().collect::<Vec<u8>>();
            std::str::from_utf8(&rev_qual_buf).unwrap_or("*")
        } else {
            std::str::from_utf8(qual).unwrap_or("*")
        }
    } else {
        "*"
    };

    let mapq = if is_secondary { 0 } else { 60 };
    let xs = if pos_strand { '+' } else { '-' };

    format!(
        "{qname_s}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq_s}\t{qual_s}\tAS:i:{}\tXS:A:{xs}\n",
        aln.score
    )
}

#[inline]
fn push_cigar_op(ops: &mut Vec<(u32, char)>, len: u32, ch: char) {
    if len == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        if last.1 == ch {
            last.0 += len;
            return;
        }
    }
    ops.push((len, ch));
}

fn build_sam_pos_and_cigar(
    aln: &ReadAlignment,
    nodes: &[crate::index::NodeRecord],
    pos_strand: bool,
) -> (u32, String) {
    // Circular back-splice alignments cannot be represented as a single monotonic linear SAM CIGAR:
    // emit the primary (longest) linear segment and soft-clip the back-spliced overhang cleanly.
    if aln.has_circ {
        let mut best_idx = 0usize;
        for i in 1..aln.segs.len() {
            if aln.segs[i].read_len > aln.segs[best_idx].read_len {
                best_idx = i;
            }
        }
        let s = aln.segs[best_idx];
        let left_clip: u32 = aln.soft_left
            + aln.segs[..best_idx].iter().map(|x| x.read_len).sum::<u32>();
        let right_clip: u32 = aln.soft_right
            + aln.segs[best_idx + 1..].iter().map(|x| x.read_len).sum::<u32>();
        let mut ops: Vec<(u32, char)> = Vec::new();
        push_cigar_op(&mut ops, left_clip, 'S');
        emit_seg_cigar(&mut ops, s);
        push_cigar_op(&mut ops, right_clip, 'S');
        if !pos_strand {
            ops.reverse();
        }
        let mut cigar = String::new();
        for (len, ch) in ops {
            let _ = write!(cigar, "{len}{ch}");
        }
        let nd = &nodes[s.node_id as usize];
        let pos = if pos_strand {
            nd.genomic_coord + s.node_offset
        } else {
            nd.genomic_coord + nd.len.saturating_sub(s.node_offset + s.ref_len)
        };
        return (pos.max(1), cigar);
    }

    let mut ops: Vec<(u32, char)> = Vec::new();
    push_cigar_op(&mut ops, aln.soft_left, 'S');

    for i in 0..aln.segs.len() {
        let s = aln.segs[i];
        emit_seg_cigar(&mut ops, s);
        if i + 1 < aln.segs.len() {
            let n_cur = &nodes[s.node_id as usize];
            let n_nxt = &nodes[aln.segs[i + 1].node_id as usize];
            let intron: i64 = if pos_strand {
                (n_nxt.genomic_coord as i64) - ((n_cur.genomic_coord + n_cur.len) as i64)
            } else {
                (n_cur.genomic_coord as i64) - ((n_nxt.genomic_coord + n_nxt.len) as i64)
            };
            if intron > 0 {
                push_cigar_op(&mut ops, intron as u32, 'N');
            }
        }
    }
    push_cigar_op(&mut ops, aln.soft_right, 'S');
    if !pos_strand {
        ops.reverse();
    }

    let mut cigar = String::new();
    for (len, ch) in ops {
        let _ = write!(cigar, "{len}{ch}");
    }

    let leftmost_pos = if pos_strand {
        let first = aln.segs[0];
        nodes[first.node_id as usize].genomic_coord + first.node_offset
    } else {
        let last = aln.segs[aln.segs.len() - 1];
        let nd = &nodes[last.node_id as usize];
        let right_end_in_node = last.node_offset + last.ref_len;
        nd.genomic_coord + nd.len.saturating_sub(right_end_in_node)
    };

    (leftmost_pos.max(1), cigar)
}

/// Emit CIGAR operations for one aligned node segment preserving exact `read_len` and `ref_len`
/// consumption (emitting `I` when `read_len > ref_len`, `D` when `ref_len > read_len`).
#[inline]
fn emit_seg_cigar(ops: &mut Vec<(u32, char)>, s: crate::align::AlignedNodeSeg) {
    if s.read_len == s.ref_len {
        push_cigar_op(ops, s.read_len, 'M');
    } else if s.read_len > s.ref_len {
        push_cigar_op(ops, s.ref_len, 'M');
        push_cigar_op(ops, s.read_len - s.ref_len, 'I');
    } else {
        push_cigar_op(ops, s.read_len, 'M');
        push_cigar_op(ops, s.ref_len - s.read_len, 'D');
    }
}
