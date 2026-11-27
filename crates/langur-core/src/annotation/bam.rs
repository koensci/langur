//! Optional BAM/BAI splice-site supplementation for `langur index`.
//! Matches Whippet2's `--bam`, `--bam-min-reads`, and `--bam-both-novel` semantics,
//! while filtering secondary/supplementary alignments, honoring STAR/HISAT2 `XS`
//! strand tags when present, and keeping `known` immutable during read counting so
//! single-read noise cannot chain past `one_known_required`.

use crate::annotation::RefSet;
use crate::error::{Result, LangurError};
use noodles::bam;
use noodles::core::{Position, Region};
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct BamSupplementOptions {
    /// Minimum read support for a novel splice site to be added (default 1).
    pub min_reads: usize,
    /// If `true`, allow junctions where neither donor nor acceptor is annotated (`--bam-both-novel`).
    /// Default `false` (requires at least one known splice site, matching Whippet2).
    pub allow_both_novel: bool,
    /// If `true`, ignore read strand when `XS` tag is absent (useful for unstranded BAMs
    /// without `XS` tags). When `XS` is present, `XS` always dictates junction strand.
    pub unstranded_fallback: bool,
}

impl Default for BamSupplementOptions {
    fn default() -> Self {
        Self {
            min_reads: 1,
            allow_both_novel: false,
            unstranded_fallback: false,
        }
    }
}

/// Extract junction strand from the `XS` tag if present (`+` -> `Some(true)`, `-` -> `Some(false)`).
fn extract_xs_strand(rec: &bam::Record) -> Option<bool> {
    let xs_tag = Tag::new(b'X', b'S');
    let data = rec.data();
    let val = data.get(&xs_tag)?.ok()?;
    match val {
        Value::Character(b'+') => Some(true),
        Value::Character(b'-') => Some(false),
        Value::String(s) => match AsRef::<[u8]>::as_ref(s) {
            b"+" => Some(true),
            b"-" => Some(false),
            _ => None,
        },
        _ => None,
    }
}

/// Supplement a `RefSet` in place with novel donor/acceptor splice sites observed in
/// an indexed coordinate-sorted BAM file (expects `<bam>.bai` beside the BAM).
/// Returns the total number of newly added splice sites.
pub fn supplement_refset_from_bam(
    refset: &mut RefSet,
    bam_path: &Path,
    opts: &BamSupplementOptions,
) -> Result<usize> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| {
            LangurError::Bam(format!(
                "failed to open indexed BAM {}: {e}",
                bam_path.display()
            ))
        })?;

    let header = reader
        .read_header()
        .map_err(|e| LangurError::Bam(format!("failed to read BAM header: {e}")))?;

    let bam_chroms: FxHashSet<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    let one_known_required = !opts.allow_both_novel;
    let mut total_novel_ss = 0usize;

    for gene in &mut refset.genes {
        if !bam_chroms.contains(&gene.chrom) {
            continue;
        }
        let (min_s, max_e) = gene.genomic_span();
        if min_s == 0 || max_e < min_s {
            continue;
        }
        let Ok(start_pos) = Position::try_from(min_s as usize) else {
            continue;
        };
        let Ok(end_pos) = Position::try_from(max_e as usize) else {
            continue;
        };
        let region = Region::new(gene.chrom.as_str(), start_pos..=end_pos);

        // Keep `known` strictly immutable during read accumulation so single-read noise
        // cannot bootstrap other unannotated splice sites before `min_reads` filtering.
        let known: FxHashSet<u32> = gene.acc.iter().chain(gene.don.iter()).copied().collect();
        let mut novel_don_counts: FxHashMap<u32, usize> = FxHashMap::default();
        let mut novel_acc_counts: FxHashMap<u32, usize> = FxHashMap::default();

        let query = reader
            .query(&header, &region)
            .map_err(|e| LangurError::Bam(format!("BAM query failed for {}: {e}", gene.gene_id)))?;

        for rec_res in query {
            let rec =
                rec_res.map_err(|e| LangurError::Bam(format!("BAM record error: {e}")))?;
            let flags = rec.flags();
            // Skip unmapped, secondary, and supplementary alignments
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }

            // Honor STAR/HISAT2 `XS` junction strand tag when present; fall back to read
            // strand only when `XS` is absent (or accept either if `unstranded_fallback`).
            let strand_ok = match extract_xs_strand(&rec) {
                Some(xs_pos) => xs_pos == gene.strand,
                None => {
                    opts.unstranded_fallback || (!flags.is_reverse_complemented() == gene.strand)
                }
            };
            if !strand_ok {
                continue;
            }

            let Some(Ok(aln_start)) = rec.alignment_start() else {
                continue;
            };
            let mut refpos = usize::from(aln_start) as u32;
            for op_res in rec.cigar().iter() {
                let op = op_res.map_err(|e| LangurError::Bam(format!("CIGAR error: {e}")))?;
                let len = op.len() as u32;
                match op.kind() {
                    CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                        refpos += len;
                    }
                    CigarKind::Deletion => {
                        refpos += len;
                    }
                    CigarKind::Skip => {
                        let donor = refpos.saturating_sub(1);
                        refpos += len;
                        let accep = refpos;
                        if one_known_required && !known.contains(&donor) && !known.contains(&accep) {
                            continue;
                        }
                        *novel_don_counts.entry(donor).or_insert(0) += 1;
                        *novel_acc_counts.entry(accep).or_insert(0) += 1;
                    }
                    _ => {}
                }
            }
        }

        let anno_don_set: FxHashSet<u32> = gene.don.iter().copied().collect();
        let anno_acc_set: FxHashSet<u32> = gene.acc.iter().copied().collect();

        let mut new_don: Vec<u32> = novel_don_counts
            .into_iter()
            .filter(|&(pos, cnt)| cnt >= opts.min_reads && pos >= min_s && pos <= max_e)
            .map(|(pos, _)| pos)
            .collect();
        let mut new_acc: Vec<u32> = novel_acc_counts
            .into_iter()
            .filter(|&(pos, cnt)| cnt >= opts.min_reads && pos >= min_s && pos <= max_e)
            .map(|(pos, _)| pos)
            .collect();

        new_don.sort_unstable();
        new_acc.sort_unstable();

        let prev_don_len = gene.don.len();
        let prev_acc_len = gene.acc.len();

        gene.don.extend_from_slice(&new_don);
        gene.don.sort_unstable();
        gene.don.dedup();
        gene.acc.extend_from_slice(&new_acc);
        gene.acc.sort_unstable();
        gene.acc.dedup();

        // Re-clean txst/txen against updated splice sites (Whippet2 refset.jl L274-275)
        let updated_acc: FxHashSet<u32> = gene.acc.iter().copied().collect();
        let updated_don: FxHashSet<u32> = gene.don.iter().copied().collect();
        gene.txst.retain(|x| !updated_acc.contains(x));
        gene.txen.retain(|x| !updated_don.contains(x));

        gene.novel_don = new_don
            .into_iter()
            .filter(|x| !anno_don_set.contains(x))
            .collect();
        gene.novel_acc = new_acc
            .into_iter()
            .filter(|x| !anno_acc_set.contains(x))
            .collect();

        total_novel_ss += (gene.don.len() - prev_don_len) + (gene.acc.len() - prev_acc_len);
    }

    Ok(total_novel_ss)
}
