//! End-to-end deterministic `langur quant` sample pipeline and `--batch`
//! cohort runner (loads `.lidx` once via `mmap` across all samples).

use crate::align::{
    align_paired_reads, align_single_read, AlignConfig, AlignScratch, PairedAlignOutcome,
    SingleAlignOutcome,
};
use crate::error::{Result, LangurError};
use crate::events::{quantify_all_events, EventConfig};
use crate::index::MmapIndex;
use crate::io::fastq::{FastqBatchReader, FastqChunk};
use crate::io::writers::{
    format_sam_line, write_gene_tpm_gz, write_isoform_tpm_gz, write_jnc_gz, write_map_gz,
    write_psi_gz, write_sam_header,
};
use crate::quant::bias::BiasModel;
use crate::quant::fld::single_end_effective_length;
use crate::quant::shard::ReadBiasObs;
use crate::quant::{
    route_multimappers_to_stage2, solve_cc_vbem, QuantShard, SampleCounts, VbemConfig,
};
use rayon::prelude::*;
use smallvec::SmallVec;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct QuantConfig {
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
    pub out_prefix: PathBuf,
    pub index_label: String,
    pub force_gz: bool,
    pub bias_correct: bool,
    pub emit_sam: bool,
    pub align: AlignConfig,
    pub vbem: VbemConfig,
    pub events: EventConfig,
    /// Number of reads (or read pairs) per sub-chunk processed in parallel (default `4096`).
    pub chunk_size: usize,
    /// Number of parallel chunks per deterministic merge wave (default `32`).
    pub chunks_per_wave: usize,
}

impl Default for QuantConfig {
    fn default() -> Self {
        Self {
            r1: PathBuf::new(),
            r2: None,
            out_prefix: PathBuf::from("./output"),
            index_label: String::from("graph.lidx"),
            force_gz: false,
            bias_correct: false,
            emit_sam: false,
            align: AlignConfig::default(),
            vbem: VbemConfig::default(),
            events: EventConfig::default(),
            chunk_size: 16384,
            chunks_per_wave: 32,
        }
    }
}

#[derive(Clone, Debug)]
pub struct QuantSummary {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub multimap_reads: u64,
    pub decoy_filtered_reads: u64,
    pub mean_read_len: f64,
    pub mean_frag_len: Option<f64>,
    pub num_events: usize,
    pub num_vbem_components: usize,
    pub vbem_pass1_iters: usize,
    pub vbem_pass2_iters: usize,
    pub align_secs: f64,
    pub vbem_secs: f64,
    pub events_secs: f64,
    pub write_secs: f64,
}

/// Run complete end-to-end short-read quantification for one sample against an open `MmapIndex`.
pub fn quantify_sample(idx: &MmapIndex, cfg: &QuantConfig) -> Result<QuantSummary> {
    let t_align0 = Instant::now();
    let is_paired = cfg.r2.is_some();
    let sample_counts;
    let mut locked_primer: Option<Vec<f64>> = None;
    let mut init_sc = SampleCounts::with_index_capacity(idx);

    let mut stdout_handle = if cfg.emit_sam {
        let mut h = std::io::BufWriter::new(std::io::stdout());
        write_sam_header(&mut h, idx)?;
        Some(h)
    } else {
        None
    };

    // 1. 3-stage overlapped pipeline:
    //    - Stage A: Concurrent flat-buffer R1/R2 `zlib-ng` decompression -> bounded `sync_channel(1)`
    //    - Stage B: `rayon` parallel alignment + deterministic binary-tree shard reduction (32 -> 1)
    //    - Stage C: Dedicated merger thread merges 1 reduced wave shard into `SampleCounts`
    if let Some(ref r2_path) = cfg.r2 {
        let r1_p = cfg.r1.clone();
        let r2_p = r2_path.clone();
        let force_gz = cfg.force_gz;
        let emit_sam = cfg.emit_sam;
        let chunk_size = cfg.chunk_size;
        let chunks_per_wave = cfg.chunks_per_wave;

        let (tx_wave, rx_wave) =
            std::sync::mpsc::sync_channel::<Result<Vec<(FastqChunk, FastqChunk)>>>(2);

        let producer = std::thread::spawn(move || {
            let (tx1, rx1) = std::sync::mpsc::sync_channel::<Result<FastqChunk>>(4);
            let (tx2, rx2) = std::sync::mpsc::sync_channel::<Result<FastqChunk>>(4);

            let t1 = std::thread::spawn(move || {
                let mut r1 = match FastqBatchReader::open(&r1_p, force_gz, emit_sam) {
                    Ok(r) => r,
                    Err(e) => {
                        let _ = tx1.send(Err(e));
                        return;
                    }
                };
                loop {
                    match r1.next_chunk(chunk_size) {
                        Ok(c) if c.is_empty() => break,
                        Ok(c) => {
                            if tx1.send(Ok(c)).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            let _ = tx1.send(Err(e));
                            break;
                        }
                    }
                }
            });

            let t2 = std::thread::spawn(move || {
                let mut r2 = match FastqBatchReader::open(&r2_p, force_gz, emit_sam) {
                    Ok(r) => r,
                    Err(e) => {
                        let _ = tx2.send(Err(e));
                        return;
                    }
                };
                loop {
                    match r2.next_chunk(chunk_size) {
                        Ok(c) if c.is_empty() => break,
                        Ok(c) => {
                            if tx2.send(Ok(c)).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            let _ = tx2.send(Err(e));
                            break;
                        }
                    }
                }
            });

            loop {
                let mut wave: Vec<(FastqChunk, FastqChunk)> = Vec::with_capacity(chunks_per_wave);
                let mut err_opt: Option<LangurError> = None;
                for _ in 0..chunks_per_wave {
                    match (rx1.recv(), rx2.recv()) {
                        (Ok(Ok(c1)), Ok(Ok(c2))) => {
                            if c1.len() != c2.len() {
                                err_opt = Some(LangurError::Fasta(
                                    "paired FASTQ files have different numbers of reads".into(),
                                ));
                                break;
                            }
                            wave.push((c1, c2));
                        }
                        (Err(_), Err(_)) => break,
                        (Ok(Err(e)), _) | (_, Ok(Err(e))) => {
                            err_opt = Some(e);
                            break;
                        }
                        _ => {
                            err_opt = Some(LangurError::Fasta(
                                "paired FASTQ files have different numbers of reads".into(),
                            ));
                            break;
                        }
                    }
                }
                if let Some(e) = err_opt {
                    let _ = tx_wave.send(Err(e));
                    break;
                }
                if wave.is_empty() {
                    break;
                }
                if tx_wave.send(Ok(wave)).is_err() {
                    break;
                }
            }
            drop(rx1);
            drop(rx2);
            let _ = t1.join();
            let _ = t2.join();
        });

        let (tx_merge, rx_merge) = std::sync::mpsc::sync_channel::<(QuantShard, Vec<String>)>(2);
        let merger = std::thread::spawn(
            move || -> Result<(SampleCounts, Option<std::io::BufWriter<std::io::Stdout>>)> {
                for (wave_shard, sam_lines) in rx_merge {
                    init_sc.merge_shard(wave_shard);
                    if let Some(ref mut out) = stdout_handle {
                        for line in sam_lines {
                            out.write_all(line.as_bytes())?;
                        }
                    }
                }
                Ok((init_sc, stdout_handle))
            },
        );

        for (wave_idx, wave_res) in rx_wave.into_iter().enumerate() {
            let wave = wave_res?;
            let pw_ref: Option<&[f64]> = locked_primer.as_deref();
            let wave_results: Vec<(QuantShard, Vec<String>)> = wave
                .into_par_iter()
                .map_init(
                    || AlignScratch::new(&cfg.align),
                    |scratch, (c1, c2)| process_paired_chunk(idx, &c1, &c2, cfg, pw_ref, scratch),
                )
                .collect();
            let mut reduced = reduce_wave_shards(wave_results);
            if wave_idx == 0 && cfg.bias_correct {
                let w0_model = reduced.0.bias.finalize_primer();
                reduced.0.collapse_wave0_hexamers(&w0_model.primer_weight);
                locked_primer = Some(w0_model.primer_weight);
            }
            if tx_merge.send(reduced).is_err() {
                break;
            }
        }
        drop(tx_merge);
        let _ = producer.join();
        let (merged_sc, out_h) = merger
            .join()
            .map_err(|_| LangurError::InvalidParam("merger thread panicked".into()))??;
        sample_counts = merged_sc;
        stdout_handle = out_h;
    } else {
        let r1_p = cfg.r1.clone();
        let force_gz = cfg.force_gz;
        let emit_sam = cfg.emit_sam;
        let chunk_size = cfg.chunk_size;
        let chunks_per_wave = cfg.chunks_per_wave;

        let (tx_wave, rx_wave) = std::sync::mpsc::sync_channel::<Result<Vec<FastqChunk>>>(2);
        let producer = std::thread::spawn(move || {
            let mut reader = match FastqBatchReader::open(&r1_p, force_gz, emit_sam) {
                Ok(r) => r,
                Err(e) => {
                    let _ = tx_wave.send(Err(e));
                    return;
                }
            };
            loop {
                let mut wave: Vec<FastqChunk> = Vec::with_capacity(chunks_per_wave);
                let mut err_opt: Option<LangurError> = None;
                for _ in 0..chunks_per_wave {
                    match reader.next_chunk(chunk_size) {
                        Ok(chunk) if chunk.is_empty() => break,
                        Ok(chunk) => wave.push(chunk),
                        Err(e) => {
                            err_opt = Some(e);
                            break;
                        }
                    }
                }
                if let Some(e) = err_opt {
                    let _ = tx_wave.send(Err(e));
                    break;
                }
                if wave.is_empty() {
                    break;
                }
                if tx_wave.send(Ok(wave)).is_err() {
                    break;
                }
            }
        });

        let (tx_merge, rx_merge) = std::sync::mpsc::sync_channel::<(QuantShard, Vec<String>)>(2);
        let merger = std::thread::spawn(
            move || -> Result<(SampleCounts, Option<std::io::BufWriter<std::io::Stdout>>)> {
                for (wave_shard, sam_lines) in rx_merge {
                    init_sc.merge_shard(wave_shard);
                    if let Some(ref mut out) = stdout_handle {
                        for line in sam_lines {
                            out.write_all(line.as_bytes())?;
                        }
                    }
                }
                Ok((init_sc, stdout_handle))
            },
        );

        for (wave_idx, wave_res) in rx_wave.into_iter().enumerate() {
            let wave = wave_res?;
            let pw_ref: Option<&[f64]> = locked_primer.as_deref();
            let wave_results: Vec<(QuantShard, Vec<String>)> = wave
                .into_par_iter()
                .map_init(
                    || AlignScratch::new(&cfg.align),
                    |scratch, chunk| process_single_chunk(idx, &chunk, cfg, pw_ref, scratch),
                )
                .collect();
            let mut reduced = reduce_wave_shards(wave_results);
            if wave_idx == 0 && cfg.bias_correct {
                let w0_model = reduced.0.bias.finalize_primer();
                reduced.0.collapse_wave0_hexamers(&w0_model.primer_weight);
                locked_primer = Some(w0_model.primer_weight);
            }
            if tx_merge.send(reduced).is_err() {
                break;
            }
        }
        drop(tx_merge);
        let _ = producer.join();
        let (merged_sc, out_h) = merger
            .join()
            .map_err(|_| LangurError::InvalidParam("merger thread panicked".into()))??;
        sample_counts = merged_sc;
        stdout_handle = out_h;
    }

    if let Some(mut out) = stdout_handle {
        out.flush()?;
    }
    let align_secs = t_align0.elapsed().as_secs_f64();
    let t_vbem0 = Instant::now();

    // 2. Finalize deterministic post-pass FragLenDist & base transcript effective lengths (Safeguard #2)
    let mean_read_len = sample_counts.stats.mean_read_len();
    let fld = if is_paired {
        Some(sample_counts.fld.finalize((mean_read_len * 2.0).max(150.0), 40.0))
    } else {
        None
    };

    let base_eff_lens: Vec<f64> = idx
        .transcripts()
        .iter()
        .map(|t| match fld {
            Some(ref f) => f.effective_length(t.seq_len),
            None => single_end_effective_length(t.seq_len, mean_read_len),
        })
        .collect();

    // 3. Finalize 5'-hexamer primer weights (reusing deterministic Wave 0 locked primer weights
    //    paired with full-sample observed GC counts `sample_counts.bias.gc_fore`)
    let mut bias_model = if cfg.bias_correct {
        if let Some(pw) = locked_primer {
            BiasModel::with_locked_primer(pw, sample_counts.bias.gc_fore)
        } else {
            sample_counts.bias.finalize_primer()
        }
    } else {
        BiasModel::default()
    };

    // 4. Stage 1 Pass 1: Connected-Component Transcript VBEM + SQUAREM (primer-weighted)
    let mut sorted_ec_keys: Vec<SmallVec<[u32; 4]>> =
        sample_counts.ec_counts.keys().cloned().collect();
    sorted_ec_keys.sort_unstable();

    let ecs_pass1: Vec<(SmallVec<[u32; 4]>, f64)> = sorted_ec_keys
        .iter()
        .map(|k| {
            let w = sample_counts.ec_counts[k].primer_adjusted(&bias_model, cfg.bias_correct);
            (k.clone(), w)
        })
        .collect();

    let vbem_pass1 = solve_cc_vbem(idx, &ecs_pass1, base_eff_lens.clone(), &cfg.vbem, None);
    let pass1_iters = vbem_pass1.total_iters;

    // 5. Stage 1 Pass 2 (when --biascorrect): fit 20-bin GC weights against tx_expected_gc
    //    and run warm-started VBEM using primer-only `ecs_pass1` + GC-adjusted `gc_eff_lens`
    //    (avoiding GC double-counting while keeping `tx_counts` on the observed-read scale!).
    let (vbem_final, pass2_iters) = if cfg.bias_correct {
        bias_model.fit_gc_from_tx_counts(idx, &vbem_pass1.tx_counts);

        let gc_eff_lens: Vec<f64> = base_eff_lens
            .iter()
            .enumerate()
            .map(|(t, &l)| (l * bias_model.transcript_gc_capture_rate(idx, t)).max(1.0))
            .collect();

        let res2 = solve_cc_vbem(
            idx,
            &ecs_pass1,
            gc_eff_lens,
            &cfg.vbem,
            Some(&vbem_pass1.tx_counts),
        );
        let it2 = res2.total_iters;
        (res2, it2)
    } else {
        (vbem_pass1, 0)
    };

    // 6. Stage 1 -> Stage 2 routing (converged transcript/gene posteriors -> per-gene nodes/edges)
    let (final_nodes, final_edges) = route_multimappers_to_stage2(
        idx,
        &sample_counts,
        &vbem_final,
        &bias_model,
        cfg.bias_correct,
    );
    let vbem_secs = t_vbem0.elapsed().as_secs_f64();

    // 7. Stage 2: Local splice-graph event PSI quantification (with PCR-stack soft-capping & Inf_Var)
    let t_ev0 = Instant::now();
    let mut ev_cfg = cfg.events.clone();
    ev_cfg.allow_circ = cfg.align.allow_circ;
    ev_cfg.mean_read_len = mean_read_len;
    ev_cfg.min_anno_overhang = cfg.align.min_anno_overhang;
    ev_cfg.min_overhang = cfg.align.min_novel_overhang;

    let events = quantify_all_events(idx, &final_nodes, &final_edges, &ev_cfg);
    let events_secs = t_ev0.elapsed().as_secs_f64();

    // 8. Write all 5 parallel-gzip output files
    let t_wr0 = Instant::now();
    let prefix = cfg.out_prefix.to_string_lossy();
    let psi_path = PathBuf::from(format!("{prefix}.psi.gz"));
    let g_tpm_path = PathBuf::from(format!("{prefix}.gene.tpm.gz"));
    let i_tpm_path = PathBuf::from(format!("{prefix}.isoform.tpm.gz"));
    let jnc_path = PathBuf::from(format!("{prefix}.jnc.gz"));
    let map_path = PathBuf::from(format!("{prefix}.map.gz"));

    write_psi_gz(&psi_path, idx, &events)?;
    write_gene_tpm_gz(&g_tpm_path, idx, &vbem_final)?;
    write_isoform_tpm_gz(&i_tpm_path, idx, &vbem_final)?;
    write_jnc_gz(&jnc_path, idx, &final_edges)?;
    write_map_gz(
        &map_path,
        &cfg.index_label,
        &sample_counts.stats,
        fld.as_ref(),
        &vbem_final,
        &final_nodes,
        &final_edges,
        idx,
    )?;
    let write_secs = t_wr0.elapsed().as_secs_f64();

    Ok(QuantSummary {
        total_reads: sample_counts.stats.total_reads,
        mapped_reads: sample_counts.stats.mapped_reads,
        multimap_reads: sample_counts.stats.multimap_reads,
        decoy_filtered_reads: sample_counts.stats.decoy_filtered_reads,
        mean_read_len,
        mean_frag_len: fld.map(|f| f.mean),
        num_events: events.len(),
        num_vbem_components: vbem_final.num_components,
        vbem_pass1_iters: pass1_iters,
        vbem_pass2_iters: pass2_iters,
        align_secs,
        vbem_secs,
        events_secs,
        write_secs,
    })
}

/// Deterministically reduce a wave's shards pairwise in parallel (`(0+1), (2+3), ...`) inside
/// `rayon` until a single merged wave shard remains, preserving fixed tree order across thread counts.
fn reduce_wave_shards(mut shards: Vec<(QuantShard, Vec<String>)>) -> (QuantShard, Vec<String>) {
    if shards.is_empty() {
        return (QuantShard::default(), Vec::new());
    }
    while shards.len() > 1 {
        let mut pairs = Vec::with_capacity(shards.len().div_ceil(2));
        let mut iter = shards.into_iter();
        while let Some(a) = iter.next() {
            pairs.push((a, iter.next()));
        }
        shards = pairs
            .into_par_iter()
            .map(|((mut a_sh, mut a_sam), b_opt)| {
                if let Some((b_sh, b_sam)) = b_opt {
                    a_sh.merge_from(b_sh);
                    if !b_sam.is_empty() {
                        a_sam.extend(b_sam);
                    }
                }
                (a_sh, a_sam)
            })
            .collect();
    }
    shards.pop().unwrap()
}

fn process_single_chunk(
    idx: &MmapIndex,
    chunk: &FastqChunk,
    cfg: &QuantConfig,
    locked_primer: Option<&[f64]>,
    scratch: &mut AlignScratch,
) -> (QuantShard, Vec<String>) {
    let mut shard = QuantShard::with_chunk_hint(locked_primer.is_some() || !cfg.bias_correct, chunk.len());
    let mut sam_lines = Vec::new();

    for i in 0..chunk.len() {
        let seq = chunk.seq(i);
        shard.stats.total_reads += 1;
        shard.stats.read_len_sum += seq.len() as u64;
        shard.stats.read_len_cnt += 1;

        match align_single_read(idx, seq, &cfg.align, scratch) {
            SingleAlignOutcome::Unmapped => {}
            SingleAlignOutcome::DecoyFiltered => {
                shard.stats.decoy_filtered_reads += 1;
            }
            SingleAlignOutcome::Mapped(alns) => {
                shard.stats.mapped_reads += 1;
                let hex_buf: SmallVec<[u16; 2]>;
                let bias_obs = if !cfg.bias_correct {
                    ReadBiasObs::None
                } else if let Some(pw) = locked_primer {
                    let (w_sum, n_hex, gc_bin) = shard.bias.observe_read_fast(seq, pw);
                    let primer_w = if n_hex > 0 {
                        Some(w_sum / (n_hex as f64))
                    } else {
                        None
                    };
                    ReadBiasObs::Fast { primer_w, gc_bin }
                } else {
                    let (h, gc_bin) = shard.bias.observe_read(seq);
                    hex_buf = h;
                    ReadBiasObs::Wave0 {
                        hexamers: &hex_buf,
                        gc_bin,
                    }
                };

                if alns.len() == 1 {
                    shard.record_single_unique(idx, &alns[0], bias_obs);
                } else {
                    shard.stats.multimap_reads += 1;
                    shard.record_single_multi(idx, &alns, bias_obs);
                }

                if cfg.emit_sam {
                    let id = chunk.id(i);
                    let qual = chunk.qual(i);
                    for (ai, a) in alns.iter().enumerate() {
                        sam_lines.push(format_sam_line(
                            idx,
                            id,
                            seq,
                            qual,
                            a,
                            false,
                            true,
                            ai > 0,
                        ));
                    }
                }
            }
        }
    }

    (shard, sam_lines)
}

fn process_paired_chunk(
    idx: &MmapIndex,
    c1: &FastqChunk,
    c2: &FastqChunk,
    cfg: &QuantConfig,
    locked_primer: Option<&[f64]>,
    scratch: &mut AlignScratch,
) -> (QuantShard, Vec<String>) {
    let mut shard = QuantShard::with_chunk_hint(locked_primer.is_some() || !cfg.bias_correct, c1.len());
    let mut sam_lines = Vec::new();

    for i in 0..c1.len() {
        let s1 = c1.seq(i);
        let s2 = c2.seq(i);
        shard.stats.total_reads += 1;
        shard.stats.read_len_sum += (s1.len() + s2.len()) as u64;
        shard.stats.read_len_cnt += 2;

        match align_paired_reads(idx, s1, s2, &cfg.align, scratch) {
            PairedAlignOutcome::Unmapped => {}
            PairedAlignOutcome::DecoyFiltered => {
                shard.stats.decoy_filtered_reads += 1;
            }
            PairedAlignOutcome::Mapped(pairs) => {
                shard.stats.mapped_reads += 1;
                let mut hex_buf: SmallVec<[u16; 4]> = SmallVec::new();
                let bias_obs = if !cfg.bias_correct {
                    ReadBiasObs::None
                } else if let Some(pw) = locked_primer {
                    let (w1, n1, g1) = shard.bias.observe_read_fast(s1, pw);
                    let (w2, n2, g2) = shard.bias.observe_read_fast(s2, pw);
                    let n_tot = n1 + n2;
                    let primer_w = if n_tot > 0 {
                        Some((w1 + w2) / (n_tot as f64))
                    } else {
                        None
                    };
                    let gc_bin = match (g1, g2) {
                        (Some(a), Some(b)) => Some((a + b) / 2),
                        (Some(a), None) => Some(a),
                        (None, Some(b)) => Some(b),
                        (None, None) => None,
                    };
                    ReadBiasObs::Fast { primer_w, gc_bin }
                } else {
                    let (h1, g1) = shard.bias.observe_read(s1);
                    let (h2, g2) = shard.bias.observe_read(s2);
                    hex_buf.extend_from_slice(&h1);
                    for x in h2 {
                        if hex_buf.len() < 4 {
                            hex_buf.push(x);
                        }
                    }
                    let gc_bin = match (g1, g2) {
                        (Some(a), Some(b)) => Some((a + b) / 2),
                        (Some(a), None) => Some(a),
                        (None, Some(b)) => Some(b),
                        (None, None) => None,
                    };
                    ReadBiasObs::Wave0 {
                        hexamers: &hex_buf,
                        gc_bin,
                    }
                };

                if pairs.len() == 1 {
                    shard.record_paired_unique(&pairs[0], bias_obs);
                } else {
                    shard.stats.multimap_reads += 1;
                    shard.record_paired_multi(&pairs, bias_obs);
                }

                if cfg.emit_sam {
                    let id1 = c1.id(i);
                    let q1 = c1.qual(i);
                    let id2 = c2.id(i);
                    let q2 = c2.qual(i);
                    for (pi, p) in pairs.iter().enumerate() {
                        let (a1, a2) = if p.r1_is_5p {
                            (&p.mate_5p, &p.mate_3p)
                        } else {
                            (&p.mate_3p, &p.mate_5p)
                        };
                        sam_lines.push(format_sam_line(
                            idx,
                            id1,
                            s1,
                            q1,
                            a1,
                            true,
                            true,
                            pi > 0,
                        ));
                        sam_lines.push(format_sam_line(
                            idx,
                            id2,
                            s2,
                            q2,
                            a2,
                            true,
                            false,
                            pi > 0,
                        ));
                    }
                }
            }
        }
    }

    (shard, sam_lines)
}

/// Quantify a multi-sample cohort from a CSV manifest (`sample,r1,r2,out_prefix` or `sample,r1,out_prefix`)
/// while holding a single zero-copy `MmapIndex` open across all samples.
pub fn quantify_batch_manifest(
    idx: &MmapIndex,
    manifest_path: &Path,
    base_cfg: &QuantConfig,
) -> Result<Vec<(String, QuantSummary)>> {
    let f = File::open(manifest_path).map_err(|source| LangurError::IoPath {
        path: manifest_path.to_path_buf(),
        source,
    })?;
    let reader = BufReader::new(f);
    let mut results = Vec::new();

    for (line_idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split(',').map(str::trim).collect();
        // Skip header row if first column is "sample" / "sample_id"
        if line_idx == 0 && matches!(cols[0].to_ascii_lowercase().as_str(), "sample" | "sample_id") {
            continue;
        }
        if cols.len() < 3 {
            return Err(LangurError::InvalidParam(format!(
                "manifest line {} needs at least 3 columns (sample,r1,out_prefix or sample,r1,r2,out_prefix)",
                line_idx + 1
            )));
        }

        let sample_name = cols[0].to_owned();
        let mut scfg = base_cfg.clone();
        scfg.r1 = PathBuf::from(cols[1]);
        if cols.len() >= 4 {
            let r2s = cols[2];
            scfg.r2 = if r2s.is_empty() || r2s == "-" || r2s.eq_ignore_ascii_case("NA") {
                None
            } else {
                Some(PathBuf::from(r2s))
            };
            scfg.out_prefix = PathBuf::from(cols[3]);
        } else {
            scfg.r2 = None;
            scfg.out_prefix = PathBuf::from(cols[2]);
        }

        let summary = quantify_sample(idx, &scfg)?;
        results.push((sample_name, summary));
    }

    Ok(results)
}
