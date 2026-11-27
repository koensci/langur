//! Lockstep streaming `.psi` / `.psi.gz` cohort loader supporting both 14-column
//! Whippet2 and 15-column Langur (`Inf_Var`) outputs, with safe `Inf_Var`
//! effective-depth scaling (`n_eff = max(1.0, raw_reads * clamp(v_nom / v_j, 0.1, 1.0))`).

use crate::delta::design::DesignSpec;
use crate::error::{Result, LangurError};
use crate::events::psi::beta_posterior_var;
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Default canonical Whippet event types included in `langur delta`.
pub const DEFAULT_WHITELIST_EVENTS: [&str; 8] = ["TS", "TE", "AF", "AL", "RI", "CE", "AA", "AD"];

#[derive(Clone, Debug)]
pub struct CohortReadOptions {
    /// Minimum raw `Total_Reads` required for an event in a sample to be included (default `5`).
    pub min_reads: u32,
    /// Minimum valid samples (per group in 2-group/contrast mode, or total in general GLM mode; default `1`).
    pub min_samples: usize,
    /// Include circular back-splicing (`BS`) events in differential testing (default `false`).
    pub include_circ: bool,
    /// Disable `Inf_Var` effective-depth scaling even when column 15 `Inf_Var` is present (default `false`).
    pub ignore_inf_var: bool,
}

impl Default for CohortReadOptions {
    fn default() -> Self {
        Self {
            min_reads: 5,
            min_samples: 1,
            include_circ: false,
            ignore_inf_var: false,
        }
    }
}

/// One testable splicing event assembled across all valid samples in the cohort.
#[derive(Clone, Debug)]
pub struct GlmEventData {
    pub gene: String,
    pub node: String,
    pub coord: String,
    pub strand: String,
    pub etype: String,
    pub complexity: usize,
    pub entropy: f64,
    /// 0-based sample indices into `DesignSpec::files` / `DesignSpec::x`.
    pub valid_idx: Vec<usize>,
    /// Observed PSI values in `[0.0, 1.0]` for `valid_idx`.
    pub psi: Vec<f64>,
    /// Effective read depth (`n_eff`) used in Beta-Binomial likelihood (`<= raw_reads`, `>= 1.0`).
    pub eff_reads: Vec<f64>,
    pub samp_a: usize,
    pub samp_b: usize,
    /// Mean raw `Total_Reads` in group A (or overall in general GLM mode).
    pub mean_reads_a: f64,
    /// Mean raw `Total_Reads` in group B (0.0 in general GLM mode).
    pub mean_reads_b: f64,
}

struct PsiStream {
    reader: Box<dyn BufRead>,
    inf_var_col: Option<usize>,
    line_buf: String,
}

impl PsiStream {
    fn open(path: &Path) -> Result<Self> {
        let mut f = File::open(path).map_err(|source| LangurError::IoPath {
            path: path.to_path_buf(),
            source,
        })?;
        let mut magic = [0u8; 2];
        let n = f.read(&mut magic)?;
        drop(f);

        let f = File::open(path).map_err(|source| LangurError::IoPath {
            path: path.to_path_buf(),
            source,
        })?;
        let mut reader: Box<dyn BufRead> = if n == 2 && magic == [0x1F, 0x8B] {
            Box::new(BufReader::with_capacity(
                128 * 1024,
                MultiGzDecoder::new(f),
            ))
        } else {
            Box::new(BufReader::with_capacity(128 * 1024, f))
        };

        // Read header line and locate optional `Inf_Var` column by header name
        let mut hdr = String::new();
        reader.read_line(&mut hdr)?;
        let inf_var_col = hdr
            .trim_end_matches(['\r', '\n'])
            .split('\t')
            .position(|c| c == "Inf_Var");

        Ok(Self {
            reader,
            inf_var_col,
            line_buf: String::with_capacity(512),
        })
    }

    /// Read next non-empty line into `self.line_buf`. Returns `Ok(false)` on EOF.
    fn next_line(&mut self) -> Result<bool> {
        loop {
            self.line_buf.clear();
            let n = self.reader.read_line(&mut self.line_buf)?;
            if n == 0 {
                return Ok(false);
            }
            if !self.line_buf.trim().is_empty() {
                return Ok(true);
            }
        }
    }
}

/// Compute `Inf_Var`-adjusted effective depth from raw `Total_Reads` and observed `Inf_Var`.
///
/// - Only shrinks depth when inferential variance `v_j` exceeds the nominal Beta variance `v_nom`
///   (e.g. complex `K_n >= 2` multi-path events or future Gibbs/bootstrap variance).
/// - Guards against `0.0` / underflowed `v_j` (`v_j > 1e-15 && v_nom > 1e-15`), clamps ratio to `[0.1, 1.0]`,
///   and floors `n_eff >= 1.0`.
#[inline]
pub fn scale_effective_reads(
    psi: f64,
    raw_reads: f64,
    inf_var: Option<f64>,
    ignore_inf_var: bool,
) -> f64 {
    if ignore_inf_var {
        return raw_reads;
    }
    let Some(v_j) = inf_var else {
        return raw_reads;
    };
    let v_nom = beta_posterior_var(psi, raw_reads);
    if v_j.is_finite() && v_j > 1e-15 && v_nom.is_finite() && v_nom > 1e-15 {
        let ratio = (v_nom / v_j).clamp(0.1, 1.0);
        (raw_reads * ratio).max(1.0)
    } else {
        raw_reads
    }
}

/// Parsed single-sample observation from a `.psi`/`.psi.gz` row that passed `--min-reads`.
struct ParsedSampleRow {
    psi: f64,
    raw_reads: f64,
    eff_reads: f64,
    complexity: usize,
    entropy: f64,
}

#[inline]
fn parse_psi_row(
    fields: &[&str],
    inf_var_col: Option<usize>,
    min_r: f64,
    ignore_inf_var: bool,
) -> Option<ParsedSampleRow> {
    let psi_str = fields[5];
    let num_str = fields[8];
    if psi_str == "NA" || num_str == "NA" {
        return None;
    }
    let psi_v = psi_str.parse::<f64>().ok()?;
    let num_v = num_str.parse::<f64>().ok()?;
    // Apply --min-reads strictly on RAW Total_Reads before Inf_Var scaling!
    if !psi_v.is_finite() || !num_v.is_finite() || psi_v < 0.0 || num_v < min_r {
        return None;
    }
    let psi_c = psi_v.clamp(0.0, 1.0);
    let inf_var = inf_var_col
        .and_then(|ci| fields.get(ci).copied())
        .filter(|&v| v != "NA")
        .and_then(|v| v.parse::<f64>().ok());
    let n_eff = scale_effective_reads(psi_c, num_v, inf_var, ignore_inf_var);

    let complexity = fields[9]
        .strip_prefix('K')
        .and_then(|rest| rest.parse::<usize>().ok())
        .unwrap_or(0);
    let entropy = if fields[10] != "NA" {
        fields[10]
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v > 0.0)
            .unwrap_or(0.0)
    } else {
        0.0
    };

    Some(ParsedSampleRow {
        psi: psi_c,
        raw_reads: num_v,
        eff_reads: n_eff,
        complexity,
        entropy,
    })
}

/// Sparse accumulator for sample-specific circular back-splicing (`BS`) events.
struct CircEventAccum {
    max_comp: usize,
    max_ent: f64,
    /// `(sample_idx, psi, raw_reads, eff_reads)` for valid samples passing `--min-reads`.
    samples: Vec<(usize, f64, f64, f64)>,
}

type CircKey = (String, String, String, String);

fn siphon_bs_row(
    circ_map: &mut FxHashMap<CircKey, CircEventAccum>,
    sample_idx: usize,
    fields: &[&str],
    inf_var_col: Option<usize>,
    min_r: f64,
    ignore_inf_var: bool,
) {
    let Some(obs) = parse_psi_row(fields, inf_var_col, min_r, ignore_inf_var) else {
        return;
    };
    let key: CircKey = (
        fields[0].to_owned(),
        fields[1].to_owned(),
        fields[2].to_owned(),
        fields[3].to_owned(),
    );
    let entry = circ_map.entry(key).or_insert_with(|| CircEventAccum {
        max_comp: 0,
        max_ent: 0.0,
        samples: Vec::new(),
    });
    if obs.complexity > entry.max_comp {
        entry.max_comp = obs.complexity;
    }
    if obs.entropy > entry.max_ent {
        entry.max_ent = obs.entropy;
    }
    if let Some(existing) = entry.samples.iter_mut().find(|x| x.0 == sample_idx) {
        *existing = (sample_idx, obs.psi, obs.raw_reads, obs.eff_reads);
    } else {
        entry
            .samples
            .push((sample_idx, obs.psi, obs.raw_reads, obs.eff_reads));
    }
}

#[allow(clippy::too_many_arguments)]
fn finalize_event_from_samples(
    design: &DesignSpec,
    opts: &CohortReadOptions,
    gene: String,
    node: String,
    coord: String,
    strand: String,
    etype: String,
    max_comp: usize,
    max_ent: f64,
    samples: &[(usize, f64, f64, f64)],
) -> Option<GlmEventData> {
    if design.is_two_group {
        let mut na = 0usize;
        let mut nb = 0usize;
        let mut sum_ra = 0.0f64;
        let mut sum_rb = 0.0f64;
        for &(s, _, raw_r, _) in samples {
            if design.group_a_mask[s] {
                na += 1;
                sum_ra += raw_r;
            } else if design.group_b_mask[s] {
                nb += 1;
                sum_rb += raw_r;
            }
        }
        if na < opts.min_samples || nb < opts.min_samples {
            return None;
        }
        let mut vidx = Vec::with_capacity(samples.len());
        let mut vpsi = Vec::with_capacity(samples.len());
        let mut veff = Vec::with_capacity(samples.len());
        for &(s, psi, _, eff_r) in samples {
            vidx.push(s);
            vpsi.push(psi);
            veff.push(eff_r);
        }
        Some(GlmEventData {
            gene,
            node,
            coord,
            strand,
            etype,
            complexity: max_comp,
            entropy: max_ent,
            valid_idx: vidx,
            psi: vpsi,
            eff_reads: veff,
            samp_a: na,
            samp_b: nb,
            mean_reads_a: sum_ra / (na as f64),
            mean_reads_b: sum_rb / (nb as f64),
        })
    } else {
        let nv = samples.len();
        if nv < opts.min_samples.max(2) {
            return None;
        }
        let mut sum_r = 0.0f64;
        let mut vidx = Vec::with_capacity(nv);
        let mut vpsi = Vec::with_capacity(nv);
        let mut veff = Vec::with_capacity(nv);
        for &(s, psi, raw_r, eff_r) in samples {
            sum_r += raw_r;
            vidx.push(s);
            vpsi.push(psi);
            veff.push(eff_r);
        }
        Some(GlmEventData {
            gene,
            node,
            coord,
            strand,
            etype,
            complexity: max_comp,
            entropy: max_ent,
            valid_idx: vidx,
            psi: vpsi,
            eff_reads: veff,
            samp_a: nv,
            samp_b: 0,
            mean_reads_a: sum_r / (nv as f64),
            mean_reads_b: 0.0,
        })
    }
}

/// Read all `.psi`/`.psi.gz` files in lockstep, filtering to whitelisted canonical events passing
/// `--min-reads` (on raw `Total_Reads`) and `--min-samples`.
///
/// Sample-specific circular back-splicing (`BS`) rows are siphoned into a side map during streaming
/// (when `opts.include_circ` is true) so lockstep alignment of the 8 canonical event types is never
/// disturbed, then appended in deterministic sorted order after EOF. Truncated files with missing
/// canonical rows raise an explicit error naming the truncated/mismatched sample.
pub fn read_psi_cohort(design: &DesignSpec, opts: &CohortReadOptions) -> Result<Vec<GlmEventData>> {
    let ns = design.files.len();
    if ns == 0 {
        return Ok(Vec::new());
    }
    let mut streams: Vec<PsiStream> = design
        .files
        .iter()
        .map(|p| PsiStream::open(p))
        .collect::<Result<Vec<_>>>()?;

    let canonical_whitelist: FxHashSet<&'static str> =
        DEFAULT_WHITELIST_EVENTS.into_iter().collect();

    let min_r = opts.min_reads as f64;
    let mut events = Vec::new();
    let mut circ_map: FxHashMap<CircKey, CircEventAccum> = FxHashMap::default();
    let mut row_samples: Vec<(usize, f64, f64, f64)> = Vec::with_capacity(ns);

    loop {
        let mut ref_gene = String::new();
        let mut ref_node = String::new();
        let mut ref_coord = String::new();
        let mut ref_strand = String::new();
        let mut ref_type = String::new();
        let mut s0_eof = false;
        let mut max_comp: usize = 0;
        let mut max_ent: f64 = 0.0;
        row_samples.clear();

        for s in 0..ns {
            let inf_var_col = streams[s].inf_var_col;
            loop {
                if !streams[s].next_line()? {
                    if s == 0 {
                        s0_eof = true;
                        break;
                    } else {
                        return Err(LangurError::InvalidParam(format!(
                            "truncated PSI file '{}': reached EOF early while '{}' still has event ({ref_gene}, {ref_node}, {ref_type})",
                            design.sample_names[s], design.sample_names[0]
                        )));
                    }
                }
                let line = streams[s].line_buf.trim_end_matches(['\r', '\n']);
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 11 {
                    continue;
                }
                let etype = fields[4];
                if etype == "BS" {
                    if opts.include_circ {
                        siphon_bs_row(
                            &mut circ_map,
                            s,
                            &fields,
                            inf_var_col,
                            min_r,
                            opts.ignore_inf_var,
                        );
                    }
                    continue;
                }
                if !canonical_whitelist.contains(etype) {
                    continue;
                }

                if s == 0 {
                    ref_gene = fields[0].to_owned();
                    ref_node = fields[1].to_owned();
                    ref_coord = fields[2].to_owned();
                    ref_strand = fields[3].to_owned();
                    ref_type = etype.to_owned();
                } else if fields[0] != ref_gene || fields[1] != ref_node || fields[4] != ref_type {
                    return Err(LangurError::InvalidParam(format!(
                        "mismatched events across PSI files: expected ({ref_gene}, {ref_node}, {ref_type}) from {}, got ({}, {}, {}) in {}",
                        design.sample_names[0],
                        fields[0],
                        fields[1],
                        fields[4],
                        design.sample_names[s]
                    )));
                }

                if let Some(obs) = parse_psi_row(&fields, inf_var_col, min_r, opts.ignore_inf_var) {
                    row_samples.push((s, obs.psi, obs.raw_reads, obs.eff_reads));
                    if obs.complexity > max_comp {
                        max_comp = obs.complexity;
                    }
                    if obs.entropy > max_ent {
                        max_ent = obs.entropy;
                    }
                }
                break;
            }
            if s0_eof {
                break;
            }
        }

        if s0_eof {
            // Sample 0 hit EOF: drain any trailing BS rows from samples 1..ns and verify no
            // remaining canonical rows exist (which would indicate sample 0 was truncated).
            for s in 1..ns {
                let inf_var_col = streams[s].inf_var_col;
                while streams[s].next_line()? {
                    let line = streams[s].line_buf.trim_end_matches(['\r', '\n']);
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 11 {
                        continue;
                    }
                    let etype = fields[4];
                    if etype == "BS" {
                        if opts.include_circ {
                            siphon_bs_row(
                                &mut circ_map,
                                s,
                                &fields,
                                inf_var_col,
                                min_r,
                                opts.ignore_inf_var,
                            );
                        }
                        continue;
                    }
                    if canonical_whitelist.contains(etype) {
                        return Err(LangurError::InvalidParam(format!(
                            "truncated PSI file '{}': reached EOF early while '{}' still has event ({}, {}, {})",
                            design.sample_names[0],
                            design.sample_names[s],
                            fields[0],
                            fields[1],
                            fields[4]
                        )));
                    }
                }
            }
            break;
        }

        if let Some(ev) = finalize_event_from_samples(
            design,
            opts,
            ref_gene,
            ref_node,
            ref_coord,
            ref_strand,
            ref_type,
            max_comp,
            max_ent,
            &row_samples,
        ) {
            events.push(ev);
        }
    }

    if opts.include_circ && !circ_map.is_empty() {
        let mut circ_entries: Vec<(CircKey, CircEventAccum)> = circ_map.into_iter().collect();
        circ_entries.sort_by(|a, b| a.0.cmp(&b.0));
        for ((gene, node, coord, strand), mut acc) in circ_entries {
            acc.samples.sort_by_key(|x| x.0);
            if let Some(ev) = finalize_event_from_samples(
                design,
                opts,
                gene,
                node,
                coord,
                strand,
                "BS".to_owned(),
                acc.max_comp,
                acc.max_ent,
                &acc.samples,
            ) {
                events.push(ev);
            }
        }
    }

    Ok(events)
}
