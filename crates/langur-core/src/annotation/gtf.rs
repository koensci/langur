//! GTF2.2 parser preserving Whippet2 CSG gene/transcript semantics while
//! gracefully supporting unsorted/interleaved transcript exon lines and gzipped GTFs.

use crate::error::{Result, LangurError};
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct GtfOptions {
    /// If true, only transcripts carrying `tag "basic"` define `txst`/`txen`/`reftx`,
    /// while non-basic transcripts still contribute splice sites (`don`, `acc`).
    pub basic_or_ss_only: bool,
    /// If true (and `!basic_or_ss_only`), skip exons whose `transcript_support_level`
    /// is neither missing (`_`/`NA`) nor `1`/`2`.
    pub suppress_low_tsl: bool,
}

impl Default for GtfOptions {
    fn default() -> Self {
        Self {
            basic_or_ss_only: false,
            suppress_low_tsl: false,
        }
    }
}

/// Annotated reference transcript within a gene.
#[derive(Clone, Debug)]
pub struct RefTx {
    pub tx_id: String,
    pub tx_start: u32,
    pub tx_end: u32,
    pub is_basic: bool,
    pub tsl: u8, // 0 = NA/unknown, 1..=5
    /// Sorted 1-based exon end coordinates (genomic donor coordinates, incl. last exon end).
    pub don: Vec<u32>,
    /// Sorted 1-based exon start coordinates (genomic acceptor coordinates, incl. first exon start).
    pub acc: Vec<u32>,
    /// Total spliced exonic length in bp.
    pub length: u32,
}

/// Aggregated reference gene ready for CSG construction (matching Whippet2's `RefGene`).
#[derive(Clone, Debug)]
pub struct RefGene {
    pub gene_id: String,
    pub gene_name: String,
    pub chrom: String,
    /// `true` for `+`, `false` for `-`.
    pub strand: bool,
    /// Sorted unique internal donor coordinates (1-based genomic exon ends, within gene range).
    pub don: Vec<u32>,
    /// Sorted unique internal acceptor coordinates (1-based genomic exon starts, within gene range).
    pub acc: Vec<u32>,
    /// Sorted unique transcript start coordinates (excluding any coordinate in `acc`).
    pub txst: Vec<u32>,
    /// Sorted unique transcript end coordinates (excluding any coordinate in `don`).
    pub txen: Vec<u32>,
    /// Sorted unique exon intervals `(start, end)` (1-based inclusive).
    pub exons: Vec<(u32, u32)>,
    /// Sorted unique novel donor coordinates added from BAM supplementation.
    pub novel_don: Vec<u32>,
    /// Sorted unique novel acceptor coordinates added from BAM supplementation.
    pub novel_acc: Vec<u32>,
    /// Mean transcript length across annotated `reftx`.
    pub mean_tx_len: f32,
    /// Annotated transcripts belonging to this gene.
    pub reftx: Vec<RefTx>,
}

impl RefGene {
    /// 1-based inclusive genomic span `(min_txst, max_txen)`.
    #[inline]
    pub fn genomic_span(&self) -> (u32, u32) {
        let min_s = self
            .txst
            .first()
            .copied()
            .or_else(|| self.exons.first().map(|e| e.0))
            .unwrap_or(0);
        let max_e = self
            .txen
            .last()
            .copied()
            .or_else(|| self.exons.iter().map(|e| e.1).max())
            .unwrap_or(0);
        (min_s, max_e)
    }
}

/// Collection of genes loaded from GTF, preserving deterministic order.
#[derive(Clone, Debug, Default)]
pub struct RefSet {
    pub genes: Vec<RefGene>,
    pub gene_index: FxHashMap<String, usize>,
}

impl RefSet {
    pub fn get(&self, gene_id: &str) -> Option<&RefGene> {
        self.gene_index.get(gene_id).map(|&i| &self.genes[i])
    }

    pub fn get_mut(&mut self, gene_id: &str) -> Option<&mut RefGene> {
        self.gene_index.get(gene_id).copied().map(|i| &mut self.genes[i])
    }

    pub fn len(&self) -> usize {
        self.genes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }
}

struct TxAccum {
    tx_id: String,
    gene_id: String,
    gene_name: String,
    chrom: String,
    strand: bool,
    is_basic: bool,
    tsl: u8,
    exons: Vec<(u32, u32)>,
}

struct GeneAccum {
    gene_id: String,
    gene_name: String,
    chrom: String,
    strand: bool,
    has_primary_tx: bool,
    don: Vec<u32>,
    acc: Vec<u32>,
    txst: Vec<u32>,
    txen: Vec<u32>,
    exons: Vec<(u32, u32)>,
    tx_len_sum: u64,
    reftx: Vec<RefTx>,
}

/// Load a GTF2.2 file (plain `.gtf` or gzip-compressed `.gtf.gz`).
pub fn load_gtf_file(path: &Path, opts: &GtfOptions) -> Result<RefSet> {
    let mut file = File::open(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    drop(file);

    let file = File::open(path).map_err(|source| LangurError::IoPath {
        path: path.to_path_buf(),
        source,
    })?;
    if n == 2 && magic == [0x1F, 0x8B] {
        let dec = MultiGzDecoder::new(file);
        load_gtf_reader(BufReader::with_capacity(256 * 1024, dec), opts)
    } else {
        load_gtf_reader(BufReader::with_capacity(256 * 1024, file), opts)
    }
}

/// Parse GTF2.2 exon records from any buffered reader.
pub fn load_gtf_reader<R: BufRead>(mut reader: R, opts: &GtfOptions) -> Result<RefSet> {
    let mut tx_order: Vec<String> = Vec::new();
    let mut tx_map: FxHashMap<String, TxAccum> = FxHashMap::default();
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_seen: FxHashSet<String> = FxHashSet::default();

    let mut line_buf = String::with_capacity(512);
    let mut line_no: usize = 0;

    loop {
        line_buf.clear();
        let bytes = reader.read_line(&mut line_buf)?;
        if bytes == 0 {
            break;
        }
        line_no += 1;
        let line = line_buf.trim_end_matches(['\r', '\n']);
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut cols = line.split('\t');
        let chrom = cols.next().ok_or_else(|| LangurError::GtfParse {
            line: line_no,
            message: "missing chrom column".into(),
        })?;
        let _source = cols.next();
        let feature = cols.next().ok_or_else(|| LangurError::GtfParse {
            line: line_no,
            message: "missing feature column".into(),
        })?;
        if feature != "exon" {
            continue;
        }
        let start_s = cols.next().ok_or_else(|| LangurError::GtfParse {
            line: line_no,
            message: "missing start column".into(),
        })?;
        let end_s = cols.next().ok_or_else(|| LangurError::GtfParse {
            line: line_no,
            message: "missing end column".into(),
        })?;
        let _score = cols.next();
        let strand_s = cols.next().ok_or_else(|| LangurError::GtfParse {
            line: line_no,
            message: "missing strand column".into(),
        })?;
        let _frame = cols.next();
        let attrs_s = cols.next().ok_or_else(|| LangurError::GtfParse {
            line: line_no,
            message: "missing attributes column".into(),
        })?;

        let start: u32 = start_s.parse().map_err(|_| LangurError::GtfParse {
            line: line_no,
            message: format!("invalid start coordinate '{start_s}'"),
        })?;
        let end: u32 = end_s.parse().map_err(|_| LangurError::GtfParse {
            line: line_no,
            message: format!("invalid end coordinate '{end_s}'"),
        })?;
        if end < start {
            return Err(LangurError::GtfParse {
                line: line_no,
                message: format!("exon end ({end}) < start ({start})"),
            });
        }
        let strand = match strand_s.as_bytes().first() {
            Some(b'+') => true,
            Some(b'-') => false,
            _ => true,
        };

        let attrs = parse_gtf_attributes(attrs_s);
        let gene_id = match attrs.gene_id {
            Some(id) if !id.is_empty() => id,
            _ => {
                return Err(LangurError::GtfParse {
                    line: line_no,
                    message: "missing gene_id attribute".into(),
                })
            }
        };
        let tx_id = match attrs.transcript_id {
            Some(id) if !id.is_empty() => id,
            _ => {
                return Err(LangurError::GtfParse {
                    line: line_no,
                    message: "missing transcript_id attribute".into(),
                })
            }
        };

        // Low-TSL filtering (matching Whippet2 load_gtf behavior)
        if !opts.basic_or_ss_only && opts.suppress_low_tsl {
            if let Some(tsl_raw) = attrs.tsl {
                let tok = tsl_raw.split_whitespace().next().unwrap_or(tsl_raw);
                if tok != "_" && tok != "1" && tok != "2" && tok != "NA" {
                    continue;
                }
            }
        }

        let tsl_num = parse_tsl_u8(attrs.tsl);
        if !gene_seen.contains(gene_id) {
            gene_seen.insert(gene_id.to_owned());
            gene_order.push(gene_id.to_owned());
        }

        if let Some(acc) = tx_map.get_mut(tx_id) {
            if attrs.is_basic {
                acc.is_basic = true;
            }
            if acc.gene_name == "_" {
                if let Some(gn) = attrs.gene_name {
                    acc.gene_name = gn.to_owned();
                }
            }
            acc.exons.push((start, end));
        } else {
            tx_order.push(tx_id.to_owned());
            tx_map.insert(
                tx_id.to_owned(),
                TxAccum {
                    tx_id: tx_id.to_owned(),
                    gene_id: gene_id.to_owned(),
                    gene_name: attrs.gene_name.unwrap_or("_").to_owned(),
                    chrom: chrom.to_owned(),
                    strand,
                    is_basic: attrs.is_basic,
                    tsl: tsl_num,
                    exons: vec![(start, end)],
                },
            );
        }
    }

    // Aggregate transcripts into genes
    let mut gene_map: FxHashMap<String, GeneAccum> = FxHashMap::default();
    for tx_id in &tx_order {
        let mut tx = tx_map.remove(tx_id).expect("tx_id present");
        if tx.exons.is_empty() {
            continue;
        }
        let mut tranacc: Vec<u32> = tx.exons.iter().map(|&(s, _)| s).collect();
        let mut trandon: Vec<u32> = tx.exons.iter().map(|&(_, e)| e).collect();
        tranacc.sort_unstable();
        trandon.sort_unstable();

        let tx_s = tranacc[0];
        let tx_e = *trandon.last().unwrap();
        let tx_len: u32 = tx.exons.iter().map(|&(s, e)| e - s + 1).sum();

        let don_internal = &trandon[..trandon.len().saturating_sub(1)];
        let acc_internal = if tranacc.len() > 1 {
            &tranacc[1..]
        } else {
            &[]
        };
        let ss_only = opts.basic_or_ss_only && !tx.is_basic;

        let g = gene_map
            .entry(tx.gene_id.clone())
            .or_insert_with(|| GeneAccum {
                gene_id: tx.gene_id.clone(),
                gene_name: tx.gene_name.clone(),
                chrom: tx.chrom.clone(),
                strand: tx.strand,
                has_primary_tx: false,
                don: Vec::new(),
                acc: Vec::new(),
                txst: Vec::new(),
                txen: Vec::new(),
                exons: Vec::new(),
                tx_len_sum: 0,
                reftx: Vec::new(),
            });

        // A gene can only reside on one chromosome (matching Whippet2)
        if g.has_primary_tx && tx.chrom != g.chrom {
            continue;
        }
        if g.gene_name == "_" && tx.gene_name != "_" {
            g.gene_name = tx.gene_name.clone();
        }

        // Exons are tracked for all transcripts of the gene (matches Whippet2 gnexons)
        g.exons.extend_from_slice(&tx.exons);
        g.don.extend_from_slice(don_internal);
        g.acc.extend_from_slice(acc_internal);

        if !ss_only {
            if !g.has_primary_tx {
                g.chrom = tx.chrom;
                g.strand = tx.strand;
                g.has_primary_tx = true;
            }
            g.txst.push(tx_s);
            g.txen.push(tx_e);
            g.tx_len_sum += tx_len as u64;
            g.reftx.push(RefTx {
                tx_id: tx.tx_id,
                tx_start: tx_s,
                tx_end: tx_e,
                is_basic: tx.is_basic,
                tsl: tx.tsl,
                don: trandon,
                acc: tranacc,
                length: tx_len,
            });
        }
        tx.exons.clear();
    }

    let mut refset = RefSet::default();
    for gene_id in gene_order {
        let Some(mut g) = gene_map.remove(&gene_id) else {
            continue;
        };
        if !g.has_primary_tx || g.txst.is_empty() || g.txen.is_empty() {
            continue;
        }

        g.don.sort_unstable();
        g.don.dedup();
        g.acc.sort_unstable();
        g.acc.dedup();
        g.txst.sort_unstable();
        g.txst.dedup();
        g.txen.sort_unstable();
        g.txen.dedup();
        g.exons.sort_unstable();
        g.exons.dedup();

        let min_s = g.txst[0];
        let max_e = *g.txen.last().unwrap();

        // Remove any txst matching an acceptor, or txen matching a donor (Whippet2 refset.jl L274-275)
        let acc_set: FxHashSet<u32> = g.acc.iter().copied().collect();
        let don_set: FxHashSet<u32> = g.don.iter().copied().collect();
        g.txst.retain(|x| !acc_set.contains(x));
        g.txen.retain(|x| !don_set.contains(x));

        // Filter splice sites to gene range [min_s, max_e]
        g.don.retain(|&x| x >= min_s && x <= max_e);
        g.acc.retain(|&x| x >= min_s && x <= max_e);

        let mean_tx_len = if g.reftx.is_empty() {
            0.0
        } else {
            (g.tx_len_sum as f64 / g.reftx.len() as f64) as f32
        };

        let idx = refset.genes.len();
        refset.gene_index.insert(g.gene_id.clone(), idx);
        refset.genes.push(RefGene {
            gene_id: g.gene_id,
            gene_name: g.gene_name,
            chrom: g.chrom,
            strand: g.strand,
            don: g.don,
            acc: g.acc,
            txst: g.txst,
            txen: g.txen,
            exons: g.exons,
            novel_don: Vec::new(),
            novel_acc: Vec::new(),
            mean_tx_len,
            reftx: g.reftx,
        });
    }

    Ok(refset)
}

struct ParsedAttrs<'a> {
    gene_id: Option<&'a str>,
    transcript_id: Option<&'a str>,
    gene_name: Option<&'a str>,
    tsl: Option<&'a str>,
    is_basic: bool,
}

fn parse_gtf_attributes(s: &str) -> ParsedAttrs<'_> {
    let mut out = ParsedAttrs {
        gene_id: None,
        transcript_id: None,
        gene_name: None,
        tsl: None,
        is_basic: false,
    };

    for field in s.split(';') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }
        // Standard GTF2.2: `key "value"`, also tolerate `key=value`
        let (key, val) = if let Some((k, v)) = field.split_once(char::is_whitespace) {
            (k.trim(), v.trim().trim_matches('"'))
        } else if let Some((k, v)) = field.split_once('=') {
            (k.trim(), v.trim().trim_matches('"'))
        } else {
            continue;
        };

        match key {
            "gene_id" => out.gene_id = Some(val),
            "transcript_id" => out.transcript_id = Some(val),
            "gene_name" => out.gene_name = Some(val),
            "transcript_support_level" => out.tsl = Some(val),
            "tag" => {
                if val == "basic" {
                    out.is_basic = true;
                }
            }
            _ => {}
        }
    }
    out
}

fn parse_tsl_u8(raw: Option<&str>) -> u8 {
    let Some(s) = raw else {
        return 0;
    };
    let tok = s.split_whitespace().next().unwrap_or("");
    match tok {
        "1" => 1,
        "2" => 2,
        "3" => 3,
        "4" => 4,
        "5" => 5,
        _ => 0,
    }
}
