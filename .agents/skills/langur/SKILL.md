---
name: langur
description: >-
  Use this skill when building splice-graph indexes (.lidx), quantifying bulk
  short-read RNA-seq alternative splicing (PSI) or gene/isoform TPMs, or running
  covariate-aware Beta-Binomial GLM differential splicing analysis with Langur
  (`langur index`, `langur inspect`, `langur quant`, `langur delta`).
---

# Langur Splice-Graph RNA-seq Workflow

Langur (`v1.0.0`) is a multithreaded, deterministic Rust engine for alignment-free Contiguous Splice Graph (CSG) splicing quantification and covariate-aware differential splicing directly from bulk short-read FASTQ files (~1.4 GB RAM, zero BAMs).

## Scope & Tool Selection

- **Use Langur for:** Bulk paired-end or single-end Illumina RNA-seq splicing (`CE`, `AA`, `AD`, `RI`, `AF`, `AL`, `TS`, `TE`, optional `BS`), gene/isoform TPMs, and multivariable differential splicing.
- **Hybrid novel/cryptic splicing:** Pass a coordinate-sorted, indexed (`.bai`) BAM from a representative genome-aligned subset (or a custom long-read GTF) to `langur index --bam` to compile novel splice sites once into `.lidx`, then quantify the full cohort alignment-free. Sample-specific unannotated junctions absent from the discovery subset still require genome alignment.
- **Do not use for:** Single-cell RNA-seq or full-length long-read (PacBio HiFi / ONT) quantification (use Isosceles, Oarfish, or IsoQuant).

## 1. Build and Verify a `.lidx` Index

Ensure chromosome naming matches between genome FASTA and GTF (e.g., both `chr1` or both `1`).

```bash
# Standard annotation index (includes selective genomic decoys by default)
langur index \
  --fasta genome.fa \
  --gtf annotations.gtf \
  -x index/my_index \
  -t 16

# Optional: augment with novel splice junctions from a sorted+indexed BAM
langur index \
  --fasta genome.fa \
  --gtf annotations.gtf \
  --bam discovery.sorted.bam \
  --bam-min-reads 3 \
  -x index/my_index_bam \
  -t 16

# Verify the index (NOTE: --index / -x flag is required; not positional!)
langur inspect --index index/my_index.lidx
```

Confirm `langur inspect` reports `CRC32` OK, `decoys=true`, `jnc_seeds=true`, expected gene/transcript counts, and `0 hash collisions`.

## 2. Quantify Samples (`langur quant`)

Always pass `--biascorrect` (unified 5'-hexamer + 20-bin GC bias) for real bulk RNA-seq.

### Strandness Selection (Critical)

| Library Protocol | Flag |
|---|---|
| Unstranded | *(omit strand flag)* |
| Reverse-stranded RF (Illumina TruSeq / dUTP; kallisto `--rf-stranded`; Salmon `-l ISR`; STAR/MAJIQ `REVERSE`) | `--stranded-reverse` |
| Forward-stranded FR (kallisto `--fr-stranded`; Salmon `-l ISF`) | `--stranded` |

If strandness is uncertain, test a small subset (e.g., 500k–1M pairs) with `--stranded-reverse` vs. `--stranded`; the wrong strand typically drops mapped rate >10×.

### Single Sample vs. Multi-Sample Batch Mode

For cohorts, prefer `--batch` so the `.lidx` index is memory-mapped once across all samples:

```bash
# Single paired-end sample
langur quant \
  -x index/my_index.lidx \
  --biascorrect \
  --stranded-reverse \
  -t 16 \
  -o results/sample1 \
  sample1_R1.fastq.gz sample1_R2.fastq.gz

# Single-end sample
langur quant \
  -x index/my_index.lidx \
  --biascorrect \
  -t 16 \
  -o results/sample1 \
  sample1.fastq.gz
```

**Batch manifest format** (CSV or TSV; header optional if first column is `sample` or `sample_id`; 4 columns for paired-end, 3 columns for single-end):

```csv
sample,r1,r2,out_prefix
S1,/path/to/S1_1.fastq.gz,/path/to/S1_2.fastq.gz,/path/to/langur_output/S1.output
S2,/path/to/S2_1.fastq.gz,/path/to/S2_2.fastq.gz,/path/to/langur_output/S2.output
```

```bash
langur quant \
  -x index/my_index.lidx \
  --batch manifest.csv \
  --biascorrect \
  --stranded-reverse \
  -t 16
```

Each sample writes `<out_prefix>.psi.gz`, `.jnc.gz`, `.gene.tpm.gz`, `.isoform.tpm.gz`, and `.map.gz`.

## 3. Differential Splicing (`langur delta`)

`langur delta` fits a Beta-Binomial GLM with Cox–Reid Empirical Bayes dispersion shrinkage, automatic Spliced vs. Tandem UTR (`TS`/`TE`) prior stratification, and column-15 `Inf_Var` weighting for complex ($K_n \ge 2$) events.

### Two-Group Comparison (`DeltaPsi = Psi_A - Psi_B`)

```bash
langur delta \
  -a treated_1.psi.gz,treated_2.psi.gz,treated_3.psi.gz \
  -b ctrl_1.psi.gz,ctrl_2.psi.gz,ctrl_3.psi.gz \
  -r 5 -s 2 \
  -t 16 \
  -o results/treated_vs_ctrl
```

### Covariate-Aware Multivariable GLM (`--design` + `--formula`)

Prepare a TSV/CSV metadata table with a `.psi.gz` path column (default name `file`, override with `--file-col`):

```tsv
sample	file	batch	donor	condition
S1	langur_output/S1.output.psi.gz	B1	D1	ctrl
S2	langur_output/S2.output.psi.gz	B1	D1	treated
S3	langur_output/S3.output.psi.gz	B2	D2	ctrl
S4	langur_output/S4.output.psi.gz	B2	D2	treated
```

```bash
langur delta \
  --design design.tsv \
  --file-col file \
  --formula "~ batch + donor + condition" \
  --test condition:treated_vs_ctrl \
  -r 5 -s 2 \
  -t 16 \
  -o results/treated_vs_ctrl_adj
```

- `--test` accepts a factor name (`condition`), explicit contrast (`condition:treated_vs_ctrl`, where positive `DeltaPsi` = higher in `treated`), interaction (`batch:condition`), or continuous predictor (`age`). Defaults to the last term in `--formula`.
- Writes `<out>.diff.gz`.

## 4. Verification & Common Pitfalls

1. **Validate outputs before deleting intermediate FASTQs:** Check non-empty `.psi.gz` with `gzip -t <file>.psi.gz` and inspect mapping rate in `<out_prefix>.map.gz`.
2. **Canonical splicing vs. tandem UTRs:** `.psi.gz` and `.diff.gz` include both canonical splice events (`CE`, `AA`, `AD`, `RI`, `AF`, `AL`) and tandem transcription start/end sites (`TS`, `TE`). Filter to `Type %in% c("CE","AA","AD","RI","AF","AL")` when analyzing canonical splice junctions only.
3. **Never mix `.psi.gz` files across different indexes or tools** (e.g., Langur vs. legacy Whippet.jl) in a single `langur delta` run—node numbering must match and `Inf_Var` is Langur-specific.
4. **macOS Bash scripting gotcha:** When wrapping `langur quant` in `#!/usr/bin/env bash` with `set -u`, expand optional flag arrays as `${STRAND_FLAGS[@]+"${STRAND_FLAGS[@]}"}` for macOS Bash 3.2 compatibility.

## Detailed Reference

For full column schemas (`.psi.gz`, `.diff.gz`, `.jnc.gz`, `.tpm.gz`), event-type definitions, and advanced CLI flags, read [references/outputs-and-flags.md](./references/outputs-and-flags.md).
