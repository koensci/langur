# Langur

**Fast, statistically calibrated Contiguous Splice Graph (CSG) quantification and differential splicing for bulk short-read RNA-seq in Rust.**

`langur` is a standalone, single-binary Rust suite (`index`, `inspect`, `quant`, `delta`) that profiles both binary ($K_1$) and arbitrarily complex ($K_n \ge 2$) alternative splicing events directly from FASTQ files—no genome alignment or BAM generation required.

- **Zero-copy `mmap` index (`.lidx`, magic `LANGUR01`)** with canonical minimizers, annotated junction-bridge seeds, and selective genomic decoys (`< 1 ms` load time).
- **Hybrid ungapped + SIMD gapped spliced aligner** (`NEON` / `AVX2`) with 1-SNP junction tolerance and strict-inequality decoy screening (`decoy > CSG`).
- **Two-stage quantification**: joint Connected-Component Variational Bayes EM (`CC-VBEM`) with unified 5'-hexamer + 20-bin GC bias correction (`--biascorrect`), microexon-safe effective lengths, overhang-diversity PCR-stack soft-capping, and targeted collapsed Gibbs inferential variance (`Inf_Var`, column 15 of `.psi.gz`).
- **Covariate-aware Beta-Binomial GLM (`delta`)**: supports both two-group (`-a` / `-b`) and DESeq2-style multivariable formulas (`--design` / `--formula` / `--test`), with constitutive-event exclusion and Spliced vs. Tandem UTR (`TS`/`TE`) stratified Cox–Reid Empirical Bayes dispersion shrinkage.
- **100% bit-for-bit deterministic** across thread counts and standalone vs. `--batch` runs.

For full algorithmic derivations, design rationale, and validation details, see **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)**.

---

## Installation & Build

Requires a recent Rust toolchain (Edition 2021, `cargo >= 1.75`). Optional: installing `igzip` (Intel ISA-L) on Unix enables automatic accelerated `.fq.gz` streaming (seamlessly falls back to built-in `zlib-ng` when unavailable).

```bash
# Standard optimized release build
cargo build --release

# Recommended: enable host-native SIMD (Apple Silicon NEON/LSE or x86_64 AVX2)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# Binary is placed at ./target/release/langur
./target/release/langur --help
```

Run the automated test suite (42 unit & integration tests):

```bash
cargo test --workspace
```

---

## Quick Start

### 1. Build a Splice-Graph + Selective Decoy Index (`.lidx`)

```bash
langur -t 16 index \
  --fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
  --gtf Homo_sapiens.GRCh37.gtf.gz \
  -x human_grch37
# Writes human_grch37.lidx, human_grch37.exons.tsv.gz, and human_grch37.paths.tsv.gz
```

Optional de novo splice-site supplementation from a coordinate-sorted, indexed BAM:

```bash
langur -t 16 index \
  --fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  --bam cohort_merged.bam --bam-min-reads 3 \
  -x human_plus_bam
```

### 2. Inspect Index Metadata & Load Latency

```bash
langur inspect -x human_grch37.lidx
```

### 3. Quantify Single Samples or Multi-Sample Cohorts

**Single paired-end sample** (unstranded; add `--stranded` for FR or `--stranded-reverse` for RF / TruSeq dUTP):

```bash
langur -t 16 quant \
  sample1_R1.fq.gz sample1_R2.fq.gz \
  -x human_grch37.lidx \
  --biascorrect \
  -o results/sample1
# Writes results/sample1.{psi,gene.tpm,isoform.tpm,jnc,map}.gz
```

**Cohort `--batch` mode** (maps `.lidx` once via `mmap` across all samples, keeping peak RAM flat at ~1.7 GB):

```csv
sample,r1,r2,out_prefix
S1,fastq/S1_1.fq.gz,fastq/S1_2.fq.gz,results/S1
S2,fastq/S2_1.fq.gz,fastq/S2_2.fq.gz,results/S2
S3,fastq/S3_1.fq.gz,fastq/S3_2.fq.gz,results/S3
S4,fastq/S4_1.fq.gz,fastq/S4_2.fq.gz,results/S4
```

```bash
langur -t 16 quant \
  --batch cohort.csv \
  -x human_grch37.lidx \
  --biascorrect
```

### 4. Differential Splicing (`delta`)

**Two-group comparison (`-a` vs. `-b`, `DeltaPsi = Psi_A - Psi_B`):**

```bash
langur -t 16 delta \
  -a S1.psi.gz,S2.psi.gz,S3.psi.gz \
  -b S4.psi.gz,S5.psi.gz,S6.psi.gz \
  -d results/ \
  -r 5 -s 2 \
  -o results/treated_vs_ctrl
# Writes results/treated_vs_ctrl.diff.gz
```

**Multivariable / covariate-aware GLM (`--design`, `--formula`, `--test`):**

```tsv
file	batch	donor	condition
S1.psi.gz	B1	D1	ctrl
S2.psi.gz	B1	D2	ctrl
S3.psi.gz	B2	D3	ctrl
S4.psi.gz	B1	D1	treated
S5.psi.gz	B1	D2	treated
S6.psi.gz	B2	D3	treated
```

```bash
langur -t 16 delta \
  --design design.tsv \
  --formula "~ batch + donor + condition" \
  --test condition:treated_vs_ctrl \
  -d results/ \
  -r 5 -s 2 \
  -o results/treated_vs_ctrl_adjusted
```

---

## Human-Scale Benchmark Highlights (`E-MTAB-3766`)

Benchmarked on the 480M-read human paired-end dataset (`E-MTAB-3766`, Soneson et al. *Genome Biology* 2016; Ensembl GRCh37.71 protein-coding: 20,410 genes, 145,342 transcripts, 1,000 truly differentially spliced genes; 16 threads, `--biascorrect`, Apple Silicon arm64) against original `Whippet.jl` v1.6.2 (Julia 1.6, 1T) and our multithreaded Julia 1.9 fork (`Whippet2-turbo`):

| Stage | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | Original `Whippet.jl` v1.6.2 *(1T)* | Langur Advantage |
|---|---:|---:|---:|---|
| **Index Build** | **22.6 s** · 4.71 GB RSS | 233.0 s · 8.74 GB RSS | 340.8 s · 8.58 GB RSS | **15.1× vs. W1 · 10.3× vs. W2** |
| **Index Load** | **0.04–0.30 ms** (`mmap`) | 6.20 s | 7.33 s | **>20,000× faster** |
| **Sample 1 `quant`** (80M reads) | **60.4 s** (`1.33M reads/s`) · **1.69 GB** | 81.1 s (`987k/s`) · 9.62 GB | 1,608.4 s (`49.7k/s`) · 7.19 GB | **26.6× vs. W1 · 1.34× vs. W2**; **5.7× less RAM** |
| **6-Sample `--batch`** (480M reads) | **373.1 s** (`1.29M reads/s`) · **1.67 GB** | 433.8 s (`1.11M/s`) · 9.83 GB | *9,650.6 s* *(6× extrap.)* · 7.19 GB | **25.9× vs. W1 · 1.16× vs. W2**; **5.9× less RAM** |
| **Gene / Isoform TPM MARD** (vs. RSEM truth) | **0.0137 / 0.0406** | 0.0443 / 0.0687 | 0.0753 / 0.1145 | **5.5× / 2.8× lower error vs. W1** |
| **`delta` GLM** (223k events, `-r 5 -s 2`) | **2.52 s** · **189 MB** | 27.38 s · 1.87 GB | 85.51 s *(`--legacy` MC)* · 609 MB | **33.9× vs. legacy · 10.9× & 9.9× RAM vs. W2** |
| **Differential Splicing** (`TPR / FDR @0.05` · AUPRC) | **0.336 / 5.4%** *(2.0% spliced)* · **0.564** *(0.567 spliced)* | 0.306 / 5.3% *(1.3% spliced)* · 0.551 *(0.549 spliced)* | 0.744 / 81.5% *(66.9% at `\|ΔΨ\|≥0.1`)* · 0.352 | **15× lower FDR vs. legacy**; **10.6× lower FDR & +0.17 AUPRC vs. `rMATS`**; **beats `Whippet2-turbo` on TPR, AUROC & AUPRC** |

*`Whippet2-turbo` timing and memory rows reflect commit `400f0f3`; TPM and stratified GLM rows reflect commit `c725028` (identical quantification accuracy). See **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)** for full tables, ablations, and mathematical specifications.*

---

## Attribution, Transparency & License

- **Credit & Lineage:** Langur is an independent Rust re-architecture building with deep gratitude on the Contiguous Splice Graph (CSG) representation and alternative-splicing motif taxonomy pioneered by **Tim Sterne-Weiler** and colleagues in **Whippet** (*Molecular Cell*, 2018; `Whippet.jl`), alongside independent implementations of published statistical concepts from `Salmon` and `edgeR`/`DESeq2`, plus overhang-diversity QC for PCR junction-stack artifacts highlighted in `MAJIQ`.
- **AI-Assisted Development:** Langur's architecture, Rust implementation, test suite, and benchmarks were developed iteratively with AI-assisted coding and review (**Google Gemini**), with all biological invariants, edge cases, and human-scale benchmarks verified through automated testing and ground-truth evaluation.
- **License:** Released under the **MIT License** — see [LICENSE](LICENSE).
