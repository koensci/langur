# Langur

**Fast, statistically calibrated Contiguous Splice Graph (CSG) quantification and differential splicing for bulk short-read RNA-seq in Rust.**

`langur` is a standalone, single-binary Rust suite (`index`, `inspect`, `quant`, `delta`) that profiles both binary ($K_1$) and arbitrarily complex ($K_n \ge 2$) alternative splicing events directly from FASTQ files—no genome alignment or BAM generation required.

- **Zero-copy `mmap` index (`.lidx`, magic `LANGUR01`)** with canonical minimizers, annotated junction-bridge seeds, and selective genomic decoys (`< 1 ms` load time).
- **Hybrid ungapped + banded Gotoh gapped spliced aligner** with two-tier repetitive minimizer screening, 1-SNP junction tolerance, and strict-inequality decoy screening (`decoy > CSG`).
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

# Recommended: enable host-native CPU optimizations (Apple Silicon arm64 or x86_64 AVX2)
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

**Cohort `--batch` mode** (maps `.lidx` once via `mmap` across all samples, keeping peak RAM flat at ~1.4 GB):

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

## Human-Scale Benchmark Highlights

### Simulated Human Cohort (`E-MTAB-3766`, GRCh37.71)

Benchmarked on the 480M-read human paired-end dataset (`E-MTAB-3766`, Soneson et al. *Genome Biology* 2016; Ensembl GRCh37.71 protein-coding: 20,410 genes, 145,342 transcripts, 1,000 truly differentially spliced genes; 16 threads, `--biascorrect`, Apple Silicon arm64) against original `Whippet.jl` v1.6.2 (Julia 1.6, 1T) and our multithreaded Julia 1.9 fork (`Whippet2-turbo`):

| Stage | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | Original `Whippet.jl` v1.6.2 *(1T)* | Langur Advantage |
|---|---:|---:|---:|---|
| **Index Build** *(uncompressed FASTA)* | **14.6 s** · 4.36 GB *(incl. 1.10 Gb decoys)* | 23.0 s · 4.21 GB *(no decoys)* | 287.1 s · 8.10 GB *(no decoys)* | **19.7× vs. W1 · 1.58× vs. W2** |
| **Index Load** | **0.03 ms** (`mmap`) | 6.20 s | 7.33 s | **>200,000× faster** |
| **Sample 1 `quant`** (80M reads) | **40.5 s** (`1.98M reads/s`) · **1.40 GB** | 81.1 s (`987k/s`) · 9.62 GB | 1,608.4 s (`49.7k/s`) · 7.19 GB | **39.7× vs. W1 · 2.00× vs. W2**; **6.9× less RAM** |
| **6-Sample `--batch`** (480M reads) | **248.9 s** (`1.93M reads/s`) · **1.42 GB** | 433.8 s (`1.11M/s`) · 9.83 GB | *9,650.6 s* *(6× extrap.)* · 7.19 GB | **38.8× vs. W1 · 1.74× vs. W2**; **6.9× less RAM** |
| **Gene / Isoform TPM MARD** (vs. RSEM truth) | **0.0137 / 0.0402** | 0.0443 / 0.0687 | 0.0753 / 0.1145 | **5.5× / 2.85× lower error vs. W1** |
| **`delta` GLM** (223k events, `-r 5 -s 2`) | **2.50 s** · **191 MB** | 27.38 s · 1.87 GB | 85.51 s *(`--legacy` MC)* · 609 MB | **34.2× vs. legacy · 11.0× & 9.8× RAM vs. W2** |
| **Differential Splicing** (`TPR / FDR @0.05` · AUPRC) | **0.336 / 5.4%** *(0.338 / 2.0% spliced)* · **0.565** *(0.567 spliced)* | 0.306 / 5.3% *(0.309 / 1.3% spliced)* · 0.551 *(0.549 spliced)* | 0.744 / 81.5% *(66.9% at `\|ΔΨ\|≥0.1`)* · 0.352 | **15× lower FDR vs. legacy**; **10.6× lower FDR & +0.17 AUPRC vs. `rMATS`**; **beats `Whippet2-turbo` on TPR, AUROC & AUPRC** |

### Real Human RNA-seq vs. PacBio HiFi Ground Truth (LRGASP WTC11, GRCh38)

Benchmarked on real human WTC11 Illumina paired-end RNA-seq (ENCODE `ENCSR673UKZ` Rep 1: `45.78M` pairs = `91.57M` reads, `--biascorrect`, 16T) against full GENCODE v38 basic (`60,649` genes, `111,832` transcripts, GRCh38) alongside genome-aligned **`STAR 2.7.11b` + `rMATS-turbo v4.3.0`** (16T), with orthogonal long-read ground-truth PSI from pooled 3-replicate matched PacBio HiFi Iso-Seq alignments (ENCODE `ENCSR507JOF`; `LR ≥ 20`, `SR ≥ 10`):

| Metric (GRCh38 / WTC11) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | Original `Whippet.jl` v1.6.2 *(1T)* | `STAR 2.7.11b` + `rMATS-turbo v4.3.0` *(16T)* | Langur Advantage |
|---|---:|---:|---:|---:|---|
| **GRCh38 Index Build** *(60,649 genes, uncompressed FASTA)* | **14.8 s** · **5.19 GB** · **1.06 GB** *(incl. 1.31 Gb decoys)* | 28.6 s · 5.33 GB · 691 MB *(no decoys)* | 753.8 s · 9.87 GB · 634 MB | 683.1 s *(11m 23s)* · 35.09 GB · **28 GiB** | **50.9× vs. W1 · 46.2× vs. STAR · 6.8× less RAM & 27× smaller index** |
| **WTC11 Rep 1 End-to-End** *(91.6M reads)* | **72.7 s** (`1.26M reads/s`) · **1.35 GB** · **0 GB BAM** | 82.3 s (`1.11M/s`) · 9.10 GB · 0 GB | 1,467.4 s (`62.4k/s`) · 7.00 GB · 0 GB | 526.1 s *(`357.0 s` STAR + `169.1 s` rMATS)* · **31.64 GB** · **4.8 GiB BAM** | **7.24× faster & 23.5× less RAM vs. STAR+rMATS**; **zero BAMs** |
| **Mapped Pairs** *(of 45.78M)* | **30.00M (65.5%)** *(382k decoy-filtered)* | 24.70M (53.9%) | 24.70M (53.9%) | 30.55M (66.7%) unique + 1.87M multi *(genome-wide)* | **+5.30M pairs vs. W2/W1**; matches STAR unique exon capture |
| **PacBio Polymorphic PSI (3-way CSG, `n = 13,359`)** | **0.9338 / 0.8715 / 0.0789** | 0.9316 / 0.8689 / 0.0821 | 0.9064 / 0.8473 / 0.1023 | *(binary `SE/A3SS/A5SS/RI` only; no `AF`/`AL`/`TS`/`TE`/TPMs)* | **27/27 clean sweep** across all 6 motifs + `K0/K1` & `K_n ≥ 2` |
| **4-Way Shared Binary (`K0/K1`, `n = 3,106`)** | **0.9441 / 0.8774 / 0.0709** | 0.9417 / 0.8753 / 0.0732 | 0.9209 / 0.8594 / 0.0925 | 0.9440 / 0.8721 / 0.0715 *(`JC`)* | **Parity with STAR+rMATS** at 7.2× speed & 23.5× less RAM |
| **4-Way Shared Complex (`K_n ≥ 2`, `n = 4,709`)** | **0.9436 / 0.8996 / 0.0672** | 0.9413 / 0.8976 / 0.0698 | 0.9095 / 0.8645 / 0.0925 | 0.8774 / 0.8009 / 0.0989 *(`JC`)* | **+0.066 $r$ & -32% MAE** on full-length complex exon inclusion |

*`Whippet2-turbo` index builds reflect commit `94a491f`; `E-MTAB-3766` `quant` timing/memory reflects commit `400f0f3` and TPM/stratified GLM reflects commit `c725028` (identical quantification accuracy). All memory figures report macOS `/usr/bin/time -l` peak memory footprint. See **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)** for full stratified tables (including exact rMATS-junction sensitivity checks), ablations, and mathematical specifications.*

---

## Attribution, Transparency & License

- **Credit & Lineage:** Langur is an independent Rust re-architecture building with deep gratitude on the Contiguous Splice Graph (CSG) representation and alternative-splicing motif taxonomy pioneered by **Tim Sterne-Weiler** and colleagues in **Whippet** (*Molecular Cell*, 2018; `Whippet.jl`), alongside independent implementations of published statistical concepts from `Salmon` and `edgeR`/`DESeq2`, plus overhang-diversity QC for PCR junction-stack artifacts highlighted in `MAJIQ`.
- **AI-Assisted Development:** Langur's architecture, Rust implementation, test suite, and benchmarks were developed iteratively with AI-assisted coding and review (**Google Gemini**), with all biological invariants, edge cases, and human-scale benchmarks verified through automated testing and ground-truth evaluation.
- **License:** Released under the **MIT License** — see [LICENSE](LICENSE).
