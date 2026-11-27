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

## Human-Scale Benchmark Highlights

### Simulated Human Cohort (`E-MTAB-3766`, GRCh37.71)

Benchmarked on the 480M-read human paired-end dataset (`E-MTAB-3766`, Soneson et al. *Genome Biology* 2016; Ensembl GRCh37.71 protein-coding: 20,410 genes, 145,342 transcripts, 1,000 truly differentially spliced genes; 16 threads, `--biascorrect`, Apple Silicon arm64) against original `Whippet.jl` v1.6.2 (Julia 1.6, 1T) and our multithreaded Julia 1.9 fork (`Whippet2-turbo`):

| Stage | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | Original `Whippet.jl` v1.6.2 *(1T)* | Langur Advantage |
|---|---:|---:|---:|---|
| **Index Build** | **22.6 s** · 4.71 GB RSS | 27.6 s · 4.90 GB RSS | 340.8 s · 8.58 GB RSS | **15.1× vs. W1 · 1.22× vs. W2** *(incl. 1.1 Gb decoys)* |
| **Index Load** | **0.04–0.30 ms** (`mmap`) | 6.20 s | 7.33 s | **>20,000× faster** |
| **Sample 1 `quant`** (80M reads) | **60.4 s** (`1.33M reads/s`) · **1.69 GB** | 81.1 s (`987k/s`) · 9.62 GB | 1,608.4 s (`49.7k/s`) · 7.19 GB | **26.6× vs. W1 · 1.34× vs. W2**; **5.7× less RAM** |
| **6-Sample `--batch`** (480M reads) | **373.1 s** (`1.29M reads/s`) · **1.67 GB** | 433.8 s (`1.11M/s`) · 9.83 GB | *9,650.6 s* *(6× extrap.)* · 7.19 GB | **25.9× vs. W1 · 1.16× vs. W2**; **5.9× less RAM** |
| **Gene / Isoform TPM MARD** (vs. RSEM truth) | **0.0137 / 0.0406** | 0.0443 / 0.0687 | 0.0753 / 0.1145 | **5.5× / 2.8× lower error vs. W1** |
| **`delta` GLM** (223k events, `-r 5 -s 2`) | **2.52 s** · **189 MB** | 27.38 s · 1.87 GB | 85.51 s *(`--legacy` MC)* · 609 MB | **33.9× vs. legacy · 10.9× & 9.9× RAM vs. W2** |
| **Differential Splicing** (`TPR / FDR @0.05` · AUPRC) | **0.336 / 5.4%** *(2.0% spliced)* · **0.564** *(0.567 spliced)* | 0.306 / 5.3% *(1.3% spliced)* · 0.551 *(0.549 spliced)* | 0.744 / 81.5% *(66.9% at `\|ΔΨ\|≥0.1`)* · 0.352 | **15× lower FDR vs. legacy**; **10.6× lower FDR & +0.17 AUPRC vs. `rMATS`**; **beats `Whippet2-turbo` on TPR, AUROC & AUPRC** |

### Real Human RNA-seq vs. PacBio HiFi Ground Truth (LRGASP WTC11, GRCh38)

Benchmarked on real human WTC11 Illumina paired-end RNA-seq (ENCODE `ENCSR673UKZ` Rep 1: `45.78M` pairs = `91.57M` reads, `--biascorrect`, 16T) against full GENCODE v38 basic (`60,649` genes, `111,832` transcripts, GRCh38) with orthogonal long-read ground-truth PSI computed from pooled 3-replicate matched PacBio HiFi Iso-Seq alignments (ENCODE `ENCSR507JOF`; `LR ≥ 20`, `SR ≥ 10` across all 3 tools):

| Metric (GRCh38 / WTC11) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | Original `Whippet.jl` v1.6.2 *(1T)* | Langur Advantage |
|---|---:|---:|---:|---|
| **GRCh38 Index Build** (60,649 genes) | **14.8 s** · **5.19 GB** *(incl. 1.31 Gb decoys)* | 28.6 s · 5.33 GB *(no decoys)* | 753.8 s *(12m 34s)* · 9.87 GB | **50.9× vs. W1 · 1.93× vs. W2** |
| **GRCh38 Index Load** | **0.11 ms** (`mmap`) | 1.83 s | 1.71 s | **>15,000× faster** |
| **WTC11 Rep 1 `quant`** (91.6M reads) | **72.7 s** (`1.26M reads/s`) · **1.35 GB** | 82.3 s (`1.11M/s`) · 9.10 GB | 1,467.4 s (`62.4k/s`) · 7.00 GB | **20.2× vs. W1 · 1.13× vs. W2**; **6.75× less RAM** |
| **Mapped Pairs** (of 45.78M) | **30.00M (65.5%)** *(382k decoy-filtered)* | 24.70M (53.9%) | 24.70M (53.9%) | **+5.30M more mapped pairs** (`+11.6 pp`) |
| **PacBio HiFi Polymorphic PSI** (`n = 13,359`, `r / ρ / MAE`) | **0.9338 / 0.8715 / 0.0789** | 0.9316 / 0.8689 / 0.0821 | 0.9064 / 0.8473 / 0.1023 | **27/27 clean sweep** across all 6 motifs + `K0/K1` & `K_n ≥ 2` |
| **Complex Events (`K_n ≥ 2`)** (`n = 7,294`, `r / ρ / MAE`) | **0.9352 / 0.8995 / 0.0746** | 0.9333 / 0.8974 / 0.0780 | 0.9018 / 0.8662 / 0.1017 | **+0.0334 $r$ & -27% MAE vs. W1** |

*`Whippet2-turbo` index builds reflect commit `94a491f`; `E-MTAB-3766` `quant` timing/memory reflects commit `400f0f3` and TPM/stratified GLM reflects commit `c725028` (identical quantification accuracy). Memory figures report macOS `/usr/bin/time -l` peak memory footprint (index build reports peak RSS on `E-MTAB-3766`). See **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)** for full stratified tables, ablations, and mathematical specifications.*

---

## Attribution, Transparency & License

- **Credit & Lineage:** Langur is an independent Rust re-architecture building with deep gratitude on the Contiguous Splice Graph (CSG) representation and alternative-splicing motif taxonomy pioneered by **Tim Sterne-Weiler** and colleagues in **Whippet** (*Molecular Cell*, 2018; `Whippet.jl`), alongside independent implementations of published statistical concepts from `Salmon` and `edgeR`/`DESeq2`, plus overhang-diversity QC for PCR junction-stack artifacts highlighted in `MAJIQ`.
- **AI-Assisted Development:** Langur's architecture, Rust implementation, test suite, and benchmarks were developed iteratively with AI-assisted coding and review (**Google Gemini**), with all biological invariants, edge cases, and human-scale benchmarks verified through automated testing and ground-truth evaluation.
- **License:** Released under the **MIT License** — see [LICENSE](LICENSE).
