# Langur

[![CI](https://github.com/koensci/langur/actions/workflows/ci.yml/badge.svg)](https://github.com/koensci/langur/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24615835.svg)](https://doi.org/10.5281/zenodo.24615835)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Fast, statistically calibrated Contiguous Splice Graph (CSG) quantification and differential splicing for bulk short-read RNA-seq in Rust.**

`langur` (`index`, `inspect`, `quant`, `delta`) profiles both binary ($K_1$) and complex ($K_n \ge 2$) alternative splicing—alongside gene and isoform expression—directly from raw FASTQ files in ~1.4 GB of memory, with zero BAM files and 100% bit-for-bit determinism.

- **Instant memory-mapped index (`.lidx`, `< 1 ms` load):** Combines canonical minimizers, annotated junction-bridge seeds, and selective genomic decoys into a compact 1.06 GB GRCh38 index.
- **Gapped, strand-aware spliced alignment:** Pairs fast ungapped matching with banded Gotoh indel recovery, 1-SNP splice-boundary tolerance, strict-inequality pseudogene filtering, and native unstranded / FR / RF (`--stranded-reverse`) library support.
- **Joint expression and local splicing quantification (`quant`):** Uses connected-component variational Bayes EM with 5'-hexamer and GC bias correction (`--biascorrect`), unbiased 3–27 bp microexon scoring, and collapsed Gibbs inferential variance (`Inf_Var`) at complex splice hubs.
- **Calibrated differential splicing (`delta`):** Fits covariate-aware Beta-Binomial GLMs (two-group `-a`/`-b` or DESeq2-style `--design`/`--formula`/`--test`) with event-class stratified Empirical Bayes dispersion shrinkage (Spliced vs. Tandem UTR) to eliminate real-cohort false positives and adjust for batch or donor effects.

For full algorithmic details, mathematical derivations, and complete benchmark tables, see **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)**.

---

## Build & Test

Requires Rust Edition 2021 (`cargo >= 1.75`). Optional: `igzip` (Intel ISA-L) accelerates `.fq.gz` streaming and automatically falls back to built-in `zlib-ng`.

```bash
# Recommended: enable host-native CPU optimizations (Apple Silicon arm64 or x86_64 AVX2)
RUSTFLAGS="-C target-cpu=native" cargo build --release
cargo test --workspace   # 42 unit & integration tests
./target/release/langur --help
```

---

## Quick Start

### 1. Build & Inspect Index (`.lidx`)

```bash
langur -t 16 index \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --gtf Homo_sapiens.GRCh38.gtf.gz \
  -x human_grch38
# Optional novel splice sites from aligned BAMs: add --bam cohort_merged.bam --bam-min-reads 3
# Writes human_grch38.{lidx,exons.tsv.gz,paths.tsv.gz}

langur inspect -x human_grch38.lidx
```

### 2. Quantify Single Sample or Cohort (`--batch`)

```bash
# Single sample (unstranded; add --stranded for FR or --stranded-reverse for RF / TruSeq dUTP)
langur -t 16 quant sample1_R1.fq.gz sample1_R2.fq.gz \
  -x human_grch38.lidx --biascorrect -o results/sample1
# Writes results/sample1.{psi,gene.tpm,isoform.tpm,jnc,map}.gz

# Cohort --batch mode (maps .lidx once via mmap; peak RAM stays ~1.4 GB)
langur -t 16 quant --batch cohort.csv -x human_grch38.lidx --biascorrect
```

<details>
<summary><b>Example <code>cohort.csv</code> &amp; <code>design.tsv</code></b></summary>

```csv
sample,r1,r2,out_prefix
S1,fastq/S1_1.fq.gz,fastq/S1_2.fq.gz,results/S1
S2,fastq/S2_1.fq.gz,fastq/S2_2.fq.gz,results/S2
S3,fastq/S3_1.fq.gz,fastq/S3_2.fq.gz,results/S3
S4,fastq/S4_1.fq.gz,fastq/S4_2.fq.gz,results/S4
```

```tsv
file	batch	donor	condition
S1.psi.gz	B1	D1	ctrl
S2.psi.gz	B1	D2	ctrl
S3.psi.gz	B2	D3	ctrl
S4.psi.gz	B1	D1	treated
S5.psi.gz	B1	D2	treated
S6.psi.gz	B2	D3	treated
```

</details>

### 3. Differential Splicing (`delta`)

```bash
# Two-group comparison (DeltaPsi = Psi_A - Psi_B)
langur -t 16 delta \
  -a S1.psi.gz,S2.psi.gz,S3.psi.gz \
  -b S4.psi.gz,S5.psi.gz,S6.psi.gz \
  -d results/ -r 5 -s 2 -o results/treated_vs_ctrl

# Multivariable covariate-adjusted GLM
langur -t 16 delta \
  --design design.tsv --formula "~ batch + donor + condition" \
  --test condition:treated_vs_ctrl \
  -d results/ -r 5 -s 2 -o results/treated_vs_ctrl_adjusted
```

---

## Benchmark Highlights

All full-human GRCh38 benchmarks (GENCODE v38 basic, 60,649 genes) were run on a 16-core Apple M4 Max (64 GB RAM; 16 threads for multithreaded tools; disk sizes in decimal GB, $1\text{ GB} = 10^9$ bytes). Exact software versions and full commands are documented in [docs/WHITEPAPER.md](docs/WHITEPAPER.md).

### 1. Real Human Long-Read Validation (LRGASP WTC11, 6-Way)

Evaluated on real human WTC11 paired-end RNA-seq (ENCODE `ENCSR673UKZ` Rep 1, 91.6M reads, `--biascorrect`) against pooled 3-replicate matched PacBio HiFi Iso-Seq alignments (ENCODE `ENCSR507JOF`; PSI accuracy reported as `r / ρ / MAE`):

| Metric (GRCh38 / WTC11) | **Langur** *(16T)* | `Whippet2-turbo` *(16T)* | `Whippet.jl` *(1T)* | `Salmon 2 + SUPPA2` *(16T)* | `STAR + MAJIQ v3` *(16T)* | `STAR + rMATS-turbo` *(16T)* |
|---|---:|---:|---:|---:|---:|---:|
| **Index time / disk** | **14.8 s / 1.06 GB** | 28.6 s / 0.69 GB | 753.8 s / 0.63 GB | 110.8 s / 8.51 GB | 700.7 s / 29.5 GB | 683.1 s / 29.5 GB |
| **End-to-end time / RAM** | **72.7 s / 1.35 GB** | 82.3 s / 9.10 GB | 1,467.4 s / 7.00 GB | 190.5 s / 12.24 GB | 370.9 s¹ / 31.64 GB *(+5.1 GB BAM)* | 526.1 s / 31.64 GB *(+5.1 GB BAM)* |
| **Mapped pairs** *(of 45.78M)* | **30.00M (65.5%)** | 24.70M (53.9%) | 24.70M (53.9%) | **30.44M (66.5%)** | 30.55M (66.7%) uniq. + 1.87M multi | 30.55M (66.7%) uniq. + 1.87M multi |
| **3-way CSG non-constitutive** (`n = 13,359`, 6 event types) | **0.9338 / 0.8715 / 0.0789** | 0.9316 / 0.8689 / 0.0821 | 0.9064 / 0.8473 / 0.1023 | — | *(LSV junctions + IR; no TPMs)* | *(no `AF`/`AL`/`TS`/`TE` or TPMs)* |
| **6-way shared non-constitutive** (`n = 3,672`) | **0.9568 / 0.9402 / 0.0653** | 0.9518 / 0.9350 / 0.0697 | 0.9188 / 0.9052 / 0.0966 | 0.6980 / 0.6802 / 0.1694 | 0.8924 / 0.8755 / 0.0948 | 0.9186 / 0.8876 / 0.0854 *(`JC`)* |
| **6-way binary (`K0/K1`, `n = 1,635`)** | **0.9540 / 0.9365 / 0.0693** | 0.9492 / 0.9330 / 0.0734 | 0.9198 / 0.9064 / 0.0976 | 0.7403 / 0.7297 / 0.1560 | 0.9431 / 0.9259 / 0.0718 | 0.9499 / 0.9266 / 0.0719 *(`JC`)* |
| **6-way complex (`K_n ≥ 2`, `n = 2,037`)** | **0.9604 / 0.9445 / 0.0621** | 0.9551 / 0.9376 / 0.0668 | 0.9181 / 0.9039 / 0.0958 | 0.6643 / 0.6390 / 0.1802 | 0.8569 / 0.8379 / 0.1133 | 0.8940 / 0.8549 / 0.0963 *(`JC`)* |
| **MAJIQ v3 native SJ LSVs** (`n = 12,899`)² | **0.9602 / 0.9445 / 0.0645** *(`.jnc`)* | 0.9458 / 0.9289 / 0.0747 *(`.jnc`)* | — | — | **0.9603 / 0.9455 / 0.0632** | — |
| **Gene / Isoform MedARD vs. `Salmon 2`** (`TPM ≥ 1`) | **0.0410 / 0.0678** *(3.2× / 2.1× closer)* | 0.1305 / 0.1399 | 0.1402 / 0.1557 | *(reference)* | *(no TPM output)* | *(no TPM output)* |

*¹ `STAR + MAJIQ v3`: `357.0 s` STAR + `13.9 s` MAJIQ v3. ² Evaluated on MAJIQ v3's own 6,194 non-constitutive annotated splice-junction LSVs (12,899 connections, including 4,343 complex `K ≥ 3`). On SUPPA2's exact binary PacBio truth (`n = 2,339`), Langur cuts MAE by 53% (`0.0734` vs. `0.1563`).*

### 2. Real Human Cohort Null Calibration & Batch Adjustment (Geuvadis, 14 Donors)

Evaluated on 14 human LCL donors (`ERP001942`, 706.9M reads, imbalanced across two sequencing centers `UNIGE` and `ICMB`) across all 10 exact 3v3 same-condition null splits and multivariable batch adjustment (`~ center + sex`):

| Metric (Geuvadis 14-Donor Cohort) | **Langur** *(16T)* | `Whippet2-turbo` *(16T)* | Legacy `Whippet.jl` *(1T)*¹ |
|---|---:|---:|---:|
| **14-sample `quant`** (706.9M reads) | **5m 47s** (`2.03M/s`) · **1.39 GB** | 7m 09s · 9.88 GB | 3h 06m serial · 6.30 GB |
| **Mapped pairs** *(of 353.4M)* | **86.05%** *(+43.8M vs. W1/W2)* | 73.66% | 73.66% |
| **10 null splits — mean / median FPs** | **0.2 / 0.0** *(0 in 9/10 splits)* | 0.8 / 0.0 *(4× more)* | **1,508.9 / 1,402.0** *(92% `TS`/`TE`)* |
| **10 null splits — FPR `@p≤0.05`** (All / Spliced) | **1.38% / 0.70%** | 1.44% / 0.83% | 1.11% / 0.12% *(call rate; 120.1 spliced FPs)* |
| **14-sample `delta` (`~ center + sex`)** | **7.0 s** · **177 MB** | 46.5 s · 2.86 GB | 104.9 s · 1.23 GB *(unadjusted only)* |
| **Batch adjustment (`~ center + sex`)** | **-15.7%** spliced $\hat\phi$ · **-35%** center bias · **2×** `chrX/Y` hits | -13.9% spliced $\hat\phi$ · 1.09× `chrX/Y` enrich. | 90 unadj. calls *(97% `TS`/`TE`, 0 spliced `chrX/Y`)* |

*¹ Legacy `Whippet.jl` crashed on raw SRA FASTQ headers and required line-3 header stripping; `Whippet2-turbo` and `Langur` parse raw SRA FASTQs natively.*

### 3. Real Reverse-Stranded PTBP1 Knockdown & eCLIP Validation (ENCODE K562, 5-Way)

Evaluated on RF-stranded ENCODE `PTBP1` shRNA knockdown vs. control in K562 (`ENCSR527IVX` vs. `ENCSR129RWD`, 2 vs. 2 replicates, 255.3M reads) using orthogonal **ENCODE K562 `PTBP1` eCLIP peaks** (`ENCFF907HNN`) and upstream CU-rich tetramer enrichment:

| Pipeline (`PTBP1` KD, 2 vs. 2) | End-to-End Time | Peak Quant RAM / Disk | Default `ΔΨ>0` CU Motif *(Fold)* | Top-100 CU Motif *(Fold)* | Top-100 eCLIP Up / Flank² | Top-100 Consensus (`n=55`) | `MAP4K4` 9 bp / `PTBP2` ex10 |
|---|---:|---:|---:|---:|---:|---:|---|
| **`Langur`** (`--stranded-reverse`) | **3m 17s** | **1.35 GB / 1.1 GB (0 BAM)** | **5.69 (2.22×)** | **4.53 (1.77×)** | **6.0% / 8.0%** | **53/55 (96.4%)** | **`#1` (`+0.74`) / `#30` (`+0.71`)** |
| `Whippet2-turbo` *(unstranded¹)* | 4m 06s | 8.80 GB / 0.7 GB (0 BAM) | 5.66 (2.21×) | 4.62 (1.81×) | 6.0% / 7.0% | 53/55 (96.4%) | `#2` / `#38` *(201 phantom genes)* |
| `STAR + MAJIQ v3` (`REVERSE`) | 23m 11s | 31.64 GB / 45.6 GB | 4.37 (1.74×) | 4.21 (1.68×) | 5.0% / 7.0% | 48/55 (87.3%) | `#9` (`+0.71`) / *filtered³* |
| `STAR + rMATS-turbo` (`fr-first`) | 24m 56s | 31.64 GB / 45.6 GB | 3.55 (1.39×) | 4.24 (1.66×) | 5.0% / 10.0% | 45/55 (81.8%) | `#7` (`+0.76`) / `#92` (`+0.68`) |
| `Salmon 2 + SUPPA2` (`-l ISR`) | 6m 44s | 12.13 GB / 8.5 GB (0 BAM) | 3.28 (1.22×) | 2.78 (1.04× ns) | 1.0% / 3.0% ns | 3/55 (5.5%) | `#272` (`+0.37`) / *missed* |

*¹ `Whippet2-turbo` lacks `--stranded-reverse`, losing 12.2M pairs and creating 201 phantom antisense-bleed genes; Langur's `--stranded-reverse` separates correct vs. wrong strand 18.7-fold. ² Background upstream eCLIP overlap is 0.3% (CSG), 0.6% (MAJIQ LSVs), and 0.4% (rMATS/SUPPA2); CU-motif background is ~2.55/100 bp across all tools. ³ `MAJIQ v3` filtered `PTBP2` exon 10 because NMD reduced control reads below its 10-read threshold; `MAJIQ heterogen` (designed for larger cohorts) called 0 events at `BH FDR ≤ 0.05` on this 2-versus-2 design.*

### Additional Production Checks
- **3–27 bp microexons:** Unbiased across standard microexons (`90/90` captured, `-0.0024` mean bias under boundary SNPs) and 7–33× cleaner than `Whippet2-turbo` / 3.4–5.0× cleaner than `STAR + rMATS-turbo` on adversarial competing `AD`/`AA` loci.
- **Single-end, trimmed FASTQs & `index --bam`:** Single-end WTC11 (`45.8M` reads) finishes in `76.8 s` at `0.97 GB` RAM (`r = 0.9234` vs. PacBio); variable-length `10–100 bp` trimmed FASTQs run with zero numerical failures; real ENCODE `index --bam` augmentation imports `615` novel splice nodes in `19.9 s`.
- **Cross-platform CI:** Verified on Linux (`x86_64`) and macOS (`arm64`) with `42/42` tests and zero Clippy warnings (`-D warnings`).

---

## AI Agent Skill (Claude, Codex, Gemini)

This repository includes a portable Agent Skill at [`.agents/skills/langur/`](.agents/skills/langur/SKILL.md) (exact CLI flags, strandness rules, `--batch` manifests, `--design` GLM formulas, and output schemas) compatible with Claude Code, OpenAI Codex, Gemini CLI / Antigravity, and Cursor:

```bash
# Auto-discovered inside this repo by Antigravity / Gemini, Codex, and Cursor.
# To install globally across all projects:
mkdir -p ~/.agents/skills ~/.claude/skills ~/.gemini/config/skills
ln -s "$(pwd)/.agents/skills/langur" ~/.agents/skills/langur
ln -s "$(pwd)/.agents/skills/langur" ~/.claude/skills/langur
ln -s "$(pwd)/.agents/skills/langur" ~/.gemini/config/skills/langur
```

---

## Attribution, Citation & License

- **Credit & Lineage:** Independent Rust re-architecture building on the Contiguous Splice Graph representation and event taxonomy pioneered by **Tim Sterne-Weiler** and colleagues in **Whippet** (*Molecular Cell*, 2018), alongside independent implementations of published statistical methods from Salmon, edgeR, DESeq2, and MAJIQ. For a modernized Julia 1.9 Whippet2 implementation (which requires a Julia 1.9 `.jls` index rather than legacy Whippet v1.6 `.jls` files), see **[Whippet2-turbo](https://github.com/koensci/Whippet2.jl)** ([DOI: 10.5281/zenodo.24519354](https://doi.org/10.5281/zenodo.24519354)).
- **AI-Assisted Development:** Developed with AI-assisted design, coding, and debugging (**Google Gemini 4 Neon via Google Antigravity**), with software behavior validated empirically via automated unit and parity tests and orthogonal long-read, eCLIP, and real-cohort benchmarks.
- **Archive / Citation:** *Langur* (`v1.0.0`), Zenodo: https://doi.org/10.5281/zenodo.24615835.
- **License:** MIT — see [LICENSE](LICENSE).
