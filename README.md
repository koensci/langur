# Langur

[![CI](https://github.com/koensci/langur/actions/workflows/ci.yml/badge.svg)](https://github.com/koensci/langur/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Fast, statistically calibrated Contiguous Splice Graph (CSG) quantification and differential splicing for bulk short-read RNA-seq in Rust.**

`langur` (`index`, `inspect`, `quant`, `delta`) profiles both binary ($K_1$) and complex ($K_n \ge 2$) alternative splicing—plus gene and isoform expression—directly from FASTQ files without genome alignment or BAM generation.

- **Zero-copy `mmap` index (`.lidx`, `< 1 ms` load):** Canonical minimizers, annotated junction-bridge seeds, and selective genomic decoys.
- **Hybrid spliced aligner:** Two-tier repetitive minimizer screening, fast ungapped path with zero-allocation stack-banded Gotoh gapped fallback, 1-SNP junction tolerance, and strict-inequality decoy filtering (`decoy > CSG`).
- **Two-stage quantification (`quant`):** Joint Connected-Component Variational Bayes EM (`CC-VBEM`) with unified 5'-hexamer + 20-bin GC bias correction (`--biascorrect`), microexon-safe effective lengths, PCR-stack soft-capping, and collapsed Gibbs inferential variance (`Inf_Var`, col. 15 of `.psi.gz`).
- **Covariate-aware Beta-Binomial GLM (`delta`):** Two-group (`-a`/`-b`) and multivariable DESeq2-style formulas (`--design`/`--formula`/`--test`) with constitutive-event exclusion and Spliced vs. Tandem UTR (`TS`/`TE`) stratified Cox–Reid Empirical Bayes dispersion shrinkage.
- **100% bit-for-bit deterministic** across thread counts and standalone vs. `--batch` runs.

Full algorithmic derivations, ablations, and benchmarks: **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)**.

---

## Build & Test

Requires Rust Edition 2021 (`cargo >= 1.75`). Optional: `igzip` (Intel ISA-L) accelerates `.fq.gz` streaming (auto-falls back to built-in `zlib-ng`).

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
# Optional de novo junctions: add --bam cohort_merged.bam --bam-min-reads 3
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

### Real Human RNA-seq vs. PacBio HiFi Ground Truth (LRGASP WTC11, GRCh38, 6-Way)

Real human WTC11 Illumina paired-end RNA-seq (ENCODE `ENCSR673UKZ` Rep 1: `45.78M` pairs = `91.57M` reads, `--biascorrect`, 16T) on full GENCODE v38 basic (`60,649` genes, `111,832` transcripts) vs. pooled 3-replicate matched PacBio HiFi Iso-Seq (ENCODE `ENCSR507JOF`; `LR ≥ 20`, `SR ≥ 10`; PSI accuracy reported as `r / ρ / MAE`):

| Metric (GRCh38 / WTC11) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(16T)* | `Whippet.jl` v1.6.2 *(1T)* | `Salmon 2.12.0` + `SUPPA2 2.4` *(16T)* | `STAR 2.7.11b` + `MAJIQ v3.0.25` *(16T)* | `STAR 2.7.11b` + `rMATS-turbo v4.3.0` *(16T)* |
|---|---:|---:|---:|---:|---:|---:|
| **Index / Setup** *(uncompressed FASTA)* | **14.8 s** · **5.19 GB** · **1.06 GB** | 28.6 s · 5.33 GB · 691 MB | 753.8 s · 9.87 GB · 634 MB | 110.8 s · 17.93 GB · **8.51 GB** | 700.7 s · **35.09 GB** · **29.5 GB** | 683.1 s · 35.09 GB · **29.5 GB** |
| **End-to-End** *(91.6M reads)* | **72.7 s** (`1.26M/s`) · **1.35 GB** · **0 GB BAM** | 82.3 s · 9.10 GB · 0 GB | 1,467.4 s · 7.00 GB · 0 GB | 190.5 s · **12.24 GB** · 0 GB | 370.9 s¹ · **31.64 GB** · **5.1 GB BAM** | 526.1 s · **31.64 GB** · **5.1 GB BAM** |
| **Mapped Pairs** *(of 45.78M)* | **30.00M (65.5%)** | 24.70M (53.9%) | 24.70M (53.9%) | **30.44M (66.5%)** | 30.55M (66.7%) unique + 1.87M multi | 30.55M (66.7%) unique + 1.87M multi |
| **3-Way CSG Polymorphic** (`n = 13,359`, 6 motifs) | **0.9338 / 0.8715 / 0.0789** | 0.9316 / 0.8689 / 0.0821 | 0.9064 / 0.8473 / 0.1023 | — | *(LSV junctions + IR; no TPMs)* | *(no `AF`/`AL`/`TS`/`TE` or TPMs)* |
| **6-Way Shared Polymorphic** (`n = 3,672`) | **0.9568 / 0.9402 / 0.0653** | 0.9518 / 0.9350 / 0.0697 | 0.9188 / 0.9052 / 0.0966 | 0.6980 / 0.6802 / 0.1694 | 0.8924 / 0.8755 / 0.0948 | 0.9186 / 0.8876 / 0.0854 *(`JC`)* |
| **6-Way Binary (`K0/K1`, `n = 1,635`)** | **0.9540 / 0.9365 / 0.0693** | 0.9492 / 0.9330 / 0.0734 | 0.9198 / 0.9064 / 0.0976 | 0.7403 / 0.7297 / 0.1560 | 0.9431 / 0.9259 / 0.0718 | 0.9499 / 0.9266 / 0.0719 *(`JC`)* |
| **6-Way Complex (`K_n ≥ 2`, `n = 2,037`)** | **0.9604 / 0.9445 / 0.0621** | 0.9551 / 0.9376 / 0.0668 | 0.9181 / 0.9039 / 0.0958 | 0.6643 / 0.6390 / 0.1802 | 0.8569 / 0.8379 / 0.1133 | 0.8940 / 0.8549 / 0.0963 *(`JC`)* |
| **MAJIQ v3 Native SJ LSVs** (`n = 12,899` conns)² | **0.9602 / 0.9445 / 0.0645** *(`.jnc`)* | 0.9458 / 0.9289 / 0.0747 *(`.jnc`)* | — | — | **0.9603 / 0.9455 / 0.0632** | — |
| **Gene / Isoform MedARD vs. `Salmon 2`** (`TPM ≥ 1`) | **0.0410 / 0.0678** *(3.2× / 2.1× closer)* | 0.1305 / 0.1399 | 0.1402 / 0.1557 | *(reference)* | *(no TPM output)* | *(no TPM output)* |

*Memory figures report macOS `/usr/bin/time -l` peak memory footprint. ¹ `STAR + MAJIQ v3.0.25`: `357.0 s` STAR + `13.9 s` MAJIQ v3 (`sj 8.88 s · build 1.39 s · psi-coverage 1.66 s · quantify 1.94 s`, `653 MB` peak, plus `~10 s` `.bai`). ² Evaluated directly on MAJIQ v3's own `6,194` polymorphic annotated splice-junction LSVs (`12,899` connections, incl. `4,343` complex `K ≥ 3`: Langur `.jnc.gz` `0.9525 / 0.0661` vs. `STAR + MAJIQ v3` `0.9538 / 0.0651`). Even under `SUPPA2`'s own binary PacBio truth (`n = 2,339`), Langur scores `0.9403 / 0.0734` vs. `Salmon 2 + SUPPA2` `0.7393 / 0.1563` (-53% MAE).*

### Real Human Cohort Null Calibration & Batch Adjustment (Geuvadis `ERP001942`, GRCh38)

14 real human CEU LCL donors (`UNIGE` `6F + 2M`, `ICMB` `2F + 4M`; `706.9M` 75 bp reads, 16T, `--biascorrect`) — all 10 exact 3v3 same-condition null splits (`6F UNIGE`, `-r 5 -s 2`) and imbalanced multivariable batch adjustment (`~ center + sex`):

| Metric (Geuvadis 14-Donor Cohort, GRCh38) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(16T)* | Original `Whippet.jl` v1.6.2 *(1T, `+`-stripped)*¹ |
|---|---:|---:|---:|
| **14-Sample `quant`** (706.9M reads) | **347.3 s** (`2.03M/s`) · **1.39 GB** | 429.2 s (`1.65M/s`) · 9.88 GB | **11,145.4 s** (`3h 06m` serial, `63.4k/s`) · 6.30 GB |
| **Mapped Pairs** *(of 353.4M)* | **86.05%** (`304.2M`; **+43.8M** vs. W1/W2) | 73.66% (`260.3M`) | 73.66% (`260.3M`, exact W2 parity) |
| **10 Null Splits — Mean / Med FPs** (`FDR ≤ 0.05` or `P ≥ 0.95`) | **0.2 / 0.0** *(9/10 splits = 0 FPs)* | 0.8 / 0.0 *(4× more FPs)* | **1,508.9 / 1,402.0** *(1,202.9 @ \|ΔΨ\| ≥ 0.1; 92% `TS`/`TE`)* |
| **10 Null Splits — FPR `@p≤0.05`** (All / Spliced) | **1.38% / 0.70%** | 1.44% / 0.83% | 1.11% / 0.12% *(call rate; `120.1` spliced FPs)* |
| **14-Sample `delta` (`~ center + sex`)** | **7.02 s** · **177 MB** | 46.51 s · 2.86 GB | 104.85 s · 1.23 GB *(unadj. only)* |
| **Batch Adjustment (`~ center + sex`)** | **-15.7%** spliced $\hat\phi$ · **-35%** center `\|ΔΨ\|` bias · **2×** `chrX/Y` hits (`1.67×` spliced enrich.) | -13.9% spliced $\hat\phi$ · 1.09× spliced `chrX/Y` enrich. | **90** unadj. calls (`87` `TS`/`TE`, `0.68×` depleted `chrX/Y`, `0` spliced `chrX/Y`) |

*¹ Original `Whippet.jl` v1.6.2 required stripping raw SRA line-3 `+ERR...` headers to bare `+` (`timbitz/Whippet.jl#135`; `Whippet2-turbo` and `Langur` parse raw SRA FASTQs natively). Full stratified tables, sensitivity checks, `Whippet2-turbo` commit details, and Geuvadis `TS/TE` corner-interaction notes: **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)**.*

### Real RF-Stranded RBP Knockdown & Production Checks

- **ENCODE `PTBP1` shRNA KD (`K562`, `ENCSR527IVX` vs. `ENCSR129RWD`, `255.3M` RF-stranded reads):** `--stranded-reverse` (`73.1%`) collapses **18.7×** under wrong-strand `--stranded` (`3.9%`) and eliminates opposite-strand antisense bleed (`TUBA1C`, `UBB`, `RBBP4`, `PIEZO1`, `FDPS`). Full 4-sample `--batch` runs in **`2m 59s` (`178.7 s`, `1.43M/s`) at `1.18 GB` RAM**; `delta` runs in **`1.39 s` (`129 MB`)**, calling **`58/58` Spliced hits (`0` spurious `TS`/`TE`)** at `FDR ≤ 0.05` while recovering `PTBP1` `3.2×` knockdown, `PTBP2` `8.2×` expression rescue + exon-10 NMD switch (`ΔΨ = +0.714`), **`MAP4K4` 9 bp microexon as the `#1` genome-wide hit (`ΔΨ = +0.742`, `FDR = 8.64e-6`)**, and canonical targets `SMAP1`, `GABBR1`, `RTN4`, `CTNND1`, `TPM1`, and `PICALM`.
- **Single-end real human WTC11 (`45.8M` R1 reads, GRCh38):** `76.8 s`, **`0.97 GB` peak RAM**, `r = 0.9234 / MAE = 0.0875` vs. PacBio HiFi (`r = 0.9799` SE-vs-PE).
- **Variable-length trimmed FASTQs & real `index --bam`:** Zero NaNs/Infs on `10–100 bp` + all-`'N'` + lowercase FASTQs (`r = 0.9929` PE); full GRCh38 `index --bam` with real ENCODE `PTBP1` BAM (`19.94 s`, `+615` BAM-derived nodes, `389` supported PSI rows, `100` `ALT_UNIQ` novel junctions).
- **CI:** Verified on Linux (`x86_64`) and macOS (`arm64`) with `42/42` tests and zero Clippy warnings (`-D warnings`).

---

## Attribution, Transparency & License

- **Credit & Lineage:** Independent Rust re-architecture building on the Contiguous Splice Graph (CSG) representation and motif taxonomy pioneered by **Tim Sterne-Weiler** and colleagues in **Whippet** (*Molecular Cell*, 2018; `Whippet.jl`), alongside independent implementations of published statistical concepts from `Salmon` and `edgeR`/`DESeq2`, plus overhang-diversity QC for PCR junction-stack artifacts highlighted in `MAJIQ`. Users needing 100% drop-in compatibility with legacy Julia 1.9 `.jls` indexes can use [Whippet2-turbo](https://github.com/koensci/Whippet2.jl).
- **AI-Assisted Development:** Developed iteratively with AI-assisted coding and review (**Google Gemini**), with all biological invariants and benchmarks verified via automated unit, parity, and full-human test suites.
- **License:** MIT — see [LICENSE](LICENSE).
