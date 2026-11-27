# Langur

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

### Simulated Human Cohort (`E-MTAB-3766`, GRCh37.71, 480M Reads)

3 vs. 3 paired-end 101 bp replicates (`20,410` protein-coding genes, `145,342` transcripts, `1,000` truly DS genes; 16T, `--biascorrect`, Apple Silicon arm64):

| Stage | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | `Whippet.jl` v1.6.2 *(1T)* | Langur Advantage |
|---|---:|---:|---:|---|
| **Index Build** *(uncompressed FASTA)* | **14.6 s** · 4.36 GB *(incl. 1.10 Gb decoys)* | 23.0 s · 4.21 GB | 287.1 s · 8.10 GB | **19.7× vs. W1 · 1.58× vs. W2** |
| **Index Load** | **0.03 ms** (`mmap`) | 6.20 s | 7.33 s | **>200,000× faster** |
| **Sample 1 `quant`** (80M reads) | **40.5 s** (`1.98M/s`) · **1.40 GB** | 81.1 s (`987k/s`) · 9.62 GB | 1,608.4 s (`49.7k/s`) · 7.19 GB | **39.7× vs. W1 · 2.00× vs. W2 · 6.9× less RAM** |
| **6-Sample `--batch`** (480M reads) | **248.9 s** (`1.93M/s`) · **1.42 GB** | 433.8 s (`1.11M/s`) · 9.83 GB | *9,650.6 s* *(6× extrap.)* · 7.19 GB | **38.8× vs. W1 · 1.74× vs. W2 · 6.9× less RAM** |
| **Gene / Isoform MedARD** *(Median ARD vs. RSEM)* | **0.0137 / 0.0402** | 0.0443 / 0.0687 | 0.0753 / 0.1145 | **5.5× / 2.85× lower error vs. W1** |
| **`delta` GLM** (223k events, `-r 5 -s 2`) | **2.50 s** · **191 MB** | 27.38 s · 1.87 GB | 85.51 s *(`--legacy`)* · 609 MB | **34.2× vs. legacy · 11.0× & 9.8× RAM vs. W2** |
| **Differential Splicing** (`TPR / FDR @0.05` · AUPRC) | **0.336 / 5.4%** *(0.338 / 2.0% spliced)* · **0.565** | 0.306 / 5.3% *(0.309 / 1.3% spliced)* · 0.551 | 0.744 / 81.5% · 0.352 | **15× lower FDR vs. legacy · 10.6× lower FDR vs. `rMATS`** |

### Real Human RNA-seq vs. PacBio HiFi Ground Truth (LRGASP WTC11, GRCh38, 5-Way)

Real human WTC11 Illumina paired-end RNA-seq (ENCODE `ENCSR673UKZ` Rep 1: `45.78M` pairs = `91.57M` reads, `--biascorrect`, 16T) on full GENCODE v38 basic (`60,649` genes, `111,832` transcripts) vs. pooled 3-replicate matched PacBio HiFi Iso-Seq (ENCODE `ENCSR507JOF`; `LR ≥ 20`, `SR ≥ 10`; PSI accuracy reported as `r / ρ / MAE`):

| Metric (GRCh38 / WTC11) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(16T)* | `Whippet.jl` v1.6.2 *(1T)* | `Salmon 2.12.0` + `SUPPA2 2.4` *(16T)* | `STAR 2.7.11b` + `rMATS-turbo v4.3.0` *(16T)* |
|---|---:|---:|---:|---:|---:|
| **Index / Setup** *(uncompressed FASTA)* | **14.8 s** · **5.19 GB** · **1.06 GB** | 28.6 s · 5.33 GB · 691 MB | 753.8 s · 9.87 GB · 634 MB | 110.8 s · 17.93 GB · **8.51 GB** | 683.1 s · 35.09 GB · **29.5 GB** |
| **End-to-End** *(91.6M reads)* | **72.7 s** (`1.26M/s`) · **1.35 GB** · **0 GB BAM** | 82.3 s · 9.10 GB · 0 GB | 1,467.4 s · 7.00 GB · 0 GB | 190.5 s · **12.24 GB** · 0 GB | 526.1 s · **31.64 GB** · **5.1 GB BAM** |
| **Mapped Pairs** *(of 45.78M)* | **30.00M (65.5%)** | 24.70M (53.9%) | 24.70M (53.9%) | **30.44M (66.5%)** | 30.55M (66.7%) unique + 1.87M multi |
| **3-Way CSG Polymorphic** (`n = 13,359`, 6 motifs) | **0.9338 / 0.8715 / 0.0789** | 0.9316 / 0.8689 / 0.0821 | 0.9064 / 0.8473 / 0.1023 | — | *(no `AF`/`AL`/`TS`/`TE` or TPMs)* |
| **5-Way Shared Polymorphic** (`n = 4,693`) | **0.9507 / 0.9350 / 0.0688** | 0.9456 / 0.9291 / 0.0728 | 0.9105 / 0.8974 / 0.1012 | 0.6832 / 0.6688 / 0.1747 | 0.9099 / 0.8822 / 0.0896 *(`JC`)* |
| **5-Way Binary (`K0/K1`, `n = 2,008`)** | **0.9467 / 0.9292 / 0.0740** | 0.9424 / 0.9245 / 0.0778 | 0.9126 / 0.8994 / 0.1032 | 0.7272 / 0.7166 / 0.1599 | 0.9414 / 0.9188 / 0.0767 *(`JC`)* |
| **5-Way Complex (`K_n ≥ 2`, `n = 2,685`)** | **0.9551 / 0.9413 / 0.0648** | 0.9493 / 0.9340 / 0.0692 | 0.9090 / 0.8959 / 0.0996 | 0.6510 / 0.6327 / 0.1859 | 0.8863 / 0.8533 / 0.0993 *(`JC`)* |
| **Gene / Isoform MedARD vs. `Salmon 2`** (`TPM ≥ 1`) | **0.0410 / 0.0678** *(3.2× / 2.1× closer)* | 0.1305 / 0.1399 | 0.1402 / 0.1557 | *(reference)* | *(no TPM output)* |

*Memory figures report macOS `/usr/bin/time -l` peak memory footprint. Even under `SUPPA2`'s own exact-event PacBio truth on binary `K0/K1` events (`n = 2,339`), Langur scores `0.9403 / 0.9232 / 0.0734` vs. `Salmon 2.12.0 + SUPPA2 2.4` `0.7393 / 0.7299 / 0.1563` (-53% MAE).*

### Real Human Cohort Null Calibration & Batch Adjustment (Geuvadis `ERP001942`, GRCh38)

14 real human CEU LCL donors (`UNIGE` `6F + 2M`, `ICMB` `2F + 4M`; `706.9M` 75 bp reads, 16T, `--biascorrect`) — all 10 exact 3v3 same-condition null splits (`6F UNIGE`, `-r 5 -s 2`) and imbalanced multivariable batch adjustment (`~ center + sex`):

| Metric (Geuvadis 14-Donor Cohort, GRCh38) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(16T)* | Legacy `Whippet.jl` Beta MC |
|---|---:|---:|---:|
| **14-Sample `quant --batch`** (706.9M reads) | **347.3 s** (`2.03M/s`) · **1.39 GB** | 429.2 s (`1.65M/s`) · 9.88 GB | *(SRA FASTQ parser issue `#135`)* |
| **Mapped Pairs** *(of 353.4M)* | **86.05%** (`304.2M`; **+43.8M** vs. W2) | 73.66% (`260.3M`) | — |
| **10 Null Splits — Mean / Med FPs** (`FDR ≤ 0.05` or `P ≥ 0.95`) | **0.2 / 0.0** *(9/10 splits = 0 FPs)* | 0.8 / 0.0 *(4× more FPs)* | **7,250.2 / 7,145.0** *(4,080 @ \|ΔΨ\| ≥ 0.1)* |
| **10 Null Splits — FPR `@p≤0.05`** (All / Spliced) | **1.38% / 0.70%** | 1.44% / 0.83% | 4.44% / 3.31% *(call rate)* |
| **14-Sample `delta` (`~ center + sex`)** | **7.02 s** · **177 MB** | 46.51 s · 2.86 GB | 86.94 s *(unadj. only; `--max-variance 0.05`)* |
| **Batch Adjustment (`~ center + sex`)** | **-15.7%** spliced $\hat\phi$ · **-35%** center `\|ΔΨ\|` bias · **2×** `chrX/Y` hits (`1.67×` spliced enrich.) | -13.9% spliced $\hat\phi$ · 1.09× spliced `chrX/Y` enrich. | **3,276** unadj. calls (`1.04×` null `chrX/Y` enrich.) |

*Full stratified tables, sensitivity checks, `Whippet2-turbo` commit details, and Geuvadis `TS/TE` corner-interaction notes: **[docs/WHITEPAPER.md](docs/WHITEPAPER.md)**.*

---

## Attribution, Transparency & License

- **Credit & Lineage:** Independent Rust re-architecture building on the Contiguous Splice Graph (CSG) representation and motif taxonomy pioneered by **Tim Sterne-Weiler** and colleagues in **Whippet** (*Molecular Cell*, 2018; `Whippet.jl`), alongside independent implementations of published statistical concepts from `Salmon` and `edgeR`/`DESeq2`, plus overhang-diversity QC for PCR junction-stack artifacts highlighted in `MAJIQ`.
- **AI-Assisted Development:** Developed iteratively with AI-assisted coding and review (**Google Gemini**), with all biological invariants and benchmarks verified via automated unit, parity, and full-human test suites.
- **License:** MIT — see [LICENSE](LICENSE).
