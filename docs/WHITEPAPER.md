# Langur: Fast, Statistically Calibrated Splice-Graph Quantification and Differential Splicing in Rust

*CLI binary: `langur` · Index format: `.lidx`*

---

## Abstract

Alignment-free quantification on **Contiguous Splice Graphs (CSGs)**, pioneered by Whippet (Sterne-Weiler et al., *Molecular Cell* 2018), profiles both binary ($K_1$) and complex ($K_n \ge 2$) alternative splicing directly from FASTQ files without genome alignment or BAM generation. However, scaling legacy CSG implementations to modern human cohorts exposes four bottlenecks: ungapped matching drops reads near indels and splice-proximal SNPs; transcriptome-only indexing allows diverged pseudogene and intronic reads to force-map onto parent genes; per-gene EM and naive Beta variance understate paralog and multi-path ambiguity; and two-group pooled Beta Monte Carlo testing (`whippet-delta`) cannot adjust for covariates and suffers severe false-discovery inflation under biological overdispersion.

**Langur** is a ground-up Rust re-architecture (`index`, `inspect`, `quant`, `delta`) combining a zero-copy `mmap` splice-graph + selective genomic decoy index (`< 1 ms` load), hybrid ungapped + zero-allocation banded Gotoh spliced alignment with 1-SNP junction tolerance, joint Connected-Component Variational Bayes EM (`CC-VBEM`) with unified 5'-hexamer/GC bias correction, targeted collapsed Gibbs inferential variance (`Inf_Var`), and a covariate-aware Beta-Binomial GLM (`delta`) with Spliced vs. Tandem UTR stratified Cox–Reid Empirical Bayes dispersion shrinkage. On the 480M-read human `E-MTAB-3766` cohort (GRCh37.71, 16 threads, `--biascorrect`), Langur builds the index in **14.6 s**, quantifies **480M reads in 4m 09s (`248.9 s`, `1.93M reads/s`)** at **1.4 GB RAM** (**38.8×** faster than `Whippet.jl`, **1.74×** faster and **6.9×** less RAM than `Whippet2-turbo`), cuts gene/isoform MedARD **5.5× / 2.85×**, and achieves **5.4% realized FDR** (`2.0%` spliced; `0.565 / 0.567` AUPRC) in **2.5 s**—a **15× FDR reduction** over legacy Whippet (`81.5%`) and **10.6× lower FDR** than published `rMATS` (`57.1%`). On real human WTC11 RNA-seq (LRGASP `ENCSR673UKZ`, GRCh38 GENCODE v38, 60,649 genes, 91.6M reads), Langur builds the index in **14.8 s** (**46.2×** faster and **28×** smaller than `STAR 2.7.11b`) and quantifies WTC11 in **72.7 s** at **1.35 GB RAM** (**2.6×** faster and **9.1×** less RAM than Rust `Salmon 2.12.0 + SUPPA2 2.4`; **7.24×** faster, **23.5×** less RAM, and zero BAMs vs. `STAR 2.7.11b + rMATS-turbo v4.3.0`), sweeping **27/27** polymorphic CSG strata against matched PacBio HiFi Iso-Seq ($r = 0.9338$, $\text{MAE} = 0.0789$), beating `Salmon 2 + SUPPA2` across every stratum even on SUPPA2's own truth (binary $K_0/K_1$: $r = 0.9403$ vs. $0.7393$, **-53% MAE**), and matching `STAR + rMATS-turbo` on binary splicing while outperforming it on full-length complex ($K_n \ge 2$) exon inclusion ($r = 0.9551$ vs. $0.8863$, **-35% MAE**). On 14 real human CEU LCL donors from Geuvadis (`ERP001942`, 706.9M reads, GRCh38), Langur quantifies the cohort in **5m 47s (`347.3 s`, `2.03M reads/s`)** at **1.39 GB RAM** (**1.24×** faster, **7.1×** less RAM, and **+43.8M** more mapped pairs than `Whippet2-turbo`), yields **0.2 mean / 0.0 median false positives** across all 10 exact 3v3 same-condition null splits (**>20,000×** fewer than legacy Beta MC and **4×** fewer than `Whippet2-turbo`), and cuts residual spliced overdispersion **15.7%** while doubling `chrX`/`chrY` signal recovery under multivariable batch adjustment (`~ center + sex`) in **7.0 s** at **177 MB RAM**.

---

## 1. Motivation, Lineage & Scope

Bulk RNA-seq splicing studies typically trade off between two paradigms:
1. **Genome-aligned callers (`STAR` + `MAJIQ`, `rMATS`, `LeafCutter`):** Strong local junction accuracy and de novo discovery, at the cost of heavy genome indices (~30 GB), high RAM (>30 GB), large coordinate-sorted BAMs, and binary-only event definitions in `rMATS` (omitting `AF`/`AL`, complex $K_n \ge 2$ hubs, and expression TPMs).
2. **Alignment-free tools (`Whippet.jl`, `Salmon` + `SUPPA2`):** Eliminate BAM generation, but `SUPPA2` derives event PSI solely from full-length isoform TPMs—propagating isoform deconvolution error into local splice sites—while legacy `Whippet.jl` suffers from ungapped matching, pseudogene force-mapping, uncalibrated pooled Beta MC testing, and single-threaded Julia overhead.

**Attribution, Scope & Transparency:**
- **Attribution:** Langur is an independent Rust re-architecture building on the CSG representation and motif taxonomy introduced by Tim Sterne-Weiler and colleagues (*Molecular Cell*, 2018), alongside independent implementations of published statistical concepts from `Salmon` (selective decoys, VBEM, Gibbs uncertainty) and `edgeR`/`DESeq2` (Cox–Reid Empirical Bayes GLMs), plus overhang-diversity QC for PCR junction-stack artifacts highlighted in `MAJIQ`.
- **Scope:** Purpose-built for **bulk short-read RNA-seq** (paired-end or single-end Illumina), with optional BAM-derived novel splice-site import at index time (`index --bam`). For single-cell or full-length long-read (PacBio HiFi / ONT) isoform analysis, use dedicated frameworks such as **Isosceles** (Kabza, Sterne-Weiler et al., *Nature Communications* 2024), **Oarfish**, or **IsoQuant**.
- **License & Transparency:** MIT License. Developed iteratively with AI-assisted coding and review (Google Gemini), with every biological invariant and benchmark audited via automated unit, parity, and full-human test suites.

---

## 2. Pipeline Architecture

```mermaid
flowchart LR
    A["Genome FASTA + GTF (+ Optional BAM)"] --> B["index: CSG + Minimizers + Selective Decoys"]
    B -->|"mmap (< 1 ms)"| C["quant (Stage 1): Hybrid Aligner + Decoy Filter + CC-VBEM"]
    D["FASTQ (.fq / .fq.gz)"] --> C
    C --> T[".gene.tpm.gz · .isoform.tpm.gz"]
    C -->|"Routed Counts"| E["quant (Stage 2): Local CSG PSI + Stack Cap + Gibbs Inf_Var"]
    E --> P[".psi.gz · .jnc.gz · .map.gz"]
    P --> F["delta: Stratified Cox-Reid EB Beta-Binomial GLM"]
    G["Design Table / Formula"] --> F
    F --> H[".diff.gz"]
```

| Bottleneck | Legacy Whippet (`Whippet.jl` / `Whippet2`) | Langur (Rust) |
|---|---|---|
| **Index loading** | Full heap deserialization every run (`1.7–7.3 s`) | Zero-copy `mmap` binary file (`0.03–0.3 ms`) |
| **Pseudogene / intronic noise** | Transcriptome-only index | Selective genomic decoys with strict-inequality filter (`decoy > CSG`) |
| **Indels & boundary SNPs** | Ungapped matching only | Fast ungapped path + zero-alloc banded Gotoh X-drop + 1-SNP splice jumps |
| **Multi-gene ambiguity** | Equal pre-splitting across genes | Joint Connected-Component Variational Bayes EM across paralog families |
| **Complex $K_n \ge 2$ uncertainty** | Closed-form Beta variance only | Targeted collapsed Gibbs posterior variance (`Inf_Var`, col. 15 of `.psi.gz`) |
| **Differential splicing** | Pooled 2-group Beta MC (no covariates, high FDR) | Covariate-aware Beta-Binomial GLM with Spliced vs. Tandem UTR EB priors |
| **Reproducibility** | Minor MC/thread drift | 100% bit-for-bit identical across thread counts and `--batch` runs |

---

## 3. Zero-Copy Splice-Graph & Selective Decoy Index (`index`)

- **CSG & Motif Parity:** Partitions each gene's overlapping exons (plus optional `--bam` novel junctions) into minimal non-overlapping sub-intervals ordered $5' \to 3'$, preserving Whippet's canonical motifs (`CE`, `AA`, `AD`, `RI`, `AF`, `AL`, optional `BS`, and tandem UTRs `TS`/`TE`).
- **Exon + Junction-Bridge Minimizers:** Indexes canonical minimizers across both exonic nodes and exact sequence bridges spanning every annotated splice junction, enabling single-lookup splice jumping even across microexons.
- **Selective Genomic Decoys:** Following Salmon's selective alignment strategy (Srivastava et al., 2020), scans non-exonic intervals in parallel and retains only genomic windows sharing seeds with exons—capturing processed pseudogenes and paralogous repeats without indexing the full 3 Gb genome.
- **Instant `mmap` Loading:** Writes a 64-byte-aligned binary container mapped directly into virtual memory (`0.03–0.3 ms` load, shared zero-copy across `--batch` cohorts).

<details>
<summary><b>Technical Details: Minimizer Table, Flank Encoding & Binary Format</b></summary>

- **Minimizers:** Canonical odd-$k$ minimizers (default $k = 17$, window $w = 8$, invertible Wang64 hash) in a 2-level radix-bucketed open table (`MinBucketOffsets`, `MinKeyLo`, `MinSlots`, `MinPostings`) with inline single-hit fast-path slots (`MinimizerLookup::Inline`).
- **Annotated junction seeds:** Extracted from the $(k-1)+(k-1)$ sequence across every annotated splice edge `(src_node, dst_node)` with exact `left_overhang` $\in [1, k-1]$.
- **Precomputed node metadata:** 2-bit 16-bp acceptor/donor boundary flanks with explicit valid base counts (`0..=16`) for single-instruction XOR + popcount 1-SNP matching; canonical splice dinucleotide flags (`GT-AG`, `GC-AG`, `AT-AC`); cumulative transcript prefix lengths (`tx_node_prefix_len`) for $O(\log N_t)$ fragment-length queries; and 20-bin expected transcript GC profiles (`tx_expected_gc`).
- **Selective decoys:** Non-CSG hits expanded by `decoy_window_pad = 250 bp`, merged, split at `N`/IUPAC runs (preventing spurious poly-`A` runs in 2-bit encoding), and indexed with `max_decoy_occ = 32`.
- **Container:** 17 `#[repr(C)]` `bytemuck::Pod` sections aligned to 64-byte boundaries with magic header (`LANGUR01`) and `CRC32` checksum.

</details>

---

## 4. Hybrid Spliced Alignment & Decoy Screening (`quant` — Alignment)

- **Two-tier seeding & zero-allocation gapped recovery:** Reads with $\ge 2$ low-copy CSG minimizers ($\le 64$ postings) automatically skip repetitive Alu/poly-A postings (falling back to all postings for high-copy paralog families). Clean reads run a fast SIMD `memcmp` / ungapped scan; genuine small-indel mismatch clusters (complexity-guarded against homopolymer tails) invoke a zero-allocation stack-banded Gotoh affine-gap X-drop aligner (`HALF_BAND = 12`, `X_DROP = 24`), paired with exact score upper-bound candidate pruning.
- **Boundary-proximal SNPs & annotated-first splice jumping:** Retains near-optimal boundary extensions so SNPs 1–3 bp from donor/acceptor sites never block splice jumping; evaluates annotated junctions first with **zero splice penalty** (`min_anno_overhang = 4 bp`), falling back to 1-SNP-tolerant novel splice-flank matching (`min_novel_overhang = 8 bp`, prioritizing canonical `GT-AG`/`GC-AG`/`AT-AC`).
- **Strict-inequality pseudogene immunity:** Perfect CSG hits skip decoy alignment completely. On imperfect hits sharing decoy seeds, decoys are scored under the identical matrix; because annotated junctions incur zero penalty, an intronless processed pseudogene can at most *tie* spliced parent-gene mRNA. Langur filters a read as `DecoyFiltered` **only when** $S_{\text{decoy}} > S_{\text{CSG}}$—discarding diverged pseudogene/intronic reads while preserving true parent-gene reads on exact ties.
- **Strand-aware paired-end concordance:** Supports unstranded, FR (`--stranded`), and RF/TruSeq dUTP (`--stranded-reverse`) libraries with exact spliced transcriptomic fragment lengths across compatible isoforms.

<details>
<summary><b>Technical Details: Scoring Matrix, Deferred Decoys & Streaming I/O</b></summary>

- **Scoring:** Match `+2`, mismatch `-4`, `'N'` `-1`, gap open `-6`, gap extend `-1`; non-canonical novel splice penalty `-2`; optional `--circ` back-splices (`dst <= src`, `min_ov = 10 bp`); boundary extension retained within `12` points of prefix peak.
- **Deferred decoy seeding:** Primary seeding sets only `has_any_decoy_posting`; decoy seeds (bucketed by `(decoy_id, strand, diag_bucket)`) are collected only if `best_csg_score < max_possible_score`, and mate 2 / gapped decoy alignment is skipped whenever exact score upper bounds cannot beat `best_csg_score`.
- **Exact transcriptomic fragment length:**

  $$
  d_t = \left(\text{prefix}_t[i_{3'}] + \text{end}_{3'}\right) - \left(\text{prefix}_t[i_{5'}] + \text{start}_{5'}\right).
  $$

- **I/O:** Auto-detects `igzip -dc` (with exit-status verification and automatic `flate2`/`zlib-ng` fallback) and writes RFC 1952 multi-member gzip outputs in parallel.

</details>

---

## 5. Two-Stage Quantification: Expression, PSI & Inferential Variance (`quant`)

Reads stream in parallel waves into sparse thread-local counters merged in strict file order, ensuring **100% bit-for-bit determinism** across thread counts at **~1.4 GB** peak RAM on human libraries.

### Stage 1: Transcript & Gene Abundance (Global Routing)
1. **Post-pass FLD:** Fits one smoothed sample-wide fragment-length curve after alignment so effective lengths are exact and order-independent.
2. **Joint CC-VBEM:** Groups transcripts sharing multi-mapping reads (paralog families and annotated pseudogenes) into connected components and solves them jointly with sparse Variational Bayes EM, routing ambiguous reads to expressed parents rather than silent paralogs.
3. **Unified 5'-hexamer + 20-bin GC bias (`--biascorrect`):** Fits 5'-hexamer weights deterministically from Wave 0 (~1.6M pairs), solves Pass 1 CC-VBEM on primer-weighted counts to fit 20-bin GC bias, then folds GC weights directly into transcript effective lengths in Pass 2 (avoiding GC double-counting).

### Stage 2: Local Splice-Graph PSI, Artifact Control & Gibbs `Inf_Var`
Routing ambiguous reads via converged Stage 1 posteriors into each gene's CSG:
- **Microexon-safe effective lengths:** Preserves full mappable overhang spans on internal splice junctions (clamping by single-node length only at true 5'/3' ends), preventing PSI inflation on short microexons (e.g., 15 bp `CE`).
- **Overhang-diversity PCR-stack soft-capping:** Tracks `Distinct_Pos` and `Max_Stack` in `.jnc.gz`; on **novel or `RI`-touching junctions only**, pileups of $\ge 10$ reads at one overhang coordinate are softly capped to `2.0` (`--stack-cap`, disable via `--no-stack-filter`).
- **Targeted collapsed Gibbs (`Inf_Var`):** Binary ($K_1$) and unambiguous events (~80–85%) use exact closed-form Beta variance; complex ($K_n \ge 2$) ambiguous events run a deterministic collapsed Gibbs sampler (`~59 µs`/event, group-balanced prior) to report posterior variance (`Inf_Var`, col. 15 of `.psi.gz`, backward-compatible with Whippet 2 cols. 1–14).

<details>
<summary><b>Mathematical Details: Effective Length, CC-VBEM, Unified GC Pass & Gibbs Sampler</b></summary>

- **Post-pass effective length:** From the 5-point Gaussian-smoothed empirical FLD PMF $p(d)$, CDF $F(L)$, and first moment $M_1(L) = \sum_{d=1}^L d\,p(d)$:

  $$
  \ell_t = (L_t + 1)\,F(\min(L_t, D_{\max})) - M_1(\min(L_t, D_{\max})).
  $$

- **CC-VBEM + SQUAREM** (default sparse Dirichlet prior $\alpha_0 = 0.01$; `--use-em` selects MLE; $\psi(\cdot)$ denotes digamma, distinct from splicing $\Psi$):

  $$
  w_t^{(k)} = \frac{\exp\!\left(\psi(\alpha_0 + c_t^{(k)})\right)}{\ell_t}, \qquad c_t^{(k+1)} = c_t^{\text{uniq}} + \sum_{e \ni t} n_e \frac{w_t^{(k)}}{\sum_{u \in e} w_u^{(k)}}.
  $$

- **Pass 2 GC effective-length update** (run on primer-only counts, with $w_{\text{GC}}[b] = \text{expected}[b] / \text{observed}[b] \in [0.15, 6.0]$):

  $$
  \ell_t^{\text{GC}} = \ell_t \cdot \sum_{b=0}^{19} \frac{E_{\text{GC}}[t, b]}{w_{\text{GC}}[b]}.
  $$

- **Targeted collapsed Gibbs (`K_n >= 2`):** Group-balanced Dirichlet prior ($\alpha_p = 1/n_{\text{inc}}$ on inclusion paths, $1/n_{\text{exc}}$ on exclusion paths $\implies$ marginal $\text{Beta}(1,1)$ on $\Psi$), deterministic `Xoshiro256++` seeded from `(seed, gene_id, node_1based)`, `100` draws (`32` burn-in, thin `2`):

  $$
  \theta_p^{(b)} \sim \text{Gamma}\!\left(\alpha_p + z_p^{(b)},\; \text{rate} = \ell_p\right), \qquad \Psi^{(b)} = \frac{\sum_{p \in \text{Inc}} \theta_p^{(b)}}{\sum_q \theta_q^{(b)}}, \qquad \text{InfVar} = \widehat{\text{Var}}(\Psi^{(1..B)}).
  $$

</details>

---

## 6. Covariate-Aware Differential Splicing (`delta`)

`delta` replaces pooled two-group Beta Monte Carlo testing with a multithreaded **Beta-Binomial GLM** featuring Cox–Reid Empirical Bayes dispersion shrinkage—analogous to `edgeR` quasi-likelihood GLMs, tailored for fractional splicing proportions ($\Psi_{ij}, N_{ij}$):

1. **Multivariable designs & collinearity resilience:** Supports two-group (`-a`/`-b`) and DESeq2-style formulas (`--design`, `--formula "~ batch + donor + condition"`, `--test`, continuous covariates, `*`/`:` interactions); per-event pivoted-QR drops only collinear nuisance columns on partial missingness.
2. **Propagating complex-event uncertainty (`Inf_Var`):** Filters on raw `Total_Reads >= --min-reads`, then smoothly down-weights complex $K_n \ge 2$ observations whose Gibbs `Inf_Var` exceeds nominal Beta variance (floor `0.1×`, $K_1$ untouched), tightening dispersion estimates and simultaneously boosting power and FDR calibration.
3. **Excluding constitutive events from dispersion fitting:** ~63% of human CSG nodes are constitutive (`0.0` or `1.0` across samples) and carry zero variance information; including them causes the Cox–Reid penalty to peak near $\phi \approx 0.64$, collapsing the EB prior. Excluding boundary events restores full power and speeds up `delta` ~3× (to **2.5 s**).
4. **Event-class stratified EB priors (Spliced vs. Tandem UTR):** Tandem UTR ends (`TS`/`TE`) exhibit **7.8× higher median overdispersion** (`phi = 0.117`) than canonical splice junctions (`phi = 0.015`). Fitting separate EB priors (`--no-stratify-eb` to opt out) boosts sensitivity **6×–25×** while cutting `TS`/`TE` false positives **>300×**.
5. **Finite-sample $F$-test LRT & sign-consistent effects:** Uses $F(\text{df},\, d_{\text{res}} + d_0)$ with BH `FDR`, and computes marginal covariate-adjusted `Psi_A`, `Psi_B`, and `DeltaPsi` sharing exact sign with `LogOddsRatio`.

<details>
<summary><b>Mathematical Details: Beta-Binomial GLM, Cox–Reid APL & Smyth Bias Correction</b></summary>

- **Effective depth & continuous Beta-Binomial log-likelihood:** With $v_j = \text{InfVar}_j$ and nominal Beta variance $v_{\text{nom}}(\Psi_j, n_j) = (\Psi_j n_j + 1)((1-\Psi_j)n_j + 1)/((n_j+2)^2(n_j+3))$,

  $$
  n_{\text{eff}, j} = \max\!\left(1.0,\; n_j \cdot \text{clamp}\!\left(\frac{v_{\text{nom}}(\Psi_j, n_j)}{v_j},\; 0.1,\; 1.0\right)\right), \qquad k_j = \Psi_j n_{\text{eff}, j}, \qquad r_j = n_{\text{eff}, j} - k_j,
  $$

  precision $M = (1-\phi)/\phi = e^\theta$, and $\mu_j = \text{logistic}(x_j^\top \beta)$:

  $$
  \ell(\beta, M) = \sum_{j=1}^{m} \left[\ln\Gamma(k_j + M\mu_j) + \ln\Gamma(r_j + M(1-\mu_j)) - \ln\Gamma(n_{\text{eff},j} + M) + \ln\Gamma(M) - \ln\Gamma(M\mu_j) - \ln\Gamma(M(1-\mu_j))\right].
  $$

  Solved via damped IRLS (`< 1e-14` `trigamma`, adaptive Cholesky damping, backtracking line search, ridge $\lambda = 10^{-4}$, intercept $10^{-6}$) for guaranteed finite convergence under complete separation.

- **Constitutive exclusion:** Events with $\sum k_j < 1$, $\sum r_j < 1$, pooled $\bar{\Psi} \le 0.01$ or $\ge 0.99$, or $\max \Psi_j - \min \Psi_j < 10^{-4}$ receive their stratum prior $\mu_\theta$ directly.

- **Cox–Reid APL & Smyth finite-df median bias correction** within each stratum (**Spliced** vs. **Tandem UTR**, fallback if $< 20$ interior events) on $\theta \in [-4.0, 11.5]$:

  $$
  \text{APL}(\theta) = \ell(\hat{\beta}(\theta), e^\theta) - \tfrac{1}{2}\ln\det\!\left(X^\top W X + 10^{-8} I\right), \qquad \mu_\theta = \text{median}(\hat{\theta}) + \ln\!\left(\frac{\text{median}(\chi^2_{d_{\text{res}}})}{d_{\text{res}}}\right), \qquad d_0 = \text{clamp}\!\left(\frac{2}{\sigma_\theta^2},\; 1,\; 30\right).
  $$

- **Pass 3 inference:** Maximizes $\text{APL}(\theta) - (\theta - \mu_\theta)^2 / (2\sigma_\theta^2)$ for shrunk $\tilde{M}_i = e^{\tilde{\theta}_i}$, computes $\text{LRT} = 2(\ell_{\text{full}} - \ell_{\text{null}})$, and evaluates $p_i = P(F_{\text{df},\; d_{\text{res}} + d_0} \ge \text{LRT}/\text{df})$.

</details>

---

## 7. Benchmarks & Validation

### 7.1 Simulated Human Cohort: `E-MTAB-3766` (GRCh37.71, 480M Reads)

Benchmarked on `E-MTAB-3766` (Soneson et al., 2016: 3 vs. 3 paired-end 101 bp replicates, **20,410 protein-coding genes**, **145,342 transcripts**, **1,000 truly DS genes**, 40M pairs/sample = **480M reads**, 16T, `--biascorrect`, Apple Silicon arm64, `/usr/bin/time -l`) against our multithreaded Julia 1.9 fork (`Whippet2-turbo`), original `Whippet.jl` v1.6.2 (Julia 1.6, 1T, Rosetta 2), and published callers.

#### A. Runtime, Throughput & Memory
| Stage (`E-MTAB-3766`, `--biascorrect`) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(Julia 1.9, 16T)* | Original `Whippet.jl` v1.6.2 *(1T)* | Langur Advantage |
|---|---:|---:|---:|---|
| **Index Build** *(uncompressed FASTA)* | **14.6 s** · **4.36 GB** *(5.06 GB RSS)* · 820 MB *(incl. 1.10 Gb decoys)* | 23.0 s · 4.21 GB *(4.37 GB RSS)* · 529 MB | 287.1 s *(4m 47s)* · 8.10 GB *(8.18 GB RSS)* · 473 MB | **19.7× vs. W1 · 1.58× vs. W2** |
| **Index Load** | **0.03 ms** (`mmap`) | 6.20 s | 7.33 s | **>200,000× faster** |
| **Sample 1 `quant`** (80M reads) | **40.5 s** (**1.98M/s**, 11.1× CPU) · **1.40 GB** | 81.1 s (`987k/s`) · 9.62 GB | 1,608.4 s *(26m 48s)* · 7.19 GB | **39.7× vs. W1 · 2.00× vs. W2 · 6.9× less RAM** |
| **6-Sample `--batch`** (480M reads) | **248.9 s** *(4m 09s, **1.93M/s**)* · **1.42 GB** | 433.8 s (`1.11M/s`) · 9.83 GB | *9,650.6 s* *(2h 41m, 6× extrap.)* · 7.19 GB | **38.8× vs. W1 · 1.74× vs. W2 · 6.9× less RAM** |
| **Mapped Pairs** (Sample 1 of 40M) | **37.94M (94.8%)** *(truth: 95.0%)* | 36.80M (92.0%) | 36.80M (92.0%) | **+1.14M more mapped pairs/sample** |
| **`delta` Testing** (223k events, `-r 5 -s 2`) | **2.50 s** · **191 MB** | 27.38 s (GLM) / 85.51 s (`--legacy`) · 1.87 GB | 85.51 s *(`--legacy`)* · 609 MB | **34.2× vs. legacy · 11.0× & 9.8× RAM vs. W2** |

*Sample 1 post-alignment stages in Langur total **0.79 s** (`vbem 0.18 s · events 0.53 s · write 0.08 s`, **15.6× faster** than `Whippet2-turbo`'s `12.3 s`). `Whippet2-turbo` commits: `1e860f9` (`TS` bugfix), `c725028` (stratified EB GLM), `400f0f3` (low-RAM `--batch`), `94a491f` (fast `whippet-index.jl`).*

#### B. Expression & PSI Accuracy (Sample 1 vs. Exact RSEM Ground Truth)
| Metric (145,342 Transcripts / 20,410 Genes) | Original `Whippet.jl` v1.6.2 | `Whippet2-turbo` (`c725028`) | **Langur (CC-VBEM)** | Langur Improvement |
|---|---:|---:|---:|---|
| **Isoform Spearman $\rho$** (all 145k tx) | 0.7954 | 0.8029 | **0.8776** | +0.082 vs. W1 · +0.075 vs. W2 |
| **Isoform MedARD** *(Median ARD, $\text{TPM} \ge 1$)* | 0.1145 | 0.0687 | **0.0402** | **2.85× vs. W1 · 1.71× vs. W2** |
| **Gene Spearman $\rho$ / $\log_2$ Pearson $r$** | 0.9797 / 0.9861 | 0.9810 / 0.9880 | **0.9956 / 0.9971** | Near-exact gene abundance |
| **Gene MedARD** ($\text{TPM} \ge 1$) | 0.0753 | 0.0443 | **0.0137** | **5.50× vs. W1 · 3.23× vs. W2** |
| **vs. `Salmon`** — Isoform / Gene MedARD ($\ge 1$) | — | 0.0602 / 0.0383 | **0.0223 / 0.0090** | Near-parity with dedicated `Salmon` |
| **`TS` PSI concordance vs. Langur** (`Reads >= 30`) | -0.1085 *(upstream reversal bug)* | **0.9666** *(fixed in `1e860f9`)* | *(reference)* | Confirmed upstream `TS` bug |

*Stage 2 PSI concordance vs. fixed `Whippet2-turbo` (`172,717` events, $\ge 30$ reads): $r = 0.9779$ ($\rho = 0.9730$, $\text{MAE} = 0.0257$), with all canonical spliced classes at $r = 0.9863\text{–}0.9956$ (`CE 0.9956`, `AA 0.9945`, `AD 0.9942`, `RI 0.9942`, `AL 0.9941`, `AF 0.9863`, `TS 0.9666`, `TE 0.8954`).*

#### C. Gene-Level Differential Splicing Accuracy (20,410 Genes, 1,000 Truly DS)
| Method (`E-MTAB-3766`, Gene-Level Simes+BH) | Event Set | TPR @0.01 | FDR @0.01 | TPR @0.05 | FDR @0.05 | TPR @0.10 | FDR @0.10 | AUROC | AUPRC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Langur Stratified + `Inf_Var` (`-r 5 -s 2`)** | All (`TS/TE` + Spliced) | **0.273** | **0.028** | **0.336** | **0.054** | **0.361** | **0.084** | **0.8302** | **0.5650** |
| **Langur Stratified + `Inf_Var` (`-r 5 -s 2`)** | Spliced (`CE/AA/AD/RI/AF/AL`) | **0.294** | **0.003** | **0.338** | **0.020** | **0.359** | **0.027** | **0.8122** | **0.5672** |
| `Whippet2-turbo` (`c725028`) Stratified GLM | All / Spliced | 0.225 / 0.252 | 0.017 / 0.008 | 0.306 / 0.309 | 0.053 / 0.013 | 0.329 / 0.332 | 0.073 / 0.018 | 0.8209 / 0.8007 | 0.5512 / 0.5486 |
| `Whippet2` `--legacy` Beta MC (`Prob ≥ 0.99/0.95/0.90`) | All | 0.585 | **0.681** | 0.744 | **0.815** | 0.867 | **0.865** | 0.8903 | **0.3517** |
| `Whippet2` `--legacy` Beta MC (`\|ΔΨ\| ≥ 0.1`) | All | 0.532 | **0.534** | 0.638 | **0.669** | 0.692 | **0.740** | 0.8649 | **0.3926** |
| Published: `rMATS_junction` (Soneson et al. 2016) | Junction | 0.336 | **0.463** | 0.379 | **0.571** | 0.409 | **0.623** | 0.7666 | **0.3937** |
| Published: `cuffdiff_ensembl` (Soneson et al. 2016) | Isoform | 0.114 | 0.000 | 0.142 | 0.000 | 0.159 | 0.000 | 0.6110 | **0.2599** |
| Published: `dexseq_splicinggraph` (Soneson et al. 2016) | Exon-bin | 0.637 | **0.172** | 0.675 | **0.252** | 0.697 | **0.294** | 0.8913 | 0.7272 |

**Key findings & ablations:**
- **Calibrated FDR & top CSG power:** Slashes realized FDR **15×** vs. legacy Beta MC (`81.5% → 5.4%` All, `2.0%` Spliced; `+0.21` AUPRC) and **10.6×** vs. `rMATS_junction` (`+0.17` AUPRC, `+0.064` AUROC). Outperforms `Whippet2-turbo` by **+30 TP genes** All (`+0.0138` AUPRC) and **+29 TP** Spliced (`+0.0186` AUPRC; `CE` alone: `163 TP / 1 FP` = `0.006` FDR, `0.3473` AUPRC) via a tighter spliced EB prior (`sd(log M) = 0.395`, `d0 = 12.83` vs. `9.77`). *(Exon-bin `DEXSeq` achieves higher recall on UTR-driven swaps but exceeds nominal 5% FDR by 5×–10×.)*
- **Ablation proof:** (1) Stratifying Spliced vs. Tandem UTR priors (`-r 5 -s 1`) boosted `TPR @0.05` **6×** All (`0.042 → 0.249`) and **25×** Spliced (`0.011 → 0.276`) while cutting UTR FPs **>300×**. (2) Enabling Gibbs `Inf_Var` (`-r 5 -s 2`) beat `--ignore-inf-var` by **+33 TP** All (`0.303 → 0.336`) and **+32 TP** Spliced (`0.306 → 0.338`), lifting AUPRC (`0.5541 → 0.5650` All; `0.5633 → 0.5672` Spliced) and tightening `d0` (`8.12 → 12.83`).

### 7.2 Automated Verification Suite & Synthetic Cohort
- **42/42 tests pass** (`cargo test --workspace`, zero warnings), covering Whippet 2 CSG/motif parity (`ENSMUSG00000038685` 18-path, `nrf2test` 6-path), bounded bridge sweeps (`LTBP2`), FR/RF/unstranded pairing, pseudogene decoy ties vs. diverged decoys, 1-SNP splice jumps, banded Gotoh indels, CC-VBEM paralog routing (`>99.5%`), 15 bp microexon PSI, PCR-stack soft-capping, Gibbs ~90% CI coverage across 200 $K_2$ replicates, constitutive exclusion + stratified EB priors, multivariable batch/interaction/continuous/collinear designs, and 100% thread determinism.
- **200-gene 3v3 synthetic cohort (600k reads, 4T):** `index` `0.032 s` (185× vs. W2), 6-sample `--batch` `0.377 s` (126× faster, 17× less RAM), `delta` `0.039 s` (122× faster, 95× less RAM); **100% TPR (100/100)** at **0.00% FDR** ($\Delta\Psi\ r = 0.9935$).

### 7.3 Real-Data Long-Read Validation: LRGASP WTC11 (GRCh38, 5-Way Illumina vs. PacBio HiFi)

Following MAJIQ V3 (2026), we benchmarked all three CSG generations alongside **`Salmon 2.12.0` + `SUPPA2 2.4`** (native arm64 Rust Salmon with selective decoys + `--gcBias --seqBias`) and **`STAR 2.7.11b` + `rMATS-turbo v4.3.0`** on real human WTC11 RNA-seq (LRGASP, Pardo-Palacios et al., 2024; ENCODE `ENCSR673UKZ` Rep 1: **45.78M pairs = 91.57M unstranded 100 bp reads**, 16T, `--biascorrect`) on full **GRCh38 + GENCODE v38 basic** (**60,649 genes**, **111,832 transcripts**). Orthogonal ground-truth PSI was tallied directly from pooled 3-replicate matched WTC11 PacBio HiFi FLNC alignments (ENCODE `ENCSR507JOF`: `6.79M` reads, `MAPQ >= 20`, `ts` strand-oriented; requiring $\text{LR} \ge 20$, $\text{SR} \ge 10$).

#### A. GRCh38 Index Build & WTC11 End-to-End Performance (5-Way)
| Stage (GRCh38 GENCODE v38) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(`94a491f`, 16T)* | `Whippet.jl` v1.6.2 *(1T)* | `Salmon 2.12.0` + `SUPPA2 2.4` *(16T)* | `STAR 2.7.11b` + `rMATS-turbo v4.3.0` *(16T)* |
|---|---:|---:|---:|---:|---:|
| **Index / Setup** *(uncompressed FASTA)* | **14.8 s** · **5.19 GB** *(6.16 GB RSS)* · **1.06 GB** *(incl. 1.31 Gb decoys)* | 28.6 s · 5.33 GB · 691 MB | 753.8 s · 9.87 GB · 634 MB | **110.8 s** *(`102.0 s` + `8.8 s`)* · **17.93 GB** *(25.56 GB RSS)* · **8.51 GB** | **683.1 s** *(11m 23s)* · **35.09 GB** *(35.02 GB RSS)* · **29.5 GB** |
| **Index Load** | **0.11 ms** (`mmap`) | 1.83 s | 1.71 s | 1.02 s | 14.0 s |
| **End-to-End** *(91.6M reads)* | **72.7 s** (`1.26M/s`) · **1.35 GB** *(2.83 GB RSS)* | 82.3 s (`1.11M/s`) · 9.10 GB | 1,467.4 s (`62.4k/s`) · 7.00 GB | **190.5 s** *(`188.7 s` + `1.8 s`)* · **12.24 GB** | **526.1 s** *(`357.0 s` + `169.1 s`)* · **31.64 GB** *(STAR)* / 493 MB *(rMATS)* |
| **BAM Disk** | **0 GB** | 0 GB | 0 GB | **0 GB** | **5.1 GB** *(sorted BAM)* |
| **Mapped Pairs** *(of 45.78M)* | **30.00M (65.5%)** *(382k decoy-filtered)* | 24.70M (53.9%) | 24.70M (53.9%) | **30.44M (66.5%)** | 30.55M (66.7%) unique + 1.87M multi *(genome)* |
| **Event & Quant Scope** | `CE/AA/AD/RI/AF/AL`, `TS/TE`, $K_n \ge 2$, TPMs | Same | Same | `SE/A3/A5/RI/AF/AL/MX` *(from TPMs)*, TPMs | Binary `SE/A3SS/A5SS/RI/MXE` only *(no `AF`/`AL`/TPMs)* |

*Langur post-alignment stages total `4.16 s` (`vbem 3.31 s · events 0.79 s · write 0.06 s`) across 394,702 events vs. `10.44 s` in W2 and `86.87 s` in W1. Minor background I/O overlapped `salmon quant` (`12.0×` avg CPU = `2,260 s` user); retired instructions directly confirm the efficiency gap (`30.4T` in `Salmon 2.12.0` vs. `10.6T` in Langur, 2.9× fewer).*

#### B. 3-Way CSG Stratified PacBio HiFi Concordance (Polymorphic, $n = 13{,}359$)
| Stratum ($0.05 \le \Psi_{\text{LR}} \le 0.95$) | $n$ | **Langur ($r$ / $\rho$ / MAE)** | `Whippet2-turbo` ($r$ / $\rho$ / MAE) | Original `Whippet.jl` ($r$ / $\rho$ / MAE) |
|---|---:|---|---|---|
| **All canonical** | 13,359 | **0.9338 / 0.8715 / 0.0789** | 0.9316 / 0.8689 / 0.0821 | 0.9064 / 0.8473 / 0.1023 |
| **CE** (Core exon) | 8,704 | **0.9251 / 0.7852 / 0.0761** | 0.9217 / 0.7822 / 0.0783 | 0.9121 / 0.7644 / 0.0870 |
| **AA** (Alt. acceptor) | 1,699 | **0.9266 / 0.9145 / 0.0758** | 0.9247 / 0.9143 / 0.0808 | 0.8440 / 0.8517 / 0.1494 |
| **AD** (Alt. donor) | 1,390 | **0.9221 / 0.9079 / 0.0821** | 0.9154 / 0.8988 / 0.0870 | 0.8567 / 0.8425 / 0.1289 |
| **RI** (Retained intron) | 815 | **0.9029 / 0.8862 / 0.0899** | 0.8996 / 0.8787 / 0.0955 | 0.8760 / 0.8537 / 0.1109 |
| **AF** (Alt. first exon) | 489 | **0.9088 / 0.8535 / 0.0965** | 0.9062 / 0.8388 / 0.0992 | 0.8936 / 0.8158 / 0.1096 |
| **AL** (Alt. last exon) | 262 | **0.8425 / 0.7417 / 0.1098** | 0.8281 / 0.7071 / 0.1183 | 0.8248 / 0.7168 / 0.1232 |
| **Binary ($K_0 / K_1$)** | 6,065 | **0.9332 / 0.8260 / 0.0841** | 0.9304 / 0.8229 / 0.0870 | 0.9107 / 0.8117 / 0.1030 |
| **Complex ($K_n \ge 2$)** | 7,294 | **0.9352 / 0.8995 / 0.0746** | 0.9333 / 0.8974 / 0.0780 | 0.9018 / 0.8662 / 0.1017 |

<details>
<summary><b>All-Matched 3-Way CSG Table (n = 117,719) &amp; Per-Tool Own-Coverage</b></summary>

| Stratum (All Matched) | $n$ | **Langur ($r$ / $\rho$ / MAE)** | `Whippet2-turbo` ($r$ / $\rho$ / MAE) | Original `Whippet.jl` ($r$ / $\rho$ / MAE) |
|---|---:|---|---|---|
| **All canonical** | 117,719 | **0.9801** / 0.6774 / **0.0162** | 0.9790 / **0.6825** / 0.0165 | 0.9738 / 0.6631 / 0.0185 |
| **CE** | 106,036 | **0.9682** / 0.5697 / **0.0128** | 0.9664 / **0.5758** / 0.0130 | 0.9649 / 0.5438 / 0.0132 |
| **AA** | 3,587 | **0.9676 / 0.9409 / 0.0491** | 0.9651 / 0.9372 / 0.0513 | 0.9352 / 0.9216 / 0.0855 |
| **AD** | 2,988 | **0.9644 / 0.9393 / 0.0527** | 0.9614 / 0.9381 / 0.0558 | 0.9394 / 0.9163 / 0.0778 |
| **RI** | 2,004 | **0.9720** / 0.9341 / 0.0533 | 0.9704 / **0.9376** / **0.0509** | 0.9599 / 0.9272 / 0.0599 |
| **AF** | 2,255 | **0.9777 / 0.8147 / 0.0280** | 0.9772 / 0.7923 / 0.0289 | 0.9735 / 0.7751 / 0.0318 |
| **AL** | 849 | **0.9500** / 0.7729 / **0.0480** | 0.9434 / **0.7730** / 0.0500 | 0.9422 / 0.7585 / 0.0519 |
| **Binary ($K_0 / K_1$)** | 91,814 | **0.9806 / 0.5518 / 0.0114** | 0.9790 / 0.5505 / 0.0118 | 0.9763 / 0.5391 / 0.0125 |
| **Complex ($K_n \ge 2$)** | 25,905 | **0.9765** / 0.8589 / 0.0330 | 0.9760 / **0.8667** / **0.0328** | 0.9666 / 0.8384 / 0.0396 |

*Per-tool own-coverage ($\text{LR} \ge 20$, tool $\text{SR} \ge 10$): Langur `120,127` events (`14,004` poly: `0.9286 / 0.8663 / 0.0821`), `Whippet2-turbo` `119,842` (`13,848` poly: `0.9289 / 0.8666 / 0.0842`), `Whippet.jl` `124,279` (`14,775` poly: `0.9003 / 0.8386 / 0.1073`).*

</details>

#### C. 5-Way Paradigm Comparison & Exact-Event Sensitivity Checks

**Strictly matched 5-way polymorphic intersection (`CE/AA/AD/RI`, $n = 4{,}693$, Node-Level PacBio Truth):**

| Stratum | $n$ | **Langur** ($r$ / $\rho$ / MAE) | `Whippet2-turbo` | `Whippet.jl` | **`Salmon 2 + SUPPA2`** | **`STAR + rMATS JC`** |
|---|---:|---|---|---|---|---|
| **All 5-way shared** | 4,693 | **0.9507 / 0.9350 / 0.0688** | 0.9456 / 0.9291 / 0.0728 | 0.9105 / 0.8974 / 0.1012 | 0.6832 / 0.6688 / 0.1747 | 0.9099 / 0.8822 / 0.0896 |
| **`CE / SE`** | 2,998 | **0.9556 / 0.9332 / 0.0686** | 0.9484 / 0.9244 / 0.0728 | 0.9388 / 0.9099 / 0.0838 | 0.6769 / 0.6619 / 0.1788 | 0.9357 / 0.9038 / 0.0824 |
| **`AA / A3`** | 977 | **0.9407 / 0.9332 / 0.0669** | 0.9362 / 0.9288 / 0.0719 | 0.8538 / 0.8673 / 0.1488 | 0.6887 / 0.6718 / 0.1692 | 0.8667 / 0.8450 / 0.0958 |
| **`AD / A5`** | 517 | 0.9564 / **0.9479 / 0.0664** | **0.9577** / 0.9463 / 0.0679 | 0.9089 / 0.9045 / 0.1115 | 0.7326 / 0.7250 / 0.1579 | 0.8896 / 0.8643 / 0.0972 |
| **`RI`** | 201 | **0.9121 / 0.9070 / 0.0865** | 0.9050 / 0.8954 / 0.0906 | 0.8874 / 0.8805 / 0.1021 | 0.6545 / 0.6345 / 0.1841 | 0.7941 / 0.7819 / 0.1477 |
| **Binary ($K_0 / K_1$)** | 2,008 | **0.9467 / 0.9292 / 0.0740** | 0.9424 / 0.9245 / 0.0778 | 0.9126 / 0.8994 / 0.1032 | 0.7272 / 0.7166 / 0.1599 | 0.9414 / 0.9188 / 0.0767 |
| **Complex ($K_n \ge 2$)** | 2,685 | **0.9551 / 0.9413 / 0.0648** | 0.9493 / 0.9340 / 0.0692 | 0.9090 / 0.8959 / 0.0996 | 0.6510 / 0.6327 / 0.1859 | 0.8863 / 0.8533 / 0.0993 |

**Exact-event sensitivity checks (`SUPPA2` and `rMATS` under their own exact flanking-junction PacBio estimands):**

| Sensitivity Check (Polymorphic, Exact Tool-Event Truth) | $n$ | **Langur ($r$ / $\rho$ / MAE)** | `Whippet2-turbo` | `Whippet.jl` | Target Baseline ($r$ / $\rho$ / MAE) |
|---|---:|---|---|---|---|
| **vs. `Salmon 2 + SUPPA2`** — Binary ($K_0/K_1$, identical estimand) | 2,339 | **0.9403 / 0.9232 / 0.0734** | 0.9349 / 0.9193 / 0.0773 | 0.8884 / 0.8778 / 0.1083 | 0.7393 / 0.7299 / 0.1563 *(`SUPPA2`)* |
| **vs. `Salmon 2 + SUPPA2`** — Complex ($K_n \ge 2$) | 2,971 | **0.8180 / 0.8062 / 0.1179** | 0.8168 / 0.8053 / 0.1207 | 0.7793 / 0.7706 / 0.1478 | 0.7262 / 0.7203 / 0.1602 *(`SUPPA2`)* |
| **vs. `Salmon 2 + SUPPA2`** — All 6 motifs (`SE/A3/A5/RI/AF/AL`) | 5,310 | **0.8718 / 0.8589 / 0.0983** | 0.8686 / 0.8565 / 0.1016 | 0.8281 / 0.8194 / 0.1304 | 0.7322 / 0.7255 / 0.1585 *(`SUPPA2`)* |
| **vs. `STAR + rMATS JC`** — Binary ($K_0/K_1$, identical estimand) | 2,747 | **0.9398 / 0.8728 / 0.0692** | 0.9356 / 0.8712 / 0.0726 | 0.9094 / 0.8514 / 0.0931 | **0.9507 / 0.8892 / 0.0665** *(`rMATS JC`)* |
| **vs. `STAR + rMATS JC`** — Complex ($K_n \ge 2$, different estimands) | 3,744 | 0.7942 / 0.7190 / 0.1216 | 0.7906 / 0.7158 / 0.1245 | 0.7616 / 0.6979 / 0.1449 | 0.8732 / 0.8110 / 0.0979 *(`rMATS JC`)* |

*Reference intersections: 4-way alignment-free across all 6 motifs (`n = 6,436` poly, node truth): Langur `0.9408 / 0.9257 / 0.0749` vs. `Salmon 2 + SUPPA2` `0.6515 / 0.6379 / 0.1863` (all `SUPPA2` `TPM >= 1` poly, `n = 9,893`: `0.5765 / 0.5625 / 0.2107`). 4-way with `STAR + rMATS JC` (`n = 7,815` poly, node truth): Langur `0.9433 / 0.8907 / 0.0687` vs. `rMATS JC` `0.9041 / 0.8285 / 0.0880` (`Binary K0/K1` `n = 3,106`: `0.9441 / 0.0709` vs. `0.9440 / 0.0715`; all `rMATS JC` poly `n = 8,531` exact truth: `0.8445 / 0.7858 / 0.1127`; `JCEC`: `0.8345 / 0.7737 / 0.1167`).*

#### D. Real-Data Gene & Isoform Expression Concordance vs. `Salmon 2.12.0`
Stage 1 abundance across `60,649` GENCODE v38 genes (`15,378` at $\text{TPM} \ge 1$) and `111,334` non-duplicate transcripts (`27,220` at $\text{TPM} \ge 1$) vs. selective-alignment `Salmon 2.12.0` (`--gcBias --seqBias`):

| Tool (WTC11 Rep 1, `--biascorrect`) | Gene $\rho$ (all / $\ge 1$) | Gene $\log_2 r$ | **Gene MedARD ($\ge 1$)** | Isoform $\rho$ (all / $\ge 1$) | Isoform $\log_2 r$ | **Isoform MedARD ($\ge 1$)** |
|---|---:|---:|---:|---:|---:|---:|
| **Langur (CC-VBEM)** | **0.9244 / 0.9684** | **0.9675** | **0.0410** *(3.2×–3.4× closer)* | **0.8989 / 0.9211** | **0.9518** | **0.0678** *(2.1×–2.3× closer)* |
| `Whippet2-turbo` | 0.9019 / 0.9607 | 0.9638 | 0.1305 | 0.8188 / 0.8991 | 0.9367 | 0.1399 |
| Original `Whippet.jl` v1.6.2 | 0.9017 / 0.9595 | 0.9624 | 0.1402 | 0.8127 / 0.8780 | 0.9266 | 0.1557 |

**Key findings:**
1. **27/27 clean sweep across CSG generations:** Langur ranks #1 across every polymorphic motif and complexity stratum, cutting `AA` MAE **-49%**, `AD` **-36%**, and complex $K_n \ge 2$ MAE **-27%** vs. `Whippet.jl` while improving every stratum over `Whippet2-turbo`.
2. **Expression parity with `Salmon 2.12.0`, decisive splicing win over `Salmon 2 + SUPPA2`:** Langur's CC-VBEM tracks `Salmon 2.12.0` **3.2×–3.4× closer** at the gene level and **2.1×–2.3× closer** at the isoform level than legacy Whippet, yet avoids the isoform-deconvolution error that leaks into `SUPPA2` PSI—beating `Salmon 2 + SUPPA2` across **27/27 strata** even under SUPPA2's own exact-event truth (**-53% MAE** on binary $K_0/K_1$) at **7.5×** faster indexing, **2.6×** faster quant, and **9.1×** lower peak RAM.
3. **Binary parity with `STAR + rMATS-turbo` at order-of-magnitude lower cost:** Matches `STAR + rMATS` on shared binary events (within `0.0027` MAE even on rMATS's exact-junction truth) while capturing full-length complex ($K_n \ge 2$) exon inclusion (**-35% MAE** vs. `rMATS`), building the index **46.2× faster** (`28×` smaller), running **7.24× faster** with **23.5× less RAM** and **zero BAM files**, and quantifying `AF`/`AL`, tandem UTRs, and TPMs in one pass.

### 7.4 Real-Cohort Null Calibration & Multivariable Batch Adjustment: Geuvadis (`ERP001942`)

Following MAJIQ V2 (Vaquero-Garcia et al., 2023), we evaluated empirical null calibration and real-world batch-confounding adjustment across **14 real human CEU LCL donors** from Geuvadis (Lappalainen et al., 2013; `ERP001942`, **353.4M PE pairs = 706.9M 75 bp reads**, GRCh38 GENCODE v38, 16T, `--biascorrect`), intentionally imbalanced across two sequencing centers (`UNIGE`: `6F + 2M`; `ICMB`: `2F + 4M`).

#### A. 14-Sample Cohort Runtime, Memory & Alignment
| Tool (14 Geuvadis Donors, GRCh38) | 14-Sample `quant --batch` | Throughput | Peak Footprint / RSS | Mapped Pairs (of 353.4M) | 6-Sample `delta` *(10-split mean)* | 14-Sample `delta` *(`~ center + sex`)* | `delta` Peak Footprint |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Langur** *(Rust, 16T)* | **347.3 s** *(5m 47s)* | **2.03M/s** | **1.39 GB** / 2.92 GB | **86.05%** *(304.2M; 74.4–93.2%)* | **1.93 s** | **7.02 s** *(4.33 s unadj.)* | **177 MB** |
| **`Whippet2-turbo`** *(Julia 1.9, 16T)* | 429.2 s *(7m 09s)* | 1.65M/s | 9.88 GB / 10.13 GB | 73.66% *(260.3M; 65.1–80.4%)* | 20.77 s *(10.8×)* | 46.51 s *(30.93 s unadj.)* | 2.86 GB *(16.2×)* |
| **Original `Whippet.jl` v1.6.2** | *Parser error*¹ | — | — | — | 61.28 s *(31.8×)* | 86.94 s² *(unadj. only)* | 668 MB |

*¹ Original `Whippet.jl` v1.6.2 `whippet-quant.jl` fails on raw SRA FASTQs (`ERR188086_1.fastq.gz`) with `Cannot encode 78` — the known upstream SRA FASTQ parser issue (`timbitz/Whippet.jl#135`, avoidable in v1.6.2 by stripping `+ERR...` line-3 headers to bare `+`, whereas `Whippet2-turbo` and `Langur` parse raw SRA FASTQs natively). Legacy Beta MC (`whippet-delta.jl --legacy`) was therefore run on `Whippet2-turbo`'s identically aligned `.psi.gz` files. ² On the 14-sample unadjusted `6M vs 8F` run, `whippet-delta.jl --legacy` also threw `DomainError("Beta: the condition α > zero(α) is not satisfied")` on discordant biological outliers (`[0.9999, 0.001, ...]`) at default `--max-variance` and required `--max-variance 0.05` to complete.*

#### B. Real-Data Same-Condition Null Calibration (10 Exact Distinct 3v3 Splits, 6 Female UNIGE Donors)
All $\binom{6}{3}/2 = 10$ unique 3v3 partitions of the 6 same-sex, same-center donors (`-r 5 -s 2`):

| Method | Event Set (mean $n$) | FPR @ $p \le 0.05$ *(mean / med)* | FPR @ $p \le 0.01$ *(mean / med)* | FP @ $\text{FDR} \le 0.05$ *(mean / med)* | Per-Split FPs (`10` splits) |
|---|---:|---:|---:|---:|---|
| **Langur (`Inf_Var`)** | All (`147,196`) | **1.38% / 1.33%** | **0.22% / 0.21%** | **0.2 / 0.0** | `[0,2,0,0,0,0,0,0,0,0]` (`2` total) |
| **Langur (`Inf_Var`, Spliced)** | Spliced (`105,189`) | **0.70% / 0.70%** | **0.04% / 0.04%** | **0.2 / 0.0** | `[0,2,0,0,0,0,0,0,0,0]` (`2` total) |
| Langur (`--ignore-inf-var`) | All / Spliced | 1.42% / 0.76% *(mean)* | 0.22% / 0.04% *(mean)* | 0.3 / 0.0 *(both)* | `[0,3,0,0,0,0,0,0,0,0]` (`3` total) |
| `Whippet2-turbo` (Strat. GLM) | All (`144,093`) | 1.44% / 1.38% | 0.22% / 0.19% | 0.8 / 0.0 | `[0,0,0,0,0,0,0,0,8,0]` (`8` total) |
| `Whippet2-turbo` (Spliced) | Spliced (`102,539`) | 0.83% / 0.84% | 0.05% / 0.06% | 0.3 / 0.0 | `[0,0,0,0,0,0,0,0,3,0]` (`3` total) |
| Legacy Beta MC (`Prob ≥ 0.95`) | All (`163,423`) | 4.44% *(call rate)* | — | **7,250.2 / 7,145.0** | `3,855.8` Spliced FPs / split |
| Legacy Beta MC (`P ≥ 0.95 & \|ΔΨ\| ≥ 0.1`) | All (`163,423`) | 2.50% *(call rate)* | — | **4,080.0 / 3,994.0** | `1,692.7` Spliced FPs / split |

#### C. Multivariable Batch Confounding (`~ center + sex`, 14 Imbalanced CEU Donors)
| Model | Spliced Median $\hat{\phi}$ | Pure Center Artifacts ($n = 87$) `med \|ΔΨ\|` / `FDR ≤ 0.05` Leaks | `chrX`/`chrY` Hits @ $p \le 0.01$ *(All / Spliced)* | `chrX`/`chrY` Enrichment vs. Autosomes *(All / Spliced)* | Total Calls (`FDR ≤ 0.05` or `Prob ≥ 0.95`) |
|---|---:|---:|---:|---:|---:|
| **Langur unadjusted (`M vs F`)** | `0.0383` | `0.080` / `0` | `11` / `4` | `1.11×` / `1.21×` | `0` (`0` Spliced) |
| **Langur adjusted (`~ center + sex`)** | **`0.0323`** *(-15.7%)* | **`0.052`** *(-35%)* / **`0`** | **`22` / `8`** *(2×)* | **`1.44×` / `1.67×`** | `13` (`0` Spliced)³ |
| `Whippet2-turbo` unadjusted (`M vs F`) | `0.0417` | `0.055` / `0` | `17` / `4` | `1.61×` / `1.08×` | `0` (`0` Spliced) |
| `Whippet2-turbo` adjusted (`~ center + sex`) | `0.0359` *(-13.9%)* | `0.045` / `0` | `27` / `6` | `1.70×` / `1.09×` | `2` (`0` Spliced) |
| Legacy Beta MC (`Prob ≥ 0.95`) | — *(no $\phi$)* | — | `108` *(1.95%)* | `1.04×` *(null)* | **3,276** *(3,166 auto.)* |
| Legacy Beta MC (`P ≥ 0.95 & \|ΔΨ\| ≥ 0.1`) | — *(no $\phi$)* | — | `57` *(1.03%)* | `1.12×` *(near-null)* | **1,593** *(1,536 auto.)* |

**Key findings:**
1. **Cohort speed, memory & alignment recovery:** Langur quantifies all 14 human donors **1.24× faster** (`5m 47s`) at **7.1× lower peak RAM** (`1.39 GB` vs. `9.88 GB`) while recovering **+43.8M** more mapped pairs (`+12.4 pp`), and runs 14-sample multivariable `delta` **6.6× faster** at **16.2× lower RAM** (`7.0 s`, `177 MB`) than `Whippet2-turbo`.
2. **>20,000× false-positive reduction on real human null splits:** Langur (`Inf_Var`) produces **zero false positives in 9/10 splits** (`0.2` mean / `0.0` median vs. `4,080` for legacy Beta MC at `Prob ≥ 0.95 & |ΔΨ| ≥ 0.1`, and `4×` fewer than `Whippet2-turbo`'s `0.8`), with a well-calibrated slightly conservative polymorphic spliced null tail (`|ΔΨ| ≥ 0.01`: `FPR@0.05 = 2.83%`, `q5 = 0.069`, `q10 = 0.101`) tightened by `Inf_Var` weighting over `--ignore-inf-var`.
3. **Batch debiasing & `chrX`/`chrY` signal recovery:** Adjusting for `~ center + sex` cuts residual spliced overdispersion $\hat\phi$ by **15.7%** (`0.0383 → 0.0323`; TS/TE `0.2638 → 0.2429`), reduces apparent sex `|ΔΨ|` bias across `87` pure concordant autosomal sequencing-center artifacts by **35%** (`0.080 → 0.052`, `0` leaks at `FDR ≤ 0.05`; e.g., center-driven `ZNF302` CE/AA/RI `ΔΨ_center = +0.671`, `FDR = 0.022` drops from unadjusted sex bias `-0.155` to `+0.067`, `FDR = 1.000`), and **doubles** `chrX`/`chrY` hits at $p \le 0.01$ (`11 → 22` All; `4 → 8` Spliced, enrichment `1.21× → 1.67×`, recovering X-linked `UBE2A`, `KDM6A`, `ENOX2`, `CXorf38`, `TRMT2B`, `ARHGEF6`, `HEPH`, `CDK16`, `HDAC6`), whereas legacy Beta MC floods the unadjusted comparison with `3,276` calls (`1.04×` near-null `chrX`/`chrY` enrichment). *³ All `13` `FDR ≤ 0.05` events in Langur `~ center + sex` (`0` Spliced) are single-stratum (2-donor) `UNIGE_M` or `ICMB_F` tandem UTR/TSS (`TS`/`TE`) corner outliers (e.g., `SLC35B2`, `BANP`, `MDM2`, `TUBB3`), where an additive main-effects logit model without a `center:sex` interaction term simultaneously attributes both `center` and `sex` effects.*

---

## 8. Key References

1. Sterne-Weiler T., Weatheritt R.J., Best A.J., Ha K.C.H., Blencowe B.J. (2018). Efficient and accurate quantitative profiling of alternative splicing patterns of any complexity on a laptop. *Molecular Cell*, 72(1):187–200.e6.
2. Srivastava A., Malik L., Sarkar H., Zakeri M., Almodaresi F., Soneson C., Love M.I., Kingsford C., Patro R. (2020). Alignment and mapping methodology influence transcript abundance estimation. *Genome Biology*, 21:239.
3. Patro R., Duggal G., Love M.I., Irizarry R.A., Kingsford C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14:417–419.
4. Soneson C., Matthes K.L., Nowicka M., Law C.W., Robinson M.D. (2016). Isoform prefiltering improves performance of count-based methods for analysis of differential transcript usage. *Genome Biology*, 17:12 (`E-MTAB-3766`).
5. Vaquero-Garcia J., Aicher J.K., Jewell S., et al. (2023). RNA splicing analysis using heterogeneous and large RNA-seq datasets. *Nature Communications*, 14:1230.
6. Aicher J.K., Slaff B., Jewell S., Gazzara M.R., Bhoj E.J.K., Barash Y. (2026). MAJIQ V3: Scalable and accurate RNA splicing analysis with splicegraph algebra and improved intron retention modeling. *NAR Genomics and Bioinformatics*, 8(4):lqag126.
7. Kabza M., Ritter A., Byrne A., Sereti K., Le D., Stephenson W., Sterne-Weiler T. (2024). Accurate long-read transcript discovery and quantification at single-cell, pseudo-bulk and bulk resolution with Isosceles. *Nature Communications*, 15:7316.
8. McCarthy D.J., Chen Y., Smyth G.K. (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Research*, 40(10):4288–4297.
9. Varadhan R., Roland C. (2008). Simple and globally convergent methods for accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*, 35(2):335–353.
10. Hansen K.D., Brenner S.E., Dudoit S. (2010). Biases in Illumina transcriptome sequencing caused by random hexamer priming. *Nucleic Acids Research*, 38(12):e131.
11. Pardo-Palacios F.J., Wang D., Reese F., et al. (2024). Systematic assessment of long-read RNA-seq methods for transcript identification and quantification (LRGASP). *Nature Methods*, 21:1349–1363.
12. Dobin A., Davis C.A., Schlesinger F., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1):15–21.
13. Wang Y., Xie Z., Kutschera E., Adams J.I., Kadash-Edmondson K.E., Xing Y. (2024). rMATS-turbo: an efficient and flexible computational tool for alternative splicing analysis of large-scale RNA-seq data. *Nature Protocols*, 19:1083–1104.
14. Trincado J.L., Entizne J.C., Hysenaj G., et al. (2018). SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. *Genome Biology*, 19:40.
15. Lappalainen T., Sammeth M., Friedländer M.R., et al. (2013). Transcriptome and genome sequencing uncovers functional variation in humans. *Nature*, 501:506–511.
