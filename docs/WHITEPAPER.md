# Langur: Fast, Statistically Calibrated Splice-Graph Quantification and Differential Splicing in Rust

*CLI binary: `langur` · Index format: `.lidx`*

---

## Abstract

Alignment-free quantification on **Contiguous Splice Graphs (CSGs)**, pioneered by Whippet (Sterne-Weiler et al., *Molecular Cell* 2018), profiles both binary ($K_1$) and complex ($K_n \ge 2$) alternative splicing directly from FASTQ files without genome alignment or BAM generation. However, scaling legacy CSG implementations to modern human cohorts exposes four bottlenecks: ungapped matching drops reads near indels and splice-proximal SNPs; transcriptome-only indexing allows diverged pseudogene and intronic reads to force-map onto parent genes; per-gene EM and naive Beta variance understate paralog and multi-path ambiguity; and two-group pooled Beta Monte Carlo testing (`whippet-delta`) cannot adjust for covariates and suffers severe false-discovery inflation under biological overdispersion.

**Langur** is a ground-up Rust re-architecture (`index`, `inspect`, `quant`, `delta`) combining a zero-copy `mmap` splice-graph + selective genomic decoy index (`< 1 ms` load), hybrid ungapped + zero-allocation banded Gotoh spliced alignment with 1-SNP junction tolerance, joint Connected-Component Variational Bayes EM (`CC-VBEM`) with unified 5'-hexamer/GC bias correction, targeted collapsed Gibbs inferential variance (`Inf_Var`), and a covariate-aware Beta-Binomial GLM (`delta`) with Spliced vs. Tandem UTR stratified Cox–Reid Empirical Bayes dispersion shrinkage. On real human WTC11 RNA-seq (LRGASP `ENCSR673UKZ`, GRCh38 GENCODE v38, 60,649 genes, 91.6M reads), Langur builds the index in **14.8 s** (**46.2×–47.3×** faster and **28×** smaller than `STAR 2.7.11b` / `STAR + MAJIQ v3`; **7.5×** faster than `Salmon 2.12.0 + SUPPA2 2.4`; **1.93×** faster than `Whippet2-turbo`) and quantifies WTC11 in **72.7 s** at **1.35 GB RAM** (**2.6×** faster and **9.1×** less RAM than Rust `Salmon 2.12.0 + SUPPA2 2.4`; **5.1×–7.24×** faster, **23.5×** less RAM, and zero BAMs vs. `STAR 2.7.11b + MAJIQ v3.0.25` and `STAR + rMATS-turbo v4.3.0`), sweeping **27/27** polymorphic CSG strata against matched PacBio HiFi Iso-Seq ($r = 0.9338$, $\text{MAE} = 0.0789$), matching `STAR + MAJIQ v3.0.25` head-to-head even on MAJIQ's own `12,899` annotated splice-junction LSV connections ($r = 0.9602$ vs. $0.9603$, $\text{MAE} = 0.0645$ vs. $0.0632$), beating `Salmon 2 + SUPPA2` across every stratum even on SUPPA2's own truth (binary $K_0/K_1$: $r = 0.9403$ vs. $0.7393$, **-53% MAE**), outperforming `STAR + rMATS-turbo` on full-length complex ($K_n \ge 2$) exon inclusion ($r = 0.9551$ vs. $0.8863$, **-35% MAE**), and tracking `Salmon 2.12.0` gene/isoform TPMs **3.2×–3.4× / 2.1×–2.3× closer** than legacy Whippet. On 14 real human CEU LCL donors from Geuvadis (`ERP001942`, 706.9M reads, GRCh38), Langur quantifies the cohort in **5m 47s (`347.3 s`, `2.03M reads/s`)** at **1.39 GB RAM** (**32.1×** faster and **4.5×** less RAM than original `Whippet.jl` v1.6.2, **1.24×** faster and **7.1×** less RAM than `Whippet2-turbo`, and **+43.8M** more mapped pairs than both), yields **0.2 mean / 0.0 median false positives** across all 10 exact 3v3 same-condition null splits (**7,545×** fewer than native `Whippet.jl` v1.6.2 `Prob ≥ 0.95` and **4×** fewer than `Whippet2-turbo`), and cuts residual spliced overdispersion **15.7%** while doubling `chrX`/`chrY` signal recovery under multivariable batch adjustment (`~ center + sex`) in **7.0 s** at **177 MB RAM**. Finally, on a 5-way real RF-stranded ENCODE `PTBP1` knockdown cohort (`K562`, `ENCSR527IVX` vs. `ENCSR129RWD`, 255.3M reads) validated against orthogonal ENCODE K562 `PTBP1` eCLIP peaks (`ENCFF907HNN`) and upstream CU-rich motifs, Langur finishes end-to-end in **3m 17s** (`3m 02s` quant+DS) at **1.35 GB RAM** with **0 BAMs** (**7.1×–7.6×** faster and **23.4×** less RAM than `STAR + MAJIQ v3.0.25` and `STAR + rMATS-turbo v4.3.0`; **2.05×** faster and **9.0×** less RAM than `Salmon 2 + SUPPA2`), calling **58/58 Spliced hits** (`0` `TS`/`TE` noise, `73.3%` `ΔΨ > 0`, `2.22×` upstream CU motif) while leading threshold-free Top-100 consensus (`96.4%`, `53/55` genes), upstream eCLIP hit rate (`6.0%`), and CU-motif enrichment (`1.77×`), ranking the `MAP4K4` 9 bp microexon **`#1` genome-wide (`+0.742`)**, and recovering the canonical `PTBP2` exon-10 NMD switch (`#30`, `+0.714`) that `MAJIQ v3` filtered out and `SUPPA2` missed.

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
3. **Excluding constitutive events from dispersion fitting:** ~63% of human CSG nodes are constitutive (`0.0` or `1.0` across samples) and carry zero variance information; including them causes the Cox–Reid penalty to peak near $\phi \approx 0.64$, collapsing the EB prior. Excluding boundary events restores full power and speeds up `delta` ~3× (to **1.4–7.0 s** on human cohorts).
4. **Event-class stratified EB priors (Spliced vs. Tandem UTR):** Tandem UTR ends (`TS`/`TE`) exhibit **4.4×–7.5× higher median overdispersion** than canonical splice junctions across real human cohorts (Geuvadis: `phi = 0.2429` vs. `0.0323` = **7.5×**; ENCODE `PTBP1` KD: `phi = 0.1305` vs. `0.0296` = **4.4×**). Fitting separate EB priors (`--no-stratify-eb` to opt out) preserves sensitivity on clean splice junctions while preventing noisy `TS`/`TE` false positives (`58/58` Spliced hits and `0` `TS`/`TE` calls in `PTBP1` KD).
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

### 7.1 Real-Data Long-Read Validation: LRGASP WTC11 (GRCh38, 6-Way Illumina vs. PacBio HiFi)

Following MAJIQ V3 (2026), we benchmarked all three CSG generations alongside **`Salmon 2.12.0` + `SUPPA2 2.4`** (native arm64 Rust Salmon with selective decoys + `--gcBias --seqBias`), **`STAR 2.7.11b` + `MAJIQ v3.0.25`** (native arm64 C++17 `majiq-v3` on official GENCODE v38 basic GFF3, `60,649` genes / `111,832` transcripts), and **`STAR 2.7.11b` + `rMATS-turbo v4.3.0`** on real human WTC11 RNA-seq (LRGASP, Pardo-Palacios et al., 2024; ENCODE `ENCSR673UKZ` Rep 1: **45.78M pairs = 91.57M unstranded 100 bp reads**, 16T, `--biascorrect`) on full **GRCh38 + GENCODE v38 basic** (**60,649 genes**, **111,832 transcripts**). Orthogonal ground-truth PSI was tallied directly from pooled 3-replicate matched WTC11 PacBio HiFi FLNC alignments (ENCODE `ENCSR507JOF`: `6.79M` reads, `MAPQ >= 20`, `ts` strand-oriented; requiring $\text{LR} \ge 20$, $\text{SR} \ge 10$).

#### A. GRCh38 Index Build & WTC11 End-to-End Performance (6-Way)
| Stage (GRCh38 GENCODE v38) | **Langur** *(Rust, 16T)* | `Whippet2-turbo` *(`94a491f`, 16T)* | `Whippet.jl` v1.6.2 *(1T)* | `Salmon 2.12.0` + `SUPPA2 2.4` *(16T)* | `STAR 2.7.11b` + `MAJIQ v3.0.25` *(16T)* | `STAR 2.7.11b` + `rMATS-turbo v4.3.0` *(16T)* |
|---|---:|---:|---:|---:|---:|---:|
| **Index / Setup** *(uncompressed FASTA)* | **14.8 s** · **5.19 GB** *(6.16 GB RSS)* · **1.06 GB** *(incl. 1.31 Gb decoys)* | 28.6 s · 5.33 GB · 691 MB | 753.8 s · 9.87 GB · 634 MB | **110.8 s** *(`102.0 s` + `8.8 s`)* · **17.93 GB** *(25.56 GB RSS)* · **8.51 GB** | **700.7 s** *(`683.1 s` STAR + `17.6 s` `gff3`)* · **35.09 GB** · **29.5 GB** | **683.1 s** *(11m 23s)* · **35.09 GB** *(35.02 GB RSS)* · **29.5 GB** |
| **Index Load** | **0.11 ms** (`mmap`) | 1.83 s | 1.71 s | 1.02 s | 14.0 s *(STAR)* | 14.0 s |
| **End-to-End** *(91.6M reads)* | **72.7 s** (`1.26M/s`) · **1.35 GB** *(2.83 GB RSS)* | 82.3 s (`1.11M/s`) · 9.10 GB | 1,467.4 s (`62.4k/s`) · 7.00 GB | **190.5 s** *(`188.7 s` + `1.8 s`)* · **12.24 GB** | **370.9 s** *(`357.0 s` STAR + `13.9 s` MAJIQ¹)* · **31.64 GB** *(STAR)* / 653 MB *(MAJIQ)* | **526.1 s** *(`357.0 s` + `169.1 s`)* · **31.64 GB** *(STAR)* / 493 MB *(rMATS)* |
| **BAM Disk** | **0 GB** | 0 GB | 0 GB | **0 GB** | **5.1 GB** *(sorted BAM + `.bai`)* | **5.1 GB** *(sorted BAM)* |
| **Mapped Pairs** *(of 45.78M)* | **30.00M (65.5%)** *(382k decoy-filtered)* | 24.70M (53.9%) | 24.70M (53.9%) | **30.44M (66.5%)** | 30.55M (66.7%) unique + 1.87M multi *(genome)* | 30.55M (66.7%) unique + 1.87M multi *(genome)* |
| **Event & Quant Scope** | `CE/AA/AD/RI/AF/AL`, `TS/TE`, $K_n \ge 2$, TPMs | Same | Same | `SE/A3/A5/RI/AF/AL/MX` *(from TPMs)*, TPMs | LSV junctions + IR + *de novo* *(no TPMs)* | Binary `SE/A3SS/A5SS/RI/MXE` only *(no `AF`/`AL`/TPMs)* |

*Langur post-alignment stages total `4.16 s` (`vbem 3.31 s · events 0.79 s · write 0.06 s`) across 394,702 events vs. `10.44 s` in W2 and `86.87 s` in W1. ¹ `majiq-v3` post-BAM stages are fast (`sj 8.88 s / 653 MB · build 1.39 s / 553 MB · psi-coverage 1.66 s / 505 MB · quantify 1.94 s / 531 MB` = `13.9 s`, plus `~10 s` `samtools index`), with runtime and memory dominated by `STAR`. Minor background I/O overlapped `salmon quant` (`12.0×` avg CPU = `2,260 s` user); retired instructions directly confirm the efficiency gap (`30.4T` in `Salmon 2.12.0` vs. `10.6T` in Langur, 2.9× fewer).*

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

#### C. 6-Way Paradigm Comparison & Exact-Event Sensitivity Checks

**Strictly matched 6-way polymorphic intersection (`CE/AA/AD/RI`, $n = 3{,}672$, Node-Level PacBio Truth):**

| Stratum | $n$ | **Langur** ($r$ / $\rho$ / MAE) | `Whippet2-turbo` | `Whippet.jl` | **`STAR + MAJIQ v3`** | **`STAR + rMATS JC`** | **`Salmon 2 + SUPPA2`** |
|---|---:|---|---|---|---|---|---|
| **All 6-way shared** | 3,672 | **0.9568 / 0.9402 / 0.0653** | 0.9518 / 0.9350 / 0.0697 | 0.9188 / 0.9052 / 0.0966 | 0.8924 / 0.8755 / 0.0948 | 0.9186 / 0.8876 / 0.0854 | 0.6980 / 0.6802 / 0.1694 |
| **`CE / SE`** | 2,473 | **0.9595 / 0.9351 / 0.0656** | 0.9532 / 0.9277 / 0.0698 | 0.9462 / 0.9152 / 0.0793 | 0.9152 / 0.8902 / 0.0886 | 0.9477 / 0.9104 / 0.0762 | 0.7015 / 0.6803 / 0.1695 |
| **`AA / A3`** | 860 | **0.9528 / 0.9459 / 0.0622** | 0.9469 / 0.9401 / 0.0683 | 0.8664 / 0.8812 / 0.1452 | 0.8534 / 0.8430 / 0.0998 | 0.8676 / 0.8450 / 0.0949 | 0.6958 / 0.6805 / 0.1672 |
| **`AD / A5`** | 175 | 0.9532 / **0.9482 / 0.0676** | **0.9584 / 0.9510** / 0.0677 | 0.8982 / 0.8952 / 0.1090 | 0.7622 / 0.7436 / 0.1581 | 0.8082 / 0.7733 / 0.1299 | 0.7290 / 0.7230 / 0.1626 |
| **`RI`** | 164 | **0.9375 / 0.9183 / 0.0742** | 0.9335 / 0.9105 / 0.0781 | 0.9160 / 0.8978 / 0.0903 | 0.8689 / 0.8334 / 0.0954 | 0.8312 / 0.8212 / 0.1270 | 0.6297 / 0.6120 / 0.1876 |
| **Binary ($K_0 / K_1$)** *(identical estimand)* | 1,635 | **0.9540 / 0.9365 / 0.0693** | 0.9492 / 0.9330 / 0.0734 | 0.9198 / 0.9064 / 0.0976 | 0.9431 / 0.9259 / 0.0718 | 0.9499 / 0.9266 / 0.0719 | 0.7403 / 0.7297 / 0.1560 |
| **Complex ($K_n \ge 2$)** *(node vs. LSV/binary)* | 2,037 | **0.9604 / 0.9445 / 0.0621** | 0.9551 / 0.9376 / 0.0668 | 0.9181 / 0.9039 / 0.0958 | 0.8569 / 0.8379 / 0.1133 | 0.8940 / 0.8549 / 0.0963 | 0.6643 / 0.6390 / 0.1802 |

**Exact-event sensitivity checks (`MAJIQ v3`, `SUPPA2`, and `rMATS` under their own exact PacBio estimands):**

| Sensitivity Check (Polymorphic, Exact Competitor-Event Truth) | $n$ | **Langur ($r$ / $\rho$ / MAE)** | `Whippet2-turbo` | `Whippet.jl` | Target Baseline ($r$ / $\rho$ / MAE) |
|---|---:|---|---|---|---|
| **vs. `STAR + MAJIQ v3`** — All annotated SJ LSVs (`6,194` LSVs, `.jnc` ratio) | 12,899 | **0.9602 / 0.9445 / 0.0645** | 0.9458 / 0.9289 / 0.0747 | — | **0.9603 / 0.9455 / 0.0632** *(`MAJIQ v3`)* |
| **vs. `STAR + MAJIQ v3`** — Binary SJ LSVs (`K = 2`, `4,278` LSVs, `.jnc` ratio) | 8,556 | **0.9633** / 0.9470 / 0.0637 | 0.9484 / 0.9302 / 0.0741 | — | 0.9629 / **0.9477 / 0.0622** *(`MAJIQ v3`)* |
| **vs. `STAR + MAJIQ v3`** — Complex SJ LSVs (`K ≥ 3`, `1,916` LSVs, `.jnc` ratio) | 4,343 | **0.9525 / 0.9350 / 0.0661** | 0.9389 / 0.9202 / 0.0761 | — | **0.9538 / 0.9365 / 0.0651** *(`MAJIQ v3`)* |
| **vs. `Salmon 2 + SUPPA2`** — Binary ($K_0/K_1$, identical estimand) | 2,339 | **0.9403 / 0.9232 / 0.0734** | 0.9349 / 0.9193 / 0.0773 | 0.8884 / 0.8778 / 0.1083 | 0.7393 / 0.7299 / 0.1563 *(`SUPPA2`)* |
| **vs. `Salmon 2 + SUPPA2`** — Complex ($K_n \ge 2$) | 2,971 | **0.8180 / 0.8062 / 0.1179** | 0.8168 / 0.8053 / 0.1207 | 0.7793 / 0.7706 / 0.1478 | 0.7262 / 0.7203 / 0.1602 *(`SUPPA2`)* |
| **vs. `Salmon 2 + SUPPA2`** — All 6 motifs (`SE/A3/A5/RI/AF/AL`) | 5,310 | **0.8718 / 0.8589 / 0.0983** | 0.8686 / 0.8565 / 0.1016 | 0.8281 / 0.8194 / 0.1304 | 0.7322 / 0.7255 / 0.1585 *(`SUPPA2`)* |
| **vs. `STAR + rMATS JC`** — Binary ($K_0/K_1$, identical estimand) | 2,747 | **0.9398 / 0.8728 / 0.0692** | 0.9356 / 0.8712 / 0.0726 | 0.9094 / 0.8514 / 0.0931 | **0.9507 / 0.8892 / 0.0665** *(`rMATS JC`)* |
| **vs. `STAR + rMATS JC`** — Complex ($K_n \ge 2$, different estimands) | 3,744 | 0.7942 / 0.7190 / 0.1216 | 0.7906 / 0.7158 / 0.1245 | 0.7616 / 0.6979 / 0.1449 | 0.8732 / 0.8110 / 0.0979 *(`rMATS JC`)* |

*Reference intersections: 5-way excluding MAJIQ (`n = 4,693` poly, node truth): Langur `0.9507 / 0.9350 / 0.0688` vs. W2 `0.9456 / 0.0728`, `STAR + rMATS JC` `0.9099 / 0.0896` (`Complex K_n >= 2`, `n = 2,685`: `0.9551 / 0.0648` vs. `0.8863 / 0.0993`, -35% MAE), `Salmon 2 + SUPPA2` `0.6832 / 0.1747`. All `16,591` polymorphic `MAJIQ v3` LSVs (`35,917` connections incl. IR & `9,254` de novo connections): `0.9491 / 0.9272 / 0.0689` (`Splice junctions`: `0.9535 / 0.0670`; `IR`: `0.8466 / 0.0867`; `De novo`: `0.8523 / 0.0702`). 4-way alignment-free across all 6 motifs (`n = 6,436` poly, node truth): Langur `0.9408 / 0.9257 / 0.0749` vs. `Salmon 2 + SUPPA2` `0.6515 / 0.6379 / 0.1863`.*

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
3. **Parity with `STAR + MAJIQ v3.0.25` and `STAR + rMATS-turbo` at order-of-magnitude lower cost:** Even on MAJIQ v3's own `6,194` polymorphic annotated splice-junction LSVs (`12,899` connections, incl. `4,343` complex `K ≥ 3`), alignment-free Langur `.jnc.gz` ratios match `STAR + MAJIQ v3` head-to-head (`r = 0.9602 / MAE = 0.0645` vs. `0.9603 / 0.0632`, `ΔMAE = 0.0013`), while beating both `STAR + MAJIQ v3` and `STAR + rMATS` on shared binary CSG events (`n = 1,635`: `0.9540 / 0.0693`) and `RI` (`0.9375 / 0.0742`). On complex hubs, both CSG nodes (Langur `r = 0.9604 / 0.0621`) and LSV junctions (MAJIQ v3 `r = 0.9538 / 0.0651`) excel on their respective estimands, with Langur building the index **46.2×–47.3× faster** (`28×` smaller), running **5.1×–7.24× faster** with **23.5× less RAM** and **zero BAM files**, and quantifying `AF`/`AL`, tandem UTRs, and TPMs in one pass (while `STAR + MAJIQ v3` remains uniquely suited for BAM-direct *de novo* junction discovery without `index --bam`).

### 7.2 Real-Cohort Null Calibration & Multivariable Batch Adjustment: Geuvadis (`ERP001942`)

Following MAJIQ V2 (Vaquero-Garcia et al., 2023), we evaluated empirical null calibration and real-world batch-confounding adjustment across **14 real human CEU LCL donors** from Geuvadis (Lappalainen et al., 2013; `ERP001942`, **353.4M PE pairs = 706.9M 75 bp reads**, GRCh38 GENCODE v38, 16T, `--biascorrect`), intentionally imbalanced across two sequencing centers (`UNIGE`: `6F + 2M`; `ICMB`: `2F + 4M`).

#### A. 14-Sample Cohort Runtime, Memory & Alignment
| Tool (14 Geuvadis Donors, GRCh38) | 14-Sample `quant` | Throughput | Peak Footprint / RSS | Mapped Pairs (of 353.4M) | 6-Sample `delta` *(10-split mean)* | 14-Sample `delta` *(`~ center + sex`)* | `delta` Peak Footprint |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Langur** *(Rust, 16T, `--batch`)* | **347.3 s** *(5m 47s)* | **2.03M/s** | **1.39 GB** / 2.92 GB | **86.05%** *(304.2M; 74.4–93.2%)* | **1.93 s** | **7.02 s** *(4.33 s unadj.)* | **177 MB** |
| **`Whippet2-turbo`** *(Julia 1.9, 16T, `--batch`)* | 429.2 s *(7m 09s)* | 1.65M/s | 9.88 GB / 10.13 GB | 73.66% *(260.3M; 65.1–80.4%)* | 20.77 s *(10.8×)* | 46.51 s *(30.93 s unadj.)* | 2.86 GB *(16.2×)* |
| **Original `Whippet.jl` v1.6.2** *(Julia 1.6, 1T)*¹ | **11,145.4 s** *(3h 06m serial)*² | **63.4k/s** *(Solo D1: 52.8k/s)* | **6.30 GB** / 6.37 GB *(Solo D1: 5.53 GB)* | **73.66%** *(260.3M; exact W2 parity)* | **57.74 s** *(29.9×)* | **104.85 s** *(unadj. only)* | **1.23 GB** *(6-samp: 660 MB)* |

*¹ Original `Whippet.jl` v1.6.2 `whippet-quant.jl` fails on raw SRA FASTQs (`ERR188086_1.fastq.gz`) with `Cannot encode 78` (`timbitz/Whippet.jl#135`); stripping line-3 `+ERR...` headers to bare `+` (without modifying any sequences or quality scores) resolved the error cleanly across all 14 samples (whereas `Whippet2-turbo` and `Langur` parse raw SRA FASTQs natively). ² Donor 1 (`ERR188086`) was timed strictly solo (`928.39 s` = `15m 28s`, `52.8k reads/s`, `5.53 GB` footprint / `5.26 GB` RSS, `71.32%` mapped); the remaining 13 donors (`685.5–907.7 s` each at 4× concurrency) summed to `11,145.44 s` (`3.10 h`) serial wall time, mapping the exact same `260,345,641 / 353,445,768` pairs as `Whippet2-turbo`.*

#### B. Real-Data Same-Condition Null Calibration (10 Exact Distinct 3v3 Splits, 6 Female UNIGE Donors)
All $\binom{6}{3}/2 = 10$ unique 3v3 partitions of the 6 same-sex, same-center donors (`-r 5 -s 2`, each tool's own native `.psi.gz` outputs):

| Method | Event Set (mean $n$) | FPR @ $p \le 0.05$ *(mean / med)* | FPR @ $p \le 0.01$ *(mean / med)* | FP @ $\text{FDR} \le 0.05$ *(mean / med)* | Per-Split FPs (`10` splits) |
|---|---:|---:|---:|---:|---|
| **Langur (`Inf_Var`)** | All (`147,196`) | **1.38% / 1.33%** | **0.22% / 0.21%** | **0.2 / 0.0** | `[0,2,0,0,0,0,0,0,0,0]` (`2` total) |
| **Langur (`Inf_Var`, Spliced)** | Spliced (`105,189`) | **0.70% / 0.70%** | **0.04% / 0.04%** | **0.2 / 0.0** | `[0,2,0,0,0,0,0,0,0,0]` (`2` total) |
| Langur (`--ignore-inf-var`) | All / Spliced | 1.42% / 0.76% *(mean)* | 0.22% / 0.04% *(mean)* | 0.3 / 0.0 *(both)* | `[0,3,0,0,0,0,0,0,0,0]` (`3` total) |
| `Whippet2-turbo` (Strat. GLM) | All (`144,093`) | 1.44% / 1.38% | 0.22% / 0.19% | 0.8 / 0.0 | `[0,0,0,0,0,0,0,0,8,0]` (`8` total) |
| `Whippet2-turbo` (Spliced) | Spliced (`102,539`) | 0.83% / 0.84% | 0.05% / 0.06% | 0.3 / 0.0 | `[0,0,0,0,0,0,0,0,3,0]` (`3` total) |
| Original `Whippet.jl` v1.6.2 (`Prob ≥ 0.95`) | All (`136,507`) | 1.11% *(call rate)* | — | **1,508.9 / 1,402.0** | `1,014–2,494` (`92.0%` `TS`/`TE`; `120.1` Spliced) |
| Original `Whippet.jl` v1.6.2 (`P ≥ 0.95 & \|ΔΨ\| ≥ 0.1`) | All (`136,507`) | 0.88% *(call rate)* | — | **1,202.9 / 1,099.5** | `851–1,996` (`63.4` Spliced) |

#### C. Multivariable Batch Confounding (`~ center + sex`, 14 Imbalanced CEU Donors)
| Model | Spliced Median $\hat{\phi}$ | Pure Center Artifacts ($n = 87$) `med \|ΔΨ\|` / `FDR ≤ 0.05` Leaks | `chrX`/`chrY` Hits @ $p \le 0.01$ *(or `Prob ≥ 0.95`)* *(All / Spliced)* | `chrX`/`chrY` Enrichment vs. Autosomes *(All / Spliced)* | Total Calls (`FDR ≤ 0.05` or `Prob ≥ 0.95`) |
|---|---:|---:|---:|---:|---:|
| **Langur unadjusted (`M vs F`)** | `0.0383` | `0.080` / `0` | `11` / `4` | `1.11×` / `1.21×` | `0` (`0` Spliced) |
| **Langur adjusted (`~ center + sex`)** | **`0.0323`** *(-15.7%)* | **`0.052`** *(-35%)* / **`0`** | **`22` / `8`** *(2×)* | **`1.44×` / `1.67×`** | `13` (`0` Spliced)³ |
| `Whippet2-turbo` unadjusted (`M vs F`) | `0.0417` | `0.055` / `0` | `17` / `4` | `1.61×` / `1.08×` | `0` (`0` Spliced) |
| `Whippet2-turbo` adjusted (`~ center + sex`) | `0.0359` *(-13.9%)* | `0.045` / `0` | `27` / `6` | `1.70×` / `1.09×` | `2` (`0` Spliced) |
| Original `Whippet.jl` v1.6.2 (`Prob ≥ 0.95`) | — *(no $\phi$)* | `0.047` / `0` | `2` / `0` *(0.04%)* | `0.68×` / `0.00×` *(depleted)* | **90** *(88 auto., 87 `TS`/`TE`, 3 Spliced)* |
| Original `Whippet.jl` v1.6.2 (`P ≥ 0.95 & \|ΔΨ\| ≥ 0.1`) | — *(no $\phi$)* | — | `2` / `0` *(0.04%)* | `0.79×` / `0.00×` *(depleted)* | **78** *(76 auto., 77 `TS`/`TE`, 1 Spliced)* |

**Key findings:**
1. **Cohort speed, memory & alignment recovery:** Langur quantifies all 14 human donors **32.1× faster** (`5m 47s` vs. `3h 06m` serial) at **4.5× lower peak RAM** (`1.39 GB` vs. `6.30 GB`) than original `Whippet.jl` v1.6.2—and **1.24× faster** at **7.1× lower RAM** than `Whippet2-turbo`—while recovering **+43.8M** more mapped pairs (`+12.4 pp`) than both, and runs 14-sample multivariable `delta` **14.9× faster** (`7.0×` less RAM) than `Whippet.jl` v1.6.2 and **6.6× faster** (`16.2×` less RAM) than `Whippet2-turbo`.
2. **>6,000× false-positive reduction on real human null splits:** Langur (`Inf_Var`) produces **zero false positives in 9/10 splits** (`0.2` mean / `0.0` median vs. `1,508.9` for native `Whippet.jl` v1.6.2 at `Prob ≥ 0.95` and `1,202.9` at `Prob ≥ 0.95 & |ΔΨ| ≥ 0.1` — a **7,545× / 6,015×** reduction overall, **600× / 317×** on Spliced events, and **4×** fewer than `Whippet2-turbo`'s `0.8`), with a well-calibrated slightly conservative polymorphic spliced null tail (`|ΔΨ| ≥ 0.01`: `FPR@0.05 = 2.83%`, `q5 = 0.069`, `q10 = 0.101`) tightened by `Inf_Var` weighting over `--ignore-inf-var`.
3. **Batch debiasing & `chrX`/`chrY` signal recovery:** Adjusting for `~ center + sex` cuts residual spliced overdispersion $\hat\phi$ by **15.7%** (`0.0383 → 0.0323`; TS/TE `0.2638 → 0.2429`), reduces apparent sex `|ΔΨ|` bias across `87` pure concordant autosomal sequencing-center artifacts by **35%** (`0.080 → 0.052`, `0` leaks at `FDR ≤ 0.05`; e.g., center-driven `ZNF302` CE/AA/RI `ΔΨ_center = +0.671`, `FDR = 0.022` drops from unadjusted sex bias `-0.155` to `+0.067`, `FDR = 1.000`), and **doubles** `chrX`/`chrY` hits at $p \le 0.01$ (`11 → 22` All; `4 → 8` Spliced, enrichment `1.21× → 1.67×`, recovering X-linked `UBE2A`, `KDM6A`, `ENOX2`, `CXorf38`, `TRMT2B`, `ARHGEF6`, `HEPH`, `CDK16`, `HDAC6`), whereas native `Whippet.jl` v1.6.2 (unable to adjust for `center`) calls `90` unadjusted events that are `96.7%` `TS`/`TE` noise (`88` autosomal, `0.68×` depleted on `chrX`/`chrY`, `0` spliced `chrX`/`chrY`). *³ All `13` `FDR ≤ 0.05` events in Langur `~ center + sex` (`0` Spliced) are single-stratum (2-donor) `UNIGE_M` or `ICMB_F` tandem UTR/TSS (`TS`/`TE`) corner outliers (e.g., `SLC35B2`, `BANP`, `MDM2`, `TUBB3`), where an additive main-effects logit model without a `center:sex` interaction term simultaneously attributes both `center` and `sex` effects.*

### 7.3 Real RF-Stranded Differential Splicing & Orthogonal eCLIP Validation (`ENCODE PTBP1` KD, 5-Way)

To benchmark real biological differential splicing, strand-specific antisense separation, and orthogonal RBP-binding validation on GRCh38 GENCODE v38 basic, we ran all five modern pipelines on an RF-stranded ENCODE K562 `PTBP1` shRNA knockdown cohort (`ENCSR527IVX` vs. matched non-target control `ENCSR129RWD`; 2 vs. 2 replicates, **127.6M PE pairs = 255.3M 100 bp reads**, 16T, `/usr/bin/time -l`), aligned to identical directionality (`ΔΨ = Ψ_KD - Ψ_Ctrl`, where `ΔΨ > 0` = de-repression of `PTBP1`-repressed exons). Orthogonal validation used official **ENCODE K562 `PTBP1` eCLIP IDR peaks** (`ENCFF907HNN`, `7,797` peaks) in canonical upstream repression zone `[-300, 0]` and flanking `±300 bp`, plus strand-aware upstream `[-250, -20]` CU-rich PTBP1 tetramer density (`TCTC`/`CTCT`/`TCTT`/`TTCT` per 100 bp). Stage 1 TPM controls confirmed `PTBP1` knockdown (`3.21×` in Langur `157.0 → 48.9 TPM`; `3.03×` W2; `3.34×` Salmon) and canonical `PTBP2` NMD escape (`8.26×` in Langur `11.5 → 94.8 TPM`; `8.07×` W2; `8.67×` Salmon).

#### A. 5-Way Runtime, Peak Memory & Disk (`PTBP1` KD `2 vs. 2`, 255.3M Reads, 16T)
| Pipeline | Strand Mode | Index / Setup | 4-Sample Quant / Align | DS Testing | End-to-End Total | Peak Quant/Align RAM (Foot / RSS) | Index / BAM Disk |
|---|---|---:|---:|---:|---:|---:|---:|
| **`Langur`** (`--batch`) | `--stranded-reverse` | **14.8 s** | **180.8 s** *(3m 01s)* | **1.41 s** | **197.0 s** *(3m 17s)* | **1.35** / 3.02 GB | 1.06 GB / **0** |
| `Whippet2-turbo` (`--batch`) | unstranded¹ | 28.6 s | 196.8 s *(3m 17s)* | 20.7 s | 246.1 s *(4m 06s)* | 8.80 / 9.08 GB | 691 MB / 0 |
| `Salmon 2.12 + SUPPA2 2.4` | `-l ISR` | 130.7 s | 255.8 s *(4m 16s)* | 17.4 s | 403.9 s *(6m 44s)* | 12.13 / 12.23 GB | 7.9 GB / 0 |
| `STAR 2.7.11b + MAJIQ v3.0.25` | `--strandness REVERSE` | 697.9 s | 645.5 s *(10m 45s)* | 47.7 s² | 1,391.1 s *(23m 11s)* | **31.64** / 31.57 GB | **28 GB / 15 GB** |
| `STAR 2.7.11b + rMATS-turbo v4.3.0` | `fr-firststrand` | 692.6 s | 645.5 s *(10m 45s)* | 158.2 s | 1,496.3 s *(24m 56s)* | **31.64** / 31.57 GB | **28 GB / 15 GB** |

*¹ `Whippet2-turbo` supports only forward `--stranded`, not `--stranded-reverse`; running unstranded drops `12.2M` mapped pairs (`63.2%` vs. Langur `72.8%`), bleeds sense genes into overlapping antisense partners (`TUBA1C`: `897 → 568 TPM`, `UBB`: `1,048 → 410 TPM`), and creates `201` phantom antisense-bleed genes (`Langur < 1 TPM` vs. `W2 > 20 TPM`). In a 1M-pair diagnostic slice, Langur `--stranded-reverse` (`73.1%`) collapses **18.7×** under wrong-strand `--stranded` (`3.9%`) and eliminates antisense bleed (`TUBA1C`/`RP11-161H23.5`: `100.9 → 1.0`, `UBB`/`RP11-138I1.4`: `83.8 → 2.0`, `RBBP4`/`SYNC`: `26.7 → 0.0`, `PIEZO1`/`RP5-1142A6.9`: `20.0 → 0.0`, `FDPS`/`RUSC1-AS1`: `18.2 → 0.0`). ² `MAJIQ v3` DS = `sj` (`31.8 s`) + `build` (`1.53 s`) + 4× per-sample `psi-coverage` (`7.38 s`) + `deltapsi` (`7.00 s`); `majiq-v3 heterogen` (`7.62 s`) called `0` events at `BH FDR ≤ 0.05`.*

#### B. Default Operating Points & Threshold-Free Matched Top-100 Orthogonal Validation
| Method & Operating Point | Sig Events *(CE/SE · `% ΔΨ>0`)* | Default `ΔΨ>0` Up CU `/100bp` *(Fold, MW $p$)* | Default `ΔΨ>0` eCLIP Up / Flank | Top-100 Up CU `/100bp` *(Fold, MW $p$)* | Top-100 eCLIP Up / Flank³ | Top-100 Consensus (`≥3/5`, `n=55`) | Default Consensus Precision |
|---|---:|---|---|---|---|---:|---:|
| **`Langur`** (`FDR ≤ 0.05`) | **58** *(45 · **73.3%**; `0` TS/TE)* | **5.69** (**2.22×**, `1.1e-13`) | **6.1% / 6.1%** *(0% `ΔΨ<0`)* | **4.53** (**1.77×**, `8.9e-18`) | **6.0% / 8.0%** *(OR `24.1`)* | **53/55 (96.4%)** | **84.0%** (`42/50`) |
| `Whippet2-turbo` (`FDR ≤ 0.05`, unstr.) | 70 *(55 · 74.5%)* | 5.66 (2.21×, `1.3e-16`) | 7.3% / 7.3% | 4.62 (1.81×, `6.5e-16`) | 6.0% / 7.0% *(OR `23.4`)* | 53/55 (96.4%) | 88.1% (`52/59`) |
| `STAR + MAJIQ v3` (`P(\|ΔΨ\|≥0.20)≥0.95`) | 359 *(149 · 75.2%)* | 4.37 (1.74×, `4.8e-17`) | 5.4% / 8.0% | 4.21 (1.68×, `2.1e-14`) | 5.0% / 7.0% *(OR `8.8`)* | 48/55 (87.3%) | 28.4% *(at `≥0.10`)* |
| `STAR + rMATS-turbo JC` (`FDR≤0.05, \|ΔΨ\|≥0.1`) | 2,130 *(1,489 · **49.7%**)* | 3.55 (1.39×, `1.8e-31`) | 2.0% / 2.7% | 4.24 (1.66×, `1.4e-10`) | 5.0% / 10.0% *(OR `11.9`)* | 45/55 (81.8%) | **10.0%** (`147/1463`) |
| `Salmon 2 + SUPPA2` (`p≤0.05, \|ΔΨ\|≥0.1`) | 1,518 *(578 · 57.6%)* | 3.28 (1.22×, `5.3e-5`) | 0.6% / 1.5% *(n.s.)* | 2.78 (1.04×, `0.92` n.s.) | 1.0% / 3.0% *(n.s.)* | **3/55 (5.5%)** | **12.5%** (`128/1024`) |

*³ Background note: upstream CU-motif baselines are nearly identical across all 5 tools (`2.51–2.68/100bp`), and `ΔΨ < 0` negative controls are `0.91×–1.07×` (n.s., $p > 0.47$) in every tool. Tool-matched expressed background eCLIP Up `[-300,0]` hit rates reflect each tool's tested exon universe (`0.3%` across `67,709` CSG exons in Langur/W2 vs. `0.6%` across `11,440` alternative LSV exons in MAJIQ v3 vs. `0.4%` in rMATS/SUPPA2); raw hit rates (`6.0%` Langur Top-100 vs. `5.0%` MAJIQ/rMATS vs. `1.0%` SUPPA2) are reported alongside ORs for full transparency. Langur `delta` dispersion split: `spliced median φ = 0.0296 (d0 = 25.0)` vs. `TS/TE median φ = 0.1305 (d0 = 1.0)` (`4.4×`), yielding `58/58` Spliced hits and `0` spurious `TS`/`TE` calls.*

#### C. Canonical RT-PCR-Validated `PTBP1` Targets & `MAP4K4` 9 bp Microexon Recovery
| Target Gene & Exon (GRCh38) | **Langur** (`ΔΨ`, rank, `FDR`) | `Whippet2-turbo` | `STAR + MAJIQ v3` (`P(\|ΔΨ\|≥0.1)`) | `STAR + rMATS JC` (`FDR`) | `Salmon + SUPPA2` (`p`) |
|---|---|---|---|---|---|
| **`MAP4K4` 9 bp microexon** (`chr2:101868029-101868037`) | **`+0.742` (`#1`, `8.6e-6`)** | `+0.720` (`#2`, `2.0e-5`) | `+0.709` (`#9`, `1.000`) | `+0.760` (`#7`, `0.0`) | `+0.370` (`#272`, `2.8e-3`) |
| **`PTBP2` ex10 NMD** (`34 bp`, `chr1:96806419-96806452`) | **`+0.714` (`#30`, `4.6e-3`)** | `+0.664` (`#38`, `1.2e-2`) | *filtered* (`passed=0`)⁴ | `+0.679` (`#92`, `1.1e-13`) | *missed as SE* |
| **`SMAP1`** (`81 bp`, `chr6:70791689-70791769`) | **`+0.873` (`#3`, `8.8e-6`)** | `+0.833` (`#6`, `1.2e-4`) | `+0.916` (`#8`, `1.000`) | `+0.903` (`#1`, `0.0`) | `+0.421` (`#237`, `2.0e-3`) |
| **`GABBR1`** (`151 bp`, `chr6:29609229-29609379`) | **`+0.778` (`#5`, `8.8e-6`)** | `+0.791` (`#4`, `2.0e-5`) | `+0.689` (`#6`, `1.000`) | `+0.778` (`#6`, `0.0`) | `+0.087` (`#3085`, `0.089` n.s.) |
| **`RTN4`** (`2400 bp`, `chr2:55025086-55027485`) | **`+0.619` (`#9`, `1.2e-4`)** | `+0.485` (`#9`, `1.4e-3`) | `+0.818` (`#23`, `1.000`) | `+0.673` (`#10`, `0.0`) | `+0.906` (`#109` `AF`, `0.0`) |
| **`CTNND1`** (`119 bp`, `chr11:57789037-57789155`) | **`+0.671` (`#12`, `7.4e-4`)** | `+0.627` (`#10`, `1.9e-3`) | `+0.646` (`#35`, `1.000`) | `+0.694` (`#8` adj., `0.0`) | `-0.565` (`#6954` `AF`, `0.17` n.s.) |
| **`TPM1`** (`76 bp`, `chr15:63061713-63061788`) | **`+0.596` (`#14`, `1.4e-3`)** | `+0.560` (`#29`, `8.4e-3`) | `+0.646` (`#17`, `1.000`) | `+0.523` (`#24`, `0.0`) | `-0.365` (`#710` `MX`, `0.016`) |
| **`PICALM`** (`129 bp`, `chr11:85990250-85990378`) | **`+0.499` (`#24`, `2.1e-3`)** | `+0.502` (`#22`, `4.7e-3`) | `+0.453` (`#43`, `1.000`) | `+0.587` (`#17`, `0.0`) | `+0.371` (`#149`, `0.0`) |

*⁴ `MAJIQ v3` filtered out `PTBP2` exon 10 (`passed = 0`) because NMD in control K562 cells reduced exon-10 inclusion reads to `2.0` and `8.0` in control replicates (`< quantify_minreads = 10`) vs. `340.0` and `382.0` in KD. Pairwise dynamic cassette-exon (`|ΔΨ| ≥ 0.1`) concordance: Langur vs. Whippet2 `r = 0.888` (`n = 3,740`), Langur vs. MAJIQ v3 `r = 0.769` (`n = 1,932`), Langur vs. rMATS `r = 0.583–0.613` (`n = 4,283–4,500`), Langur vs. SUPPA2 `r = 0.471` (`n = 1,832`).*

**Additional production robustness checks:**
- **Single-end real human WTC11 (`ENCFF766OAK_R1.fastq.gz`, `45.8M` reads, GRCh38, `-t 16 --biascorrect`):** Completes in **`76.8 s`** at **`0.97 GB` peak RAM** (`67.4%` mapped, `717k` decoy-filtered), achieving polymorphic PacBio HiFi concordance **`r = 0.9234 / ρ = 0.8618 / MAE = 0.0875`** (`n = 12,670`, `LR ≥ 20, SR ≥ 10`; all-matched `n = 112,048`: `r = 0.9757 / MAE = 0.0181`) and **`r = 0.9799`** direct SE-vs-PE concordance (`n = 12,600`).
- **Variable-length messy FASTQ stress test:** `200k` uncompressed post-`fastp` pairs mixing `10/16/17/25/50/75/100 bp` reads, all-`'N'` reads, lowercase bases, and custom `+header` lines complete in `1.25 s` (PE) / `1.22 s` (SE) with **`0` NaNs / `0` Infs** across all `.psi.gz`/`.tpm.gz` tables and high concordance vs. full WTC11 PE (`r = 0.9929` PE / `0.9779` SE).
- **Real ENCODE `index --bam` supplementation:** Full GRCh38 `index --bam` (`ENCFF436JMY` `chr21`/`chr22` slice, `--bam-min-reads 3`, incl. `ERCC-*` contigs and secondary/supplementary records) builds in **`19.94 s`** with zero panics (`+615` CSG nodes: `417,213 → 417,828`), passes `inspect` `mmap + CRC32` in **`0.038 ms`** (`bam_novel=true`), and yields **`389` supported new/updated PSI rows** (of `799`) plus **`100` `ALT_UNIQ` novel junctions** (`117.8` reads) in `500k` pairs (`3.34 s` 4-sample `quant`, `0.58 s` `delta`).
- **Cross-platform CI:** Automated GitHub Actions workflow (`.github/workflows/ci.yml`) runs `cargo test --workspace` (`42/42` tests), `cargo clippy --workspace --all-targets -- -D warnings` (`0` warnings), and `cargo build --release` on both `ubuntu-latest` (`x86_64`) and `macos-latest` (`arm64`).

### 7.4 Automated Verification Suite & Synthetic Cohort
- **42/42 tests pass** (`cargo test --workspace`, zero warnings), covering Whippet 2 CSG/motif parity (`ENSMUSG00000038685` 18-path, `nrf2test` 6-path), bounded bridge sweeps (`LTBP2`), FR/RF/unstranded pairing, pseudogene decoy ties vs. diverged decoys, 1-SNP splice jumps, banded Gotoh indels, CC-VBEM paralog routing (`>99.5%`), 15 bp microexon PSI, PCR-stack soft-capping, Gibbs ~90% CI coverage across 200 $K_2$ replicates, constitutive exclusion + stratified EB priors, multivariable batch/interaction/continuous/collinear designs, and 100% thread determinism.
- **200-gene 3v3 synthetic cohort (600k reads, 4T):** `index` `0.032 s` (185× vs. W2), 6-sample `--batch` `0.377 s` (126× faster, 17× less RAM), `delta` `0.039 s` (122× faster, 95× less RAM); **100% TPR (100/100)** at **0.00% FDR** ($\Delta\Psi\ r = 0.9935$).

---

## 8. Key References

1. Sterne-Weiler T., Weatheritt R.J., Best A.J., Ha K.C.H., Blencowe B.J. (2018). Efficient and accurate quantitative profiling of alternative splicing patterns of any complexity on a laptop. *Molecular Cell*, 72(1):187–200.e6.
2. Srivastava A., Malik L., Sarkar H., Zakeri M., Almodaresi F., Soneson C., Love M.I., Kingsford C., Patro R. (2020). Alignment and mapping methodology influence transcript abundance estimation. *Genome Biology*, 21:239.
3. Patro R., Duggal G., Love M.I., Irizarry R.A., Kingsford C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14:417–419.
4. Vaquero-Garcia J., Aicher J.K., Jewell S., et al. (2023). RNA splicing analysis using heterogeneous and large RNA-seq datasets. *Nature Communications*, 14:1230.
5. Aicher J.K., Slaff B., Jewell S., Gazzara M.R., Bhoj E.J.K., Barash Y. (2026). MAJIQ V3: Scalable and accurate RNA splicing analysis with splicegraph algebra and improved intron retention modeling. *NAR Genomics and Bioinformatics*, 8(4):lqag126.
6. Kabza M., Ritter A., Byrne A., Sereti K., Le D., Stephenson W., Sterne-Weiler T. (2024). Accurate long-read transcript discovery and quantification at single-cell, pseudo-bulk and bulk resolution with Isosceles. *Nature Communications*, 15:7316.
7. McCarthy D.J., Chen Y., Smyth G.K. (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Research*, 40(10):4288–4297.
8. Varadhan R., Roland C. (2008). Simple and globally convergent methods for accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*, 35(2):335–353.
9. Hansen K.D., Brenner S.E., Dudoit S. (2010). Biases in Illumina transcriptome sequencing caused by random hexamer priming. *Nucleic Acids Research*, 38(12):e131.
10. Pardo-Palacios F.J., Wang D., Reese F., et al. (2024). Systematic assessment of long-read RNA-seq methods for transcript identification and quantification (LRGASP). *Nature Methods*, 21:1349–1363.
11. Dobin A., Davis C.A., Schlesinger F., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1):15–21.
12. Wang Y., Xie Z., Kutschera E., Adams J.I., Kadash-Edmondson K.E., Xing Y. (2024). rMATS-turbo: an efficient and flexible computational tool for alternative splicing analysis of large-scale RNA-seq data. *Nature Protocols*, 19:1083–1104.
13. Trincado J.L., Entizne J.C., Hysenaj G., et al. (2018). SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. *Genome Biology*, 19:40.
14. Lappalainen T., Sammeth M., Friedländer M.R., et al. (2013). Transcriptome and genome sequencing uncovers functional variation in humans. *Nature*, 501:506–511.
