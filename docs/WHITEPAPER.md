# Langur: Fast, Statistically Calibrated Splice-Graph Quantification and Differential Splicing in Rust

*CLI binary: `langur` · Index format: `.lidx`*

---

## Abstract

Alignment-free quantification on **Contiguous Splice Graphs (CSGs)**, pioneered by Whippet (Sterne-Weiler et al., *Molecular Cell* 2018), made it possible to profile both binary ($K_1$) and arbitrarily complex ($K_n \ge 2$) alternative splicing events directly from FASTQ files without genome alignment or BAM generation. However, scaling CSG analysis to modern human cohorts exposes four key limitations in legacy implementations: ungapped matching drops reads near indels and splice-proximal SNPs; transcriptome-only indexing permits reads from diverged processed pseudogenes and intronic loci to force-map onto annotated genes; per-gene EM and naive Beta intervals understate ambiguity in paralog families and complex splicing hubs; and legacy two-group pooled Beta Monte Carlo testing (`whippet-delta`) cannot adjust for covariates and suffers severe false-discovery inflation under biological overdispersion.

**Langur** is a ground-up Rust re-architecture of the CSG workflow (`index`, `inspect`, `quant`, `delta`) designed for bulk short-read RNA-seq. It combines a zero-copy memory-mapped splice-graph and selective genomic decoy index (`< 1 ms` load time), hybrid ungapped + gapped spliced alignment with 1-SNP junction tolerance, joint connected-component Variational Bayes isoform quantification with unified 5'-hexamer/GC bias correction, targeted Gibbs inferential variance (`Inf_Var`) for complex $K_n \ge 2$ events, and a covariate-aware Beta-Binomial generalized linear model (`delta`) with event-class stratified Empirical Bayes dispersion shrinkage. On the 480M-read human `E-MTAB-3766` benchmark (Ensembl GRCh37.71, 20,410 genes, 16 threads, `--biascorrect`), Langur builds the human index in **22.6 s** (**15.1× faster** than original `Whippet.jl`), quantifies **80M reads in 60.4 s** (**26.6× faster** than original `Whippet.jl`, **1.34× faster** than `Whippet2-turbo`, our multithreaded Julia 1.9 fork) and the full **6-sample cohort (480M reads) in 6 min 13 s (`373.1 s`, `1.29M reads/s`)** using **~1.7 GB RAM** (**4.3×–5.9× less memory** than `Whippet.jl` / `Whippet2-turbo`), reduces gene and isoform TPM error by **5.5× / 2.8×** vs. original `Whippet.jl` (**3.2× / 1.7×** vs. `Whippet2-turbo`), and completes cohort-wide differential splicing across 223k events (`-r 5 -s 2`) in **2.5 s** while achieving **5.4% realized FDR** across all events (**2.0%** across spliced events; `0.5642 / 0.5672` AUPRC, outperforming `Whippet2-turbo` on TPR, AUROC, and AUPRC) at a 5% nominal target—a **15× reduction in false discoveries** over legacy Whippet (`81.5%` FDR) and **10.6× lower FDR** than published `rMATS` (`57.1%`).

---

## 1. Motivation, Lineage & Scope

### Limitations of Existing Splicing Workflows
In practice, RNA splicing studies typically choose between two imperfect workflows:
1. **Genome-aligned callers (`STAR` + `MAJIQ`, `rMATS`, or `LeafCutter`):** Strong splicing accuracy and de novo junction discovery, but they require heavy genome alignment, coordinate-sorted BAMs, and substantial disk/compute overhead across large cohorts.
2. **Legacy alignment-free CSG quantification (`Whippet.jl`):** Eliminates BAM generation and quantifies complex multi-path events ($K_n \ge 2$) directly from FASTQ, but carries four methodological bottlenecks on full human cohorts:
   - **Mapping artifacts:** Pure ungapped matching loses sensitivity around small indels and SNPs adjacent to donor/acceptor sites, and indexing only annotated exons allows reads originating from diverged processed pseudogenes or unannotated intronic/intergenic regions to force-map onto parent genes.
   - **Paralog & bias decoupling:** Splitting multi-gene reads evenly prior to per-gene EM leaks signal from expressed parent genes into silent paralogs/pseudogenes, and applying GC bias after isoform EM leaves transcript TPMs uncorrected.
   - **Complex-event uncertainty & statistical calibration:** Simple Beta variance treats ambiguous multi-path ($K_n \ge 2$) events as if every read uniquely identified a single path. Downstream, legacy `whippet-delta` pools reads across replicates rather than modeling sample-to-sample biological overdispersion, leading to extreme false-positive rates and no support for confounders (batch, donor, strandedness, or continuous covariates).
   - **Runtime & memory overhead:** JIT startup, index deserialization, and garbage-collected heap footprints inflate runtime and RAM usage in multi-sample cohorts.

### Attribution, Scope & Development Transparency
- **Attribution:** Langur is an independent Rust project building on the CSG representation and motif taxonomy introduced by Tim Sterne-Weiler and colleagues (*Molecular Cell*, 2018), alongside independent implementations of published statistical concepts from `Salmon` (selective genomic decoys, Variational Bayes EM, and Gibbs inferential uncertainty) and `edgeR`/`DESeq2` (Cox–Reid Empirical Bayes dispersion shrinkage in GLMs), as well as overhang-diversity QC for PCR junction-stack artifacts highlighted in the splicing literature (e.g., `MAJIQ`).
- **Intended scope:** Purpose-built for **bulk short-read RNA-seq** (paired-end or single-end Illumina), with optional BAM-derived novel splice-site supplementation at index time (`index --bam`). For single-cell or full-length long-read (PacBio HiFi / ONT) isoform analysis, we recommend dedicated frameworks such as **Isosceles** (Kabza, Sterne-Weiler et al., *Nature Communications* 2024), **Oarfish**, or **IsoQuant**.
- **License & transparency:** Released under the **MIT License**. Its architecture, Rust codebase, tests, and benchmarks were developed iteratively with AI-assisted coding and review (Google Gemini), with every biological invariant, edge case, and benchmark audited via automated unit, parity, and full-human test suites.

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
| **Index loading** | Full deserialization into heap RAM every run (`6.2–7.3 s`) | Memory-mapped (`mmap`) zero-copy binary file (`0.04–0.3 ms`) |
| **Pseudogene/intronic noise** | Transcriptome-only index | Selective genomic decoys with strict-inequality filtering (`decoy > CSG`) |
| **Indels & boundary SNPs** | Ungapped matching only | Fast ungapped path + automatic gapped SIMD fallback + 1-SNP junction jumping |
| **Multi-gene ambiguity** | Equal pre-splitting across genes | Joint Connected-Component Variational Bayes EM across linked paralog families |
| **Complex $K_n \ge 2$ uncertainty** | Closed-form Beta variance only | Targeted collapsed Gibbs posterior variance (`Inf_Var`, col. 15 of `.psi.gz`) |
| **Differential splicing** | Pooled 2-group Beta MC (no covariates, high FDR) | Covariate-aware Beta-Binomial GLM with Spliced vs. Tandem UTR stratified EB priors |
| **Reproducibility** | Minor MC/thread drift across runs | 100% bit-for-bit identical across thread counts and standalone vs. `--batch` |

---

## 3. Zero-Copy Splice-Graph & Selective Decoy Index (`index`)

### Contiguous Splice Graph & Motif Parity
From a reference FASTA and GTF2.2 annotation (plus optional novel splice junctions imported from a strand-aware BAM via `--bam` / `--bam-min-reads`), each gene's overlapping exons are partitioned into minimal non-overlapping sub-intervals ordered $5' \to 3'$ along the gene's transcriptional strand. Every node preserves Whippet's canonical motif taxonomy: spliced events (`CE`, `AA`, `AD`, `RI`, `AF`, `AL`, optional circular `BS`) and tandem UTR ends (`TS`, `TE`).

### Minimizer Indexing Across Exons and Annotated Junctions
Instead of an FM-index, Langur indexes short sequence fingerprints (**canonical minimizers**) across both exonic nodes and the exact sequence bridges spanning every annotated splice junction. Indexing annotated junction bridges directly allows reads spanning a splice site—including reads crossing very short exons—to jump straight to the exact splice coordinate in a single lookup.

### Selective Genomic Decoys
Following Salmon's selective alignment strategy (Srivastava et al., *Genome Biology* 2020), `index` scans non-exonic genomic intervals (introns and intergenic regions) in parallel and extracts only regions that share sequence seeds with annotated exons—capturing processed pseudogenes and paralogous repeats without indexing the full 3 Gb genome.

### Instant Zero-Copy Loading (`mmap`)
Rather than parsing or deserializing the index into RAM at startup, Langur writes the index as a flat, cache-aligned binary file mapped directly into virtual memory by the OS (`mmap`). Opening a full human index takes **0.3 ms**, and `--batch` cohort runs share the same read-only mapping with zero duplication.

<details>
<summary><b>Technical Details: Minimizer Table, Flank Encoding & Binary Format</b></summary>

- **Minimizers:** Canonical odd-$k$ minimizers (default $k = 17$, window $w = 8$, invertible Wang64 hash) stored in a 2-level radix-bucketed open table (`MinBucketOffsets`, `MinKeyLo`, `MinSlots`, `MinPostings`) with inline single-hit fast-path slots (`MinimizerLookup::Inline`).
- **Annotated junction seeds:** Extracted from the $(k-1)+(k-1)$ sequence across every annotated splice edge `(src_node, dst_node)` with exact `left_overhang` $\in [1, k-1]$.
- **Precomputed node metadata:** 2-bit 16-bp acceptor and donor boundary flanks with explicit valid base counts (`0..=16`) for single-instruction XOR + popcount 1-SNP junction matching; canonical splice dinucleotide flags (`GT-AG`, `GC-AG`, `AT-AC`); cumulative transcript prefix lengths (`tx_node_prefix_len`) for $O(\log N_t)$ fragment-length queries across the $N_t$ nodes of transcript $t$; and 20-bin expected transcript GC profiles (`tx_expected_gc`).
- **Selective decoys:** Non-CSG genomic hits are expanded by `decoy_window_pad = 250 bp`, merged, split at `N`/IUPAC runs (preventing spurious poly-`A` homopolymers in 2-bit encoding), and indexed with `max_decoy_occ = 32`.
- **Container format:** 17 `#[repr(C)]` `bytemuck::Pod` sections aligned to 64-byte boundaries with magic header and `CRC32` checksum.

</details>

---

## 4. Hybrid Spliced Alignment & Decoy Screening (`quant` — Alignment)

- **Fast ungapped path with automatic gapped recovery:** Most Illumina reads contain no indels and are aligned by a fast ungapped scan. When a read exhibits a mismatch cluster inside an exon, the aligner automatically invokes a hardware-accelerated (NEON / AVX2) banded gapped aligner to recover insertions and deletions cleanly without slowing down clean reads.
- **Boundary-proximal SNPs & annotated-first splice jumping:** Simple prefix-max aligners often stop 1–3 bp short of a donor or acceptor if a SNP sits adjacent to the splice site. Langur retains near-optimal boundary extensions so splice jumps cross cleanly, tests annotated splice junctions first with **zero splice penalty** (`min_anno_overhang = 4 bp`), and falls back to 1-SNP-tolerant novel splice-flank matching (`min_novel_overhang = 8 bp`, prioritizing canonical `GT-AG`/`GC-AG`/`AT-AC` pairs).
- **Strict-inequality pseudogene immunity & zero-cost perfect hits:**
  - If a read already achieves the theoretical maximum score on the CSG, no decoy can beat it—so decoy seeding and alignment are skipped completely.
  - On imperfect CSG hits that also share seeds with a decoy, candidate decoys are aligned under the exact same scoring matrix. Because annotated splice junctions incur **zero penalty**, an intronless processed pseudogene can at most *tie* a true spliced mRNA read from its parent gene. Langur therefore filters a read as `DecoyFiltered` **only when** $S_{\text{decoy}} > S_{\text{CSG}}$: diverged pseudogene or intronic reads (carrying decoy-specific variants) are discarded, while genuine spliced parent-gene reads survive exact ties intact.
- **Strand-aware paired-end concordance:** Supports unstranded, FR (`--stranded`), and RF/TruSeq dUTP (`--stranded-reverse`) libraries, evaluating exact spliced transcriptomic fragment lengths across all compatible annotated isoforms for each concordant pair.

<details>
<summary><b>Technical Details: Scoring Matrix, Deferred Decoys & Streaming I/O</b></summary>

- **Scoring:** Match `+2`, mismatch `-4`, `'N'` `-1`, gap open `-6`, gap extend `-1`; non-canonical novel splice penalty `-2`; optional `--circ` back-splices (`dst <= src`, `min_ov = 10 bp`); boundary extension retained within `12` points of prefix peak.
- **Deferred decoy seeding:** Primary seeding sets only `has_any_decoy_posting`; decoy seeds (bucketed by `(decoy_id, strand, diag_bucket)` to prevent spurious early seeds on long contigs from shadowing true hits) are collected only if `best_csg_score < max_possible_score`.
- **Exact transcriptomic fragment length:**

  $$
  d_t = \left(\text{prefix}_t[i_{3'}] + \text{end}_{3'}\right) - \left(\text{prefix}_t[i_{5'}] + \text{start}_{5'}\right).
  $$
- **I/O:** Auto-detects `igzip -dc` (with exit-status verification and automatic `flate2`/`zlib-ng` fallback) for FASTQ decompression, and writes RFC 1952 multi-member gzip outputs in parallel.

</details>

---

## 5. Two-Stage Quantification: Expression, PSI & Inferential Variance (`quant`)

Reads are processed in parallel waves using compact, sparse thread-local counters merged in strict file order, guaranteeing **100% bit-for-bit identical outputs** across any thread count and keeping peak RAM around **~1.7 GB** on full human libraries.

### Stage 1: Transcript & Gene Abundance (Global Routing)
1. **Post-pass fragment length distribution (FLD):** After all reads are aligned, Langur fits one smoothed sample-wide fragment-length curve so transcript effective lengths are exact and independent of chunk processing order.
2. **Joint Connected-Component Variational Bayes EM (CC-VBEM):** Transcripts that share multi-mapping reads—including paralog families and annotated pseudogenes—are grouped into connected components and solved jointly using sparse Variational Bayes EM (similar to Salmon's VBEM). Because expressed parent genes accumulate unique reads in divergent exons/UTRs while silent paralogs/pseudogenes do not, the sparse prior cleanly routes ambiguous reads back to the expressed gene.
3. **Unified 5'-hexamer and 20-bin GC bias correction (`--biascorrect`):**
   - 5'-random-hexamer priming weights are fitted deterministically from Wave 0 (~1.6M read pairs), immediately freeing per-feature hexamer lists for the rest of the library.
   - Pass 1 runs CC-VBEM on primer-weighted counts to estimate initial transcript abundances, which are then used to fit 20-bin expected-vs-observed GC bias. Pass 2 folds GC bias directly into transcript effective lengths (preventing GC double-counting) and re-solves CC-VBEM—improving both `.isoform.tpm.gz` and `.gene.tpm.gz`.

### Stage 2: Local Splice-Graph PSI, Artifact Control & Gibbs `Inf_Var`
Using converged Stage 1 posteriors to route ambiguous reads into each gene's splice graph:
- **Microexon-safe effective lengths:** Internal splice junctions retain their full mappable overhang span (clamping by single-node length only at true 5' or 3' transcript ends), preventing artificial PSI inflation on short internal microexons (e.g., 15 bp `CE`).
- **Overhang-diversity PCR-stack soft-capping:** Each junction tracks its distinct overhang start positions (`Distinct_Pos`) and maximum pileup (`Max_Stack`) in `.jnc.gz`. On **novel or `RI`-touching junctions only** (never annotated junctions), if $\ge 10$ reads collapse onto a single start coordinate despite room for multiple overhangs, the count is softly capped to `2.0` (`--stack-cap`, disable via `--no-stack-filter`).
- **Targeted collapsed Gibbs inferential variance (`Inf_Var`):**
  - Binary ($K_1$) and unambiguous events (~80–85% of events) retain exact closed-form Beta credible intervals and variance.
  - Complex ($K_n \ge 2$) events with shared ambiguous edges run a fast, deterministic collapsed Gibbs sampler (`~59 µs` per event, using a group-balanced prior so inclusion vs. exclusion stays neutral regardless of path count) to quantify true posterior variance (**`Inf_Var`**, column 15 of `.psi.gz`, preserving backward compatibility with Whippet 2's columns 1–14) and realistic 90% credible intervals.

<details>
<summary><b>Mathematical Details: Effective Length, CC-VBEM, Unified GC Pass & Gibbs Sampler</b></summary>

- **Post-pass effective length:** From the 5-point Gaussian-smoothed empirical FLD PMF $p(d)$, CDF $F(L)$, and first moment $M_1(L) = \sum_{d=1}^L d\,p(d)$:

  $$
  \ell_t = (L_t + 1)\,F(\min(L_t, D_{\max})) - M_1(\min(L_t, D_{\max})).
  $$

- **CC-VBEM + SQUAREM** (default sparse Dirichlet prior $\alpha_0 = 0.01$; `--use-em` selects MLE; $\psi(\cdot)$ denotes the digamma function, distinct from splicing inclusion $\Psi$):

  $$
  w_t^{(k)} = \frac{\exp\!\left(\psi(\alpha_0 + c_t^{(k)})\right)}{\ell_t}, \qquad c_t^{(k+1)} = c_t^{\text{uniq}} + \sum_{e \ni t} n_e \frac{w_t^{(k)}}{\sum_{u \in e} w_u^{(k)}}.
  $$

- **Pass 2 GC effective-length update** (run on primer-only counts to avoid double-counting GC, where $E_{\text{GC}}[t, b]$ is transcript $t$'s expected GC fraction in bin $b$ and $w_{\text{GC}}[b] = \text{expected}[b] / \text{observed}[b] \in [0.15, 6.0]$ so $1 / w_{\text{GC}}[b]$ scales $\ell_t$ by each bin's relative observed capture rate):

  $$
  \ell_t^{\text{GC}} = \ell_t \cdot \sum_{b=0}^{19} \frac{E_{\text{GC}}[t, b]}{w_{\text{GC}}[b]}.
  $$

- **Targeted collapsed Gibbs (`K_n >= 2`):** Group-balanced Dirichlet prior ($\alpha_p = 1/n_{\text{inc}}$ on inclusion paths, $1/n_{\text{exc}}$ on exclusion paths, yielding marginal $\text{Beta}(1,1)$ on $\Psi$), deterministic `Xoshiro256++` seeded from `(seed, gene_id, node_1based)`, `100` draws (`32` burn-in, thinning `2`), alternating $O(|G_g|)$ multinomial edge allocation with exact rate-parameterized Gamma draws:

  $$
  \theta_p^{(b)} \sim \text{Gamma}\!\left(\alpha_p + z_p^{(b)},\; \text{rate} = \ell_p\right), \qquad \Psi^{(b)} = \frac{\sum_{p \in \text{Inc}} \theta_p^{(b)}}{\sum_q \theta_q^{(b)}}, \qquad \text{InfVar} = \widehat{\text{Var}}(\Psi^{(1..B)}).
  $$

</details>

---

## 6. Covariate-Aware Differential Splicing (`delta`)

`delta` replaces Whippet's pooled two-group Beta Monte Carlo test with a multithreaded **Beta-Binomial Generalized Linear Model (GLM)** featuring Empirical Bayes overdispersion shrinkage—analogous in design philosophy to `edgeR` quasi-likelihood GLMs, adapted specifically for fractional splicing proportions ($\Psi_{ij}, N_{ij}$).

### Key Modeling Decisions & Why They Matter
1. **Arbitrary multivariable designs & graceful missingness:** Supports both classic two-group runs (`-a` / `-b`) and DESeq2-style metadata formulas (`--design`, `--formula "~ batch + donor + condition"`, `--test condition:treated_vs_ctrl`, continuous covariates, and `*`/`:` interactions). When an event lacks coverage in one batch level, per-event pivoted-QR automatically drops only the collinear nuisance column rather than failing or discarding the event.
2. **Propagating complex-event uncertainty (`Inf_Var`):** Samples are filtered on raw `Total_Reads >= --min-reads`, then complex $K_n \ge 2$ observations whose Gibbs `Inf_Var` exceeds nominal Beta variance have their effective depth smoothly down-weighted (down to `0.1×`, leaving $K_1$ events untouched). On `E-MTAB-3766`, propagating `Inf_Var` **simultaneously increased sensitivity and lowered realized FDR** by reducing dispersion noise across complex events.
3. **Excluding constitutive events from dispersion estimation (critical real-genome finding):** On a full human transcriptome, ~63% of testable CSG nodes are constitutive (`PSI = 0.0` or `1.0` in every sample). Because boundary events carry no information about sample-to-sample variance, including them causes the Cox–Reid adjustment term alone to produce a spurious peak near $\phi \approx 0.64$, collapsing the global dispersion prior and wiping out power on 3v3 cohorts. Excluding constitutive/boundary events from dispersion fitting fixes the prior collapse and speeds up `delta` nearly 3× (to **2.5 s**).
4. **Separate Empirical Bayes priors for Spliced vs. Tandem UTR (`TS`/`TE`) events:** Benchmarking on human RNA-seq revealed two distinct dispersion regimes: exon-body tandem UTR ends (`TS`/`TE`) exhibit **7.8× higher median overdispersion** (`median phi = 0.117`) than canonical splice-junction events (`CE/AA/AD/RI/AF/AL`, `median phi = 0.015`). Pooling both classes into one prior over-penalizes true spliced changes while under-penalizing noisy UTR ends. Stratifying the EB prior by event class (**Spliced** vs. **Tandem UTR**, opt out via `--no-stratify-eb`) boosts sensitivity **6×–25×** while slashing `TS`/`TE` false positives by **>300×**.
5. **Finite-sample $F$-test LRT & consistent effect sizes:** Significance uses a finite-sample $F(\text{df},\, d_{\text{res}} + d_0)$ likelihood-ratio test (capturing uncertainty in the shrunk dispersion estimate, similar to `edgeR::glmQLFTest`) with Benjamini–Hochberg `FDR`. Marginal covariate-adjusted `Psi_A`, `Psi_B`, and `DeltaPsi` are computed over the observed covariate distribution so `DeltaPsi` and `LogOddsRatio` always share the same sign.

<details>
<summary><b>Mathematical Details: Beta-Binomial GLM, Cox–Reid APL & Smyth Bias Correction</b></summary>

- **Effective depth & continuous Beta-Binomial log-likelihood:** With empirical posterior variance $v_j = \text{InfVar}_j$ (`Inf_Var` in `.psi.gz`) and nominal Beta variance $v_{\text{nom}}(\Psi_j, n_j) = (\Psi_j n_j + 1)((1-\Psi_j)n_j + 1)/((n_j+2)^2(n_j+3))$,

  $$
  n_{\text{eff}, j} = \max\!\left(1.0,\; n_j \cdot \text{clamp}\!\left(\frac{v_{\text{nom}}(\Psi_j, n_j)}{v_j},\; 0.1,\; 1.0\right)\right), \qquad k_j = \Psi_j n_{\text{eff}, j}, \qquad r_j = n_{\text{eff}, j} - k_j,
  $$

  precision $M = (1-\phi)/\phi = e^\theta$, and $\mu_j = \text{logistic}(x_j^\top \beta)$:

  $$
  \ell(\beta, M) = \sum_{j=1}^{m} \left[\ln\Gamma(k_j + M\mu_j) + \ln\Gamma(r_j + M(1-\mu_j)) - \ln\Gamma(n_{\text{eff},j} + M) + \ln\Gamma(M) - \ln\Gamma(M\mu_j) - \ln\Gamma(M(1-\mu_j))\right].
  $$

  Solved by damped Fisher scoring (IRLS) with `< 1e-14` `trigamma`, adaptive Cholesky damping, backtracking line search, and light ridge regularization ($\lambda = 10^{-4}$, intercept $10^{-6}$) for guaranteed finite convergence under complete separation (`0.0` vs. `1.0`).

- **Constitutive exclusion criterion:** Events with $\sum k_j < 1$, $\sum r_j < 1$, pooled $\bar{\Psi} \le 0.01$ or $\ge 0.99$, or $\max \Psi_j - \min \Psi_j < 10^{-4}$ skip Pass 2 raw estimation and Pass 3 MAP optimization, receiving their stratum prior $\mu_\theta$ directly.

- **Cox–Reid APL & Smyth finite-df median bias correction:** Within each stratum (**Spliced** vs. **Tandem UTR**, falling back to combined if $< 20$ interior events), Pass 2 maximizes

  $$
  \text{APL}(\theta) = \ell(\hat{\beta}(\theta), e^\theta) - \tfrac{1}{2}\ln\det\!\left(X^\top W X + 10^{-8} I\right)
  $$

  on $\theta \in [-4.0, 11.5]$. The robust median/MAD normal prior $\mathcal{N}(\mu_\theta, \sigma_\theta^2)$ corrects finite-df bias (downward dispersion bias / upward bias in $\hat{\theta} = \ln\hat{M}$) via the negative shift

  $$
  \mu_\theta = \text{median}(\hat{\theta}) + \ln\!\left(\frac{\text{median}(\chi^2_{d_{\text{res}}})}{d_{\text{res}}}\right), \qquad d_0 = \text{clamp}\!\left(\frac{2}{\sigma_\theta^2},\; 1,\; 30\right).
  $$
- **Pass 3 inference:** Maximizes $\text{APL}(\theta) - (\theta - \mu_\theta)^2 / (2\sigma_\theta^2)$ for shrunk $\tilde{M}_i = e^{\tilde{\theta}_i}$, computes $\text{LRT} = 2(\ell_{\text{full}} - \ell_{\text{null}})$, and evaluates $p_i = P(F_{\text{df},\; d_{\text{res}} + d_0} \ge \text{LRT}/\text{df})$ (reducing to $\chi^2_{\text{df}}$ under `--fixed-dispersion`).

</details>

---

## 7. Benchmarks & Validation

### 7.1 Full Human-Scale Benchmark: `E-MTAB-3766` (Ensembl GRCh37.71, 480M Reads)

We benchmarked Langur against our optimized native Apple Silicon Julia 1.9 fork (**`Whippet2-turbo`**), **Original `Whippet.jl`** (Julia 1.6, single-threaded under Rosetta 2), and published baselines from Soneson et al. (*Genome Biology* 2016) on **`E-MTAB-3766`**: 3 vs. 3 paired-end 101 bp replicates (**20,410 protein-coding genes**, **145,342 transcripts**, **1,000 truly differentially spliced genes**, 40M pairs = 80M reads/sample, **480M total reads**, 16 threads, `--biascorrect`, `/usr/bin/time -l`). *(In Tables A–C, `Whippet2-turbo` commits reflect progressive upgrades ported from this benchmark: `1e860f9` fixed the upstream `TS` PSI reversal bug; `c725028` ported constitutive-event exclusion and Spliced/UTR stratified EB priors into `whippet-delta.jl`; and `400f0f3` reduced `--batch` peak RAM with identical quantification accuracy.)*

#### A. Runtime, Throughput & Memory
| Stage (`E-MTAB-3766`, `--biascorrect`) | Langur (Rust, 16T, arm64) | `Whippet2-turbo` (`400f0f3`, Julia 1.9, 16T, arm64) | Original `Whippet.jl` v1.6.2 (Julia 1.6, 1T, Rosetta 2) | Langur Advantage |
|---|---:|---:|---:|---:|
| **Index Build** (20,410 genes, 145k tx) | **22.6 s** · **4.71 GB** RSS · 820 MB *(incl. 1.1 Gb decoys)* | 233.0 s · 8.74 GB RSS · 529 MB *(no decoys)* | 340.8 s · 8.58 GB RSS · 473 MB *(no decoys)* | **15.1× vs. W1 · 10.3× vs. W2** |
| **Index Load** | **0.04–0.30 ms** (`mmap`) | 6.20 s | 7.33 s | **>20,000× faster** |
| **Sample 1 `quant`** (80M reads) | **60.4 s** (**1.33M reads/s**, 13.3× CPU) · **1.69 GB** footprint | 81.1 s (987k reads/s) · 9.62 GB footprint | 1,608.4 s *(26m 48s, 49.7k reads/s)* · 7.19 GB footprint | **26.6× vs. W1 · 1.34× vs. W2**; **4.3×–5.7× lower RAM** |
| **6-Sample Cohort (`--batch`, 480M reads)** | **373.1 s** *(6m 13s, **1.29M reads/s**)* · **1.67 GB** footprint | 433.8 s (1.11M reads/s) · 9.83 GB footprint | *9,650.6 s* *(2h 41m, 6× extrap.)* · 7.19 GB footprint | **25.9× vs. W1 · 1.16× vs. W2**; **5.9× lower RAM** |
| **Mapped Pairs (Sample 1 of 40M)** | **37.91M (94.8%)** *(truth: 95.0%)* | 36.80M (92.0%) | 36.80M (92.0%) | **+1.11M more mapped pairs/sample** |
| **`delta` Testing** (222k–223k events, `-r 5 -s 2`) | **2.52 s** · **189 MB** footprint | 27.38 s (GLM) / 85.51 s (`--legacy`) · 1.87 GB | 85.51 s *(`--legacy` MC)* · 609 MB | **33.9× vs. legacy · 10.9× & 9.9× RAM vs. W2** |

*Note: Original `Whippet.jl` v1.6.2 requires Julia 1.6 (single-threaded, x86_64 under Rosetta 2 on Apple Silicon). In Langur, all post-alignment stages combined on Sample 1 (2-pass GC CC-VBEM: `0.28 s`; 481k events + 105k Gibbs samplers: `0.55 s`; 5 parallel gzip writers: `0.07 s`) take just **0.90 s** total (**13.7× faster** than `Whippet2-turbo`'s `12.3 s`). On the 6-sample `--batch` run, Langur finishes 61 s faster than `Whippet2-turbo` (`373.1 s` vs. `433.8 s`) while simultaneously performing selective genomic decoy screening, SIMD gapped indel alignment, CC-VBEM, and 105k targeted Gibbs runs per sample at 5.9× lower peak RAM.*

#### B. Expression & PSI Quantification Accuracy (Sample 1 vs. Exact RSEM Ground Truth)
| Metric (Sample 1, 145,342 Transcripts / 20,410 Genes) | Original `Whippet.jl` v1.6.2 | `Whippet2-turbo` (`c725028`) | **Langur (CC-VBEM)** | Langur Improvement |
|---|---:|---:|---:|---|
| **Isoform Spearman $\rho$** (all 145k tx) | 0.7954 | 0.8029 | **0.8772** | +0.082 vs. W1 · +0.074 vs. W2 |
| **Isoform MARD** ($\text{TPM} \ge 1$) | 0.1145 | 0.0687 | **0.0406** | **2.82× vs. W1 · 1.69× vs. W2** |
| **Gene Spearman $\rho$ / $\log_2$ Pearson $r$** | 0.9797 / 0.9861 | 0.9810 / 0.9880 | **0.9953 / 0.9970** | Near-exact gene abundance |
| **Gene MARD** ($\text{TPM} \ge 1$) | 0.0753 | 0.0443 | **0.0137** | **5.50× vs. W1 · 3.23× vs. W2** |
| **vs. `Salmon`** — Isoform / Gene MARD ($\text{TPM} \ge 1$) | — | 0.0602 / 0.0383 | **0.0226 / 0.0091** | Near-parity with dedicated `Salmon` |
| **`TS` PSI concordance vs. Langur** (`Reads >= 30`) | -0.1085 *(upstream reversal bug)* | **0.9663** *(fixed in `1e860f9`)* | *(reference)* | Confirmed upstream `TS` bug |

- **Stage 2 PSI Concordance vs. Fixed `Whippet2-turbo` (`1e860f9`):** Across 172,725 high-depth events ($\ge 30$ reads), overall concordance is **$r = 0.9778$** ($\rho = 0.9729$, $\text{MAE} = 0.0258$), with every canonical spliced class at **$r = 0.9863\text{–}0.9956$** (`CE = 0.9956`, `AA = 0.9945`, `AD = 0.9942`, `RI = 0.9941`, `AL = 0.9940`, `AF = 0.9863`, `TS = 0.9663`, `TE = 0.8950`). All top-10 discrepancies were 1–22 bp terminal `TE`/`TS` sub-nodes where Langur matched the dominant RSEM ground-truth transcript end while `Whippet2-turbo` dropped reads on the short sub-node.

#### C. Gene-Level Differential Splicing Accuracy (20,410 Genes, 1,000 Truly DS)

| Method (`E-MTAB-3766`, Gene-Level Simes+BH) | Event Set | TPR @0.01 | FDR @0.01 | TPR @0.05 | FDR @0.05 | TPR @0.10 | FDR @0.10 | AUROC | AUPRC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Langur Stratified + `Inf_Var` (default, `-r 5 -s 2`)** | All (`TS/TE` + Spliced) | **0.270** | **0.029** | **0.336** | **0.054** | **0.360** | **0.084** | **0.8294** | **0.5642** |
| **Langur Stratified + `Inf_Var` (default, `-r 5 -s 2`)** | Spliced (`CE/AA/AD/RI/AF/AL`) | **0.294** | **0.003** | **0.336** | **0.020** | **0.357** | **0.030** | **0.8119** | **0.5672** |
| **`Whippet2-turbo` (`c725028`) Stratified GLM (`-r 5 -s 2`)** | All / Spliced | 0.225 / 0.252 | 0.017 / 0.008 | 0.306 / 0.309 | 0.053 / 0.013 | 0.329 / 0.332 | 0.073 / 0.018 | 0.8209 / 0.8007 | 0.5512 / 0.5486 |
| `Whippet2` **`--legacy` Beta MC** (`Prob ≥ 0.99/0.95/0.90`) | All | 0.585 | **0.681** | 0.744 | **0.815** | 0.867 | **0.865** | 0.8903 | **0.3517** |
| `Whippet2` **`--legacy` Beta MC** (`\|ΔΨ\| ≥ 0.1`) | All | 0.532 | **0.534** | 0.638 | **0.669** | 0.692 | **0.740** | 0.8649 | **0.3926** |
| Published: `rMATS_junction` (Soneson et al. 2016) | Junction | 0.336 | **0.463** | 0.379 | **0.571** | 0.409 | **0.623** | 0.7666 | **0.3937** |
| Published: `cuffdiff_ensembl` (Soneson et al. 2016) | Isoform | 0.114 | 0.000 | 0.142 | 0.000 | 0.159 | 0.000 | 0.6110 | **0.2599** |
| Published: `dexseq_splicinggraph` (Soneson et al. 2016) | Exon-bin | 0.637 | **0.172** | 0.675 | **0.252** | 0.697 | **0.294** | 0.8913 | 0.7272 |

**Key findings:**
1. **15× FDR reduction over legacy Whippet & clear margin over `Whippet2-turbo`:** Original `--legacy` pooled Beta Monte Carlo exhibits **68%–87% realized FDR** (`53%–74%` even with `|ΔΨ| ≥ 0.1`) because it ignores biological overdispersion and multiple testing. Langur's Stratified EB Beta-Binomial GLM brings realized FDR directly onto target at nominal 0.05 (**5.4%** across all events, **2.0%** across spliced events, and **0.3%** on spliced events at nominal 0.01) while lifting AUPRC by **+0.21** (`0.3517 → 0.5642`). Langur outperforms `Whippet2-turbo` (`c725028`) across both **All events** (**+30 TP genes** at `0.05`: `0.336` vs. `0.306`, **+0.0130 AUPRC**, **+0.0085 AUROC** at matching `5.4%` vs. `5.3%` FDR) and **Spliced events** (**+27 TP genes**: `0.336` vs. `0.309`, **+0.0186 AUPRC**, **+0.0112 AUROC**, with `CE` alone reaching `163 TP / 1 FP` = **0.006 FDR** and `0.3467` AUPRC vs. `141 / 1` and `0.3368` in `Whippet2-turbo`), driven by Langur's tighter Spliced EB prior (`median phi = 0.0150`, `sd(log M) = 0.399`, **`d0 = 12.54`** vs. `9.77`).
2. **Strong margin over published splicing callers:** Outperforms `rMATS_junction` by **+0.17 AUPRC** and **+0.063 AUROC** with **10.6× lower realized FDR on all events** (**28.6× lower on spliced events**), and outperforms `cuffdiff_ensembl` by **+0.30 AUPRC**. *(Whole-gene exon-bin tools like `DEXSeq` attain higher raw recall on `E-MTAB-3766` because Soneson et al. simulated isoform swaps that frequently differ across UTR body coverage rather than a single high-depth splice junction—yet even `DEXSeq` exceeds nominal 5% FDR by 5×–10×.)*
3. **Ablation proof for both new statistical components:**
   - **Stratifying Spliced vs. Tandem UTR priors (`-r 5 -s 1`):** Increased `TPR @0.05` **6×** (`0.042 → 0.249` all) and **25×** (`0.011 → 0.276` spliced) while cutting `TE`/`TS` false positives **>300×**.
   - **Enabling Gibbs `Inf_Var` (`-r 5 -s 2`):** Beat `--ignore-inf-var` by **+33 TP genes** (`TPR @0.05`: `0.303 → 0.336` all; `0.305 → 0.336` spliced), **+0.0116 All AUPRC** (`0.5526 → 0.5642`), and **+0.0046 Spliced AUPRC** (`0.5626 → 0.5672`) at matching All realized FDR (`0.053 → 0.054`), tightening the spliced EB prior (`sd(log M): 0.503 → 0.399`, `d0: 7.92 → 12.54`).

### 7.2 Automated Verification Suite & Synthetic Cohort
- **42/42 workspace unit and integration tests pass** (`cargo test --workspace`, zero warnings), verifying exact Whippet 2 CSG/motif parity (`ENSMUSG00000038685` 18-path, `nrf2test` 6-path), bounded outward bridge-edge sweeps on multi-exon skips (`LTBP2` topology), FR/RF/unstranded pairing, strict-inequality pseudogene decoy ties vs. diverged decoys, 1-SNP splice jumps, SIMD indel recovery, CC-VBEM paralog routing (`>99.5%` to expressed parent), 15 bp microexon PSI unbiasedness, PCR-stack soft-capping, Gibbs ~90% CI coverage across 200 $K_2$ replicates, 80%-constitutive + stratified EB dispersion recovery, batch confounding (`~ batch + condition`), interactions, continuous covariates, pivoted-QR collinearity reduction, and 100% bit-for-bit determinism across 1 vs. 4 threads.
- **200-gene 3v3 synthetic cohort (600k reads, 4 threads):** `index` ran in **0.032 s** (185× faster than `Whippet2-turbo`), 6-sample `--batch` `quant` in **0.377 s** (126× faster, 17× lower RAM), and `delta` in **0.039 s** (122× faster, 95× lower RAM), achieving **100% sensitivity (100/100)** at **0.00% realized FDR** ($\Delta\Psi$ $r = 0.9935$ vs. truth).

### 7.3 Planned Real-Data Validations (Optional Next Steps)
- **LRGASP short-read vs. long-read ground truth:** Stratified PSI Pearson $r$ and MAE against matched PacBio HiFi / ONT long reads across `CE`, `AA`/`AD`, `RI`, and complex $K_n \ge 2$ events (mirroring MAJIQ V3, 2026).
- **Real multi-sample cohort null splits:** Empirical FPR at $\alpha = 0.05$ on random same-condition splits of real human replicates (e.g., GTEx) and multivariable batch/donor adjustment (mirroring MAJIQ V2, 2023).

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
