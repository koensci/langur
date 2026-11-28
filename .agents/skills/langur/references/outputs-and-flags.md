# Langur Output Schemas & Advanced Flags

## Output File Schemas

### `.psi.gz` (15 tab-separated columns; Whippet2 cols 1–14 + `Inf_Var`)

```text
Gene  Node  Coord  Strand  Type  Psi  CI_Width  CI_Lo,Hi  Total_Reads  Complexity  Entropy  Inc_Paths  Exc_Paths  Edges  Inf_Var
```

- **`Type`:** `CE` (core/cassette exon), `AA` (alt. 3' acceptor), `AD` (alt. 5' donor), `RI` (retained intron), `AF` (alt. first exon), `AL` (alt. last exon), `TS` (tandem TSS), `TE` (tandem APA/3' end), `BS` (back-splice/circRNA, requires `--circ`).
- **`Complexity`:** `K0`/`K1` (binary, ≤ 2 paths) vs. `K2+` (complex multi-path hub; $K_n = \lceil\log_2(\text{paths})\rceil$).
- **`Inf_Var`:** Posterior variance of $\Psi$ (exact Beta variance for binary/unambiguous nodes; targeted collapsed Gibbs variance for ambiguous $K_n \ge 2$ nodes). Automatically propagated by `langur delta`.

### `.diff.gz` (from `langur delta`)

**Two-group (`-a`/`-b`) or explicit 2-level contrast (`--test cond:A_vs_B`, 20 columns):**
```text
Gene  Node  Coord  Strand  Type  Complexity  Entropy  SampleCnt_A  SampleCnt_B  MeanReads_A  MeanReads_B  Psi_A  Psi_B  DeltaPsi  LogOddsRatio  SE  Dispersion  LRT  Pvalue  FDR
```
- `DeltaPsi = Psi_A - Psi_B` (covariate-adjusted marginal $\Psi$ in `--design` mode, guaranteed to share the exact sign of `LogOddsRatio`).
- `Dispersion` is the shrunk Beta-Binomial overdispersion $\tilde{\phi}$.

**General GLM mode (continuous covariate or multi-df omnibus test, 18 columns):**
```text
Gene  Node  Coord  Strand  Type  Complexity  Entropy  SampleCnt  MeanReads  MeanPsi  DeltaPsi  LogOddsRatio  SE  Dispersion  LRT  DF  Pvalue  FDR
```

### `.jnc.gz` (8 tab-separated columns)
```text
CHROM  DONOR_COORD  ACCEPTOR_COORD  GENEID:NODES:TYPE  READCOUNT  STRAND  Distinct_Pos  Max_Stack
```
- `Distinct_Pos` and `Max_Stack` report overhang-position diversity for PCR-stack QC.

### `.gene.tpm.gz` & `.isoform.tpm.gz` (3 tab-separated columns)
```text
Gene     TpM  Read_Counts
Isoform  TpM  Read_Counts
```

## Key Optional CLI Flags

### `langur index`
- `-k <17>`, `-w <8>`: Minimizer k-mer and window size (`k <= 27` is 54-bit collision-free).
- `--basic-or-ss-only`: Restrict transcripts to GENCODE `basic` while using other transcripts for splice sites.
- `--suppress-low-tsl`: Drop low-quality TSL 3+ transcripts.
- `--bam <sorted.bam>` + `--bam-min-reads <N>` + `--bam-both-novel`: Supplement index with BAM-supported novel splice sites.
- `--no-decoys` / `--decoys <fa>` / `--decoy-pad <120>`: Customize selective genomic decoy screening.
- `--skip-aux-tsv`: Skip writing `.exons.tsv.gz` and `.paths.tsv.gz`.

### `langur quant`
- `--biascorrect`: Unified 5'-hexamer + 20-bin GC bias correction (recommended).
- `--stranded` (FR) / `--stranded-reverse` (RF) / `--pair-same-strand`: Library strand orientation.
- `--circ`: Enable circular back-splicing (`BS`) detection.
- `--sam`: Stream SAM alignments to stdout (cannot be combined with `--batch`).
- `--vb-prior <0.01>` / `--use-em`: Stage 1 transcript VBEM prior or MLE EM.
- `--gibbs-samples <100>` / `--seed <42>`: Complex-node Gibbs posterior draws (`0` disables).
- `--no-stack-filter` / `--stack-cap <2>`: Control PCR-stack soft-capping on novel/`RI` junctions.
- `--max-frag-len <2000>`, `-S / --score-min <0.65>`: Paired fragment-length and alignment-score thresholds.

### `langur delta`
- `-r / --min-reads <5>`, `-s / --min-samples <1>`: Minimum raw `Total_Reads` and valid sample counts (recommend `-s 2` or higher for replicated studies).
- `--directory <dir>`: Base directory for relative paths or `-a`/`-b` glob patterns.
- `--include-circ`: Include `BS` events.
- `--no-shrink`, `--no-stratify-eb`, `--fixed-dispersion <phi>`, `--ignore-inf-var`: Statistical ablations (keep defaults for production runs).
