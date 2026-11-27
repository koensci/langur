#!/usr/bin/env python3
"""Synthetic cohort generator and concordance evaluator for Langur vs Whippet2-turbo benchmarks.

Subcommands:
  generate --out-dir DIR [--genes 200] [--pairs-per-sample 50000] [--samples-per-group 3] [--seed 42]
  compare  --truth TRUTH.tsv --w3-dir W3_DIR --w2-dir W2_DIR --w3-diff W3.diff.gz --w2-diff W2.diff.gz
"""

import argparse
import gzip
import math
import os
import random
from pathlib import Path

RC_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choices("ACGT", k=length))


def make_intron(rng: random.Random, length: int, strand: str) -> str:
    body = random_dna(rng, length - 4)
    if strand == "+":
        return "GT" + body + "AG"
    else:
        # Reverse-complement of GT...AG is CT...AC on the '+' genome strand
        return "CT" + body + "AC"


def splice_isoform(chrom_seq: str, exons_1based: list[tuple[int, int]], strand: str) -> str:
    parts = [chrom_seq[s - 1 : e] for s, e in exons_1based]
    seq = "".join(parts)
    if strand == "-":
        seq = revcomp(seq)
    return seq


def mutate_seq(rng: random.Random, seq: str, err_rate: float) -> str:
    if err_rate <= 0.0:
        return seq
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < err_rate:
            b = out[i]
            out[i] = rng.choice([c for c in "ACGT" if c != b])
    return "".join(out)


def generate_benchmark(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    n_genes = args.genes
    n_rep = args.samples_per_group
    pairs_per_sample = args.pairs_per_sample
    read_len = 100

    # Build 2 chromosomes with realistic GT-AG introns and diverse event architectures:
    #   - 60% binary CE cassette exons (K_1): half null (dPSI=0), half DE (dPSI=+/-0.40)
    #   - 20% alternative 3'/5' splice sites (AA/AD): half null, half DE (dPSI=+/-0.35)
    #   - 20% complex 4-isoform tandem cassette exons (K_2): half null, half DE (dPSI=+/-0.35)
    chrom_chunks: dict[str, list[str]] = {"chr1": [], "chr2": []}
    chrom_pos: dict[str, int] = {"chr1": 0, "chr2": 0}

    gtf_lines: list[str] = []
    # Gene catalog: list of dicts with isoform sequences, true group PSIs, event coord/type
    genes = []

    for g in range(n_genes):
        chrom = "chr1" if g < n_genes // 2 else "chr2"
        strand = "+" if g % 2 == 0 else "-"
        gid = f"GENE_{g:04d}"

        # Intergenic spacer (800 bp)
        spacer = random_dna(rng, 800)
        chrom_chunks[chrom].append(spacer)
        chrom_pos[chrom] += len(spacer)
        cur = chrom_pos[chrom]

        kind_mod = g % 10
        is_de = (g % 2 == 1)

        if kind_mod < 6:
            # Binary CE cassette exon: E1(220) - I1(250) - E2(140) - I2(250) - E3(220)
            e1_s, e1_e = cur + 1, cur + 220
            i1 = make_intron(rng, 250, strand)
            e2_s, e2_e = e1_e + 250 + 1, e1_e + 250 + 140
            i2 = make_intron(rng, 250, strand)
            e3_s, e3_e = e2_e + 250 + 1, e2_e + 250 + 220
            seg = random_dna(rng, 220) + i1 + random_dna(rng, 140) + i2 + random_dna(rng, 220)
            chrom_chunks[chrom].append(seg)
            chrom_pos[chrom] += len(seg)

            tx_exons = [
                (f"{gid}_T1", [(e1_s, e1_e), (e2_s, e2_e), (e3_s, e3_e)]),
                (f"{gid}_T2", [(e1_s, e1_e), (e3_s, e3_e)]),
            ]
            if is_de:
                psi_a, psi_b = (0.75, 0.25) if (g // 2) % 2 == 0 else (0.25, 0.75)
            else:
                psi_a, psi_b = 0.50, 0.50
            event_coord = f"{chrom}:{e2_s}-{e2_e}"
            event_type = "CE"
            event_class = "K1_CE"
            iso_weights_fn = lambda psi: [psi, 1.0 - psi]

        elif kind_mod < 8:
            # Alternative 3' or 5' splice site (AA/AD):
            # E1(220) - I1(250) - [E2_ext(60) + E2_core(180)] - I2(250) - E3(220)
            # On '+' strand: shorter intron uses E2_ext+E2_core (AA node is E2_ext)
            e1_s, e1_e = cur + 1, cur + 220
            i1 = make_intron(rng, 250, strand)
            e2_ext_s, e2_ext_e = e1_e + 250 + 1, e1_e + 250 + 60
            e2_core_s, e2_core_e = e2_ext_e + 1, e2_ext_e + 180
            i2 = make_intron(rng, 250, strand)
            e3_s, e3_e = e2_core_e + 250 + 1, e2_core_e + 250 + 220
            # Ensure internal acceptor AG (for +) or donor CT (for -) at boundary of e2_ext
            e2_ext_seq = random_dna(rng, 58) + ("AG" if strand == "+" else "AC")
            seg = random_dna(rng, 220) + i1 + e2_ext_seq + random_dna(rng, 180) + i2 + random_dna(rng, 220)
            chrom_chunks[chrom].append(seg)
            chrom_pos[chrom] += len(seg)

            tx_exons = [
                (f"{gid}_T1", [(e1_s, e1_e), (e2_ext_s, e2_core_e), (e3_s, e3_e)]),
                (f"{gid}_T2", [(e1_s, e1_e), (e2_core_s, e2_core_e), (e3_s, e3_e)]),
            ]
            if is_de:
                psi_a, psi_b = (0.72, 0.28) if (g // 2) % 2 == 0 else (0.28, 0.72)
            else:
                psi_a, psi_b = 0.50, 0.50
            event_coord = f"{chrom}:{e2_ext_s}-{e2_ext_e}"
            event_type = "AA" if strand == "+" else "AD"
            event_class = "K1_ALTSS"
            iso_weights_fn = lambda psi: [psi, 1.0 - psi]

        else:
            # Complex 4-isoform tandem cassette (K_2):
            # E1(220) - I1(220) - E2(130) - I2(220) - E3(130) - I3(220) - E4(220)
            e1_s, e1_e = cur + 1, cur + 220
            i1 = make_intron(rng, 220, strand)
            e2_s, e2_e = e1_e + 220 + 1, e1_e + 220 + 130
            i2 = make_intron(rng, 220, strand)
            e3_s, e3_e = e2_e + 220 + 1, e2_e + 220 + 130
            i3 = make_intron(rng, 220, strand)
            e4_s, e4_e = e3_e + 220 + 1, e3_e + 220 + 220
            seg = (
                random_dna(rng, 220)
                + i1
                + random_dna(rng, 130)
                + i2
                + random_dna(rng, 130)
                + i3
                + random_dna(rng, 220)
            )
            chrom_chunks[chrom].append(seg)
            chrom_pos[chrom] += len(seg)

            tx_exons = [
                (f"{gid}_T1", [(e1_s, e1_e), (e2_s, e2_e), (e3_s, e3_e), (e4_s, e4_e)]),
                (f"{gid}_T2", [(e1_s, e1_e), (e2_s, e2_e), (e4_s, e4_e)]),
                (f"{gid}_T3", [(e1_s, e1_e), (e3_s, e3_e), (e4_s, e4_e)]),
                (f"{gid}_T4", [(e1_s, e1_e), (e4_s, e4_e)]),
            ]
            if is_de:
                psi_a, psi_b = (0.72, 0.28) if (g // 2) % 2 == 0 else (0.28, 0.72)
            else:
                psi_a, psi_b = 0.50, 0.50
            # Track E2 inclusion: T1 + T2 have E2; T3 + T4 exclude E2
            event_coord = f"{chrom}:{e2_s}-{e2_e}"
            event_type = "CE"
            event_class = "K2_COMPLEX"
            iso_weights_fn = lambda psi: [0.5 * psi, 0.5 * psi, 0.5 * (1.0 - psi), 0.5 * (1.0 - psi)]

        for tid, exs in tx_exons:
            for es, ee in exs:
                gtf_lines.append(
                    f'{chrom}\tSIM\texon\t{es}\t{ee}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n'
                )

        genes.append(
            {
                "gid": gid,
                "chrom": chrom,
                "strand": strand,
                "tx_exons": tx_exons,
                "psi_a": psi_a,
                "psi_b": psi_b,
                "is_de": is_de,
                "event_coord": event_coord,
                "event_type": event_type,
                "event_class": event_class,
                "iso_weights_fn": iso_weights_fn,
            }
        )

    # Trailing 1000 bp per chromosome
    chrom_seqs: dict[str, str] = {}
    for ch in ("chr1", "chr2"):
        chrom_chunks[ch].append(random_dna(rng, 1000))
        chrom_seqs[ch] = "".join(chrom_chunks[ch])

    # Pre-splice isoform sequences
    for gene in genes:
        ch_seq = chrom_seqs[gene["chrom"]]
        gene["iso_seqs"] = [
            splice_isoform(ch_seq, exs, gene["strand"]) for _, exs in gene["tx_exons"]
        ]

    # Write genome.fa and genes.gtf
    fa_path = out_dir / "genome.fa"
    with open(fa_path, "w") as f:
        for ch in ("chr1", "chr2"):
            f.write(f">{ch}\n")
            s = chrom_seqs[ch]
            for i in range(0, len(s), 80):
                f.write(s[i : i + 80] + "\n")

    gtf_path = out_dir / "genes.gtf"
    with open(gtf_path, "w") as f:
        f.writelines(gtf_lines)

    # Generate 2*n_rep samples (A1..An, B1..Bn) with mild biological Beta-Binomial jitter
    samples = [f"A{i+1}" for i in range(n_rep)] + [f"B{i+1}" for i in range(n_rep)]
    truth_rows = []
    qual = "I" * read_len

    # Assign sample-specific true PSIs
    sample_psis: dict[str, list[float]] = {s: [] for s in samples}
    for gene in genes:
        for s in samples:
            base = gene["psi_a"] if s.startswith("A") else gene["psi_b"]
            # Realistic biological replicate dispersion (SD ~ 0.03)
            p_rep = min(0.95, max(0.05, rng.gauss(base, 0.03)))
            sample_psis[s].append(p_rep)
        truth_rows.append(
            (
                gene["gid"],
                gene["event_coord"],
                gene["event_type"],
                gene["event_class"],
                1 if gene["is_de"] else 0,
                gene["psi_a"],
                gene["psi_b"],
                gene["psi_a"] - gene["psi_b"],
                sample_psis["A1"][-1],
            )
        )

    with open(out_dir / "truth.tsv", "w") as f:
        f.write("Gene\tCoord\tType\tClass\tIs_DE\tTrue_Psi_A\tTrue_Psi_B\tTrue_DeltaPsi\tTrue_Psi_A1\n")
        for row in truth_rows:
            f.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t"
                f"{row[5]:.5f}\t{row[6]:.5f}\t{row[7]:.5f}\t{row[8]:.5f}\n"
            )

    # Write paired-end FASTQ.gz files (FR paired-end, mean fragment length 240 +/- 25 bp, 0.2% error)
    manifest_rows = ["sample,r1,r2,out_prefix"]
    pairs_per_gene = max(1, pairs_per_sample // n_genes)

    for s in samples:
        r1_path = out_dir / f"{s}_R1.fq.gz"
        r2_path = out_dir / f"{s}_R2.fq.gz"
        manifest_rows.append(f"{s},{r1_path},{r2_path},{out_dir / 'w3_quant' / s}")

        with gzip.open(r1_path, "wt", compresslevel=1) as f1, gzip.open(
            r2_path, "wt", compresslevel=1
        ) as f2:
            rid = 0
            for gi, gene in enumerate(genes):
                psi_s = sample_psis[s][gi]
                weights = gene["iso_weights_fn"](psi_s)
                seqs = gene["iso_seqs"]
                # Weight sampling by (molecule_weight * effective_length)
                eff_w = [w * max(1, len(sq) - 200) for w, sq in zip(weights, seqs)]
                tot_w = sum(eff_w)
                probs = [x / tot_w for x in eff_w]
                cum = []
                acc = 0.0
                for p in probs:
                    acc += p
                    cum.append(acc)

                for _ in range(pairs_per_gene):
                    u = rng.random()
                    ti = 0
                    while ti < len(cum) - 1 and u > cum[ti]:
                        ti += 1
                    tx = seqs[ti]
                    tlen = len(tx)
                    flen = int(round(rng.gauss(240.0, 25.0)))
                    flen = max(read_len + 10, min(tlen, flen))
                    start = rng.randint(0, tlen - flen)
                    frag = tx[start : start + flen]
                    r1_seq = mutate_seq(rng, frag[:read_len], 0.002)
                    r2_seq = mutate_seq(rng, revcomp(frag[-read_len:]), 0.002)
                    f1.write(f"@{s}_{rid}/1\n{r1_seq}\n+\n{qual}\n")
                    f2.write(f"@{s}_{rid}/2\n{r2_seq}\n+\n{qual}\n")
                    rid += 1

    with open(out_dir / "w3_batch.csv", "w") as f:
        f.write("\n".join(manifest_rows) + "\n")

    print(
        f"Generated benchmark cohort in {out_dir}: {n_genes} genes, {len(samples)} samples, "
        f"{pairs_per_gene * n_genes} pairs/sample ({read_len} bp PE)."
    )


def pearson_r(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


def mae(xs: list[float], ys: list[float]) -> float:
    return sum(abs(x - y) for x, y in zip(xs, ys)) / max(1, len(xs))


def rmse(xs: list[float], ys: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(xs, ys)) / max(1, len(xs)))


def load_psi_by_coord(path: Path) -> dict[tuple[str, str, str], tuple[float, float]]:
    out = {}
    with gzip.open(path, "rt") as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            gene, coord, etype, psi_s, reads_s = parts[0], parts[2], parts[4], parts[5], parts[8]
            if psi_s == "NA" or reads_s == "NA":
                continue
            try:
                out[(gene, coord, etype)] = (float(psi_s), float(reads_s))
            except ValueError:
                continue
    return out


def load_diff_by_coord(path: Path) -> dict[tuple[str, str, str], dict[str, float]]:
    out = {}
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(hdr)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            key = (parts[0], parts[2], parts[4])
            try:
                dpsi = float(parts[idx["DeltaPsi"]])
                lor = float(parts[idx["LogOddsRatio"]])
                pval = float(parts[idx["Pvalue"]])
                fdr = float(parts[idx["FDR"]])
            except (ValueError, KeyError):
                continue
            if all(math.isfinite(v) for v in (dpsi, lor, pval, fdr)):
                out[key] = {"dpsi": dpsi, "lor": lor, "pval": pval, "fdr": fdr}
    return out


def compare_benchmark(args: argparse.Namespace) -> None:
    truth = {}
    with open(args.truth) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            key = (parts[0], parts[1], parts[2])
            truth[key] = {
                "cls": parts[3],
                "is_de": int(parts[4]),
                "psi_a": float(parts[5]),
                "psi_b": float(parts[6]),
                "dpsi": float(parts[7]),
                "psi_a1": float(parts[8]),
            }

    w3_a1 = load_psi_by_coord(Path(args.w3_dir) / "A1.psi.gz")
    w2_a1 = load_psi_by_coord(Path(args.w2_dir) / "A1.psi.gz")

    # 1. Sample A1 PSI accuracy vs Ground Truth
    t_w3_x, t_w3_y = [], []
    t_w2_x, t_w2_y = [], []
    for k, tr in truth.items():
        if k in w3_a1:
            t_w3_x.append(tr["psi_a1"])
            t_w3_y.append(w3_a1[k][0])
        if k in w2_a1:
            t_w2_x.append(tr["psi_a1"])
            t_w2_y.append(w2_a1[k][0])

    # 2. Head-to-head PSI concordance across all matched events in A1
    shared_psi = sorted(set(w3_a1.keys()) & set(w2_a1.keys()))
    h_w3 = [w3_a1[k][0] for k in shared_psi]
    h_w2 = [w2_a1[k][0] for k in shared_psi]

    # 3. Differential splicing GLM concordance & DE accuracy
    w3_diff = load_diff_by_coord(Path(args.w3_diff))
    w2_diff = load_diff_by_coord(Path(args.w2_diff))

    shared_diff = sorted(set(w3_diff.keys()) & set(w2_diff.keys()))
    d_w3_dpsi = [w3_diff[k]["dpsi"] for k in shared_diff]
    d_w2_dpsi = [w2_diff[k]["dpsi"] for k in shared_diff]
    d_w3_lor = [w3_diff[k]["lor"] for k in shared_diff]
    d_w2_lor = [w2_diff[k]["lor"] for k in shared_diff]

    def eval_de(diff_map):
        tp = fp = fn = tn = 0
        xs, ys = [], []
        for k, tr in truth.items():
            if k not in diff_map:
                continue
            r = diff_map[k]
            xs.append(tr["dpsi"])
            ys.append(r["dpsi"])
            called = r["fdr"] < 0.05 and abs(r["dpsi"]) >= 0.10
            if tr["is_de"]:
                if called:
                    tp += 1
                else:
                    fn += 1
            else:
                if called:
                    fp += 1
                else:
                    tn += 1
        sens = tp / max(1, tp + fn)
        fdr_emp = fp / max(1, tp + fp)
        return len(xs), pearson_r(xs, ys), mae(xs, ys), sens, fdr_emp, tp, fp, fn, tn

    w3_de = eval_de(w3_diff)
    w2_de = eval_de(w2_diff)

    print("### Accuracy & Concordance Summary\n")
    print("| Metric | Langur (Rust) | Whippet2-turbo (Julia) |")
    print("|---|---:|---:|")
    print(
        f"| A1 PSI vs Truth (n={len(t_w3_x)}/{len(t_w2_x)}) Pearson r | {pearson_r(t_w3_x, t_w3_y):.5f} | {pearson_r(t_w2_x, t_w2_y):.5f} |"
    )
    print(
        f"| A1 PSI vs Truth MAE | {mae(t_w3_x, t_w3_y):.5f} | {mae(t_w2_x, t_w2_y):.5f} |"
    )
    print(
        f"| A1 PSI vs Truth RMSE | {rmse(t_w3_x, t_w3_y):.5f} | {rmse(t_w2_x, t_w2_y):.5f} |"
    )
    print(
        f"| Delta dPSI vs Truth (n={w3_de[0]}/{w2_de[0]}) Pearson r | {w3_de[1]:.5f} | {w2_de[1]:.5f} |"
    )
    print(f"| Delta dPSI vs Truth MAE | {w3_de[2]:.5f} | {w2_de[2]:.5f} |")
    print(
        f"| DE Call Sensitivity (FDR<0.05, |dPSI|>=0.1) | {w3_de[3]*100:.1f}% ({w3_de[5]}/{w3_de[5]+w3_de[7]}) | {w2_de[3]*100:.1f}% ({w2_de[5]}/{w2_de[5]+w2_de[7]}) |"
    )
    print(
        f"| DE Empirical FDR (FP / Called) | {w3_de[4]*100:.2f}% ({w3_de[6]}/{w3_de[5]+w3_de[6]}) | {w2_de[4]*100:.2f}% ({w2_de[6]}/{w2_de[5]+w2_de[6]}) |"
    )
    print()
    print("### Head-to-Head Concordance (Langur vs Whippet2-turbo)\n")
    print(f"- **A1 PSI matched events ({len(shared_psi)})**: Pearson r = **{pearson_r(h_w3, h_w2):.5f}**, MAE = **{mae(h_w3, h_w2):.5f}**")
    print(
        f"- **Delta GLM matched events ({len(shared_diff)})**: DeltaPsi Pearson r = **{pearson_r(d_w3_dpsi, d_w2_dpsi):.5f}** (MAE = **{mae(d_w3_dpsi, d_w2_dpsi):.5f}**), LogOddsRatio Pearson r = **{pearson_r(d_w3_lor, d_w2_lor):.5f}**"
    )


def main() -> None:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("generate")
    g.add_argument("--out-dir", required=True)
    g.add_argument("--genes", type=int, default=200)
    g.add_argument("--pairs-per-sample", type=int, default=50000)
    g.add_argument("--samples-per-group", type=int, default=3)
    g.add_argument("--seed", type=int, default=42)

    c = sub.add_parser("compare")
    c.add_argument("--truth", required=True)
    c.add_argument("--w3-dir", required=True)
    c.add_argument("--w2-dir", required=True)
    c.add_argument("--w3-diff", required=True)
    c.add_argument("--w2-diff", required=True)

    args = p.parse_args()
    if args.cmd == "generate":
        generate_benchmark(args)
    else:
        compare_benchmark(args)


if __name__ == "__main__":
    main()
