#!/usr/bin/env bash
# Head-to-head end-to-end benchmark: Langur (Rust) vs Whippet2-turbo (Julia 1.9)
# Measures wall time, index load time, reads/sec, peak RSS (/usr/bin/time -l),
# index disk size, PSI accuracy/concordance, and Beta-Binomial GLM delta concordance.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WHIPPET2_DIR="${WHIPPET2_DIR:-/Users/david/Bioinformatics/Whippet2.jl}"
JULIA_BIN="${JULIA_BIN:-julia}"
if [[ -z "${CARGO_BIN:-}" ]]; then
  if command -v cargo >/dev/null 2>&1; then
    CARGO_BIN="cargo"
  else
    CARGO_BIN="$HOME/.cargo/bin/cargo"
  fi
fi
WORK_DIR="${WORK_DIR:-$ROOT_DIR/target/bench_vs_whippet2}"
GENES="${GENES:-200}"
PAIRS_PER_SAMPLE="${PAIRS_PER_SAMPLE:-50000}"
SAMPLES_PER_GROUP="${SAMPLES_PER_GROUP:-3}"
THREADS="${THREADS:-4}"

mkdir -p "$WORK_DIR/w3_quant" "$WORK_DIR/w2_quant"

W3_BIN="$ROOT_DIR/target/release/langur"
if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  echo "==> Building langur release binary..."
  "$CARGO_BIN" build --offline --release --manifest-path "$ROOT_DIR/Cargo.toml"
fi

if [[ "${SKIP_GEN:-0}" != "1" ]]; then
  echo "==> Generating synthetic multi-gene cohort ($GENES genes, ${SAMPLES_PER_GROUP}v${SAMPLES_PER_GROUP}, $PAIRS_PER_SAMPLE pairs/sample)..."
  python3 "$ROOT_DIR/scripts/bench_helper.py" generate \
    --out-dir "$WORK_DIR" \
    --genes "$GENES" \
    --pairs-per-sample "$PAIRS_PER_SAMPLE" \
    --samples-per-group "$SAMPLES_PER_GROUP" \
    --seed 42
fi

# Helper to measure wall time (sec) and Peak RSS (MB) on macOS (/usr/bin/time -l) or Linux (-v)
measure_cmd() {
  local label="$1"
  local logfile="$2"
  shift 2
  local t0 t1 wall rss_bytes rss_mb
  t0=$(python3 -c 'import time; print(time.perf_counter())')
  if [[ "$(uname -s)" == "Darwin" ]]; then
    /usr/bin/time -l "$@" >"$logfile.out" 2>"$logfile.err"
    rss_bytes=$(awk '/maximum resident set size/ {print $1; exit}' "$logfile.err")
    rss_mb=$(python3 -c "print(f'{int(${rss_bytes:-0})/1048576.0:.1f}')")
  else
    /usr/bin/time -v "$@" >"$logfile.out" 2>"$logfile.err"
    rss_kb=$(awk -F: '/Maximum resident set size/ {gsub(/ /,"",$2); print $2; exit}' "$logfile.err")
    rss_mb=$(python3 -c "print(f'{int(${rss_kb:-0})/1024.0:.1f}')")
  fi
  t1=$(python3 -c 'import time; print(time.perf_counter())')
  wall=$(python3 -c "print(f'{$t1 - $t0:.3f}')")
  echo "$wall $rss_mb"
}

echo "==> [1/6] Indexing with Langur..."
read -r W3_IDX_WALL W3_IDX_RSS < <(
  measure_cmd "w3_index" "$WORK_DIR/w3_index" \
    "$W3_BIN" -t "$THREADS" index \
      --fasta "$WORK_DIR/genome.fa" \
      --gtf "$WORK_DIR/genes.gtf" \
      -x "$WORK_DIR/w3_graph"
)
W3_IDX_KB=$(du -k "$WORK_DIR/w3_graph.lidx" | awk '{print $1}')

echo "==> [2/6] Indexing with Whippet2-turbo..."
read -r W2_IDX_WALL W2_IDX_RSS < <(
  measure_cmd "w2_index" "$WORK_DIR/w2_index" \
    "$JULIA_BIN" --project="$WHIPPET2_DIR" "$WHIPPET2_DIR/bin/whippet-index.jl" \
      --fasta "$WORK_DIR/genome.fa" \
      --gtf "$WORK_DIR/genes.gtf" \
      -x "$WORK_DIR/w2_graph"
)
W2_IDX_KB=$(du -k "$WORK_DIR/w2_graph.jls" | awk '{print $1}')

echo "==> [3/6] Quantifying sample A1 (single-sample timing & index load)..."
read -r W3_Q1_WALL W3_Q1_RSS < <(
  measure_cmd "w3_q1" "$WORK_DIR/w3_q1" \
    "$W3_BIN" -t "$THREADS" quant \
      "$WORK_DIR/A1_R1.fq.gz" "$WORK_DIR/A1_R2.fq.gz" \
      -x "$WORK_DIR/w3_graph.lidx" \
      -o "$WORK_DIR/w3_quant/A1"
)
W3_LOAD_MS=$(grep -oE 'loaded .* in [0-9.]+ ms' "$WORK_DIR/w3_q1.err" | awk '{print $(NF-1)}' || echo "NA")

read -r W2_Q1_WALL W2_Q1_RSS < <(
  measure_cmd "w2_q1" "$WORK_DIR/w2_q1" \
    "$JULIA_BIN" --project="$WHIPPET2_DIR" -t "$THREADS" "$WHIPPET2_DIR/bin/whippet-quant.jl" \
      "$WORK_DIR/A1_R1.fq.gz" "$WORK_DIR/A1_R2.fq.gz" \
      -x "$WORK_DIR/w2_graph.jls" \
      -o "$WORK_DIR/w2_quant/A1"
)
W2_LOAD_S=$(grep -m1 -oE '[0-9.]+ seconds' "$WORK_DIR/w2_q1.err" | awk '{print $1}' || echo "NA")

echo "==> [4/6] Quantifying full $((SAMPLES_PER_GROUP * 2))-sample cohort..."
read -r W3_QALL_WALL W3_QALL_RSS < <(
  measure_cmd "w3_qall" "$WORK_DIR/w3_qall" \
    "$W3_BIN" -t "$THREADS" quant \
      --batch "$WORK_DIR/w3_batch.csv" \
      -x "$WORK_DIR/w3_graph.lidx"
)

W2_QALL_T0=$(python3 -c 'import time; print(time.perf_counter())')
W2_QALL_RSS="$W2_Q1_RSS"
A_LIST="A1.psi.gz"
B_LIST=""
for ((i=2; i<=SAMPLES_PER_GROUP; i++)); do
  A_LIST="${A_LIST},A${i}.psi.gz"
  read -r _s_wall _s_rss < <(
    measure_cmd "w2_q_A${i}" "$WORK_DIR/w2_q_A${i}" \
      "$JULIA_BIN" --project="$WHIPPET2_DIR" -t "$THREADS" "$WHIPPET2_DIR/bin/whippet-quant.jl" \
        "$WORK_DIR/A${i}_R1.fq.gz" "$WORK_DIR/A${i}_R2.fq.gz" \
        -x "$WORK_DIR/w2_graph.jls" \
        -o "$WORK_DIR/w2_quant/A${i}"
  )
  W2_QALL_RSS=$(python3 -c "print(f'{max(float($W2_QALL_RSS), float($_s_rss)):.1f}')")
done
for ((i=1; i<=SAMPLES_PER_GROUP; i++)); do
  [[ -n "$B_LIST" ]] && B_LIST="${B_LIST},"
  B_LIST="${B_LIST}B${i}.psi.gz"
  read -r _s_wall _s_rss < <(
    measure_cmd "w2_q_B${i}" "$WORK_DIR/w2_q_B${i}" \
      "$JULIA_BIN" --project="$WHIPPET2_DIR" -t "$THREADS" "$WHIPPET2_DIR/bin/whippet-quant.jl" \
        "$WORK_DIR/B${i}_R1.fq.gz" "$WORK_DIR/B${i}_R2.fq.gz" \
        -x "$WORK_DIR/w2_graph.jls" \
        -o "$WORK_DIR/w2_quant/B${i}"
  )
  W2_QALL_RSS=$(python3 -c "print(f'{max(float($W2_QALL_RSS), float($_s_rss)):.1f}')")
done
W2_QALL_T1=$(python3 -c 'import time; print(time.perf_counter())')
W2_QALL_WALL=$(python3 -c "print(f'{($W2_QALL_T1 - $W2_QALL_T0) + $W2_Q1_WALL:.3f}')")

echo "==> [5/6] Running differential splicing Beta-Binomial GLM (delta)..."
read -r W3_D_WALL W3_D_RSS < <(
  measure_cmd "w3_delta" "$WORK_DIR/w3_delta" \
    "$W3_BIN" -t "$THREADS" delta \
      -a "$A_LIST" \
      -b "$B_LIST" \
      -d "$WORK_DIR/w3_quant" \
      -o "$WORK_DIR/w3_diff"
)

read -r W2_D_WALL W2_D_RSS < <(
  measure_cmd "w2_delta" "$WORK_DIR/w2_delta" \
    "$JULIA_BIN" --project="$WHIPPET2_DIR" -t "$THREADS" "$WHIPPET2_DIR/bin/whippet-delta.jl" \
      -a "$A_LIST" \
      -b "$B_LIST" \
      -d "$WORK_DIR/w2_quant" \
      -o "$WORK_DIR/w2_diff"
)

TOTAL_READS_1=$((PAIRS_PER_SAMPLE * 2))
TOTAL_READS_ALL=$((TOTAL_READS_1 * SAMPLES_PER_GROUP * 2))
W3_RPS_1=$(python3 -c "print(f'{$TOTAL_READS_1 / max(float($W3_Q1_WALL), 1e-6):,.0f}')")
W2_RPS_1=$(python3 -c "print(f'{$TOTAL_READS_1 / max(float($W2_Q1_WALL), 1e-6):,.0f}')")
W3_RPS_ALL=$(python3 -c "print(f'{$TOTAL_READS_ALL / max(float($W3_QALL_WALL), 1e-6):,.0f}')")
W2_RPS_ALL=$(python3 -c "print(f'{$TOTAL_READS_ALL / max(float($W2_QALL_WALL), 1e-6):,.0f}')")

echo ""
echo "==================================================================="
echo "  Langur (Rust) vs Whippet2-turbo (Julia) — Benchmark Report"
echo "  Cohort: $GENES genes, ${SAMPLES_PER_GROUP}v${SAMPLES_PER_GROUP} PE samples ($PAIRS_PER_SAMPLE pairs/sample), $THREADS threads"
echo "==================================================================="
echo ""
echo "### Runtime, Throughput & Memory"
echo ""
echo "| Stage | Langur (Rust) | Whippet2-turbo (Julia) | Speedup |"
echo "|---|---:|---:|---:|"
python3 - <<PY
def sp(w2, w3):
    return f"{float(w2)/max(float(w3), 1e-6):.1f}x"
print(f"| index wall time | ${W3_IDX_WALL} s | ${W2_IDX_WALL} s | {sp(${W2_IDX_WALL}, ${W3_IDX_WALL})} |")
print(f"| index peak RSS | ${W3_IDX_RSS} MB | ${W2_IDX_RSS} MB | {sp(${W2_IDX_RSS}, ${W3_IDX_RSS})} lower |")
print(f"| index size on disk | ${W3_IDX_KB} KB | ${W2_IDX_KB} KB | — |")
print(f"| index load time (quant A1) | ${W3_LOAD_MS} ms | ${W2_LOAD_S} s (startup+load) | — |")
print(f"| quant A1 wall time (1 sample) | ${W3_Q1_WALL} s | ${W2_Q1_WALL} s | {sp(${W2_Q1_WALL}, ${W3_Q1_WALL})} |")
print(f"| quant A1 throughput (reads/s) | ${W3_RPS_1} | ${W2_RPS_1} | {sp(${W2_Q1_WALL}, ${W3_Q1_WALL})} |")
print(f"| quant A1 peak RSS | ${W3_Q1_RSS} MB | ${W2_Q1_RSS} MB | {sp(${W2_Q1_RSS}, ${W3_Q1_RSS})} lower |")
print(f"| quant 6-sample cohort wall time | ${W3_QALL_WALL} s (--batch) | ${W2_QALL_WALL} s (6x CLI) | {sp(${W2_QALL_WALL}, ${W3_QALL_WALL})} |")
print(f"| quant 6-sample throughput (reads/s) | ${W3_RPS_ALL} | ${W2_RPS_ALL} | {sp(${W2_QALL_WALL}, ${W3_QALL_WALL})} |")
print(f"| quant 6-sample peak RSS | ${W3_QALL_RSS} MB | ${W2_QALL_RSS} MB | {sp(${W2_QALL_RSS}, ${W3_QALL_RSS})} lower |")
print(f"| delta GLM wall time | ${W3_D_WALL} s | ${W2_D_WALL} s | {sp(${W2_D_WALL}, ${W3_D_WALL})} |")
print(f"| delta GLM peak RSS | ${W3_D_RSS} MB | ${W2_D_RSS} MB | {sp(${W2_D_RSS}, ${W3_D_RSS})} lower |")
PY
echo ""
echo "==> [6/6] Evaluating PSI & Delta GLM concordance..."
python3 "$ROOT_DIR/scripts/bench_helper.py" compare \
  --truth "$WORK_DIR/truth.tsv" \
  --w3-dir "$WORK_DIR/w3_quant" \
  --w2-dir "$WORK_DIR/w2_quant" \
  --w3-diff "$WORK_DIR/w3_diff.diff.gz" \
  --w2-diff "$WORK_DIR/w2_diff.diff.gz"
