#!/usr/bin/env bash
# Run bench_cert_step over every .txt case and store a RESULT summary.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CASES="${1:?usage: run_bench.sh <case-dir> [binary] [tau] [out.tsv]}"
BIN="${2:-$ROOT/build/bench_cert_step}"
TAU="${3:-0.02}"
OUT="${4:-$ROOT/results/bench_summary.tsv}"

mkdir -p "$(dirname "$OUT")"
: > "$OUT"
for f in "$CASES"/*.txt; do
  [ -e "$f" ] || continue
  res="$("$BIN" "$f" "$TAU" 2>/dev/null | grep 'RESULT' || true)"
  printf "%s\t%s\t%s\n" "$(basename "$CASES")" "$(basename "$f")" "$res" >> "$OUT"
done
echo "wrote $OUT"
