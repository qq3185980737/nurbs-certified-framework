#!/bin/bash
# Compile & run the SISL SSI baseline harness.
# Requires libsisl.a built by scripts/build_sisl.sh
set -uo pipefail
ROOT=/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines
SISLSRC=/mnt/d/sisl_dl/src/SISL-master
LIB=$ROOT/sisl_build/libsisl.a
OUT=$ROOT/sisl_build
mkdir -p "$OUT" "$ROOT/results"

if [ ! -f "$LIB" ]; then
  echo "ERROR: $LIB missing, run scripts/build_sisl.sh first" >&2
  exit 1
fi

g++ -O2 -std=c++17 -Wall -I"$SISLSRC/include" \
    "$ROOT/src/sisl_ssi_harness.cpp" "$LIB" -lm -o "$OUT/sisl_ssi_harness" \
    2>"$ROOT/results/sisl_ssi_compile.err" || {
      echo "COMPILE FAILED:"; cat "$ROOT/results/sisl_ssi_compile.err"; exit 1; }
echo "compiled: $OUT/sisl_ssi_harness (warnings: $(wc -l < "$ROOT/results/sisl_ssi_compile.err"))"

echo "=== selftest ==="
"$OUT/sisl_ssi_harness" --selftest || { echo "SELFTEST FAILED"; exit 1; }

echo "=== fixtures A-E (surface/surface s1859+s1310) + F (surface/curve s1858) ==="
timeout 600 "$OUT/sisl_ssi_harness" > "$ROOT/results/sisl_ssi_raw.txt" 2>&1
echo "exit=$?"
# SISL prints debug noise ("Zone equality") from s1310 on degenerate input; keep
# the CSV rows only.
grep -E '^(fixture|[A-F],)' "$ROOT/results/sisl_ssi_raw.txt" > "$ROOT/results/sisl_ssi.csv"
cat "$ROOT/results/sisl_ssi.csv"

echo
echo "=== detail (guide points / traced curve endpoints / residuals), epsge=1e-6 ==="
timeout 600 "$OUT/sisl_ssi_harness" --detail > "$ROOT/results/sisl_ssi_detail.txt" 2>&1
echo "exit=$?"
grep -v 'Zone equality' "$ROOT/results/sisl_ssi_detail.txt" | head -80
