#!/bin/bash
# Build and run the SISL failure-mode probes.
#   part 1: multi-root curve/curve (s1857) + self-intersection probe
#   part 2: tolerance behaviour on the degenerate tangential fixture D (bounded)
set -uo pipefail
ROOT=/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines
SRC=/mnt/d/sisl_dl/src/SISL-master
OUT=$ROOT/sisl_build
RES=$ROOT/results
mkdir -p "$OUT" "$RES"

g++ -O2 -std=c++17 -I"$SRC/include" "$ROOT/src/sisl_failure_probe.cpp" \
    "$OUT/libsisl.a" -lm -o "$OUT/sisl_failure_probe" \
    2>"$RES/sisl_failure_compile.err" || {
  echo "COMPILE FAILED"; cat "$RES/sisl_failure_compile.err"; exit 1; }

{
  echo "###### part 1: multi-root (s1857) + self-intersection ######"
  timeout 120 "$OUT/sisl_failure_probe" multi 2>&1 | grep -v 'Zone equality'
  timeout 120 "$OUT/sisl_failure_probe" self 2>&1 | grep -v 'Zone equality'
  echo
  echo "###### part 2: fixture D (degenerate tangential) vs epsge ######"
  for e in 1e-4 1e-6; do
    timeout 600 "$OUT/sisl_failure_probe" D "$e" 2>&1 | grep -v 'Zone equality'
    echo "  [exit=$?]"
  done
} | tee "$RES/sisl_failure_probe.txt"
