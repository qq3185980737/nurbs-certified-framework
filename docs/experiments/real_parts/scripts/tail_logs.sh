#!/bin/bash
# tail_logs.sh — 快速查看已完成 log 的 RESULT / PAIRSTAT / META
RP=/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/real_parts
cd "$RP/results/batch" || exit 1
for g in A H B; do
  echo "================ group $g ================"
  for f in "$g"/*.log; do
    [ -e "$f" ] || continue
    echo "--- $(basename "$f")"
    grep -E "^(parsed|RESULT|PAIRSTAT|META)" "$f" 2>/dev/null | sed 's/  */ /g' || echo "  (no result: timed out early)"
  done
done
