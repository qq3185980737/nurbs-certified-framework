#!/bin/bash
# status.sh — 查看批量任务进度
RP=/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/real_parts
cd "$RP/results/batch" 2>/dev/null || exit 1
for g in A H B C; do
  n=$(ls "$g"/*.log 2>/dev/null | wc -l)
  echo "group $g: logs=$n"
done
echo "--- running bench processes ---"
pgrep -af bench_real | sed 's/.*bench_real/bench_real/' | head -20
echo "--- load ---"
uptime
