#!/bin/bash
# count_compare.sh — 对比"原 bench AABB（负坐标有反转 bug）"与"修正 AABB"的
#   face-pair / patch-pair / bbox 预筛数量。只做计数，不调用 classify（npos=nneg=0）。
RP=/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/real_parts
OUT="$RP/results/bbox_bug_counts.txt"
: > "$OUT"
for name in 228 387 107 411 136 10 38 248 115; do
  echo "### part $name" >> "$OUT"
  /usr/bin/timeout 300 "$RP/build/probe_pairs_v2" "$RP/data/$name.txt" /tmp/empty_pairs.csv \
      0.02 4000000 0 0 7 2>&1 | grep '^PROBE' >> "$OUT"
done
cat "$OUT"
