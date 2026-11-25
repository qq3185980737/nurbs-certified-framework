#!/bin/bash
# progress.sh — 看 /dev/shm 里正在跑的 detail 行数（= 已处理 face-pair 数）
for f in /dev/shm/rp_*.detail; do
  [ -e "$f" ] || continue
  n=$(wc -l < "$f" 2>/dev/null)
  echo "$(basename "$f") lines=$n"
done
df -h /dev/shm | tail -1
