#!/bin/bash
# run_verify_extract.sh — 快速核验：对 2 个零件重新抽取并与现存 data/*.txt 逐项比对。
# 用法: wsl -e bash /mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/real_parts/scripts/run_verify_extract.sh
set -u
ROOT=/mnt/d/qqdownload/nurbs-certified-framework
RP=$ROOT/docs/experiments/real_parts
PY=/root/miniconda3/envs/occ/bin/python
OUT=$RP/results/extract_verify
mkdir -p "$OUT"
for name in "$@"; do
  echo "=========== $name ==========="
  echo "--- variant A: NurbsConvert ON ---"
  timeout 300 "$PY" "$RP/scripts/extract_simjeb.py" \
      --stp "$RP/stp/$name.stp" --out "$OUT/${name}_A.txt" --meta "$OUT/${name}_meta.csv"
  echo "rc=$?"
  echo "--- variant B: NurbsConvert OFF ---"
  timeout 300 "$PY" "$RP/scripts/extract_simjeb.py" \
      --stp "$RP/stp/$name.stp" --out "$OUT/${name}_B.txt" --no-nurbs-convert
  echo "rc=$?"
  echo "--- compare A vs data/$name.txt ---"
  timeout 300 python3 "$RP/scripts/compare_extract.py" "$RP/data/$name.txt" "$OUT/${name}_A.txt" || true
  echo "--- compare B vs data/$name.txt ---"
  timeout 300 python3 "$RP/scripts/compare_extract.py" "$RP/data/$name.txt" "$OUT/${name}_B.txt" || true
done
