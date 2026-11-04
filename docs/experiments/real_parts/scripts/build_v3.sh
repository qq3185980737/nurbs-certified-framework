#!/bin/bash
# build_v3.sh — 修正 AABB 预筛 bug 后重新编译插桩基准 / 探针 / 调试工具
set -e
ROOT=/mnt/d/qqdownload/nurbs-certified-framework
RP=$ROOT/docs/experiments/real_parts
cd "$ROOT/cpp"
g++ -O2 -std=c++17 -Iinclude "$RP/scripts/bench_real_batch_probe.cpp" -o "$RP/build/bench_real_batch_probe_v3"
echo OK bench_real_batch_probe_v3
g++ -O2 -std=c++17 -Iinclude "$RP/scripts/probe_pairs.cpp" -o "$RP/build/probe_pairs_v2"
echo OK probe_pairs_v2
g++ -O2 -std=c++17 -Iinclude "$RP/scripts/debug_pair.cpp" -o "$RP/build/debug_pair"
echo OK debug_pair
"$RP/build/debug_pair" "$RP/data/113.txt" 96 0 117 0
