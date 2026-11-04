#!/bin/bash
# build.sh — 编译真实零件基准所需二进制（全部输出到 docs/experiments/real_parts/build/）
#   bench_real_batch_linux  : 原始基准（工作区 cert_ssi.hpp，与历史 logs 同一程序）
#   bench_real_batch_probe  : 插桩副本（PAIRSTAT：中位/最大单对耗时、node_cap 触顶、pair 级 type 分布）
#   probe_pairs             : patch-pair 抽样探针（供 OCCT 交叉校验）
#   bench_real_batch_head   : HEAD 基线版（git HEAD 的 cert_ssi.hpp，3 参 classify_pair）
set -e
ROOT=/mnt/d/qqdownload/nurbs-certified-framework
RP=$ROOT/docs/experiments/real_parts
mkdir -p "$RP/build" "$RP/results/head_baseline/include/nurbs"
cd "$ROOT/cpp"

g++ -O2 -std=c++17 -Iinclude tests/bench_real_batch.cpp -o build/bench_real_batch_linux
echo "OK bench_real_batch_linux"

g++ -O2 -std=c++17 -Iinclude "$RP/scripts/bench_real_batch_probe.cpp" -o "$RP/build/bench_real_batch_probe"
echo "OK bench_real_batch_probe"

g++ -O2 -std=c++17 -Iinclude "$RP/scripts/probe_pairs.cpp" -o "$RP/build/probe_pairs"
echo "OK probe_pairs"

git -C "$ROOT" show HEAD:cpp/include/nurbs/cert_ssi.hpp > "$RP/results/head_baseline/include/nurbs/cert_ssi.hpp"
g++ -O2 -std=c++17 -I"$RP/results/head_baseline/include" -Iinclude \
    "$RP/scripts/bench_real_batch_head.cpp" -o "$RP/build/bench_real_batch_head"
echo "OK bench_real_batch_head (= git HEAD cert_ssi.hpp)"
echo BUILD_OK
