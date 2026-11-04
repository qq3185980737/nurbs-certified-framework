#!/bin/bash
# run_one.sh NAME TAU NODE_CAP TIMEOUT_S [DETAIL]
# 单件有界运行 bench_real_batch，记录 RESULT / 墙钟 / detail 进度。
set -u
ROOT=/mnt/d/qqdownload/nurbs-certified-framework
RP=$ROOT/docs/experiments/real_parts
BIN=$ROOT/cpp/build/bench_real_batch_linux
NAME=$1; TAU=${2:-0.02}; CAP=${3:-4000000}; TO=${4:-180}; DET=${5:-/dev/shm/realparts_$NAME.detail}
mkdir -p "$RP/results/logs"
LOG=$RP/results/logs/${NAME}_tau${TAU}_cap${CAP}.log
rm -f "$DET"
T0=$(date +%s%N)
timeout "$TO" "$BIN" "$RP/data/$NAME.txt" "$TAU" "$CAP" "$DET" > "$LOG" 2>&1
RC=$?
T1=$(date +%s%N)
MS=$(( (T1 - T0) / 1000000 ))
DETN=0; [ -f "$DET" ] && DETN=$(wc -l < "$DET")
printf 'META name=%s rc=%d wall_ms=%d timeout_s=%d tau=%s node_cap=%s detail_lines=%d\n' \
  "$NAME" "$RC" "$MS" "$TO" "$TAU" "$CAP" "$DETN" >> "$LOG"
echo "DONE name=$NAME rc=$RC wall_ms=$MS detail_lines=$DETN log=$LOG"
