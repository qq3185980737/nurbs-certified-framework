#!/bin/bash
# Build SISL static lib with gcc in our own build dir
set -x
SRC=/mnt/d/sisl_dl/src/SISL-master
BLD=/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines/sisl_build
mkdir -p "$BLD"
cd "$BLD" || exit 1
if command -v ninja >/dev/null 2>&1; then
  GEN="-G Ninja"
else
  GEN=""
fi
cmake -S "$SRC" -B "$BLD" $GEN \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_FLAGS="-O2 -fPIC" \
  -Dsisl_COMPILE_APPS=OFF -Dsisl_COMPILE_EXAMPLES=OFF -Dsisl_COMPILE_VIEWER=OFF \
  2>&1 | tail -20
time cmake --build "$BLD" -j 8 2>&1 | tail -30
ls -la "$BLD"/libsisl.a 2>/dev/null
