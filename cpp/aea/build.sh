#!/usr/bin/env bash
# Build the AEA benchmark and the C++/GMP certificate verifier from the canonical
# sources in cpp/include and cpp/tests.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p build
CXX="${CXX:-g++}"
$CXX -std=c++17 -O2 -I ../include ../tests/bench_cert_step.cpp -o build/bench_cert_step
$CXX -std=c++17 -O2 ../tests/cert_verify.cpp -lgmp -o build/cert_verify
echo "built: $(pwd)/build/bench_cert_step"
echo "built: $(pwd)/build/cert_verify"
