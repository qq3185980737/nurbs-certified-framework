#!/bin/bash
set -u
cd /mnt/d/sisl_dl
echo "=== build dir ==="
ls -la build | head -20
find build -maxdepth 3 -name "*.a" -o -maxdepth 3 -name "*.so" 2>/dev/null | head -20
echo "=== src top ==="
ls src/SISL-master | head -20
echo "=== CMakeLists ==="
head -50 src/SISL-master/CMakeLists.txt 2>/dev/null
echo "=== src subdirs ==="
ls src/SISL-master/src 2>/dev/null | head; ls src/SISL-master/src/* 2>/dev/null | head -20
