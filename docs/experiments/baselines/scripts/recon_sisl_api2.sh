#!/bin/bash
set -u
SRC=/mnt/d/sisl_dl/src/SISL-master
echo "=== doc/manual ==="
ls -la $SRC/doc/manual | head -30
echo "=== grep intersect-ish names in sisl.h ==="
grep -n -i "intersect\|int_" $SRC/include/sisl.h | head -60
echo "=== all SISL functions with 'int' in purpose (s1xxx) ==="
grep -n "^extern.*SISLint\|SISLint.*(" $SRC/include/sisl.h | head -60
echo "=== src files named int* ==="
ls $SRC/src | grep -i "^int\|surf\|sf" | head -40
echo "=== count src ==="
ls $SRC/src | wc -l
