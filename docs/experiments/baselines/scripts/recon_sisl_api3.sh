#!/bin/bash
set -u
SRC=/mnt/d/sisl_dl/src/SISL-master
echo "=== build state ==="
ls -la /mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines/
ls /mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines/sisl_build 2>/dev/null | head
echo "=== manual func dir sample ==="
ls $SRC/doc/manual/func | head -20
ls $SRC/doc/manual/func | wc -l
echo "=== func docs mentioning surface intersection ==="
grep -ril "surface.*intersection\|intersection.*surface" $SRC/doc/manual/func | head -30
echo "=== sislP.h decls containing 'int' word ==="
grep -n "SISLint\|s9\|s6" $SRC/include/sislP.h | head -40
echo "=== app/examples dirs ==="
ls $SRC/app | head; echo ---; ls $SRC/examples | head -30
