#!/bin/bash
set -u
SRC=/mnt/d/sisl_dl/src/SISL-master
echo "=== version info ==="
head -20 $SRC/ChangeLog
grep -m5 -i "version" $SRC/README
echo "=== doc dir ==="
ls $SRC/doc
echo "=== header size ==="
wc -l $SRC/include/*.h
echo "=== SISL surface-surface / intersection decls ==="
grep -n "SISLint" $SRC/include/sisl.h | head -80
