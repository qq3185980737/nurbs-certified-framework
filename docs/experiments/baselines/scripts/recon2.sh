#!/bin/bash
set -u
echo "=== /mnt/d/sisl_dl ==="
ls -la /mnt/d/sisl_dl 2>/dev/null | head -30
du -sh /mnt/d/sisl_dl 2>/dev/null
find /mnt/d/sisl_dl -maxdepth 2 -type d 2>/dev/null | head -40

echo "=== occ env numpy ==="
/root/miniconda3/envs/occ/bin/python - <<'PY'
import numpy
print("numpy", numpy.__version__)
try:
    import scipy; print("scipy", scipy.__version__)
except Exception as e:
    print("scipy MISSING:", e)
PY

echo "=== apt / git ==="
which apt-get curl wget unzip tar 2>/dev/null
id -u

echo "=== net ==="
timeout 15 curl -sI https://github.com 2>&1 | head -3
echo "curl-exit=$?"
timeout 20 bash -c 'curl -s https://pypi.org/simple/ -o /dev/null -w "pypi:%{http_code}\n"' 2>&1
