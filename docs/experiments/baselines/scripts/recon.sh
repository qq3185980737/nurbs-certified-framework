#!/bin/bash
# Recon for SISL / Busé-Mourrain baselines
set -u
cd /mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/baselines

echo "=== OS / toolchain ==="
cmake --version | head -1
g++ --version | head -1
nproc
git --version

echo "=== pythons ==="
python3 -V
python3 - <<'PY'
try:
    import numpy
    print("system numpy", numpy.__version__)
except Exception as e:
    print("system numpy MISSING", e)
PY
/root/miniconda3/envs/occ/bin/python - <<'PY'
import numpy, scipy
print("occ-env numpy", numpy.__version__, "scipy", scipy.__version__)
import OCC
print("OCC", OCC.VERSION)
PY

echo "=== net check github ==="
timeout 25 git ls-remote --heads https://github.com/SINTEF/SISL 2>&1 | head -5
echo "git-ls-remote-exit=$?"

echo "=== existing sisl? ==="
ls -d sisl_src/* 2>/dev/null | head
find / -maxdepth 4 -iname "*sisl*" -type d 2>/dev/null | head
