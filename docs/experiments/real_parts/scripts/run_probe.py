#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_probe.py — 对 20 个零件跑 patch-pair 抽样探针 + OCCT GeomAPI_IntSS 交叉校验。

每个零件：
  1) probe_pairs：count face-pair/patch-pair/bbox 统计 + 抽样 50 个 bbox 相交 patch-pair
     （正样本，classify_pair 逐对计时/分类）+ 50 个 bbox 不相交 patch-pair（阴性对照）
  2) occt_check.py：对同样 100 对做 GeomAPI_IntSS，输出交线数 + 交线采样点
输出：results/probe/<part>_pairs.csv / <part>_occt.csv / <part>_points.csv / <part>_probe.log
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

ROOT = "/mnt/d/qqdownload/nurbs-certified-framework"
RP = os.path.join(ROOT, "docs/experiments/real_parts")
PROBE_BIN = os.path.join(RP, "build/probe_pairs_v2")
OCC_PY = os.path.join(RP, "scripts/occ_py.sh")
OCC_SCRIPT = os.path.join(RP, "scripts/occt_check.py")
OUTDIR = os.path.join(RP, "results/probe")

PARTS = ["248", "115", "69", "387", "228", "411", "506", "70", "538", "570",
         "113", "205", "546", "107", "534", "200", "614", "38", "10", "136"]

TAU = 0.02
CAP = 4000000
NPOS = 50
NNEG = 50


def one(part):
    os.makedirs(OUTDIR, exist_ok=True)
    pairs = os.path.join(OUTDIR, part + "_pairs.csv")
    occt = os.path.join(OUTDIR, part + "_occt.csv")
    points = os.path.join(OUTDIR, part + "_points.csv")
    log = os.path.join(OUTDIR, part + "_probe.log")
    data = os.path.join(RP, "data", part + ".txt")
    seed = 1000 + int(part)
    t0 = time.time()
    p = subprocess.run(["timeout", "240", PROBE_BIN, data, pairs, str(TAU), str(CAP),
                        str(NPOS), str(NNEG), str(seed)], capture_output=True, text=True)
    with open(log, "w", encoding="utf-8") as fh:
        fh.write(p.stdout + p.stderr)
    n_pairs = 0
    if os.path.exists(pairs):
        with open(pairs, "r", encoding="utf-8") as fh:
            n_pairs = sum(1 for _ in fh) - 1
    t_probe = time.time() - t0
    t1 = time.time()
    q = subprocess.run([OCC_PY, OCC_SCRIPT, "--data", data, "--pairs", pairs,
                        "--out", occt, "--points", points, "--part", part],
                       capture_output=True, text=True)
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(q.stdout + q.stderr)
    print("probe %-4s rc=%d pairs=%d occt_rc=%d probe_s=%.1f occt_s=%.1f"
          % (part, p.returncode, n_pairs, q.returncode, t_probe, time.time() - t1), flush=True)
    return part, p.returncode, n_pairs, q.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--parts", default="")
    args = ap.parse_args()
    parts = args.parts.split(",") if args.parts else PARTS
    os.makedirs(OUTDIR, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        list(ex.map(one, parts))


if __name__ == "__main__":
    sys.exit(main())
