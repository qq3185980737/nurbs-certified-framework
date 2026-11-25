#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_cover.py — 需求 4：覆盖检验。对 OCCT 判定"有交线"的抽样 patch-pair，
用 cover_check 检查 OCCT 交线采样点是否落在认证 cover 盒的三维像内。

输出 results/cover/<part>_cover.csv 与 results/cover/<part>_pairs_in.csv
"""
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

ROOT = "/mnt/d/qqdownload/nurbs-certified-framework"
RP = os.path.join(ROOT, "docs/experiments/real_parts")
RES = os.path.join(RP, "results")
COVER_BIN = os.path.join(RP, "build/cover_check_v2")
PARTS = ["228", "387", "107", "10", "38", "411"]
TAU = 0.02
CAP = 4000000


def one(part):
    outdir = os.path.join(RES, "cover")
    os.makedirs(outdir, exist_ok=True)
    occ = os.path.join(RES, "probe", part + "_occt.csv")
    pts = os.path.join(RES, "probe", part + "_points.csv")
    if not (os.path.exists(occ) and os.path.exists(pts)):
        print("skip", part, "(no probe data)")
        return
    rows = [r for r in csv.DictReader(open(occ, encoding="utf-8")) if int(r["occt_lines"]) >= 1]
    if not rows:
        print("skip", part, "(no OCCT intersection sample)")
        return
    inp = os.path.join(outdir, part + "_pairs_in.csv")
    with open(inp, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["part", "fidA", "piA", "fidB", "piB", "src", "bbox_overlap"])
        for r in rows:
            w.writerow([part, r["fidA"], r["piA"], r["fidB"], r["piB"], r["src"], r["bbox_overlap"]])
    out = os.path.join(outdir, part + "_cover.csv")
    p = subprocess.run(["timeout", "300", COVER_BIN, os.path.join(RP, "data", part + ".txt"),
                        inp, pts, out, str(TAU), str(CAP)], capture_output=True, text=True)
    print("cover %-4s rc=%d %s" % (part, p.returncode, (p.stdout + p.stderr).strip()[:150]), flush=True)


def main():
    parts = sys.argv[1].split(",") if len(sys.argv) > 1 else PARTS
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(one, parts))


if __name__ == "__main__":
    sys.exit(main())
