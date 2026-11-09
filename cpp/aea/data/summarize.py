#!/usr/bin/env python3
"""Aggregate RESULT lines produced by bench_cert_step into per-family totals."""
import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

KEYS = ("face-pairs", "disjoint", "1D_certified", "2D_coincidence",
        "coplanar_contact", "tangent_candidate", "cover_limit",
        "patch_boundary", "undecided", "adjacent_skipped")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=Path("results/e11_bench_summary.tsv"))
    args = ap.parse_args()
    agg = defaultdict(Counter)
    models = Counter()
    for line in args.input.read_text(encoding="utf-8-sig", errors="ignore").splitlines():
        if "RESULT" not in line:
            continue
        family = line.split("\t", 1)[0].replace("_data", "")
        models[family] += 1
        for k, v in re.findall(r"([A-Za-z_0-9-]+)=(-?\d+)", line):
            if k in KEYS:
                agg[family][k] += int(v)
    for family, c in agg.items():
        print(f"{family}: models={models[family]} {dict(c)}")


if __name__ == "__main__":
    main()
