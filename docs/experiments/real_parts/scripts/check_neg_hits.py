#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_neg_hits.py — 复核"bbox 不相交但 classify_pair 报 1D/有交"的可疑样本。

对给定 part + pairs csv 中的可疑行，重新解析 data/<part>.txt，用**直角坐标控制点**
直接计算保守 AABB（正权时凸包性质成立），检查 bbox 是否真的不相交；
同时检查 face 内 SURF 记录的 pi 是否与解析顺序一致（probe 用 (fid, 索引) 定位 patch，
occt_check 用 (fid, pi) 定位，两者必须一致）。
"""
import csv
import sys

RP = "/mnt/d/qqdownload/nurbs-certified-framework/docs/experiments/real_parts"


def parse(path):
    toks = open(path, "r", encoding="utf-8", errors="replace").read().split()
    faces = {}
    order_bad = []
    i, n = 0, len(toks)
    while i < n:
        if toks[i] == "NSURF":
            i += 2
        elif toks[i] == "SURF":
            fid, pi, du, dv, nu, nv = (int(toks[i + k]) for k in range(1, 7))
            i += 7
            cnt = nu * nv
            vals = list(map(float, toks[i:i + 4 * cnt]))
            i += 4 * cnt
            if i < n and toks[i] == "ENDSURF":
                i += 1
            pts = [(vals[4 * k], vals[4 * k + 1], vals[4 * k + 2], vals[4 * k + 3]) for k in range(cnt)]
            lst = faces.setdefault(fid, [])
            if pi != len(lst):
                order_bad.append((fid, pi, len(lst)))
            lst.append((pi, nu, nv, pts))
        else:
            i += 1
    return faces, order_bad


def bb(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; zs = [p[2] for p in pts]
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def overlap(a, b):
    return a[1] >= b[0] and b[1] >= a[0] and a[3] >= b[2] and b[3] >= a[2] and a[5] >= b[4] and b[5] >= a[4]


def main(part):
    data = "%s/data/%s.txt" % (RP, part)
    faces, bad = parse(data)
    print("part", part, "faces", len(faces), "pi_order_violations", len(bad), bad[:5])
    minw = 1e9; maxw = -1e9
    for fid, lst in faces.items():
        for (pi, nu, nv, pts) in lst:
            for p in pts:
                minw = min(minw, p[3]); maxw = max(maxw, p[3])
    print("  weight range: [%.6g, %.6g]" % (minw, maxw))
    pairs = "%s/results/probe/%s_pairs.csv" % (RP, part)
    nflag = 0
    for r in csv.DictReader(open(pairs, encoding="utf-8")):
        if r["src"] != "neg":
            continue
        a = faces[int(r["fidA"])][int(r["piA"])]
        b = faces[int(r["fidB"])][int(r["piB"])]
        ba, bb2 = bb(a[3]), bb(b[3])
        ov = overlap(ba, bb2)
        if ov or int(r["type"]) != 0:
            nflag += 1
            print("  NEG flag: pa=%s pb=%s type=%s bbox_overlap_cart=%s | A[%s] [%.4f,%.4f]x[%.4f,%.4f]x[%.4f,%.4f] B[%s] [%.4f,%.4f]x[%.4f,%.4f]x[%.4f,%.4f]"
                  % (r["piA"], r["piB"], r["type"], ov, r["fidA"], ba[0], ba[1], ba[2], ba[3], ba[4], ba[5],
                     r["fidB"], bb2[0], bb2[1], bb2[2], bb2[3], bb2[4], bb2[5]))
    print("  flagged", nflag)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        main(p)
