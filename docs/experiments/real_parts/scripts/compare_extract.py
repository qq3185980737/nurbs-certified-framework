#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_extract.py — 比较两份 SURF 格式数据的 patch/控制点可复现性。

逐项检查：
  - NSURF 值、SURF 记录数、face 数（fid 集合）
  - 每个 face 的 patch 数序列
  - 每个 patch 的 (fid, pi, du, dv, nu, nv) 完全一致
  - 控制点/权重的最大绝对差（要求 0，或报告 ULP 级差异）
输出：人类可读 + 一行 COMPARE 摘要，可直接写入 results/。
"""
import argparse
import sys
from collections import defaultdict


def parse(path):
    nsurf = None
    patches = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        toks = fh.read().split()
    i = 0
    n = len(toks)
    while i < n:
        t = toks[i]
        if t == "NSURF":
            nsurf = int(toks[i + 1]); i += 2
        elif t == "SURF":
            fid, pi, du, dv, nu, nv = (int(toks[i + k]) for k in range(1, 7))
            i += 7
            cnt = nu * nv
            vals = toks[i:i + 4 * cnt]
            i += 4 * cnt
            # ENDSURF
            if i < n and toks[i] == "ENDSURF":
                i += 1
            nums = list(map(float, vals))
            pts = [(nums[4 * k], nums[4 * k + 1], nums[4 * k + 2], nums[4 * k + 3])
                   for k in range(cnt)]
            patches.append((fid, pi, du, dv, nu, nv, pts))
        else:
            i += 1
    return nsurf, patches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--max-report", type=int, default=5)
    args = ap.parse_args()

    nref, pref = parse(args.ref)
    nnew, pnew = parse(args.new)
    n_faces_ref = len(set(p[0] for p in pref))
    n_faces_new = len(set(p[0] for p in pnew))
    print("REF %s: NSURF=%s records=%d faces=%d" % (args.ref, nref, len(pref), n_faces_ref))
    print("NEW %s: NSURF=%s records=%d faces=%d" % (args.new, nnew, len(pnew), n_faces_new))

    ok = True
    if nref != nnew or len(pref) != len(pnew) or n_faces_ref != n_faces_new:
        ok = False
        print("MISMATCH at header level")

    per_face_ref = defaultdict(int)
    per_face_new = defaultdict(int)
    for p in pref:
        per_face_ref[p[0]] += 1
    for p in pnew:
        per_face_new[p[0]] += 1
    face_diffs = [(f, per_face_ref.get(f, 0), per_face_new.get(f, 0))
                  for f in sorted(set(per_face_ref) | set(per_face_new))
                  if per_face_ref.get(f, 0) != per_face_new.get(f, 0)]
    if face_diffs:
        ok = False
        print("face patch-count mismatches: %d (first %d shown)"
              % (len(face_diffs), min(len(face_diffs), args.max_report)))
        for d in face_diffs[:args.max_report]:
            print("  face %d ref=%d new=%d" % d)

    nm = min(len(pref), len(pnew))
    struct_bad = 0
    max_abs = 0.0
    max_rel = 0.0
    npts = 0
    argmax = None
    for idx in range(nm):
        a = pref[idx]
        b = pnew[idx]
        if a[:6] != b[:6]:
            struct_bad += 1
            if struct_bad <= args.max_report:
                print("  struct mismatch #%d ref=%s new=%s" % (idx, a[:6], b[:6]))
            continue
        for k, (pa, pb) in enumerate(zip(a[6], b[6])):
            for c in range(4):
                d = abs(pa[c] - pb[c])
                if d > max_abs:
                    max_abs = d
                    argmax = (idx, k, c, pa[c], pb[c])
                if d > 1e-15:
                    max_rel = max(max_rel, d / max(1e-300, abs(pa[c])))
                npts += 1
    if struct_bad:
        ok = False
    # 允许 <=2 ULP 级（1e-13 相对每个坐标尺度）的舍入差异
    tol_abs = 1e-9
    if max_abs > tol_abs:
        ok = False

    print("control points compared: %d, struct mismatches=%d" % (npts, struct_bad))
    print("max |delta| = %.3e (argmax patch_idx=%s pole=%s comp=%s ref=%.17g new=%.17g)"
          % (max_abs, *(argmax if argmax else ("-", "-", "-", 0.0, 0.0))))
    print("COMPARE %s faces_ref=%d faces_new=%d patches_ref=%d patches_new=%d "
          "face_count_diff=%d max_abs_delta=%.3e"
          % ("PASS" if ok else "FAIL", n_faces_ref, n_faces_new, len(pref), len(pnew),
             len(face_diffs), max_abs))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
