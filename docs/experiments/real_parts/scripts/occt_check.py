#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""occt_check.py — 用 OCCT GeomAPI_IntSS 对抽样的 patch-pair 做交叉校验。

输入：data/<part>.txt（与 C++ 认证器完全相同的 patch 数据）+ probe_pairs 输出的 pairs CSV
输出：
  occt_pairs.csv   每对一行：part,fidA,piA,fidB,piB,src,bbox_overlap,type_bench,
                   cop_overlap,coincidence,cover,certified,undecided,ms,occt_lines,occt_error
  occt_points.csv  交线采样点（供 cover 覆盖检验）：part,fidA,piA,fidB,piB,line,k,x,y,z

仅当 IntSS 成功时才记录 occt_lines；异常/未完成记为 occt_error 与 occt_lines=-1。
用法（WSL，OCCT 环境）：
  occ_py.sh occt_check.py --data data/228.txt --pairs results/probe/228_pairs.csv \
      --out results/probe/228_occt.csv --points results/probe/228_points.csv
"""
import argparse
import csv
import sys

from OCC.Core.gp import gp_Pnt
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.Geom import Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntSS
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve


def parse_data(path):
    """返回 dict[(fid,pi)] = (nu, nv, [(x,y,z,w),...])"""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        toks = fh.read().split()
    out = {}
    i, n = 0, len(toks)
    while i < n:
        t = toks[i]
        if t == "NSURF":
            i += 2
        elif t == "SURF":
            fid, pi, du, dv, nu, nv = (int(toks[i + k]) for k in range(1, 7))
            i += 7
            cnt = nu * nv
            vals = list(map(float, toks[i:i + 4 * cnt]))
            i += 4 * cnt
            if i < n and toks[i] == "ENDSURF":
                i += 1
            pts = [(vals[4 * k], vals[4 * k + 1], vals[4 * k + 2], vals[4 * k + 3])
                   for k in range(cnt)]
            out[(fid, pi)] = (nu, nv, pts)
        else:
            i += 1
    return out


def make_surface(nu, nv, pts):
    arr = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            x, y, z, w = pts[i * nv + j]
            arr.SetValue(i + 1, j + 1, gp_Pnt(x, y, z))
    s = Geom_BezierSurface(arr)
    rational = False
    for i in range(nu):
        for j in range(nv):
            w = pts[i * nv + j][3]
            if abs(w - 1.0) > 1e-15:
                s.SetWeight(i + 1, j + 1, w)
                rational = True
    return s, rational


def sample_line_points(res, li, max_pts=40, max_lines=5):
    """从 IntSS 结果第 li 条交线采样三维点。"""
    pts = []
    try:
        line = res.Line(li)
        mk = BRepBuilderAPI_MakeEdge(line)
        if not mk.IsDone():
            return pts
        ad = BRepAdaptor_Curve(mk.Edge())
        t0, t1 = ad.FirstParameter(), ad.LastParameter()
        for k in range(max_pts):
            t = t0 + (t1 - t0) * k / max(1, max_pts - 1)
            p = ad.Value(t)
            pts.append((p.X(), p.Y(), p.Z()))
    except Exception:
        pass
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--points", default=None)
    ap.add_argument("--part", default=None)
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()
    part = args.part
    if part is None:
        import os
        part = os.path.basename(args.data).replace(".txt", "")

    data = parse_data(args.data)
    rows = list(csv.DictReader(open(args.pairs, "r", encoding="utf-8")))
    fh = open(args.out, "w", newline="", encoding="utf-8")
    w = csv.writer(fh)
    w.writerow(["part", "fidA", "piA", "fidB", "piB", "src", "bbox_overlap",
                "type_bench", "cop_overlap", "coincidence", "cover", "certified",
                "undecided", "ms", "occt_lines", "occt_error"])
    pfh = None
    pw = None
    if args.points:
        pfh = open(args.points, "w", newline="", encoding="utf-8")
        pw = csv.writer(pfh)
        pw.writerow(["part", "fidA", "piA", "fidB", "piB", "line", "k", "x", "y", "z"])

    n_ok = n_err = 0
    for r in rows:
        ka = (int(r["fidA"]), int(r["piA"]))
        kb = (int(r["fidB"]), int(r["piB"]))
        base = [part, ka[0], ka[1], kb[0], kb[1], r["src"], r["bbox_overlap"],
                r["type"], r.get("cop_overlap", ""), r.get("coincidence", ""),
                r.get("cover", ""), r.get("certified", ""), r.get("undecided", ""),
                r.get("ms", "")]
        lines = -1
        err = ""
        try:
            nuA, nvA, pA = data[ka]
            nuB, nvB, pB = data[kb]
            sA, _ = make_surface(nuA, nvA, pA)
            sB, _ = make_surface(nuB, nvB, pB)
            res = GeomAPI_IntSS(sA, sB, args.tol)
            if not res.IsDone():
                lines = -1
                err = "NotDone"
            else:
                lines = res.NbLines()
                if pw is not None and lines > 0:
                    for li in range(1, min(lines, 5) + 1):
                        pts = sample_line_points(res, li)
                        for k, (x, y, z) in enumerate(pts):
                            pw.writerow([part, ka[0], ka[1], kb[0], kb[1], li, k, x, y, z])
                n_ok += 1
        except Exception as e:
            lines = -1
            err = type(e).__name__ + ":" + str(e)[:120]
            n_err += 1
        w.writerow(base + [lines, err])
    fh.close()
    if pfh:
        pfh.close()
    print("OCCT part=%s pairs=%d ok=%d err=%d out=%s" % (part, len(rows), n_ok, n_err, args.out), flush=True)


if __name__ == "__main__":
    sys.exit(main())
