# -*- coding: utf-8 -*-
"""P0-2：支架（STEP 零件）+ 横杆（手工构造平面片）装配干涉 → 认证边-面求交数据。"""
import re
from OCC.Core.Geom import Geom_BezierSurface, Geom_BezierCurve
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt
from OCC.Core.TColStd import TColStd_Array2OfReal, TColStd_Array1OfReal
from OCC.Core.gp import gp_Pnt

BASE = "/mnt/c/Users/15632/Documents/deepseek-harness/default-workspace/cad_review12/step_case"

# ---------- 读取支架数据（已导出的 case_data.txt）----------
lines = open(BASE + "/case_data.txt", "r", encoding="utf-8").read().split()
i = 0
surfs = {}   # fid -> list of (du,dv,nu,nv,rows)
curves = {}  # eid -> list of (deg,rows)
pairs = []
def num(t): return float(t)
while i < len(lines):
    t = lines[i]
    if t == "NSURF":
        i += 2
    elif t == "SURF":
        fid, pi, du, dv, nu, nv = int(lines[i+1]), int(lines[i+2]), int(lines[i+3]), int(lines[i+4]), int(lines[i+5]), int(lines[i+6])
        i += 7
        rows = []
        for _ in range(nu * nv):
            rows.append((num(lines[i]), num(lines[i+1]), num(lines[i+2]), num(lines[i+3]))); i += 4
        surfs.setdefault(fid, []).append((du, dv, nu, nv, rows))
        assert lines[i] == "ENDSURF"; i += 1
    elif t == "NCURVE":
        i += 2
    elif t == "CURVE":
        eid, si, deg, n = int(lines[i+1]), int(lines[i+2]), int(lines[i+3]), int(lines[i+4])
        i += 5
        rows = []
        for _ in range(n):
            rows.append((num(lines[i]), num(lines[i+1]), num(lines[i+2]), num(lines[i+3]))); i += 4
        curves.setdefault(eid, []).append((deg, rows))
        assert lines[i] == "ENDCURVE"; i += 1
    elif t == "NPAIR":
        i += 2
    elif t == "PAIR":
        pi_ = len(pairs); e, s, m = int(lines[i+1]), int(lines[i+2]), int(lines[i+3]); i += 4
        pts = []
        for _ in range(m):
            pts.append((num(lines[i]), num(lines[i+1]), num(lines[i+2]))); i += 3
        pairs.append({"e": e, "s": s, "pts": pts})
        assert lines[i] == "ENDPAIR"; i += 1
    else:
        i += 1
print("支架数据: surfs=%d curves=%d pairs=%d" % (len(surfs), len(curves), len(pairs)), flush=True)

# ---------- 横杆：长方体 [5,75]x[22,38]x[16,26] 的 6 个平面片 ----------
X0, X1, Y0, Y1, Z0, Z1 = 5.0, 75.0, 22.0, 38.0, 16.0, 26.0
# 双线性 Bezier 片：控制点按 i(u) 主序，顺序 [(u0,v0),(u0,v1),(u1,v0),(u1,v1)]
BAR_FACES = [
    [(X0, Y0, Z0), (X0, Y1, Z0), (X1, Y0, Z0), (X1, Y1, Z0)],   # z=Z0
    [(X0, Y0, Z1), (X0, Y1, Z1), (X1, Y0, Z1), (X1, Y1, Z1)],   # z=Z1
    [(X0, Y0, Z0), (X0, Y0, Z1), (X1, Y0, Z0), (X1, Y0, Z1)],   # y=Y0
    [(X0, Y1, Z0), (X0, Y1, Z1), (X1, Y1, Z0), (X1, Y1, Z1)],   # y=Y1
    [(X0, Y0, Z0), (X0, Y0, Z1), (X0, Y1, Z0), (X0, Y1, Z1)],   # x=X0
    [(X1, Y0, Z0), (X1, Y0, Z1), (X1, Y1, Z0), (X1, Y1, Z1)],   # x=X1
]
BAR_EDGES = [
    ((X0, Y0, Z0), (X1, Y0, Z0)), ((X1, Y0, Z0), (X1, Y1, Z0)),
    ((X1, Y1, Z0), (X0, Y1, Z0)), ((X0, Y1, Z0), (X0, Y0, Z0)),
    ((X0, Y0, Z1), (X1, Y0, Z1)), ((X1, Y0, Z1), (X1, Y1, Z1)),
    ((X1, Y1, Z1), (X0, Y1, Z1)), ((X0, Y1, Z1), (X0, Y0, Z1)),
    ((X0, Y0, Z0), (X0, Y0, Z1)), ((X1, Y0, Z0), (X1, Y0, Z1)),
    ((X1, Y1, Z0), (X1, Y1, Z1)), ((X0, Y1, Z0), (X0, Y1, Z1)),
]

SURF_BASE = max(surfs.keys()) + 1
CURVE_BASE = max(curves.keys()) + 1
for k, c in enumerate(BAR_FACES):
    surfs[SURF_BASE + k] = [(1, 1, 2, 2, [(p[0], p[1], p[2], 1.0) for p in c])]
for k, (a, b) in enumerate(BAR_EDGES):
    curves[CURVE_BASE + k] = [(1, [(a[0], a[1], a[2], 1.0), (b[0], b[1], b[2], 1.0)])]
print("横杆加入: 6 面 (ids %d..%d), 12 边 (ids %d..%d)" %
      (SURF_BASE, SURF_BASE + 5, CURVE_BASE, CURVE_BASE + 11), flush=True)


def bez_surf(du, dv, nu, nv, rows):
    grid = TColgp_Array2OfPnt(1, nu, 1, nv)
    wts = TColStd_Array2OfReal(1, nu, 1, nv)
    rational = False
    for iu in range(nu):
        for jv in range(nv):
            x, y, z, w = rows[iu * nv + jv]
            grid.SetValue(iu + 1, jv + 1, gp_Pnt(x, y, z))
            wts.SetValue(iu + 1, jv + 1, w)
            if abs(w - 1.0) > 1e-12:
                rational = True
    if rational:
        return Geom_BezierSurface(grid, wts)
    return Geom_BezierSurface(grid)


def bez_curve(deg, rows):
    arr = TColgp_Array1OfPnt(1, len(rows))
    wts = TColStd_Array1OfReal(1, len(rows))
    rational = False
    for i2, (x, y, z, w) in enumerate(rows):
        arr.SetValue(i2 + 1, gp_Pnt(x, y, z))
        wts.SetValue(i2 + 1, w)
        if abs(w - 1.0) > 1e-12:
            rational = True
    if rational:
        return Geom_BezierCurve(arr, wts)
    return Geom_BezierCurve(arr)


# ---------- OCCT 参考交点：支架边 × 横杆面，横杆边 × 支架面 ----------
new_pairs = []
def add_pair(eid, sid):
    pts = []
    for (deg, rows) in curves[eid]:
        c = bez_curve(deg, rows)
        for (du, dv, nu, nv, srows) in surfs[sid]:
            s = bez_surf(du, dv, nu, nv, srows)
            try:
                res = GeomAPI_IntCS(c, s)
            except Exception:
                continue
            if not res.IsDone():
                continue
            for k in range(1, res.NbPoints() + 1):
                p = res.Point(k)
                pts.append((p.X(), p.Y(), p.Z()))
    if pts:
        new_pairs.append({"e": eid, "s": sid, "pts": pts})


for eid in list(curves.keys()):
    if eid >= CURVE_BASE:
        continue                     # 支架边 × 横杆面
    for sid in range(SURF_BASE, SURF_BASE + 6):
        add_pair(eid, sid)
for sid in list(surfs.keys()):
    if sid >= SURF_BASE:
        continue                     # 横杆边 × 支架面
    for eid in range(CURVE_BASE, CURVE_BASE + 12):
        add_pair(eid, sid)
print("装配跨零件边-面对:", len(new_pairs), "| 交点:", sum(len(p["pts"]) for p in new_pairs), flush=True)

# ---------- 写出增强数据 ----------
with open(BASE + "/assembly_case.txt", "w", encoding="utf-8") as fh:
    fh.write("NSURF %d\n" % (max(surfs.keys()) + 1))
    for fid in sorted(surfs.keys()):
        for pi_, (du, dv, nu, nv, rows) in enumerate(surfs[fid]):
            fh.write("SURF %d %d %d %d %d %d\n" % (fid, pi_, du, dv, nu, nv))
            for (x, y, z, w) in rows:
                fh.write("%.17g %.17g %.17g %.17g\n" % (x, y, z, w))
            fh.write("ENDSURF\n")
    fh.write("NCURVE %d\n" % (max(curves.keys()) + 1))
    for eid in sorted(curves.keys()):
        for si_, (deg, rows) in enumerate(curves[eid]):
            fh.write("CURVE %d %d %d %d\n" % (eid, si_, deg, len(rows)))
            for (x, y, z, w) in rows:
                fh.write("%.17g %.17g %.17g %.17g\n" % (x, y, z, w))
            fh.write("ENDCURVE\n")
    fh.write("NPAIR %d\n" % len(new_pairs))
    for p in new_pairs:
        fh.write("PAIR %d %d %d\n" % (p["e"], p["s"], len(p["pts"])))
        for (x, y, z) in p["pts"]:
            fh.write("%.17g %.17g %.17g\n" % (x, y, z))
        fh.write("ENDPAIR\n")
print("已写出", BASE + "/assembly_case.txt", flush=True)
