# -*- coding: utf-8 -*-
"""扩展类别 OCCT 对比（论文 Table 11）：3D 曲线-曲线、曲线-曲面、双三次 SSI、
弯曲 2D 重叠、有理圆柱面（横截 / 曲线贴合）。

环境：conda env `occ`（pythonocc-core，OCCT 8.0.1），tolerance 1e-9。
运行：/root/miniconda3/envs/occ/bin/python compare_ext.py
输出：JSON 到 stdout，供与框架结果对比。
"""
import json
import math
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColStd import TColStd_Array2OfReal
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierCurve, Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_ExtremaCurveCurve, GeomAPI_IntCS, GeomAPI_IntSS


def bez(points):
    a = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, 1):
        a.SetValue(i, gp_Pnt(*p))
    return Geom_BezierCurve(a)


def surf(rows):
    nu, nv = len(rows), len(rows[0])
    a = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            a.SetValue(i + 1, j + 1, gp_Pnt(*rows[i][j]))
    return Geom_BezierSurface(a)


out = {}

# ---------------- A) 3D 曲线-曲线 ----------------
c1 = bez([(0.0, 0.2, 0.2), (1.0, 0.2, 0.2)])       # 沿 x 的直线
c2 = bez([(0.2, 0.0, 0.2), (0.2, 1.0, 0.2)])       # 沿 y 的直线，交于 (0.2,0.2,0.2)
e = GeomAPI_ExtremaCurveCurve(c1, c2)
out["A_3d_lines_transversal"] = {
    "NbExtrema": e.NbExtrema(),
    "zero": sum(1 for i in range(1, e.NbExtrema() + 1) if e.Distance(i) < 1e-9),
    "expect": 1,
}

c3 = bez([(0.0, 0.0, 0.0), (0.5, 0.8, 0.5), (1.0, 0.0, 0.0)])   # 二次曲线，t=0.5 → (0.5,0.4,0.25)
c4 = bez([(0.5, 0.0, 0.25), (0.5, 1.0, 0.25)])                  # 过该点的直线
e2 = GeomAPI_ExtremaCurveCurve(c3, c4)
out["A2_3d_quad_vs_line"] = {
    "NbExtrema": e2.NbExtrema(),
    "zero": sum(1 for i in range(1, e2.NbExtrema() + 1) if e2.Distance(i) < 1e-9),
    "expect": 1,
}

# ---------------- B) 3D 曲线-曲面 ----------------
patch = surf([[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
              [(1.0, 0.0, 0.0), (1.0, 1.0, 0.0)]])              # 平面片 z=0
line = bez([(0.5, 0.5, -1.0), (0.5, 0.5, 1.0)])                 # 竖直横截
r = GeomAPI_IntCS(line, patch)
out["B_curve_surface_transversal"] = {"NbPoints": r.NbPoints(), "NbSegments": r.NbSegments(), "expect": 1}

tang = bez([(0.0, 0.0, 1.0), (0.5, 0.5, -1.0), (1.0, 0.0, 1.0)])  # z=(1-2t)^2，t=0.5 处触及 z=0
r2 = GeomAPI_IntCS(tang, patch)
out["B2_curve_surface_tangent"] = {"NbPoints": r2.NbPoints(), "NbSegments": r2.NbSegments(), "expect": 1}

# ---------------- C) 双三次 × 双三次 ----------------
S1 = surf([[(i / 3.0, j / 3.0, 0.0) for j in range(4)] for i in range(4)])
S2 = surf([
    [(0.0, 0.0, -0.4), (0.0, 0.333, -0.4), (0.0, 0.667, -0.4), (0.0, 1.0, -0.4)],
    [(0.333, 0.0, -0.13), (0.333, 0.333, -0.13), (0.333, 0.667, -0.13), (0.333, 1.0, -0.13)],
    [(0.667, 0.0, 0.13), (0.667, 0.333, 0.13), (0.667, 0.667, 0.13), (0.667, 1.0, 0.13)],
    [(1.0, 0.0, 0.4), (1.0, 0.333, 0.4), (1.0, 0.667, 0.4), (1.0, 1.0, 0.4)],
])
iss = GeomAPI_IntSS(S1, S2, 1e-9)
out["C_bicubic_ssi"] = {"NbLines": iss.NbLines(), "expect": 1}

# ---------------- D) 弯曲 2D 重叠 ----------------
def bowl(x, y):
    return x * x + y * y

S1b = surf([[(i / 3.0, j / 3.0, bowl(i / 3.0, j / 3.0)) for j in range(4)] for i in range(4)])
S2b = surf([[(0.25 + 0.5 * i / 3.0, 0.25 + 0.5 * j / 3.0,
              bowl(0.25 + 0.5 * i / 3.0, 0.25 + 0.5 * j / 3.0)) for j in range(4)] for i in range(4)])
iss2 = GeomAPI_IntSS(S1b, S2b, 1e-9)
out["D_curved_overlap_2d"] = {"NbLines": iss2.NbLines(), "expect": "2D overlap (region, not a line)"}

# ---------------- E) 有理四分之一圆柱面 ----------------
sq = math.sqrt(2) / 2
poles = TColgp_Array2OfPnt(1, 3, 1, 2)
w = TColStd_Array2OfReal(1, 3, 1, 2)
pts = [[(1.0, 0.0, 0.0), (1.0, 0.0, 1.0)],
       [(1.0, 1.0, 0.0), (1.0, 1.0, 1.0)],
       [(0.0, 1.0, 0.0), (0.0, 1.0, 1.0)]]
ws = [[1.0, 1.0], [sq, sq], [1.0, 1.0]]
for i in range(3):
    for j in range(2):
        poles.SetValue(i + 1, j + 1, gp_Pnt(*pts[i][j]))
        w.SetValue(i + 1, j + 1, ws[i][j])
cyl = Geom_BezierSurface(poles, w)

# E1) 水平横截线（y=0.5, z=0.5）→ 弧上 1 点
r3 = GeomAPI_IntCS(bez([(-0.5, 0.5, 0.5), (1.2, 0.5, 0.5)]), cyl)
uniq = set()
for i in range(1, r3.NbPoints() + 1):
    q = r3.Point(i)
    uniq.add((round(q.X(), 6), round(q.Y(), 6), round(q.Z(), 6)))
out["E1_rational_surface_transversal"] = {"NbPoints": r3.NbPoints(), "unique": len(uniq), "expect": 1}

# E2) 曲线贴合曲面：竖直线沿 45° 母线 → 交集为 1D 段
r4 = GeomAPI_IntCS(bez([(sq, sq, -0.5), (sq, sq, 1.5)]), cyl)
out["E2_curve_on_rational_surface"] = {"NbPoints": r4.NbPoints(), "NbSegments": r4.NbSegments(),
                                       "expect": "1D segment"}

print(json.dumps(out, indent=1, ensure_ascii=False))
