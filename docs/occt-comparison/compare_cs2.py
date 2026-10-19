# -*- coding: utf-8 -*-
"""曲线-曲面三交点：先用采样确认几何真解，再用 BRepAlgoAPI_Section(edge,face) 复核。"""
import json
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierCurve, Geom_BezierSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX


def bezier_curve(points):
    arr = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, 1):
        arr.SetValue(i, gp_Pnt(p[0], p[1], p[2]))
    return Geom_BezierCurve(arr)


def bezier_graph(fz, du=2, dv=2):
    n, m = du + 1, dv + 1
    arr = TColgp_Array2OfPnt(1, n, 1, m)
    for i in range(n):
        u = i / du
        for j in range(m):
            v = j / dv
            arr.SetValue(i + 1, j + 1, gp_Pnt(u, v, fz(u, v)))
    return Geom_BezierSurface(arr)


paraboloid = bezier_graph(lambda u, v: u * u + v * v)
a, b, c = 0.2, 0.5, 0.8
k = 0.6
c0 = 0.25 - k * a * b * c
c1 = k * (a * b + a * c + b * c)
c2 = 1.0 - k * (a + b + c)
c3 = k
P0 = c0
P3 = c0 + c1 + c2 + c3
P1 = (c1 + 3 * P0) / 3.0
P2 = (c2 - 3 * P0 + 6 * P1) / 3.0
C_C = bezier_curve([[0.0, 0.5, P0], [1 / 3.0, 0.5, P1], [2 / 3.0, 0.5, P2], [1.0, 0.5, P3]])

# 1) geometric confirmation by sampling residual z - (x^2+y^2)
geom = []
for t in [0.2, 0.5, 0.8]:
    p = C_C.Value(t)
    geom.append({"t": t, "resid": p.Z() - (p.X() ** 2 + p.Y() ** 2),
                 "pt": [round(p.X(), 6), round(p.Y(), 6), round(p.Z(), 6)]})
print(json.dumps({"geometricResiduals": geom}, ensure_ascii=False))

# 2) BRepAlgoAPI_Section(edge, face)
e = BRepBuilderAPI_MakeEdge(C_C).Edge()
f = BRepBuilderAPI_MakeFace(paraboloid, 0.0, 1.0, 0.0, 1.0, 1e-7).Face()
sec = BRepAlgoAPI_Section(e, f)
sec.Build()
expl = TopExp_Explorer(sec.Shape(), TopAbs_VERTEX)
n = 0
while expl.More():
    n += 1
    expl.Next()
print(json.dumps({"sectionVertices": n, "expected": 3}, ensure_ascii=False))
