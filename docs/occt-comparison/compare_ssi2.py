# -*- coding: utf-8 -*-
"""SSI 退化切触面：S2 z=(s-t)^2（在 s=t 处 z=0 且参数法向退化），S1 z=0。
真实交线为 s=t（退化切触）。看 OCCT GeomAPI_IntSS 是否返回。"""
import json
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntSS


def bezier_graph(fz):
    arr = TColgp_Array2OfPnt(1, 3, 1, 3)
    us = [0.0, 0.5, 1.0]
    vs = [0.0, 0.5, 1.0]
    for i, ui in enumerate(us, 1):
        for j, vj in enumerate(vs, 1):
            arr.SetValue(i, j, gp_Pnt(ui, vj, fz(ui, vj)))
    return Geom_BezierSurface(arr)


def fD(u, v):
    return (u - v) ** 2  # z=(s-t)^2, degenerate tangency along s=t


S1 = bezier_graph(lambda u, v: 0.0)
S2 = bezier_graph(fD)
for tol in [1e-3, 1e-6, 1e-9]:
    ss = GeomAPI_IntSS(S1, S2, tol)
    types = []
    for k in range(1, ss.NbLines() + 1):
        gc = ss.Line(k)
        types.append(gc.DynamicType().Name())
    print(json.dumps({"tol": tol, "done": bool(ss.IsDone()),
                      "nbLines": ss.NbLines(), "curveTypes": types}, ensure_ascii=False))
print("DONE")
