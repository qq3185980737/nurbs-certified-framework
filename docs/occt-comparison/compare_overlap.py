# -*- coding: utf-8 -*-
"""共面重叠对比：两共面图曲面 S1 z=0 over [0,1]^2, S2 z=0 over [0.25,0.75]^2。
真实重叠 = [0.25,0.75]^2。看 OCCT GeomAPI_IntSS 与 BRepAlgoAPI_Section 是否识别。"""
import json
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntSS
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX


def bezier_graph_rect(u0, u1, v0, v1, fz):
    arr = TColgp_Array2OfPnt(1, 3, 1, 3)
    for i in range(3):
        u = u0 + i * (u1 - u0) / 2.0
        for j in range(3):
            v = v0 + j * (v1 - v0) / 2.0
            arr.SetValue(i + 1, j + 1, gp_Pnt(u, v, fz(u, v)))
    return Geom_BezierSurface(arr)


S1 = bezier_graph_rect(0.0, 1.0, 0.0, 1.0, lambda u, v: 0.0)
S2 = bezier_graph_rect(0.25, 0.75, 0.25, 0.75, lambda u, v: 0.0)

res = {}
# GeomAPI_IntSS on coplanar overlap
for tol in [1e-3, 1e-6, 1e-9]:
    ss = GeomAPI_IntSS(S1, S2, tol)
    res["IntSS_tol" + str(tol)] = {
        "done": bool(ss.IsDone()),
        "nbLines": ss.NbLines(),
    }

# BRepAlgoAPI_Section(face1, face2) coplanar
f1 = BRepBuilderAPI_MakeFace(S1, 0.0, 1.0, 0.0, 1.0, 1e-7).Face()
f2 = BRepBuilderAPI_MakeFace(S2, 0.0, 1.0, 0.0, 1.0, 1e-7).Face()
sec = BRepAlgoAPI_Section(f1, f2)
sec.Build()
counts = {"faces": 0, "edges": 0, "vertices": 0}
for k, t in [(TopAbs_FACE, "faces"), (TopAbs_EDGE, "edges"), (TopAbs_VERTEX, "vertices")]:
    ex = TopExp_Explorer(sec.Shape(), k)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    counts[t] = n
res["Section"] = {"done": bool(sec.IsDone()), "counts": counts}
print(json.dumps(res, ensure_ascii=False, indent=1))
