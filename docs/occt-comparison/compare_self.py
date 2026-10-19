# -*- coding: utf-8 -*-
"""曲面自交检测对比：S(u,v)=(u^2-v^2, 2uv, u^2+v^2), (u,v) in [-1,1]^2.
(u,v) 与 (-u,-v) 折叠同点 -> 自交沿 u=0 折叠线 (x=-z, y=0)。
看 OCCT BRepCheck_Analyzer 与 GeomAPI_IntSS(S,S) 是否检测。"""
import json
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntSS
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepCheck import BRepCheck_Analyzer

arr = TColgp_Array2OfPnt(1, 3, 1, 3)
us = [-1.0, 0.0, 1.0]
vs = [-1.0, 0.0, 1.0]
for i, u in enumerate(us, 1):
    for j, v in enumerate(vs, 1):
        arr.SetValue(i, j, gp_Pnt(u * u - v * v, 2.0 * u * v, u * u + v * v))
S = Geom_BezierSurface(arr)

res = {}
# 1) self-intersection via IntSS(S,S)
ss = GeomAPI_IntSS(S, S, 1e-6)
res["IntSS_self"] = {"done": bool(ss.IsDone()), "nbLines": ss.NbLines()}

# 2) BRepCheck_Analyzer on the face
face = BRepBuilderAPI_MakeFace(S, -1.0, 1.0, -1.0, 1.0, 1e-7).Face()
ana = BRepCheck_Analyzer(face)
res["BRepCheck_Analyzer"] = {"valid": bool(ana.IsValid())}
print(json.dumps(res, ensure_ascii=False, indent=1))
