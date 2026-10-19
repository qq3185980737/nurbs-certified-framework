# -*- coding: utf-8 -*-
"""用 OCCT 正规曲线求交 API (BRepAlgoAPI_Section) 复核：多根 / 单根 / 内部相切。"""
import json
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierCurve
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX


def bezier(points):
    arr = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, 1):
        arr.SetValue(i, gp_Pnt(p[0], p[1], 0.0))
    return Geom_BezierCurve(arr)


def section_vertices(c1, c2):
    e1 = BRepBuilderAPI_MakeEdge(c1).Edge()
    e2 = BRepBuilderAPI_MakeEdge(c2).Edge()
    sec = BRepAlgoAPI_Section(e1, e2)
    sec.Build()
    expl = TopExp_Explorer(sec.Shape(), TopAbs_VERTEX)
    n = 0
    while expl.More():
        n += 1
        expl.Next()
    return n


roots = [0.08, 0.25, 0.45, 0.7, 0.92]
k = 400.0


def P(x):
    p = 1.0
    for r in roots:
        p *= (x - r)
    return k * p


multi = bezier([[i / 5.0, P(i / 5.0)] for i in range(6)])  # degree5
line_y0 = bezier([[0, 0], [1, 0]])

res = {}
res["multi_5roots"] = {"expected": 5, "sectionVertices": section_vertices(multi, line_y0)}
res["single_transverse"] = {"expected": 1,
                            "sectionVertices": section_vertices(
                                bezier([[0, 0.3], [1, 0.3]]), bezier([[0, 0], [1, 1]]))}
res["tangent_interior"] = {"expected": 1,
                           "sectionVertices": section_vertices(
                               bezier([[0, 1], [0.5, -1], [1, 1]]), line_y0)}
print(json.dumps(res, ensure_ascii=False, indent=1))
