# -*- coding: utf-8 -*-
"""NURBS 有理曲线对比：有理 Bézier 1/4 圆（weights 1,sqrt2/2,1）与整圆 vs 直线求交。
OCCT Geom_BezierCurve(poles,weights) + BRepAlgoAPI_Section。"""
import json
import math
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.TColStd import TColStd_Array1OfReal
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_Circ
from OCC.Core.Geom import Geom_BezierCurve, Geom_Circle, Geom_Line
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX


def rational_bezier(points, weights):
    pa = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, 1):
        pa.SetValue(i, gp_Pnt(p[0], p[1], 0.0))
    wa = TColStd_Array1OfReal(1, len(weights))
    for i, w in enumerate(weights, 1):
        wa.SetValue(i, w)
    return Geom_BezierCurve(pa, wa)


def section_vertices(c1, c2):
    e1 = BRepBuilderAPI_MakeEdge(c1).Edge()
    e2 = BRepBuilderAPI_MakeEdge(c2).Edge()
    sec = BRepAlgoAPI_Section(e1, e2)
    sec.Build()
    ex = TopExp_Explorer(sec.Shape(), TopAbs_VERTEX)
    n = 0
    while ex.More():
        n += 1
        ex.Next()
    return n


w2 = math.sqrt(2.0) / 2.0
# rational 1/4 circle (1,0)-(0,1), weights 1, w2, 1
quarter = rational_bezier([(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], [1.0, w2, 1.0])
# full circle as rational Bezier (4 arcs) is tedious; use Geom_Circle for the full-circle check
full = Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.0)
# lines
line_y03 = rational_bezier([(0.0, 0.3), (1.2, 0.3)], [1.0, 1.0])  # straight (linear, weights trivial)
line_x2 = rational_bezier([(0.0, -1.2), (0.0, 1.2)], [1.0, 1.0])  # x=0

res = {}
res["rational_quarter_circle_y03"] = {"expected": 1,
                                      "sectionVertices": section_vertices(quarter, line_y03)}
# full circle via Geom_Circle (exact rational conic) with y=0.3 -> 2 roots
line_e = BRepBuilderAPI_MakeEdge(line_y03).Edge()
cir_e = BRepBuilderAPI_MakeEdge(full).Edge()
sec = BRepAlgoAPI_Section(cir_e, line_e)
sec.Build()
ex = TopExp_Explorer(sec.Shape(), TopAbs_VERTEX)
n = 0
while ex.More():
    n += 1
    ex.Next()
res["full_circle_y03"] = {"expected": 2, "sectionVertices": n}
print(json.dumps(res, ensure_ascii=False, indent=1))
