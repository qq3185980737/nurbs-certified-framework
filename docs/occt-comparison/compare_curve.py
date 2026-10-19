# -*- coding: utf-8 -*-
"""曲线-曲线求交充分对比：OCCT GeomAPI_ExtremaCurveCurve vs 本框架。
用例：横截单交 / 双横截交点 / 端点相切 / 内部相切 / 高阶切触。
OCCT Extrema 求最近点极值，Distance≈0 即交点。输出每个极值的距离与参数。"""
import json
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierCurve
from OCC.Core.GeomAPI import GeomAPI_ExtremaCurveCurve


def bezier(points):
    arr = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, 1):
        arr.SetValue(i, gp_Pnt(p[0], p[1], 0.0))
    return Geom_BezierCurve(arr)


def line(p0, p1):
    return bezier([p0, p1])


CASES = {
    # C1, C2, expected description
    "transverse": (line([0.0, 0.3], [1.0, 0.3]), line([0.0, 0.0], [1.0, 1.0]),
                   "single transverse x=0.3,y=0.3"),
    "two_cross": (bezier([[0, 0], [0.5, 0], [1, 1]]),  # y=x^2
                  line([0, 0], [1, 0.5]),             # y=0.5x
                  "two transverse roots x=0, x=0.5"),
    "tangent_endpoint": (bezier([[0, 0], [0.5, 0], [1, 1]]),  # y=x^2
                         line([0, 0], [1, 0]),                # y=0
                         "tangent at endpoint (0,0)"),
    "tangent_interior": (bezier([[0, 1], [0.5, -1], [1, 1]]),  # y=4(x-.5)^2 vertex (0.5,0)
                         line([0, 0], [1, 0]),                # y=0
                         "tangent at interior (0.5,0)"),
    "osculate": (bezier([[0, 0], [0.25, 0], [0.5, 1 / 6], [0.75, 0.5], [1, 2]]),  # y=x^2+x^4
                 bezier([[0, 0], [0.5, 0], [1, 1]]),                              # y=x^2
                 "second-order osculating at (0,0)"),
}

res = {}
for name, (c1, c2, desc) in CASES.items():
    e = GeomAPI_ExtremaCurveCurve(c1, c2)
    n = e.NbExtrema()
    dists = []
    for i in range(1, n + 1):
        try:
            d = e.Distance(i)
            (u1, u2) = e.Parameters(i)
        except Exception:
            d = -1
            u1 = u2 = None
        dists.append({"d": d, "u1": u1, "u2": u2})
    res[name] = {"desc": desc, "nbExtrema": n, "extrema": dists,
                 "numIntersectionsDetected": sum(1 for x in dists if x["d"] is not None and x["d"] < 1e-6)}
    print(json.dumps(res[name], ensure_ascii=False))
print("DONE")
