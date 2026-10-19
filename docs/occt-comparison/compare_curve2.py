# -*- coding: utf-8 -*-
"""曲线求交：多横截根（degree5 Bezier × y=0，真实5交点），看 OCCT ExtremaCurveCurve 是否全找。"""
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


# degree5: y = k * prod(x - r_i), 5 transverse roots, none at Bezier knots
roots = [0.08, 0.25, 0.45, 0.7, 0.92]
k = 400.0


def P(x):
    p = 1.0
    for r in roots:
        p *= (x - r)
    return k * p


n = 6  # degree5 -> 6 ctrl pts
cpts = [[i / (n - 1), P(i / (n - 1))] for i in range(n)]
c1 = bezier(cpts)
c2 = bezier([[0, 0], [1, 0]])  # y=0

e = GeomAPI_ExtremaCurveCurve(c1, c2)
dists = []
for i in range(1, e.NbExtrema() + 1):
    try:
        d = e.Distance(i)
        (u1, u2) = e.Parameters(i)
    except Exception:
        d, u1, u2 = -1, None, None
    dists.append({"d": d, "u1": u1})
found = sum(1 for x in dists if x["d"] is not None and x["d"] < 1e-6)
print(json.dumps({
    "roots": roots,
    "nbExtrema": e.NbExtrema(),
    "foundIntersections": found,
    "expected": 5,
    "allFound": found == 5,
    "dists": dists,
}, ensure_ascii=False, indent=1))
