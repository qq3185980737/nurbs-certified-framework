# -*- coding: utf-8 -*-
"""曲线-曲面求交对比：OCCT GeomAPI_IntCS vs 本框架。
用例：平面 z=0 横截穿入 / 曲线切平面(相切) / 曲线三横截穿越图曲面。"""
import json
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierCurve, Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntCS


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


plane_z0 = bezier_graph(lambda u, v: 0.0)
paraboloid = bezier_graph(lambda u, v: u * u + v * v)  # z=x^2+y^2 over [0,1]^2

# --- cases ---
# A: transverse through plane z=0 at (0.5,0.5,0)
C_A = bezier_curve([[0.5, 0.5, -0.5], [0.5, 0.5, 0.5]])
# B: curve tangent to plane z=0 at x=0.5 -> (0.5,0.5,0), dz/dx=0
C_B = bezier_curve([[0.0, 0.5, 0.25], [0.5, 0.5, -0.25], [1.0, 0.5, 0.25]])  # z=(x-0.5)^2 tangent at x=0.5
# C: three transverse roots: curve (t,0.5,z(t)) crosses paraboloid slice z=x^2+0.25 at x=a,b,c
a, b, c = 0.2, 0.5, 0.8
k = 0.6
# zC(t) = t^2+0.25 + k*(t-a)(t-b)(t-c) = c0+c1*t+c2*t^2+c3*t^3
c0 = 0.25 - k * a * b * c
c1 = k * (a * b + a * c + b * c)
c2 = 1.0 - k * (a + b + c)
c3 = k
# Bernstein control points for degree-3 polynomial f=c0+c1 t+c2 t^2+c3 t^3
P0 = c0
P3 = c0 + c1 + c2 + c3
P1 = (c1 + 3 * P0) / 3.0
P2 = (c2 - 3 * P0 + 6 * P1) / 3.0
C_C = bezier_curve([[0.0, 0.5, P0], [1 / 3.0, 0.5, P1], [2 / 3.0, 0.5, P2], [1.0, 0.5, P3]])

res = {}
for name, C, S, expected in [
    ("transverse_plane", C_A, plane_z0, 1),
    ("tangent_plane", C_B, plane_z0, 1),
    ("three_transverse", C_C, paraboloid, 3),
]:
    ss = GeomAPI_IntCS(C, S)
    n = ss.NbPoints()
    pts = []
    for i in range(1, n + 1):
        p = ss.Point(i)
        pts.append([round(p.X(), 6), round(p.Y(), 6), round(p.Z(), 6)])
    res[name] = {"expected": expected, "nbPoints": n, "points": pts,
                 "allFound": n == expected}
    print(json.dumps(res[name], ensure_ascii=False))
print("DONE")
