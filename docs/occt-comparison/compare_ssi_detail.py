# -*- coding: utf-8 -*-
"""采样 OCCT 对 fixture A (s^2-t^2=0) 返回的那条交线，确认是 s=t 还是 s=-t，验证是否漏一条。"""
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


def fA(u, v):
    return u * u - v * v


S1 = bezier_graph(lambda u, v: 0.0)
S2 = bezier_graph(fA)
ss = GeomAPI_IntSS(S1, S2, 1e-9)
print("IsDone:", ss.IsDone(), "nbLines:", ss.NbLines())
for k in range(1, ss.NbLines() + 1):
    gc = ss.Line(k)
    print(f"--- Line {k}: {gc.DynamicType().Name()} ---")
    print("  FirstParam:", gc.FirstParameter(), "LastParam:", gc.LastParameter())
    for t in [gc.FirstParameter(),
              (gc.FirstParameter() + gc.LastParameter()) / 2,
              gc.LastParameter()]:
        p = gc.Value(t)
        x, y, z = p.X(), p.Y(), p.Z()
        tag = "s=+t" if abs(x - y) < 1e-6 else ("s=-t" if abs(x + y) < 1e-6 else "?")
        print(f"  t={t:.3f} -> ({x:.4f},{y:.4f},{z:.4f})  ({tag})")
