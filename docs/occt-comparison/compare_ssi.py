# -*- coding: utf-8 -*-
"""OCCT (pythonocc) 曲面-曲面求交对比：复现论文 §6.3 三个 graph-patch SSI fixture。
S1=(s,t,0) 平面，S2=(s,t,f)，交线为 f(s,t)=0。
fixture A: f=s^2-t^2 (两条对角线 s=±t)
fixture B: f=0.3 (恒正，全域不交)
fixture C: f=s+0.5t (在 [0,1]^2 仅角点触碰)
输出 JSON 到 stdout，供与框架覆盖带对比。"""
import json
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntSS


def bezier_graph(fz):
    """S(u,v)=(u,v,fz(u,v)) 用 degree2x2 Bezier 精确表示（fz 次数<=2）。"""
    arr = TColgp_Array2OfPnt(1, 3, 1, 3)
    us = [0.0, 0.5, 1.0]
    vs = [0.0, 0.5, 1.0]
    for i, ui in enumerate(us, 1):
        for j, vj in enumerate(vs, 1):
            arr.SetValue(i, j, gp_Pnt(ui, vj, fz(ui, vj)))
    return Geom_BezierSurface(arr)


def fA(u, v):
    return u * u - v * v


def fB(u, v):
    return 0.3


def fC(u, v):
    return u + 0.5 * v


res = {}
for name, f in [("A_s2mt2", fA), ("B_const", fB), ("C_s0p5t", fC)]:
    S1 = bezier_graph(lambda u, v: 0.0)
    S2 = bezier_graph(f)
    for tol in [1e-3, 1e-6, 1e-9]:
        ss = GeomAPI_IntSS(S1, S2, tol)
        ncur = ss.NbLines()
        types = []
        for k in range(1, ncur + 1):
            gc = ss.Line(k)
            types.append(gc.DynamicType().Name() if gc else "?")
        res.setdefault(name, {})[tol] = {
            "done": bool(ss.IsDone()),
            "nbLines": ncur,
            "curveTypes": types,
        }
print(json.dumps(res, ensure_ascii=False, indent=1))
