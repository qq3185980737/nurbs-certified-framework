#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""OCCT IntPatch-based baseline on the SAME fixtures A-F used for the SISL
baseline (see src/sisl_ssi_harness.cpp).  GeomAPI_IntSS wraps IntPatch for
surface/surface (and GeomAPI_IntCS for surface/curve).  Tolerances are the
same as the SISL runs: 1e-6 and 1e-4.

Residual definition (same spirit as the SISL harness): the returned curve is
densely sampled and every sample is projected onto BOTH input surfaces with
GeomAPI_ProjectPointOnSurf; the reported residual is max over samples of the
larger of the two projection distances.  For surface/curve fixtures the
distance is measured to the surface and to the curve.
"""
import time
import math
import sys

from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.TColStd import TColStd_Array2OfReal
from OCC.Core.gp import gp_Pnt
from OCC.Core.Geom import Geom_BezierSurface
from OCC.Core.GeomAPI import GeomAPI_IntSS, GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve

NS = 40  # samples per intersection curve


def binom(n, k):
    if k < 0 or k > n:
        return 0.0
    r = 1.0
    for i in range(k):
        r = r * (n - i) / (i + 1)
    return r


def mono_to_bezier(a, n):
    """a[i][j] = coefficient of u^i v^j (monomial) -> n x n Bernstein net."""
    z = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = 0.0
            for k in range(i + 1):
                for l in range(j + 1):
                    ak = a[k][l] if (k < len(a) and l < len(a[k])) else 0.0
                    s += ak * binom(i, k) / binom(n - 1, k) * binom(j, l) / binom(n - 1, l)
            z[i][j] = s
    return z


def graph_surface(a, n):
    z = mono_to_bezier(a, n)
    arr = TColgp_Array2OfPnt(1, n, 1, n)
    for i in range(n):
        for j in range(n):
            arr.SetValue(i + 1, j + 1, gp_Pnt(i / (n - 1), j / (n - 1), z[i][j]))
    return Geom_BezierSurface(arr)


def quarter_cylinder():
    """x=cos(t), y=sin(t), z=s, exact rational biquadratic x linear Bezier."""
    w = [1.0, math.sqrt(2.0) / 2.0, 1.0]
    px = [1.0, 1.0, 0.0]
    py = [0.0, 1.0, 1.0]
    poles = TColgp_Array2OfPnt(1, 3, 1, 2)
    weights = TColStd_Array2OfReal(1, 3, 1, 2)
    for i in range(3):
        for j in range(2):
            poles.SetValue(i + 1, j + 1, gp_Pnt(px[i], py[i], float(j)))
            weights.SetValue(i + 1, j + 1, w[i])
    return Geom_BezierSurface(poles, weights)


def bezier_line(p0, p1):
    from OCC.Core.Geom import Geom_BezierCurve
    from OCC.Core.TColgp import TColgp_Array1OfPnt
    arr = TColgp_Array1OfPnt(1, 2)
    arr.SetValue(1, gp_Pnt(*p0))
    arr.SetValue(2, gp_Pnt(*p1))
    return Geom_BezierCurve(arr)


def surf_residual(curve, surf):
    """max over samples of distance from curve point to surf."""
    try:
        u0 = curve.FirstParameter()
        u1 = curve.LastParameter()
    except Exception:
        return float('nan'), 0
    worst = 0.0
    n = 0
    for k in range(NS + 1):
        u = u0 + (u1 - u0) * k / NS
        try:
            p = curve.Value(u)
        except Exception:
            continue
        proj = GeomAPI_ProjectPointOnSurf(p, surf)
        if proj.NbPoints() > 0:
            d = proj.LowerDistance()
            worst = max(worst, d)
        n += 1
    return worst, n


def curve_point_residual(p, surf):
    proj = GeomAPI_ProjectPointOnSurf(p, surf)
    if proj.NbPoints() == 0:
        return float('nan')
    return proj.LowerDistance()


def curve_curve_residual(p, curve):
    proj = GeomAPI_ProjectPointOnCurve(p, curve)
    if proj.NbPoints() == 0:
        return float('nan')
    return proj.LowerDistance()


FIXTURES = {
    # name: (monomial matrix for S1, n1, S2 z-matrix, n2, truth)
    "A": ([[0, 0, -1], [0, 0, 0], [1, 0, 0]], 3, [[0.0]], 2,
          "1 curve: diagonal s=t"),
    "B": ([[0.3]], 2, [[0.0]], 2, "0 intersections"),
    "C": ([[0.0, 0.5], [1.0, 0.0]], 2, [[0.0]], 2,
          "1 point: corner (0,0) tangency"),
    "D": ([[0, 0, 1], [0, -2, 0], [1, 0, 0]], 3, [[0.0]], 2,
          "1 tangential/degenerate curve s=t"),
    "E": (None, 4, None, 4, "1 interior curve: u^3+v^3 = 0.1"),
}


def run_ss_fixture(name, tol):
    """surface/surface via GeomAPI_IntSS (IntPatch)."""
    if name == "E":
        a1 = [[0.0] * 4 for _ in range(4)]
        a1[3][0] = 1.0
        a2 = [[0.0] * 4 for _ in range(4)]
        a2[0][0] = 0.1
        a2[0][3] = -1.0
        s1 = graph_surface(a1, 4)
        s2 = graph_surface(a2, 4)
    else:
        a1, n1, a2, n2, _ = FIXTURES[name]
        s1 = graph_surface(a1, n1)
        s2 = graph_surface(a2, n2)

    t0 = time.perf_counter()
    it = GeomAPI_IntSS(s1, s2, tol)
    dt = (time.perf_counter() - t0) * 1e3
    done = bool(it.IsDone())
    nb_lines = it.NbLines() if done else 0
    # OCCT 8.0.1 GeomAPI_IntSS exposes only intersection *curves*; isolated
    # intersection points are not part of the public result of this API.
    nb_points = 0
    worst = 0.0
    samples = 0
    for k in range(1, nb_lines + 1):
        c = it.Line(k)
        if c is None:
            continue
        r1, n = surf_residual(c, s1)
        r2, _ = surf_residual(c, s2)
        d = max(r1, r2)
        if not math.isnan(d):
            worst = max(worst, d)
        samples += n
    for k in range(1, nb_points + 1):
        p = it.Point(k)
        d = max(curve_point_residual(p, s1), curve_point_residual(p, s2))
        if not math.isnan(d):
            worst = max(worst, d)
    return done, nb_lines, nb_points, samples, worst, dt


def run_sc_fixture(tol):
    """surface/curve via GeomAPI_IntCS (IntPatch)."""
    cyl = quarter_cylinder()
    line = bezier_line((0.5, 0.5, 0.5), (1.1, 0.5, 0.5))
    t0 = time.perf_counter()
    # NOTE: GeomAPI_IntCS has no tolerance argument in OCCT 8.0 (neither C++ nor
    # the python binding); the internal IntPatch tolerance is used.
    it = GeomAPI_IntCS(line, cyl)
    dt = (time.perf_counter() - t0) * 1e3
    done = bool(it.IsDone())
    nb_points = it.NbPoints() if done else 0
    nb_seg = it.NbSegments() if done else 0
    worst = 0.0
    samples = 0
    for k in range(1, nb_points + 1):
        p = it.Point(k)
        d = max(curve_point_residual(p, cyl), curve_curve_residual(p, line))
        if not math.isnan(d):
            worst = max(worst, d)
        samples += 1
    for k in range(1, nb_seg + 1):
        c = it.Segment(k)
        if c is None:
            continue
        r1, n = surf_residual(c, cyl)
        r2 = 0.0
        u0, u1 = c.FirstParameter(), c.LastParameter()
        for j in range(NS + 1):
            u = u0 + (u1 - u0) * j / NS
            d = curve_curve_residual(c.Value(u), line)
            if not math.isnan(d):
                r2 = max(r2, d)
        worst = max(worst, r1, r2)
        samples += n
    return done, nb_points, nb_seg, samples, worst, dt


def main():
    print("OCCT_VERSION_INFO:", end=" ")
    try:
        import OCC
        print(getattr(OCC, "VERSION", "?"))
    except Exception as e:
        print("?", e)

    out = ["fixture,method,tol,is_done,nb_curves,nb_points,samples,"
           "max_residual,time_ms,truth"]
    for name in ["A", "B", "C", "D", "E"]:
        truth = FIXTURES[name][4]
        for tol in (1e-6, 1e-4):
            done, nl, np_, ns, worst, dt = run_ss_fixture(name, tol)
            out.append(f"{name},GeomAPI_IntSS,{tol:.0e},{int(done)},{nl},{np_},"
                       f"{ns},{worst:.3e},{dt:.3f},\"{truth}\"")
    truth = "1 point: (0.866025,0.5,0.5)"
    for tol in (1e-6, 1e-4):
        done, np_, nseg, ns, worst, dt = run_sc_fixture(tol)
        out.append(f"F,GeomAPI_IntCS,{tol:.0e},{int(done)},{nseg},{np_},"
                   f"{ns},{worst:.3e},{dt:.3f},\"{truth}\"")
    print("\n".join(out))


if __name__ == "__main__":
    main()
