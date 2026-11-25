#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Matrix-representation baseline (minimal reproduction).

WHAT THIS IS
------------
The paper's certified framework solves f1(s,t)-f2(s,t)=0 (graph-patch SSI) by
interval/Krawczyk style certified subdivision.  The matrix-representation line
of work (Thang-Buse-Mourrain 2009 for curve/surface; Jia-Li-Cheng, CAD 150
(2022) 103303 for surface/surface) instead turns the intersection into an
eigenvalue problem: the intersection parameters are eigenvalues of a matrix
whose entries are polynomials ("matrix representation"), so all roots are
obtained in one shot without subdivision.

This script implements the *matrix -> eigenvalue* step in the simplest possible
honest setting:

  (1) curve x bicubic graph surface  (the 2009 curve/surface setting, reduced):
      a 3D Bezier curve C(t)=(x(t),y(t),z(t)) intersects the graph patch
      S(u,v)=(u,v,f(u,v)) with f bicubic.  The intersection equation is the
      univariate polynomial
          g(t) = z(t) - f(x(t), y(t))                (polynomial curve)
          g(t) = Z(t) W^(p+q-1) - sum a_ij X^i Y^j W^(...)   (rational curve)
      whose coefficients are the matrix representation; the roots are computed
      as eigenvalues of the companion matrix (numpy.linalg.eig).

  (2) 2x2 graph-patch surface/surface (Jia 2022, minimal):  for two graph
      patches the intersection curve is {g(s,t)=0} with g = f1-f2 bidegree
      (2,2) (or (3,3)).  For each t the s-roots are eigenvalues of the 2x2 (or
      3x3) companion matrix of g(.,t) -- the "2x2 graph-patch" eigenvalue
      version requested as a minimal reproduction of Jia et al. 2022.

WHAT THIS IS NOT
----------------
This is *not* a full implementation of Jia-Li-Cheng 2022 (or of Thang-Buse-
Mourrain 2009): the general construction builds the matrix representation from
moving surfaces / Dixon resultants of two bidegree-(2,2) rational
parameterizations, which needs multivariate resultant machinery we did not
implement.  The graph-patch specialization lets us skip implicitization (the
surface is already a height field), so only the eigenvalue step is reproduced.
All numbers below are produced by this script and must be cited as
"minimal reproduction", not as the published algorithm's performance.

Outputs CSV to stdout (see scripts/run_matrix_baseline.sh).
"""
import math
import time
import numpy as np

TOL = 1e-9


# ---------------------------------------------------------------- polynomials
# bivariate polynomial: dict {(i,j): coeff} of s^i t^j
def p2_add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0.0) + v
    return {k: v for k, v in out.items() if v != 0.0}


def p2_mul(a, b):
    out = {}
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            k = (i1 + i2, j1 + j2)
            out[k] = out.get(k, 0.0) + c1 * c2
    return {k: v for k, v in out.items() if v != 0.0}


def p2_pow(a, n):
    r = {(0, 0): 1.0}
    for _ in range(n):
        r = p2_mul(r, a)
    return r


def bern_mono(n):
    """B_i^n(u) in the monomial basis: dict i -> {k: coeff}, i,k = 0..n."""
    # B_i^n(u) = C(n,i) u^i (1-u)^{n-i} = C(n,i) sum_k (-1)^k C(n-i,k) u^{i+k}
    from math import comb
    out = []
    for i in range(n + 1):
        poly = {}
        for k in range(n - i + 1):
            poly[i + k] = comb(n, i) * ((-1) ** k) * comb(n - i, k)
        out.append(poly)
    return out


def bezier_net_to_power(z):
    """z[i][j] Bernstein net on [0,1]^2 -> dict {(i,j): coeff} monomial."""
    n = len(z) - 1
    bi = bern_mono(n)
    out = {}
    for i in range(n + 1):
        for j in range(n + 1):
            c = z[i][j]
            if c == 0.0:
                continue
            for k1, c1 in bi[i].items():
                for k2, c2 in bi[j].items():
                    key = (k1, k2)
                    out[key] = out.get(key, 0.0) + c * c1 * c2
    return {k: v for k, v in out.items() if abs(v) > 1e-15}


def poly1_mul(a, b):
    """1D polynomial (numpy array, ascending) multiplication."""
    return np.convolve(a, b)


def poly1_pow(a, n):
    r = np.array([1.0])
    for _ in range(n):
        r = np.convolve(r, a)
    return r


# ------------------------------------------------- Jacobi/Durand-Kerner roots
def poly_roots_eig(coeffs):
    """Roots of a univariate polynomial as eigenvalues of its companion matrix.

    coeffs: ascending order [a0, a1, ..., an].  Trailing (numerically) zero
    coefficients are stripped first (they correspond to a lower degree).
    """
    a = np.asarray(coeffs, dtype=float)
    nz = np.nonzero(np.abs(a) > 1e-14)[0]
    if nz.size == 0:
        return np.array([])
    a = a[: nz[-1] + 1]
    an = a[-1]
    if abs(an) < 1e-300:
        raise ValueError("leading coefficient ~ 0")
    c = -a[:-1] / an              # companion (first) row: -a_{n-1}/a_n ...
    n = len(c)
    C = np.zeros((n, n))
    C[0, :] = c[::-1]
    if n > 1:
        C[1:, :-1] = np.eye(n - 1)
    return np.linalg.eigvals(C)


def peval(a, x):
    """Evaluate ascending-coefficient polynomial at x."""
    return np.polyval(np.asarray(a, dtype=float)[::-1], x)


def pder(a):
    """Derivative of ascending-coefficient polynomial (ascending)."""
    return np.polyder(np.asarray(a, dtype=float)[::-1])[::-1]


def newton_refine(coeffs, x, tol=1e-13, itmax=60):
    a = np.asarray(coeffs, dtype=float)
    d = pder(a)
    for _ in range(itmax):
        fx = peval(a, x)
        dfx = peval(d, x) if len(d) else 0.0
        if dfx == 0:
            break
        dx = fx / dfx
        x -= dx
        if abs(dx) < tol:
            break
    return x


# ---------------------------------------------------------------- fixture data
def graph_fixture(name):
    """Returns (dict monomial poly f, n) for graph patch z=f(u,v)."""
    if name == "A":
        return {(2, 0): 1.0, (0, 2): -1.0}, 3
    if name == "B":
        return {(0, 0): 0.3}, 2
    if name == "C":
        return {(1, 0): 1.0, (0, 1): 0.5}, 2
    if name == "D":
        # (u-v)^2
        return {(2, 0): 1.0, (1, 1): -2.0, (0, 2): 1.0}, 3
    if name == "E":
        return {(3, 0): 1.0, (0, 3): 1.0, (0, 0): -0.1}, 4
    raise KeyError(name)


def bih(f, name):
    return f[name]


# ============================================================== (1) curve/surf
def _padd(a, n):
    out = np.zeros(n)
    a = np.asarray(a, dtype=float)
    out[: len(a)] = a
    return out


def _padd_sum(a, b):
    n = max(len(a), len(b))
    return _padd(a, n) + _padd(b, n)


def curve_surface_polynomial(cx, cy, cz, cw, f):
    """g(t) = Z W^(p+q-1) - sum a_ij X^i Y^j W^(p+q-i-j), ascending coeffs.

    cx,cy,cz,cw: 1D numpy coefficient arrays of the homogeneous Bezier curve
    (ascending powers of t).  f: dict monomial poly of the graph height field.
    The curve is rational; for a polynomial curve set cw=[1].
    """
    p = max(i for (i, j) in f.keys())
    q = max(j for (i, j) in f.keys())
    d = p + q
    g = poly1_mul(cz, poly1_pow(cw, d - 1)) if d >= 1 else np.asarray(cz, float)
    for (i, j), a in f.items():
        term = poly1_pow(cx, i)
        term = poly1_mul(term, poly1_pow(cy, j))
        term = poly1_mul(term, poly1_pow(cw, d - i - j))
        g = _padd_sum(g, -a * term)
    return g


def cluster_count(vals, tol=1e-6):
    """Number of distinct clusters among sorted real values (multiplicity-aware
    deduplication: a double eigenvalue gives one cluster)."""
    if not vals:
        return 0
    v = sorted(vals)
    n = 1
    for i in range(1, len(v)):
        if v[i] - v[i - 1] > tol:
            n += 1
    return n


def curve_surface_case(name, cx, cy, cz, cw, f, truth_roots, tol=TOL):
    g = curve_surface_polynomial(cx, cy, cz, cw, f)
    t0 = time.perf_counter()
    roots = poly_roots_eig(g)
    dt = (time.perf_counter() - t0) * 1e3
    found = []
    for r in roots:
        if abs(r.imag) > 1e-7:
            continue
        tr = newton_refine(g, r.real)
        if not (0.0 - 1e-9 <= tr <= 1.0 + 1e-9):
            continue
        found.append(tr)
    found.sort()
    distinct = cluster_count(found)
    # residual: |C(t) - S(x(t),y(t))| with rational evaluation
    def evalc(t):
        num = [peval(cx, t), peval(cy, t), peval(cz, t)]
        den = peval(cw, t)
        return np.array(num) / den
    worst = 0.0
    for t in found:
        P = evalc(t)
        zs = sum(a * P[0] ** i * P[1] ** j for (i, j), a in f.items())
        worst = max(worst, abs(P[2] - zs))
    return name, len(found), distinct, len(roots), truth_roots, worst, dt, found


def make_bezier_curve_1d(control, degree):
    """control: list of coefficients (in the monomial basis, ascending) of the
    coordinate as a function of t; converts to Bezier by evaluating at the
    Bernstein basis?  We keep everything in the monomial basis, so this is just
    a helper returning np.array(control)."""
    return np.asarray(control, dtype=float)


def cs_cases():
    """4 curve x bicubic-graph-surface fixtures, defined by explicit polynomials.

    f(u,v) = 0.5 u^3 - 0.4 v^3 + 0.15 u v + 0.05  (bicubic graph patch)

    Each case returns (name, cx, cy, cz, cw, truth_root_count) where the
    coordinates are *monomial* coefficient arrays (ascending powers of t) of the
    homogeneous Bezier curve.
    """
    out = []
    # (a) transversal single root: x=0.2+0.6t, y=0.3+0.4t, z=0.05+0.4t
    out.append(("CS1_transverse",
                np.array([0.2, 0.6]), np.array([0.3, 0.4]),
                np.array([0.05, 0.4]), np.array([1.0]), 1))
    # (b) 3 roots: x=t, y=0.5, z = f(t,0.5) + 1.5*(t-0.2)(t-0.5)(t-0.8)
    #     f(t,0.5) = 0.5 t^3 + 0.075 t
    #     g(t) = 1.5(t-0.2)(t-0.5)(t-0.8) -> roots 0.2, 0.5, 0.8
    prod = np.array([1.0])
    for r in (0.2, 0.5, 0.8):
        prod = np.convolve(prod, np.array([-r, 1.0]))
    zt = np.array([0.0, 0.075, 0.0, 0.5]) + 1.5 * prod
    out.append(("CS3_three_roots", np.array([0.0, 1.0]), np.array([0.5]),
                zt, np.array([1.0]), 3))
    # (c) tangential double root: x=t, y=0.5, z = f(t,0.5) + (t-0.5)^2
    zt = np.array([0.0, 0.075, 0.0, 0.5]) + np.array([0.25, -1.0, 1.0, 0.0])
    out.append(("CS1_tangential_double", np.array([0.0, 1.0]), np.array([0.5]),
                zt, np.array([1.0]), 1))
    # (d) rational quarter arc: x=(1-t^2)/(1+t^2), z=2t/(1+t^2), y=0.5
    out.append(("CS1_rational_arc",
                np.array([1.0, 0.0, -1.0]), np.array([0.5]),
                np.array([0.0, 2.0]), np.array([1.0, 0.0, 1.0]), 1))
    return out


def run_curve_surface():
    # bicubic graph: f(u,v) = 0.5 u^3 - 0.4 v^3 + 0.15 u v + 0.05
    f = {(3, 0): 0.5, (0, 3): -0.4, (1, 1): 0.15, (0, 0): 0.05}
    out = []
    for name, cx_, cy_, cz_, cw_, truth in cs_cases():
        out.append(curve_surface_case(name, cx_, cy_, cz_, cw_, f, truth))
    return out


# ============================================================ (2) surf/surf 2x2
def ss_poly(name):
    f1 = graph_fixture(name)[0]
    if name in ("A", "B", "C", "D"):
        f2 = {}
    else:
        f2 = {}
    return f1  # S2 is z=0 for A-D; E's g already includes -0.1


def ss_reduce_in_s(g):
    """g(s,t) -> coefficients a_i(t) as 1D numpy polys (ascending in t),
    deg_s = max i."""
    deg_s = max(i for (i, j) in g.keys())
    deg_t = max(j for (i, j) in g.keys())
    A = [np.zeros(deg_t + 1) for _ in range(deg_s + 1)]
    for (i, j), c in g.items():
        A[i][j] += c
    return A


def companion_roots_at(A, tval):
    """Eigenvalues of the companion matrix of sum_i a_i(t) s^i at t=tval."""
    coeffs = np.array([np.polyval(A[i][::-1], tval) for i in range(len(A))])
    # strip (numerically) zero leading coefficients
    while len(coeffs) > 1 and abs(coeffs[-1]) < 1e-14:
        coeffs = coeffs[:-1]
    if len(coeffs) == 1:
        return np.array([])
    if len(coeffs) == 2:  # explicit 1x1 companion = -a0/a1
        return np.array([-coeffs[0] / coeffs[1]])
    return poly_roots_eig(coeffs)


def ss_case(name, nt=201, tol=1e-9):
    g = ss_poly(name)
    A = ss_reduce_in_s(g)
    deg_t = max(j for (i, j) in g.keys())
    t0 = time.perf_counter()
    per_t = []
    roots_total = 0
    real_curve_pts = 0
    complex_pts = 0
    worst = 0.0
    samples = 0
    for k in range(nt):
        tv = k / (nt - 1)
        ev = companion_roots_at(A, tv)
        nreal = 0
        for r in ev:
            samples += 1
            if abs(r.imag) > 1e-7:
                complex_pts += 1
                continue
            sv = newton_refine(np.array([np.polyval(A[i][::-1], tv) for i in range(len(A))]), r.real)
            if not (-1e-9 <= sv <= 1 + 1e-9):
                continue
            nreal += 1
            # residual |g(s,t)|
            val = sum(c * sv ** i * tv ** j for (i, j), c in g.items())
            worst = max(worst, abs(val))
            real_curve_pts += 1
        roots_total += nreal
        per_t.append(nreal)
    dt = (time.perf_counter() - t0) * 1e3
    # the intersection curve is 1D: for each t at most deg_s roots lie in [0,1]
    return name, nt, roots_total, real_curve_pts, complex_pts, worst, dt, per_t


# =============================================================== certified 1D
def interval_eval_1d(coeffs, lo, hi):
    """Rigorous interval Horner evaluation of a 1D polynomial on [lo,hi]."""
    if len(coeffs) == 0:
        return 0.0, 0.0
    rlo = rhi = coeffs[-1]
    for c in coeffs[-2::-1]:
        xlo, xhi = lo if rlo >= 0 else hi, hi if rhi >= 0 else lo
        # multiply interval [rlo,rhi] by [lo,hi]
        p = [rlo * lo, rlo * hi, rhi * lo, rhi * hi]
        pmin, pmax = min(p), max(p)
        rlo = np.nextafter(pmin + c, -np.inf)
        rhi = np.nextafter(pmax + c, np.inf)
    return rlo, rhi


def certified_1d_count(coeffs, tol=1e-12):
    """Certified root count/cover for 1D: recursive interval exclusion.
    Returns (cover_boxes, excluded_boxes, distinct_components).  A box with 0
    not in the rigorous interval image is excluded; a box narrower than tol that
    still contains 0 is part of the cover.  Components are the maximal runs of
    adjacent cover boxes (a rigorous count of root *clusters*)."""
    stack = [(0.0, 1.0)]
    cover = 0
    excl = 0
    boxes = []
    while stack:
        lo, hi = stack.pop()
        ilo, ihi = interval_eval_1d(coeffs, lo, hi)
        if not (ilo <= 0.0 <= ihi):
            excl += 1
            continue
        if hi - lo <= tol:
            cover += 1
            boxes.append((lo, hi))
            continue
        mid = 0.5 * (lo + hi)
        stack.append((lo, mid))
        stack.append((mid, hi))
    boxes.sort()
    comps = 0
    last_hi = -1.0
    for lo, hi in boxes:
        if lo > last_hi + 1e-15:
            comps += 1
        last_hi = hi
    return cover, excl, comps


def run_certified_curve_surface():
    """Certified root count for the same 4 curve/surface fixtures."""
    f = {(3, 0): 0.5, (0, 3): -0.4, (1, 1): 0.15, (0, 0): 0.05}
    out = []
    for name, cx, cy, cz, cw, truth in cs_cases():
        g = curve_surface_polynomial(cx, cy, cz, cw, f)
        t0 = time.perf_counter()
        cover, excl, comps = certified_1d_count(g, 1e-10)
        dt = (time.perf_counter() - t0) * 1e3
        out.append((name, cover, excl, comps, dt))
    return out


def main():
    print("# matrix-representation baseline -- MINIMAL REPRODUCTION")
    print("# not a full implementation of Thang-Buse-Mourrain 2009 / Jia et al. 2022")
    print("# numpy", np.__version__)
    print()
    print("=== (1) curve x bicubic graph surface: eigenvalues of companion matrix ===")
    print("case,n_real_roots_in_[0,1],n_distinct_roots,n_eigenvalues_total,"
          "n_roots_truth,max_residual,time_ms,roots")
    cs = run_curve_surface()
    for name, found, distinct, nev, truth, worst, dt, roots in cs:
        rs = " ".join(f"{r:.6f}" for r in roots)
        print(f'{name},{found},{distinct},{nev},{truth},{worst:.3e},{dt:.3f},"{rs}"')
    print()
    print("=== (1b) certified interval exclusion on the same reduced polynomial ===")
    print("case,cover_boxes,excluded_boxes,certified_root_clusters,time_ms")
    for name, cover, excl, comps, dt in run_certified_curve_surface():
        print(f"{name},{cover},{excl},{comps},{dt:.3f}")
    print()
    print("=== (2) 2x2/3x3 graph-patch surface/surface (Jia 2022 minimal) ===")
    print("fixture,nt,roots_total,curve_points,complex_pts,max_residual,time_ms,"
          "roots_per_t_min,roots_per_t_max")
    for name in ["A", "B", "C", "D", "E"]:
        n, nt, total, pts, cplx, worst, dt, per_t = ss_case(name)
        print(f"{name},{nt},{total},{pts},{cplx},{worst:.3e},{dt:.3f},"
              f"{min(per_t)},{max(per_t)}")


if __name__ == "__main__":
    main()
