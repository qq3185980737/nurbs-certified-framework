// ============================================================================
// nurbs_surface.hpp -- Unified NURBS / Bezier surface storage structure.
//   A single data representation shared by all surface workloads
//   (self-intersection, overlap extraction, curve x surface, surface x surface):
//     control point grid + u/v knot vectors + u/v degrees + weights
//   Evaluation / partial derivatives: tensor-product de Boor
//     (de Boor along u per control-row, then de Boor along v).
//   AABB: control-point convex hull (NURBS surface lies inside it).
//   Note: all current use cases have unit weights (non-rational); the weight
//     field is kept for future rational extension.
// ============================================================================
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

namespace nurbs {

struct SurfacePoint { double x = 0, y = 0, z = 0; };

// Single-variable B-spline de Boor evaluation.
// P: control points (n+1), U: knot vector (length n+p+2, clamped), p: degree.
inline SurfacePoint deBoor(const std::vector<SurfacePoint>& P, const std::vector<double>& U, int p, double u) {
    const int n = (int)P.size() - 1;
    if (u <= U[0]) return P[0];
    if (u >= U[U.size() - 1]) return P[n];
    int k = p;
    for (int i = p; i < n; ++i) { if (u < U[i + 1]) { k = i; break; } }
    std::vector<SurfacePoint> d(P);
    for (int r = 1; r <= p; ++r) {
        for (int i = k; i >= k - p + r; --i) {   // standard de Boor: descending i
            const double den = U[i + p - r + 1] - U[i];
            const double a = (den > 0) ? (u - U[i]) / den : 0.0;
            const double b = 1.0 - a;
            d[i].x = b * d[i-1].x + a * d[i].x;
            d[i].y = b * d[i-1].y + a * d[i].y;
            d[i].z = b * d[i-1].z + a * d[i].z;
        }
    }
    return d[k];
}

// Single-variable B-spline first derivative (de Boor on the derivative curve).
inline SurfacePoint deBoorDeriv(const std::vector<SurfacePoint>& P, const std::vector<double>& U, int p, double u) {
    const int n = (int)P.size() - 1;
    std::vector<SurfacePoint> Q(n);
    for (int i = 0; i < n; ++i) {
        const double den = U[i + p + 1] - U[i + 1];
        const double s = (den > 0) ? (double)p / den : 0.0;
        Q[i].x = s * (P[i+1].x - P[i].x);
        Q[i].y = s * (P[i+1].y - P[i].y);
        Q[i].z = s * (P[i+1].z - P[i].z);
    }
    std::vector<double> W(U.size() - 2);
    for (int i = 0; i < (int)W.size(); ++i) W[i] = U[i + 1];
    return deBoor(Q, W, p - 1, u);
}

// Unified NURBS / Bezier surface.
struct NurbsSurface {
    int nu = 0, nv = 0;                 // control grid size (nu x nv)
    int pu = 0, pv = 0;                 // u/v degrees
    std::vector<double> uKnots, vKnots; // length nu+pu+1, nv+pv+1
    std::vector<double> w;              // nu x nv weights (currently 1, non-rational)
    std::vector<SurfacePoint> C;        // nu x nv control points (row-major, C[j*nu+i], j=v / i=u)

    // S(u,v): de Boor along u per row -> nv middle points -> de Boor along v.
    SurfacePoint eval(double u, double v) const {
        std::vector<SurfacePoint> mid(nv);
        for (int j = 0; j < nv; ++j) {
            std::vector<SurfacePoint> row(nu);
            for (int i = 0; i < nu; ++i) row[i] = C[j * nu + i];
            mid[j] = deBoor(row, uKnots, pu, u);
        }
        return deBoor(mid, vKnots, pv, v);
    }
    // dS/du: derivative along u per row -> de Boor along v.
    SurfacePoint du(double u, double v) const {
        std::vector<SurfacePoint> mid(nv);
        for (int j = 0; j < nv; ++j) {
            std::vector<SurfacePoint> row(nu);
            for (int i = 0; i < nu; ++i) row[i] = C[j * nu + i];
            mid[j] = deBoorDeriv(row, uKnots, pu, u);
        }
        return deBoor(mid, vKnots, pv, v);
    }
    // dS/dv: de Boor along u per row -> middle points (a curve along v) -> derivative along v.
    SurfacePoint dv(double u, double v) const {
        std::vector<SurfacePoint> mid(nv);
        for (int j = 0; j < nv; ++j) {
            std::vector<SurfacePoint> row(nu);
            for (int i = 0; i < nu; ++i) row[i] = C[j * nu + i];
            mid[j] = deBoor(row, uKnots, pu, u);
        }
        return deBoorDeriv(mid, vKnots, pv, v);
    }
    // AABB (control-point convex hull).
    SurfacePoint aabbMin() const {
        SurfacePoint m{1e300, 1e300, 1e300};
        for (const auto& p : C) { m.x = std::min(m.x, p.x); m.y = std::min(m.y, p.y); m.z = std::min(m.z, p.z); }
        return m;
    }
    SurfacePoint aabbMax() const {
        SurfacePoint M{-1e300, -1e300, -1e300};
        for (const auto& p : C) { M.x = std::max(M.x, p.x); M.y = std::max(M.y, p.y); M.z = std::max(M.z, p.z); }
        return M;
    }

    // Biquadratic Bezier (3x3, pu=pv=2, clamped [0,0,0,1,1,1], unit weights).
    static NurbsSurface makeBezierBiquadratic(const SurfacePoint* cp9) {
        NurbsSurface s; s.nu = 3; s.nv = 3; s.pu = 2; s.pv = 2;
        s.uKnots = {0,0,0,1,1,1}; s.vKnots = {0,0,0,1,1,1};
        s.C.assign(cp9, cp9 + 9); s.w.assign(9, 1.0); return s;
    }
    // Bilinear patch (2x2, pu=pv=1, clamped [0,0,1,1], unit weights).
    static NurbsSurface makeBilinear(const SurfacePoint* cp4) {
        NurbsSurface s; s.nu = 2; s.nv = 2; s.pu = 1; s.pv = 1;
        s.uKnots = {0,0,1,1}; s.vKnots = {0,0,1,1};
        s.C.assign(cp4, cp4 + 4); s.w.assign(4, 1.0); return s;
    }
};

} // namespace nurbs
