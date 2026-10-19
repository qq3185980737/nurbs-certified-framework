// ============================================================================
// nurbs_surface.hpp -- Unified NURBS / Bezier surface storage structure.
//   A single data representation shared by all surface workloads
//   (self-intersection, overlap extraction, curve x surface, surface x surface):
//     control point grid + u/v knot vectors + u/v degrees + weights
//   Evaluation / partial derivatives: tensor-product de Boor
//     (de Boor along u per control-row, then de Boor along v).
//   AABB: control-point convex hull (NURBS surface lies inside it).
//   Weight support:
//     * eval/du/dv : non-rational (unit weights assumed), plain tensor de Boor.
//     * rationalEval/rationalDu/rationalDv : rational (homogeneous 4D) de Boor,
//       supports arbitrary weights and non-uniform knot vectors.
//   Closed geometry (e.g. full cylinder side): model as spliced clamped rational
//   arcs (makeCylinder) so that S(0,v)==S(1,v); evaluate with rationalEval.
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
    for (int i = p; i <= n; ++i) { if (u < U[i + 1]) { k = i; break; } }
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

// ---- Homogeneous 4D point (x*w, y*w, z*w, w) for rational de Boor ----
struct H4 { double x = 0, y = 0, z = 0, w = 0; };

// Rational de Boor (homogeneous 4D) on a clamped knot vector.
inline H4 deBoorR(const std::vector<H4>& P, const std::vector<double>& U, int p, double u) {
    const int n = (int)P.size() - 1;
    if (u <= U[0]) return P[0];
    if (u >= U[U.size() - 1]) return P[n];
    int k = p;
    for (int i = p; i <= n; ++i) { if (u < U[i + 1]) { k = i; break; } }
    std::vector<H4> d(P);
    for (int r = 1; r <= p; ++r) {
        for (int i = k; i >= k - p + r; --i) {
            const double den = U[i + p - r + 1] - U[i];
            const double a = (den > 0) ? (u - U[i]) / den : 0.0;
            const double b = 1.0 - a;
            d[i].x = b * d[i-1].x + a * d[i].x;
            d[i].y = b * d[i-1].y + a * d[i].y;
            d[i].z = b * d[i-1].z + a * d[i].z;
            d[i].w = b * d[i-1].w + a * d[i].w;
        }
    }
    return d[k];
}

// Rational de Boor first derivative (homogeneous 4D) on a clamped knot vector.
inline H4 deBoorRDeriv(const std::vector<H4>& P, const std::vector<double>& U, int p, double u) {
    const int n = (int)P.size() - 1;
    std::vector<H4> Q(n);
    for (int i = 0; i < n; ++i) {
        const double den = U[i + p + 1] - U[i + 1];
        const double s = (den > 0) ? (double)p / den : 0.0;
        Q[i].x = s * (P[i+1].x - P[i].x);
        Q[i].y = s * (P[i+1].y - P[i].y);
        Q[i].z = s * (P[i+1].z - P[i].z);
        Q[i].w = s * (P[i+1].w - P[i].w);
    }
    std::vector<double> W(U.size() - 2);
    for (int i = 0; i < (int)W.size(); ++i) W[i] = U[i + 1];
    return deBoorR(Q, W, p - 1, u);
}

// Unified NURBS / Bezier surface.
struct NurbsSurface {
    int nu = 0, nv = 0;                 // control grid size (nu x nv)
    int pu = 0, pv = 0;                 // u/v degrees
    std::vector<double> uKnots, vKnots; // length nu+pu+1, nv+pv+1 (clamped)
    std::vector<double> w;              // nu x nv weights
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

    // ---- Rational evaluation (arbitrary weights, non-uniform knots) ----
    // Build homogeneous 4D control points for one row j.
    std::vector<H4> rowH(int j) const {
        std::vector<H4> row(nu);
        for (int i = 0; i < nu; ++i) {
            const double wi = w[j * nu + i];
            row[i] = { C[j*nu+i].x * wi, C[j*nu+i].y * wi, C[j*nu+i].z * wi, wi };
        }
        return row;
    }
    SurfacePoint rationalEval(double u, double v) const {
        std::vector<H4> mid(nv);
        for (int j = 0; j < nv; ++j) mid[j] = deBoorR(rowH(j), uKnots, pu, u);
        const H4 H = deBoorR(mid, vKnots, pv, v);
        return { H.x / H.w, H.y / H.w, H.z / H.w };
    }
    // dS/du (rational): S = (X,W) -> Su = (Xu*W - X*Wu) / W^2.
    SurfacePoint rationalDu(double u, double v) const {
        std::vector<H4> midU(nv), midH(nv);
        for (int j = 0; j < nv; ++j) {
            const std::vector<H4> row = rowH(j);
            midU[j] = deBoorRDeriv(row, uKnots, pu, u);
            midH[j] = deBoorR(row, uKnots, pu, u);
        }
        const H4 Hu = deBoorR(midU, vKnots, pv, v);
        const H4 H  = deBoorR(midH, vKnots, pv, v);
        const double W = H.w, Wu = Hu.w, d = W * W;
        return { (Hu.x*W - H.x*Wu)/d, (Hu.y*W - H.y*Wu)/d, (Hu.z*W - H.z*Wu)/d };
    }
    // dS/dv (rational).
    SurfacePoint rationalDv(double u, double v) const {
        std::vector<H4> midH(nv);
        for (int j = 0; j < nv; ++j) midH[j] = deBoorR(rowH(j), uKnots, pu, u);
        const H4 Hv = deBoorRDeriv(midH, vKnots, pv, v);
        const H4 H  = deBoorR(midH, vKnots, pv, v);
        const double W = H.w, Wv = Hv.w, d = W * W;
        return { (Hv.x*W - H.x*Wv)/d, (Hv.y*W - H.y*Wv)/d, (Hv.z*W - H.z*Wv)/d };
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
    // Rational quarter-cylinder side: u = 90-degree arc (w = [1,sqrt(0.5),1]),
    // v = linear height. nu=3, nv=2, pu=2, pv=1.
    static NurbsSurface makeRationalQuarterCylinder(double R, double z0, double H) {
        NurbsSurface s; s.nu = 3; s.nv = 2; s.pu = 2; s.pv = 1;
        s.uKnots = {0,0,0,1,1,1}; s.vKnots = {0,0,1,1};
        const double r2 = std::sqrt(0.5);
        // arc control points in xy: (R,0), (R,R), (0,R)
        const double wu[3] = {1.0, r2, 1.0};
        s.C.resize(6); s.w.resize(6);
        for (int j = 0; j < 2; ++j) {
            const double z = z0 + H * (double)j;
            const double ax[3] = { R, R, 0.0 };
            const double ay[3] = { 0.0, R, R };
            for (int i = 0; i < 3; ++i) {
                s.C[j*3+i] = { ax[i], ay[i], z };
                s.w[j*3+i] = wu[i];
            }
        }
        return s;
    }
    // Closed cylinder side (full 360): u = 4 spliced 90-degree rational quarter arcs
    // (9 control points, w = [1,sqrt(0.5),1] per arc, clamped non-uniform knots),
    // v = linear height. Geometry is closed (S(0,v) == S(1,v)).
    static NurbsSurface makeCylinder(double R, double z0, double H) {
        NurbsSurface s; s.nu = 9; s.nv = 2; s.pu = 2; s.pv = 1;
        s.uKnots = {0,0,0, 0.25,0.25, 0.5,0.5, 0.75,0.75, 1,1,1};
        s.vKnots = {0,0,1,1};
        const double r2 = std::sqrt(0.5);
        // 9 control points: axis points + square corners (w = sqrt(0.5) at corners)
        const double cx[9] = { R, R, 0, -R, -R, -R, 0, R, R };
        const double cy[9] = { 0, R, R, R, 0, -R, -R, -R, 0 };
        const double wu[9] = {1.0, r2, 1.0, r2, 1.0, r2, 1.0, r2, 1.0};
        s.C.resize(18); s.w.resize(18);
        for (int i = 0; i < 9; ++i) for (int j = 0; j < 2; ++j) {
            const double z = z0 + H * (double)j;
            s.C[j*9+i] = { cx[i], cy[i], z };
            s.w[j*9+i] = wu[i];
        }
        return s;
    }

    // Upper unit sphere (radius R): surface of revolution. u = full circle
    // (9 spliced 90-degree rational arcs, closed S(0,v)==S(1,v)), v = 90-degree
    // meridian arc from equator (z=0) to north pole (z=R). Biquadratic rational,
    // weight = wu[i]*wv[j] (tensor product of a revolution body).
    static NurbsSurface makeSphere(double R) {
        NurbsSurface s; s.nu = 9; s.nv = 3; s.pu = 2; s.pv = 2;
        s.uKnots = {0,0,0, 0.25,0.25, 0.5,0.5, 0.75,0.75, 1,1,1};
        s.vKnots = {0,0,0, 1,1,1};
        const double r2 = std::sqrt(0.5);
        const double wu[9] = {1.0, r2, 1.0, r2, 1.0, r2, 1.0, r2, 1.0};   // u full circle
        const double wv[3] = {1.0, r2, 1.0};                              // v 90-deg meridian
        const double cx[9] = { R, R, 0, -R, -R, -R, 0, R, R };
        const double cy[9] = { 0, R, R, R, 0, -R, -R, -R, 0 };
        const double rz[3][2] = { {R,0}, {R,R}, {0,R} };                  // meridian ctrl pts (r,z)
        s.C.resize(27); s.w.resize(27);
        for (int j = 0; j < 3; ++j) for (int i = 0; i < 9; ++i) {
            const double fr = rz[j][0] / R;                              // r_j/R
            s.C[j*9+i] = { cx[i] * fr, cy[i] * fr, rz[j][1] };
            s.w[j*9+i] = wu[i] * wv[j];
        }
        return s;
    }
};

} // namespace nurbs
