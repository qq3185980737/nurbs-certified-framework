// ============================================================================
// cert_ssi.hpp — 一般 NURBS 曲面求交的"认证覆盖"内核
// ----------------------------------------------------------------------------
// M3：有理 Bezier 面片（齐次区间控制网）、精确子片细分、凸包包络、4D 覆盖
// M4：横截认证 —— 3x4 区间 Jacobian 的广义逆 Krawczyk + Gershgorin 秩判定
// ============================================================================
#ifndef NURBS_CERT_SSI_HPP
#define NURBS_CERT_SSI_HPP

#include <array>
#include <vector>
#include <limits>
#include <algorithm>
#include <cstddef>

#include "nurbs/interval_arith.hpp"

namespace nurbs {
namespace cert {

using HP = std::array<Interval, 4>;   // 齐次控制点 (wx, wy, wz, w)

// ---------- 有理 Bezier 面片（控制网行优先：i*nv + j）----------
struct RBezierSurface {
    int nu = 0, nv = 0;
    double wmin = 0.0, wmax = 1.0;   // 控制权重范围（正权假设）
    std::vector<HP> P;

    const HP& at(int i, int j) const { return P[(std::size_t)i * nv + j]; }
    void update_weight_range() {
        wmin = 1e300; wmax = -1e300;
        for (const HP& p : P) {
            wmin = std::min(wmin, p[3].lo);
            wmax = std::max(wmax, p[3].hi);
        }
    }
};

inline RBezierSurface make_surface(int nu, int nv,
                                   const std::vector<std::array<double, 4>>& poles) {
    RBezierSurface s;
    s.nu = nu; s.nv = nv;
    s.P.resize((std::size_t)nu * nv);
    for (std::size_t k = 0; k < s.P.size(); ++k)
        for (int c = 0; c < 4; ++c)
            s.P[k][c] = Interval::point(poles[k][c]);
    s.update_weight_range();
    return s;
}

// ---------- 精确子片细分（区间 de Casteljau）----------
inline void split_1d(const std::vector<HP>& pts, std::vector<HP>& L, std::vector<HP>& R) {
    const int n = (int)pts.size();
    const Interval h(0.5, 0.5);
    std::vector<std::vector<HP>> tri;
    tri.push_back(pts);
    for (int r = 1; r < n; ++r) {
        std::vector<HP> row(n - r);
        for (int i = 0; i < n - r; ++i)
            for (int c = 0; c < 4; ++c)
                row[i][c] = iadd(imul(tri[r - 1][i][c], h), imul(tri[r - 1][i + 1][c], h));
        tri.push_back(std::move(row));
    }
    L.resize(n); R.resize(n);
    for (int r = 0; r < n; ++r) { L[r] = tri[r][0]; R[r] = tri[n - 1 - r][r]; }
}

inline void split_u(const RBezierSurface& s, RBezierSurface& L, RBezierSurface& R) {
    L.nu = s.nu; L.nv = s.nv; R.nu = s.nu; R.nv = s.nv;
    L.P.resize(s.P.size()); R.P.resize(s.P.size());
    for (int j = 0; j < s.nv; ++j) {
        std::vector<HP> col(s.nu), l, r;
        for (int i = 0; i < s.nu; ++i) col[i] = s.at(i, j);
        split_1d(col, l, r);
        for (int i = 0; i < s.nu; ++i) {
            L.P[(std::size_t)i * s.nv + j] = l[i];
            R.P[(std::size_t)i * s.nv + j] = r[i];
        }
    }
    L.update_weight_range();
    R.update_weight_range();
}

inline void split_v(const RBezierSurface& s, RBezierSurface& L, RBezierSurface& R) {
    L.nu = s.nu; L.nv = s.nv; R.nu = s.nu; R.nv = s.nv;
    L.P.resize(s.P.size()); R.P.resize(s.P.size());
    for (int i = 0; i < s.nu; ++i) {
        std::vector<HP> row(s.nv), l, r;
        for (int j = 0; j < s.nv; ++j) row[j] = s.at(i, j);
        split_1d(row, l, r);
        for (int j = 0; j < s.nv; ++j) {
            L.P[(std::size_t)i * s.nv + j] = l[j];
            R.P[(std::size_t)i * s.nv + j] = r[j];
        }
    }
    L.update_weight_range();
    R.update_weight_range();
}

// ============================================================================
// P7：NURBS 面片（齐次控制网 + 节点向量）→ 精确 Bézier 分片
//   方法：Boehm 节点插入，把每个内部节点重数升到次数，再按节点区间切片。
//   要求：两端夹紧（clamped）；权重为正（与全局假设一致）。
//   返回的 RBezierSurface 均为局部 [0,1]² 参数。
// ============================================================================
struct NurbsPatch {
    int du = 0, dv = 0, nu = 0, nv = 0;
    std::vector<HP> P;                 // 行优先 i*nv+j（齐次）
    std::vector<double> U, V;          // 节点向量
    const HP& at(int i, int j) const { return P[(std::size_t)i * nv + j]; }
};

inline void insert_knot_u(const std::vector<HP>& Pin, int nu, int nv, const std::vector<double>& U,
                          int du, double t, std::vector<HP>& Pout, std::vector<double>& Uout) {
    int k = du;
    for (int i = du; i < nu; ++i) if (U[i] <= t && t < U[i + 1]) { k = i; break; }
    Pout.assign((std::size_t)(nu + 1) * nv, HP{});
    Uout = U;
    Uout.insert(Uout.begin() + k + 1, t);
    for (int j = 0; j < nv; ++j) {
        for (int i = 0; i <= k - du; ++i) Pout[(std::size_t)i * nv + j] = Pin[(std::size_t)i * nv + j];
        for (int i = k - du + 1; i <= k; ++i) {
            double den = U[i + du] - U[i];
            double alpha = (den > 0.0) ? (t - U[i]) / den : 0.0;
            for (int c = 0; c < 4; ++c)
                Pout[(std::size_t)i * nv + j][c] =
                    iadd(imul(Interval::point(alpha), Pin[(std::size_t)i * nv + j][c]),
                         imul(Interval::point(1.0 - alpha), Pin[(std::size_t)(i - 1) * nv + j][c]));
        }
        for (int i = k + 1; i <= nu; ++i) Pout[(std::size_t)i * nv + j] = Pin[(std::size_t)(i - 1) * nv + j];
    }
}

inline void insert_knot_v(const std::vector<HP>& Pin, int nu, int nv, const std::vector<double>& V,
                          int dv, double t, std::vector<HP>& Pout, std::vector<double>& Vout) {
    int k = dv;
    for (int i = dv; i < nv; ++i) if (V[i] <= t && t < V[i + 1]) { k = i; break; }
    Pout.assign((std::size_t)nu * (nv + 1), HP{});
    Vout = V;
    Vout.insert(Vout.begin() + k + 1, t);
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j <= k - dv; ++j) Pout[(std::size_t)i * (nv + 1) + j] = Pin[(std::size_t)i * nv + j];
        for (int j = k - dv + 1; j <= k; ++j) {
            double den = V[j + dv] - V[j];
            double alpha = (den > 0.0) ? (t - V[j]) / den : 0.0;
            for (int c = 0; c < 4; ++c)
                Pout[(std::size_t)i * (nv + 1) + j][c] =
                    iadd(imul(Interval::point(alpha), Pin[(std::size_t)i * nv + j][c]),
                         imul(Interval::point(1.0 - alpha), Pin[(std::size_t)i * nv + j - 1][c]));
        }
        for (int j = k + 1; j <= nv; ++j) Pout[(std::size_t)i * (nv + 1) + j] = Pin[(std::size_t)i * nv + j - 1];
    }
}

// 齐次 de Boor 求值（验证用）
inline void nurbs_eval_h(const NurbsPatch& S, double u, double v, HP& out) {
    auto find_span = [](const std::vector<double>& K, int n, int p, double t) {
        if (t >= K[n]) return n - 1;
        int lo = p, hi = n;
        while (hi - lo > 1) { int mid = (lo + hi) / 2; if (t < K[mid]) hi = mid; else lo = mid; }
        return lo;
    };
    std::vector<HP> col(S.nu);
    int spanV = find_span(S.V, S.nv, S.dv, v);
    for (int i = 0; i < S.nu; ++i) {
        std::vector<HP> d(S.dv + 1);
        for (int k = 0; k <= S.dv; ++k) d[k] = S.at(i, spanV - S.dv + k);
        for (int r = 1; r <= S.dv; ++r)
            for (int k = S.dv; k >= r; --k) {
                int idx = spanV - S.dv + k;
                double den = S.V[idx + S.dv - r + 1] - S.V[idx];
                double a = (den > 0.0) ? (v - S.V[idx]) / den : 0.0;
                for (int c = 0; c < 4; ++c)
                    d[k][c] = iadd(imul(Interval::point(a), d[k][c]),
                                   imul(Interval::point(1.0 - a), d[k - 1][c]));
            }
        col[i] = d[S.dv];
    }
    int spanU = find_span(S.U, S.nu, S.du, u);
    std::vector<HP> d(S.du + 1);
    for (int k = 0; k <= S.du; ++k) d[k] = col[spanU - S.du + k];
    for (int r = 1; r <= S.du; ++r)
        for (int k = S.du; k >= r; --k) {
            int idx = spanU - S.du + k;
            double den = S.U[idx + S.du - r + 1] - S.U[idx];
            double a = (den > 0.0) ? (u - S.U[idx]) / den : 0.0;
            for (int c = 0; c < 4; ++c)
                d[k][c] = iadd(imul(Interval::point(a), d[k][c]),
                               imul(Interval::point(1.0 - a), d[k - 1][c]));
        }
    out = d[S.du];
}

// 转 Bézier 分片
inline std::vector<RBezierSurface> to_bezier_patches(const NurbsPatch& S) {
    std::vector<RBezierSurface> out;
    std::vector<HP> P = S.P;
    std::vector<double> U = S.U, V = S.V;
    int nu = S.nu, nv = S.nv;
    std::vector<double> distinctU;
    for (std::size_t i = 1; i + 1 < U.size(); ++i)
        if (U[i] > U[0] && U[i] < U[U.size() - 1]) distinctU.push_back(U[i]);
    for (double t : distinctU) {
        int mult = 0;
        for (double x : U) if (x == t) ++mult;
        while (mult < S.du) {
            std::vector<HP> Pn; std::vector<double> Un;
            insert_knot_u(P, nu, nv, U, S.du, t, Pn, Un);
            P = Pn; U = Un; ++nu; ++mult;
        }
    }
    std::vector<double> distinctV;
    for (std::size_t i = 1; i + 1 < V.size(); ++i)
        if (V[i] > V[0] && V[i] < V[V.size() - 1]) distinctV.push_back(V[i]);
    for (double t : distinctV) {
        int mult = 0;
        for (double x : V) if (x == t) ++mult;
        while (mult < S.dv) {
            std::vector<HP> Pn; std::vector<double> Vn;
            insert_knot_v(P, nu, nv, V, S.dv, t, Pn, Vn);
            P = Pn; V = Vn; ++nv; ++mult;
        }
    }
    std::vector<double> uk, vk;
    for (std::size_t i = 0; i + 1 < U.size(); ++i) if (U[i] < U[i + 1]) uk.push_back(U[i]);
    uk.push_back(U[U.size() - 1]);
    for (std::size_t i = 0; i + 1 < V.size(); ++i) if (V[i] < V[i + 1]) vk.push_back(V[i]);
    vk.push_back(V[V.size() - 1]);
    int segU = (int)uk.size() - 1, segV = (int)vk.size() - 1;
    for (int su = 0; su < segU; ++su) {
        for (int sv = 0; sv < segV; ++sv) {
            RBezierSurface patch;
            patch.nu = S.du + 1; patch.nv = S.dv + 1;
            patch.P.resize((std::size_t)patch.nu * patch.nv);
            for (int i = 0; i <= S.du; ++i)
                for (int j = 0; j <= S.dv; ++j)
                    patch.P[(std::size_t)i * patch.nv + j] = P[(std::size_t)(su * S.du + i) * nv + sv * S.dv + j];
            patch.update_weight_range();
            out.push_back(patch);
        }
    }
    return out;
}


// ============================================================================
// P9：退化片分类（点片 / 共线片 / 非正权）
//   - 非正权：hull_aabb 会回退全域 AABB（cover 仍是无遗漏的超集，但不排除任何盒）；
//   - 点片：像为单点 → 交集退化为"点是否在另一片上"（点包含可判）；
//   - 共线片：像为线段 → 交集退化为曲线-曲面（可退化为 1D 问题或点包含）。
// ============================================================================
struct Degeneracy {
    bool nonpositive_weight = false;
    bool point = false;
    bool collinear = false;
};

inline Degeneracy classify_degeneracy(const RBezierSurface& S) {
    Degeneracy d;
    std::vector<std::array<double, 3>> pts;
    pts.reserve(S.P.size());
    for (const HP& p : S.P) {
        if (p[3].lo <= 0.0 || p[3].hi <= 0.0) d.nonpositive_weight = true;
        if (!(p[3].mid() > 0.0)) continue;
        pts.push_back({ p[0].mid() / p[3].mid(),
                        p[1].mid() / p[3].mid(),
                        p[2].mid() / p[3].mid() });
    }
    if (pts.size() < 2) { d.point = true; return d; }
    // 直径（最大点对距离）
    double scale = 0.0;
    for (const auto& q : pts)
        scale = std::max(scale, std::fabs(q[0]) + std::fabs(q[1]) + std::fabs(q[2]));
    double diam = 0.0;
    std::size_t ia = 0, ib = 0;
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t j = i + 1; j < pts.size(); ++j) {
            double dx = pts[i][0] - pts[j][0], dy = pts[i][1] - pts[j][1], dz = pts[i][2] - pts[j][2];
            double dd = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (dd > diam) { diam = dd; ia = i; ib = j; }
        }
    if (diam <= 1e-12 * std::max(1.0, scale)) { d.point = true; return d; }
    // 共线：所有点到 ia-ib 直线的距离
    double dx = pts[ib][0] - pts[ia][0], dy = pts[ib][1] - pts[ia][1], dz = pts[ib][2] - pts[ia][2];
    double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    double maxdev = 0.0;
    for (const auto& q : pts) {
        double ex = q[0] - pts[ia][0], ey = q[1] - pts[ia][1], ez = q[2] - pts[ia][2];
        double cx = ey * dz - ez * dy, cy = ez * dx - ex * dz, cz = ex * dy - ey * dx;
        maxdev = std::max(maxdev, std::sqrt(cx * cx + cy * cy + cz * cz) / len);
    }
    if (maxdev <= 1e-12 * std::max(1.0, scale)) d.collinear = true;
    return d;
}

// ---------- 3D AABB 包络 ----------
struct AABB3 { Interval x, y, z; };

inline AABB3 hull_aabb(const RBezierSurface& s) {
    if (s.P.empty()) {
        const double inf = std::numeric_limits<double>::infinity();
        return { Interval(-inf, inf), Interval(-inf, inf), Interval(-inf, inf) };
    }
    AABB3 a{ Interval::empty(), Interval::empty(), Interval::empty() };
    for (const HP& p : s.P) {
        Interval w = p[3];
        Interval cx = idiv(p[0], w), cy = idiv(p[1], w), cz = idiv(p[2], w);
        if (cx.isEmpty() || cy.isEmpty() || cz.isEmpty()) {
            const double inf = std::numeric_limits<double>::infinity();
            return { Interval(-inf, inf), Interval(-inf, inf), Interval(-inf, inf) };
        }
        if (a.x.isEmpty()) { a.x = cx; a.y = cy; a.z = cz; }
        else {
            a.x.lo = std::min(a.x.lo, cx.lo); a.x.hi = std::max(a.x.hi, cx.hi);
            a.y.lo = std::min(a.y.lo, cy.lo); a.y.hi = std::max(a.y.hi, cy.hi);
            a.z.lo = std::min(a.z.lo, cz.lo); a.z.hi = std::max(a.z.hi, cz.hi);
        }
    }
    return a;
}

inline bool disjoint(const AABB3& a, const AABB3& b) {
    return (a.x.hi < b.x.lo) || (b.x.hi < a.x.lo) ||
           (a.y.hi < b.y.lo) || (b.y.hi < a.y.lo) ||
           (a.z.hi < b.z.lo) || (b.z.hi < a.z.lo);
}

// ---------- 4D 参数盒与覆盖 ----------
struct Box4 { Interval u, v, s, t; };

struct CoverResult {
    long excl = 0, cover = 0, sub = 0, nodes = 0;
    bool cap_hit = false;
    bool record_excluded = false;          // 测试用：记录被排除盒
    std::vector<Box4> boxes;
    std::vector<Box4> excluded;
};

namespace detail {
inline void rec(const RBezierSurface& A, const RBezierSurface& B, const Box4& box,
                int depth, double tau, int max_depth, long node_cap, CoverResult& out) {
    if (out.nodes > node_cap) {
        // 资源上限：不得丢盒——把当前盒保守并入覆盖（保持 Σ ⊆ ⋃ 覆盖盒）
        out.cap_hit = true;
        ++out.cover;
        out.boxes.push_back(box);
        return;
    }
    ++out.nodes;
    if (disjoint(hull_aabb(A), hull_aabb(B))) {
        ++out.excl;
        if (out.record_excluded) out.excluded.push_back(box);
        return;
    }
    const double wu = box.u.width(), wv = box.v.width();
    const double ws = box.s.width(), wt = box.t.width();
    const double wmax = std::max(std::max(wu, wv), std::max(ws, wt));
    if (wmax <= tau || depth >= max_depth) { ++out.cover; out.boxes.push_back(box); return; }
    ++out.sub;
    if (wu >= wv && wu >= ws && wu >= wt) {
        double mid = 0.5 * (box.u.lo + box.u.hi);
        RBezierSurface L, R; split_u(A, L, R);
        Box4 b1 = box, b2 = box;
        b1.u = Interval(box.u.lo, mid); b2.u = Interval(mid, box.u.hi);
        rec(L, B, b1, depth + 1, tau, max_depth, node_cap, out);
        rec(R, B, b2, depth + 1, tau, max_depth, node_cap, out);
    } else if (wv >= ws && wv >= wt) {
        double mid = 0.5 * (box.v.lo + box.v.hi);
        RBezierSurface L, R; split_v(A, L, R);
        Box4 b1 = box, b2 = box;
        b1.v = Interval(box.v.lo, mid); b2.v = Interval(mid, box.v.hi);
        rec(L, B, b1, depth + 1, tau, max_depth, node_cap, out);
        rec(R, B, b2, depth + 1, tau, max_depth, node_cap, out);
    } else if (ws >= wt) {
        double mid = 0.5 * (box.s.lo + box.s.hi);
        RBezierSurface L, R; split_u(B, L, R);
        Box4 b1 = box, b2 = box;
        b1.s = Interval(box.s.lo, mid); b2.s = Interval(mid, box.s.hi);
        rec(A, L, b1, depth + 1, tau, max_depth, node_cap, out);
        rec(A, R, b2, depth + 1, tau, max_depth, node_cap, out);
    } else {
        double mid = 0.5 * (box.t.lo + box.t.hi);
        RBezierSurface L, R; split_v(B, L, R);
        Box4 b1 = box, b2 = box;
        b1.t = Interval(box.t.lo, mid); b2.t = Interval(mid, box.t.hi);
        rec(A, L, b1, depth + 1, tau, max_depth, node_cap, out);
        rec(A, R, b2, depth + 1, tau, max_depth, node_cap, out);
    }
}
} // namespace detail

inline CoverResult cover(const RBezierSurface& A, const RBezierSurface& B,
                         double tau, int max_depth = 60, long node_cap = 4000000,
                         bool record_excluded = false) {
    CoverResult out;
    out.record_excluded = record_excluded;
    Box4 root{ Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0) };
    detail::rec(A, B, root, 0, tau, max_depth, node_cap, out);
    return out;
}

inline bool covered(const std::vector<Box4>& boxes, double u, double v, double s, double t, double tol) {
    for (const Box4& b : boxes)
        if (u >= b.u.lo - tol && u <= b.u.hi + tol &&
            v >= b.v.lo - tol && v <= b.v.hi + tol &&
            s >= b.s.lo - tol && s <= b.s.hi + tol &&
            t >= b.t.lo - tol && t <= b.t.hi + tol)
            return true;
    return false;
}

// ============================================================================
// M4：横截认证
// ============================================================================

// 参数相关的线性插值包络：真值 (1-t)A + tB，t∈T 时端点落在
// (1-T.lo)A + T.lo B 与 (1-T.hi)A + T.hi B 之间（端点取值的包络）。
// 比"独立区间" imul(A,1-T)+imul(B,T) 紧得多（常数控制点保持常数、权重不膨胀）。
inline Interval lerp_env(const Interval& A, const Interval& B, const Interval& T) {
    Interval e1 = iadd(imul(A, Interval::point(1.0 - T.lo)), imul(B, Interval::point(T.lo)));
    Interval e2 = iadd(imul(A, Interval::point(1.0 - T.hi)), imul(B, Interval::point(T.hi)));
    return Interval(std::min(e1.lo, e2.lo), std::max(e1.hi, e2.hi));
}

// 通用齐次 de Casteljau 求值（区间，参数相关包络）
inline void eval_net(const std::vector<HP>& P, int nu, int nv,
                     const Interval& U, const Interval& V, HP& out) {
    std::vector<HP> col(nu);
    for (int i = 0; i < nu; ++i) {
        std::vector<HP> row(nv);
        for (int j = 0; j < nv; ++j) row[j] = P[(std::size_t)i * nv + j];
        for (int r = 1; r < nv; ++r)
            for (int j = 0; j + r < nv; ++j)
                for (int c = 0; c < 4; ++c)
                    row[j][c] = lerp_env(row[j][c], row[j + 1][c], V);
        col[i] = row[0];
    }
    for (int r = 1; r < nu; ++r)
        for (int i = 0; i + r < nu; ++i)
            for (int c = 0; c < 4; ++c)
                col[i][c] = lerp_env(col[i][c], col[i + 1][c], U);
    out = col[0];
}

// 有理点 + 两个偏导（区间）
// ============================================================================
// P6：裁剪域（参数多边形 D ⊂ [0,1]²）
//   把覆盖盒按 (u,v) 投影与 D 的关系分类为 inside / outside / partial：
//   - outside 仅当"盒完全在 D 外"被**证明**（所有盒角在 D 外、所有 D 顶点在盒外、
//     且没有 D 的边穿过盒）时才排除——保证无遗漏；
//   - partial / inside 的盒保留（partial 记为 undecided）。
//   对裁剪域 D 交集的结论：Σ|_D ⊆ 保留盒，且每个 outside 盒都有参数域排除证书。
// ============================================================================
struct Poly2 {
    std::vector<std::array<double, 2>> V;
};

inline bool point_in_poly(const Poly2& P, double x, double y) {
    bool in = false;
    std::size_t n = P.V.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        double xi = P.V[i][0], yi = P.V[i][1];
        double xj = P.V[j][0], yj = P.V[j][1];
        if (((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) in = !in;
    }
    return in;
}

// 线段与轴对齐盒相交（Liang-Barsky，保守：返回 true 表示可能相交）
inline bool seg_aabb_overlap(double ax, double ay, double bx, double by,
                             const Interval& u, const Interval& v) {
    double dx = bx - ax, dy = by - ay;
    double t0 = 0.0, t1 = 1.0;
    double p[4] = { -dx, dx, -dy, dy };
    double q[4] = { ax - u.lo, u.hi - ax, ay - v.lo, v.hi - ay };
    for (int i = 0; i < 4; ++i) {
        if (p[i] == 0.0) { if (q[i] < 0.0) return false; }
        else {
            double t = q[i] / p[i];
            if (p[i] < 0.0) { if (t > t1) return false; if (t > t0) t0 = t; }
            else { if (t < t0) return false; if (t < t1) t1 = t; }
        }
    }
    return t0 <= t1;
}

enum TrimState { TRIM_OUT = 0, TRIM_IN = 1, TRIM_PARTIAL = 2 };

inline TrimState poly_box_state(const Poly2& P, const Interval& u, const Interval& v) {
    // 1) 任一 D 顶点在盒内 ⇒ 不能证明 outside
    bool any_vertex_in_box = false;
    for (const auto& p : P.V)
        if (p[0] >= u.lo && p[0] <= u.hi && p[1] >= v.lo && p[1] <= v.hi) { any_vertex_in_box = true; break; }
    // 2) 任一盒角在 D 内 ⇒ 不能证明 outside
    bool any_corner_in_poly = false;
    const double cs[4][2] = { { u.lo, v.lo }, { u.hi, v.lo }, { u.lo, v.hi }, { u.hi, v.hi } };
    for (const auto& c : cs)
        if (point_in_poly(P, c[0], c[1])) { any_corner_in_poly = true; break; }
    if (!any_vertex_in_box && !any_corner_in_poly) {
        // 3) 任一边穿过盒 ⇒ 不能证明 outside
        std::size_t n = P.V.size();
        for (std::size_t i = 0, j = n - 1; i < n; j = i++)
            if (seg_aabb_overlap(P.V[i][0], P.V[i][1], P.V[j][0], P.V[j][1], u, v)) return TRIM_PARTIAL;
        return TRIM_OUT;   // 已证明：盒完全在 D 外
    }
    // 盒完全在 D 内？ 4 角在 D 内且 D 顶点不在盒内且无边穿盒
    bool all_corners_in = true;
    for (const auto& c : cs) if (!point_in_poly(P, c[0], c[1])) { all_corners_in = false; break; }
    if (all_corners_in && !any_vertex_in_box) {
        std::size_t n = P.V.size();
        bool edge_cross = false;
        for (std::size_t i = 0, j = n - 1; i < n; j = i++)
            if (seg_aabb_overlap(P.V[i][0], P.V[i][1], P.V[j][0], P.V[j][1], u, v)) { edge_cross = true; break; }
        if (!edge_cross) return TRIM_IN;
    }
    return TRIM_PARTIAL;
}

struct TrimCoverStats {
    long cover = 0;        // 原覆盖盒
    long kept = 0;         // 裁剪后保留
    long outside = 0;      // 被正确排除（参数域证书）
    long partial = 0;      // undecided
};

// 裁剪后的覆盖：对 (u,v) 用 polyA 分类、对 (s,t) 用 polyB 分类
inline TrimCoverStats trimmed_cover(const RBezierSurface& A, const RBezierSurface& B,
                                    double tau, const Poly2& polyA, const Poly2& polyB) {
    TrimCoverStats st;
    auto r = cover(A, B, tau);
    st.cover = (long)r.boxes.size();
    for (const Box4& X : r.boxes) {
        TrimState a = poly_box_state(polyA, X.u, X.v);
        TrimState b = poly_box_state(polyB, X.s, X.t);
        if (a == TRIM_OUT || b == TRIM_OUT) { ++st.outside; continue; }
        ++st.kept;
        if (a == TRIM_PARTIAL || b == TRIM_PARTIAL) ++st.partial;
    }
    return st;
}

inline bool eval_rational_duv(const RBezierSurface& s, const Interval& U, const Interval& V,
                              Interval S[3], Interval Su[3], Interval Sv[3]) {
    HP h; eval_net(s.P, s.nu, s.nv, U, V, h);
    // 严格性关键：区间 de Casteljau 的权重包络可能被放大（含 0）；
    // 真权重 ∈ 控制权重凸包 ⊆ [wmin, wmax]，与已知范围求交得到更紧且有效的包络。
    Interval w = Interval(std::max(h[3].lo, s.wmin), std::min(h[3].hi, s.wmax));
    if (w.isEmpty() || w.lo <= 0.0) return false;
    for (int c = 0; c < 3; ++c) {
        S[c] = idiv(h[c], w);
        if (S[c].isEmpty()) return false;
    }
    // u 方向导数网
    int du = s.nu - 1;
    if (du >= 1) {
        std::vector<HP> D((std::size_t)du * s.nv);
        Interval k((double)du, (double)du);
        for (int i = 0; i < du; ++i)
            for (int j = 0; j < s.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    D[(std::size_t)i * s.nv + j][c] = imul(k, isub(s.at(i + 1, j)[c], s.at(i, j)[c]));
        HP hu; eval_net(D, du, s.nv, U, V, hu);
        for (int c = 0; c < 3; ++c) {
            Interval num = isub(hu[c], imul(S[c], hu[3]));
            Su[c] = idiv(num, w);
            if (Su[c].isEmpty()) return false;
        }
    } else {
        for (int c = 0; c < 3; ++c) Su[c] = Interval::point(0.0);
    }
    // v 方向导数网
    int dv = s.nv - 1;
    if (dv >= 1) {
        std::vector<HP> D((std::size_t)s.nu * dv);
        Interval k((double)dv, (double)dv);
        for (int i = 0; i < s.nu; ++i)
            for (int j = 0; j < dv; ++j)
                for (int c = 0; c < 4; ++c)
                    D[(std::size_t)i * dv + j][c] = imul(k, isub(s.at(i, j + 1)[c], s.at(i, j)[c]));
        HP hv; eval_net(D, s.nu, dv, U, V, hv);
        for (int c = 0; c < 3; ++c) {
            Interval num = isub(hv[c], imul(S[c], hv[3]));
            Sv[c] = idiv(num, w);
            if (Sv[c].isEmpty()) return false;
        }
    } else {
        for (int c = 0; c < 3; ++c) Sv[c] = Interval::point(0.0);
    }
    return true;
}

// Gershgorin 秩判定：M = J Jᵀ（3x3 区间），若所有 Gershgorin 下界 > 0 ⇒ rank(J)=3
inline bool rank3_certified(const Interval J[3][4], double* margin = nullptr) {
    Interval G[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            Interval acc = Interval::point(0.0);
            for (int k = 0; k < 4; ++k) acc = iadd(acc, imul(J[i][k], J[j][k]));
            G[i][j] = acc;
        }
    double m = 1e300;
    for (int i = 0; i < 3; ++i) {
        double rad = 0.0;
        for (int j = 0; j < 3; ++j) if (j != i)
            rad += std::max(std::fabs(G[i][j].lo), std::fabs(G[i][j].hi));
        m = std::min(m, G[i][i].lo - rad);
    }
    if (margin) *margin = m;
    return m > 0.0;
}

// ---------- 3x3 切片认证（沿零空间方向固定一个参数）----------
// 欠定系统 3x4 的 Krawczyk 在零空间方向无法收缩（(I-YJ) 在零方向特征值为 1），
// 因此改为：固定零空间主导坐标 k*，在 3 维切片上做方系统 Krawczyk：
//   K(X3) ⊂ int(X3) ⇒ 该切片内存在唯一解点。

inline bool invert3(const double A[3][3], double Ainv[3][3]) {
    double M[3][6];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) M[i][j] = A[i][j];
        for (int j = 0; j < 3; ++j) M[i][3 + j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < 3; ++col) {
        int piv = col; double best = std::fabs(M[col][col]);
        for (int r = col + 1; r < 3; ++r) if (std::fabs(M[r][col]) > best) { best = std::fabs(M[r][col]); piv = r; }
        if (best < 1e-14) return false;
        if (piv != col) for (int j = 0; j < 6; ++j) std::swap(M[piv][j], M[col][j]);
        double d = M[col][col];
        for (int j = 0; j < 6; ++j) M[col][j] /= d;
        for (int r = 0; r < 3; ++r) if (r != col) {
            double f = M[r][col];
            for (int j = 0; j < 6; ++j) M[r][j] -= f * M[col][j];
        }
    }
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Ainv[i][j] = M[i][3 + j];
    return true;
}

// 零空间主导坐标：y = (I - Y J)；幂迭代求其主特征向量，取分量最大的坐标
inline int null_coord(const Interval J[3][4]) {
    double Jm[3][4];
    for (int c = 0; c < 3; ++c) for (int d = 0; d < 4; ++d) Jm[c][d] = J[c][d].mid();
    double A[3][3], Ainv[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double acc = 0.0;
            for (int k = 0; k < 4; ++k) acc += Jm[i][k] * Jm[j][k];
            A[i][j] = acc;
        }
    if (!invert3(A, Ainv)) return -1;
    double Y[4][3];
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 3; ++i) {
            double acc = 0.0;
            for (int j = 0; j < 3; ++j) acc += Jm[j][k] * Ainv[j][i];
            Y[k][i] = acc;
        }
    double N[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            double acc = 0.0;
            for (int i = 0; i < 3; ++i) acc += Y[a][i] * Jm[i][b];
            N[a][b] = ((a == b) ? 1.0 : 0.0) - acc;
        }
    double v[4] = { 0.3, 0.5, 0.7, 0.9 };
    for (int it = 0; it < 64; ++it) {
        double w[4] = { 0, 0, 0, 0 };
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) w[a] += N[a][b] * v[b];
        double nrm = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        if (nrm < 1e-14) return -1;
        for (int a = 0; a < 4; ++a) v[a] = w[a] / nrm;
    }
    int k = 0; double best = std::fabs(v[0]);
    for (int a = 1; a < 4; ++a) if (std::fabs(v[a]) > best) { best = std::fabs(v[a]); k = a; }
    return k;
}

// 切片方系统 Krawczyk（固定 k*，其余三参数），带 Krawczyk 迭代收窄：
// 每轮重算区间 J（在收窄后的盒上），K ⊂ int(X3) 即认证；否则 X3 ← K ∩ X3 继续。
// 返回码：0=认证 1=零空间坐标失败 2=求逆失败 3=收缩/区间秩失败 4=不包含且无法收窄
inline int krawczyk_slice_code(const RBezierSurface& A, const RBezierSurface& B,
                               const Box4& X0, int force_k = -1) {
    Box4 X = X0;
    const int max_iter = 12;
    int kfix = -1;
    for (int iter = 0; iter < max_iter; ++iter) {
        Interval S1[3], S1u[3], S1v[3], S2[3], S2u[3], S2v[3];
        if (!eval_rational_duv(A, X.u, X.v, S1, S1u, S1v)) return 2;
        if (!eval_rational_duv(B, X.s, X.t, S2, S2u, S2v)) return 2;
        Interval J[3][4];
        for (int c = 0; c < 3; ++c) {
            J[c][0] = S1u[c]; J[c][1] = S1v[c];
            J[c][2] = Interval(-S2u[c].hi, -S2u[c].lo);
            J[c][3] = Interval(-S2v[c].hi, -S2v[c].lo);
        }
        if (kfix < 0) {
            kfix = (force_k >= 0) ? force_k : null_coord(J);
            if (kfix < 0) return 1;
        }
        int idx[3]; int t = 0;
        for (int c = 0; c < 4; ++c) if (c != kfix) idx[t++] = c;
        Interval X3[3]; double m3[3];
        for (int i = 0; i < 3; ++i) {
            const Interval* Xc = (idx[i] == 0) ? &X.u : (idx[i] == 1) ? &X.v : (idx[i] == 2) ? &X.s : &X.t;
            X3[i] = *Xc; m3[i] = Xc->mid();
        }
        Interval Fm[3];
        {
            Interval mm[4];
            double m[4] = { X.u.mid(), X.v.mid(), X.s.mid(), X.t.mid() };
            for (int i = 0; i < 4; ++i) mm[i] = Interval::point(m[i]);
            Interval a1[3], a2[3], a3[3], b1[3], b2[3], b3[3];
            if (!eval_rational_duv(A, mm[0], mm[1], a1, a2, a3)) return 2;
            if (!eval_rational_duv(B, mm[2], mm[3], b1, b2, b3)) return 2;
            for (int c = 0; c < 3; ++c) Fm[c] = isub(a1[c], b1[c]);
        }
        double J3m[3][3], J3inv[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) J3m[r][c] = J[r][idx[c]].mid();
        if (!invert3(J3m, J3inv)) return 2;
        Interval J3[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) J3[r][c] = J[r][idx[c]];

        // 收缩判据：W = I - J3inv * J3(X)，||W||_inf < 1
        Interval W[3][3];
        double wnorm = 0.0;
        for (int r = 0; r < 3; ++r) {
            for (int d = 0; d < 3; ++d) {
                Interval acc = (r == d) ? Interval::point(1.0) : Interval::point(0.0);
                for (int c = 0; c < 3; ++c)
                    acc = isub(acc, imul(Interval::point(J3inv[r][c]), J3[c][d]));
                W[r][d] = acc;
            }
            double row = 0.0;
            for (int d = 0; d < 3; ++d) row += std::max(std::fabs(W[r][d].lo), std::fabs(W[r][d].hi));
            wnorm = std::max(wnorm, row);
        }
        if (!(wnorm < 1.0)) return 3;

        // K = m3 - J3inv*Fm + W*(X3 - m3)
        Interval K[3];
        for (int r = 0; r < 3; ++r) {
            Interval kacc = Interval(m3[r], m3[r]);
            for (int c = 0; c < 3; ++c)
                kacc = isub(kacc, imul(Interval::point(J3inv[r][c]), Fm[c]));
            for (int d = 0; d < 3; ++d)
                kacc = iadd(kacc, imul(W[r][d], isub(X3[d], Interval::point(m3[d]))));
            K[r] = kacc;
        }
        bool inside = true, nonempty = true, narrower = true;
        for (int r = 0; r < 3; ++r) {
            if (!(X3[r].lo <= K[r].lo && K[r].hi <= X3[r].hi)) inside = false;
            if (!(K[r].width() <= 0.9 * X3[r].width() || X3[r].width() < 1e-12)) inside = false;
            double lo = std::max(X3[r].lo, K[r].lo), hi = std::min(X3[r].hi, K[r].hi);
            if (lo > hi) nonempty = false;
            else {
                if (!(lo > X3[r].lo || hi < X3[r].hi)) narrower = false;
                X3[r] = Interval(lo, hi);
            }
        }
        if (inside) return 0;
        if (!nonempty || !narrower) return 4;
        // 写回收窄后的切片盒
        for (int i = 0; i < 3; ++i) {
            Interval* Xc = (idx[i] == 0) ? &X.u : (idx[i] == 1) ? &X.v : (idx[i] == 2) ? &X.s : &X.t;
            *Xc = X3[i];
        }
    }
    return 4;
}

inline bool krawczyk_slice(const RBezierSurface& A, const RBezierSurface& B, const Box4& X) {
    if (krawczyk_slice_code(A, B, X) == 0) return true;
    for (int k = 0; k < 4; ++k) if (krawczyk_slice_code(A, B, X, k) == 0) return true;
    return false;
}

// 探针：以给定解点为中心取半宽 hv 的盒，返回认证码（用于单点诊断）
inline int probe_box(const RBezierSurface& A, const RBezierSurface& B,
                     double u, double v, double s, double t, double hv) {
    Box4 X;
    X.u = Interval(u - hv, u + hv); X.v = Interval(v - hv, v + hv);
    X.s = Interval(s - hv, s + hv); X.t = Interval(t - hv, t + hv);
    return krawczyk_slice_code(A, B, X);
}

// 前向声明（M6.5 定义）：3x4 -> 4x3 广义逆
inline bool pinv_y(const double J[3][4], double Y[4][3]);

// ============================================================================
// 4D Krawczyk 认证（欠定系统，J 满秩 3）：
//   若 J(X) 区间满秩（rank3_certified）且 K(X) ⊆ X（闭包含），
//   则 g(x)=x-YF(x) 把 X 映入 X，Brouwer 给不动点 x*，YF(x*)=0；
//   Y 满列秩(3) ⇒ YF=0 ⇒ F=0，故 x*∈Σ。
//   闭包含即可（允许沿零方向的等号），比切片法更强。
// ============================================================================
inline int krawczyk_4d_code(const RBezierSurface& A, const RBezierSurface& B,
                            const Box4& X, bool allow_inflate, Box4* outBox,
                            double* outQ) {
    Interval S1[3], S1u[3], S1v[3], S2[3], S2u[3], S2v[3];
    if (!eval_rational_duv(A, X.u, X.v, S1, S1u, S1v)) return 2;
    if (!eval_rational_duv(B, X.s, X.t, S2, S2u, S2v)) return 2;
    Interval J[3][4];
    for (int c = 0; c < 3; ++c) {
        J[c][0] = S1u[c]; J[c][1] = S1v[c];
        J[c][2] = Interval(-S2u[c].hi, -S2u[c].lo);
        J[c][3] = Interval(-S2v[c].hi, -S2v[c].lo);
    }
    if (!rank3_certified(J)) return 3;

    double m[4] = { X.u.mid(), X.v.mid(), X.s.mid(), X.t.mid() };
    Interval Fm[3];
    {
        Interval mm[4];
        for (int i = 0; i < 4; ++i) mm[i] = Interval::point(m[i]);
        Interval a1[3], a2[3], a3[3], b1[3], b2[3], b3[3];
        if (!eval_rational_duv(A, mm[0], mm[1], a1, a2, a3)) return 2;
        if (!eval_rational_duv(B, mm[2], mm[3], b1, b2, b3)) return 2;
        for (int c = 0; c < 3; ++c) Fm[c] = isub(a1[c], b1[c]);
    }
    double Jm[3][4];
    for (int c = 0; c < 3; ++c)
        for (int d = 0; d < 4; ++d) Jm[c][d] = J[c][d].mid();
    double Y[4][3];
    if (!pinv_y(Jm, Y)) return 2;

    // K = m - Y F(m) + (I - Y J(X))(X - m)
    Interval W[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Interval acc = (a == b) ? Interval::point(1.0) : Interval::point(0.0);
            for (int c = 0; c < 3; ++c)
                acc = isub(acc, imul(Interval::point(Y[a][c]), J[c][b]));
            W[a][b] = acc;
        }
    // 压缩因子 q = ||W||_inf（定理 2c 的唯一性条件）
    double q = 0.0;
    for (int a = 0; a < 4; ++a) {
        double row = 0.0;
        for (int b = 0; b < 4; ++b)
            row += std::max(std::fabs(W[a][b].lo), std::fabs(W[a][b].hi));
        q = std::max(q, row);
    }
    if (outQ) *outQ = q;

    Interval Xc[4] = { X.u, X.v, X.s, X.t };
    Interval K[4];
    for (int a = 0; a < 4; ++a) {
        Interval acc(m[a], m[a]);
        for (int c = 0; c < 3; ++c)
            acc = isub(acc, imul(Interval::point(Y[a][c]), Fm[c]));
        for (int b = 0; b < 4; ++b)
            acc = iadd(acc, imul(W[a][b], isub(Xc[b], Interval::point(m[b]))));
        K[a] = acc;
    }
    auto contained = [&](double inflate) {
        for (int a = 0; a < 4; ++a) {
            if (std::isnan(K[a].lo) || std::isnan(K[a].hi)) return false;   // NaN 守卫
            if (!std::isfinite(K[a].lo) && !std::isinf(K[a].lo)) return false;
            if (!std::isfinite(K[a].hi) && !std::isinf(K[a].hi)) return false;
            double e = inflate * Xc[a].width();
            double lo = std::max(0.0, Xc[a].lo - e), hi = std::min(1.0, Xc[a].hi + e);
            if (K[a].lo < lo - 1e-12 || K[a].hi > hi + 1e-12) return false;
        }
        return true;
    };
    if (contained(0.0)) {
        if (outBox) *outBox = X;
        return 0;
    }
    if (allow_inflate && contained(0.25)) {
        if (outBox) {
            Box4 Yb;
            Yb.u = Interval(std::max(0.0, X.u.lo - 0.25 * X.u.width()), std::min(1.0, X.u.hi + 0.25 * X.u.width()));
            Yb.v = Interval(std::max(0.0, X.v.lo - 0.25 * X.v.width()), std::min(1.0, X.v.hi + 0.25 * X.v.width()));
            Yb.s = Interval(std::max(0.0, X.s.lo - 0.25 * X.s.width()), std::min(1.0, X.s.hi + 0.25 * X.s.width()));
            Yb.t = Interval(std::max(0.0, X.t.lo - 0.25 * X.t.width()), std::min(1.0, X.t.hi + 0.25 * X.t.width()));
            *outBox = Yb;
        }
        return 0;
    }
    return 4;
}

struct CertStats {
    long boxes = 0, rank3 = 0, certified = 0, unique = 0, tangent = 0, undecided = 0;
    long fail[5] = { 0, 0, 0, 0, 0 };
    std::vector<char> box_certified;
};

// 深认证：直接切片失败时，把盒 16 等分（4 维各取中位）再逐个尝试
inline bool certify_box_deep(const RBezierSurface& A, const RBezierSurface& B, const Box4& X) {
    if (krawczyk_slice(A, B, X)) return true;
    const double u0 = X.u.mid(), v0 = X.v.mid(), s0 = X.s.mid(), t0 = X.t.mid();
    for (int iu = 0; iu < 2; ++iu)
        for (int iv = 0; iv < 2; ++iv)
            for (int is = 0; is < 2; ++is)
                for (int it = 0; it < 2; ++it) {
                    Box4 Y;
                    Y.u = iu ? Interval(u0, X.u.hi) : Interval(X.u.lo, u0);
                    Y.v = iv ? Interval(v0, X.v.hi) : Interval(X.v.lo, v0);
                    Y.s = is ? Interval(s0, X.s.hi) : Interval(X.s.lo, s0);
                    Y.t = it ? Interval(t0, X.t.hi) : Interval(X.t.lo, t0);
                    if (krawczyk_slice(A, B, Y)) return true;
                }
    return false;
}

// 前向声明（M5 相切认证定义在文件末尾）
inline int certify_tangency_box(const RBezierSurface& A, const RBezierSurface& B, const Box4& X);

// 对覆盖盒逐一认证（横截：零方向切片 + 3x3 Krawczyk，失败则 16 等分重试）
inline CertStats certify_boxes(const RBezierSurface& A, const RBezierSurface& B,
                               const std::vector<Box4>& boxes) {
    CertStats st;
    st.boxes = (long)boxes.size();
    st.box_certified.assign(boxes.size(), 0);
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        const Box4& X = boxes[i];
        Interval S1[3], S1u[3], S1v[3], S2[3], S2u[3], S2v[3];
        bool ok = eval_rational_duv(A, X.u, X.v, S1, S1u, S1v) &&
                  eval_rational_duv(B, X.s, X.t, S2, S2u, S2v);
        if (ok) {
            Interval J[3][4];
            for (int c = 0; c < 3; ++c) {
                J[c][0] = S1u[c]; J[c][1] = S1v[c];
                J[c][2] = Interval(-S2u[c].hi, -S2u[c].lo);
                J[c][3] = Interval(-S2v[c].hi, -S2v[c].lo);
            }
            if (rank3_certified(J)) ++st.rank3;
        }
        int code = ok ? krawczyk_slice_code(A, B, X) : 2;
        if (code != 0 && ok) {
            for (int k = 0; k < 4 && code != 0; ++k) code = krawczyk_slice_code(A, B, X, k);
        }
        bool certified = (code == 0);
        if (code == 0) ++st.unique;   // 定理 2c'：切片认证 ⇒ 该切片值下解唯一
        if (!certified && ok) certified = certify_box_deep(A, B, X);
        if (!certified && ok) {
            Box4 cb;
            double q = 1e300;
            if (krawczyk_4d_code(A, B, X, true, &cb, &q) == 0) {
                certified = true;
                if (q < 1.0) ++st.unique;   // 定理 2c：盒内解唯一
            }
        }
        if (certified) { ++st.certified; st.box_certified[i] = 1; }
        else {
            // 横截认证失败 → 相切/秩亏候选：接触函数认证
            int tcode = ok ? certify_tangency_box(A, B, X) : 2;
            if (tcode == 0) { ++st.tangent; st.box_certified[i] = 1; }
            else { ++st.undecided; if (code >= 1 && code <= 4) ++st.fail[code]; }
        }
    }
    return st;
}
// ============================================================================
// M5：相切认证 —— 接触函数 h = n·(S - P0) 的临界点 ∂h/∂u = ∂h/∂v = 0
//   点相切：2x2 区间 Newton（Hessian 非奇异）
//   线相切：固定一个参数做 1D 区间 Newton（沿横截方向）
//   面片伙伴为平面时（H 恒定），退化为"高度为零"检查
// ============================================================================

// 齐次导数网：沿 u（nu-1 次）
inline std::vector<HP> deriv_net_u(const std::vector<HP>& P, int nu, int nv) {
    if (nu < 2) return {};
    std::vector<HP> D((std::size_t)(nu - 1) * nv);
    const Interval k((double)(nu - 1), (double)(nu - 1));
    for (int i = 0; i < nu - 1; ++i)
        for (int j = 0; j < nv; ++j)
            for (int c = 0; c < 4; ++c)
                D[(std::size_t)i * nv + j][c] = imul(k, isub(P[(std::size_t)(i + 1) * nv + j][c], P[(std::size_t)i * nv + j][c]));
    return D;
}

// 齐次导数网：沿 v
inline std::vector<HP> deriv_net_v(const std::vector<HP>& P, int nu, int nv) {
    if (nv < 2) return {};
    std::vector<HP> D((std::size_t)nu * (nv - 1));
    const Interval k((double)(nv - 1), (double)(nv - 1));
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv - 1; ++j)
            for (int c = 0; c < 4; ++c)
                D[(std::size_t)i * (nv - 1) + j][c] = imul(k, isub(P[(std::size_t)i * nv + j + 1][c], P[(std::size_t)i * nv + j][c]));
    return D;
}

// 有理点 + 一阶/二阶偏导（区间）
inline bool eval_rational_d2(const RBezierSurface& s, const Interval& U, const Interval& V,
                             Interval S[3], Interval Su[3], Interval Sv[3],
                             Interval Suu[3], Interval Suv[3], Interval Svv[3]) {
    HP h; eval_net(s.P, s.nu, s.nv, U, V, h);
    // 严格性关键：区间 de Casteljau 的权重包络可能被放大（含 0）；
    // 真权重 ∈ 控制权重凸包 ⊆ [wmin, wmax]，故与已知范围求交得到更紧且有效的包络。
    Interval w = Interval(std::max(h[3].lo, s.wmin), std::min(h[3].hi, s.wmax));
    if (w.isEmpty() || w.lo <= 0.0) return false;
    for (int c = 0; c < 3; ++c) { S[c] = idiv(h[c], w); if (S[c].isEmpty()) return false; }

    std::vector<HP> Du, Dv;
    if (s.nu >= 2) Du = deriv_net_u(s.P, s.nu, s.nv);
    if (s.nv >= 2) Dv = deriv_net_v(s.P, s.nu, s.nv);
    HP hu, hv;
    if (s.nu >= 2) eval_net(Du, s.nu - 1, s.nv, U, V, hu); else for (int c = 0; c < 4; ++c) hu[c] = Interval::point(0.0);
    if (s.nv >= 2) eval_net(Dv, s.nu, s.nv - 1, U, V, hv); else for (int c = 0; c < 4; ++c) hv[c] = Interval::point(0.0);
    for (int c = 0; c < 3; ++c) {
        Su[c] = idiv(isub(hu[c], imul(S[c], hu[3])), w);
        Sv[c] = idiv(isub(hv[c], imul(S[c], hv[3])), w);
        if (Su[c].isEmpty() || Sv[c].isEmpty()) return false;
    }
    // 二阶
    HP huu, huv, hvv;
    for (int c = 0; c < 4; ++c) { huu[c] = huv[c] = hvv[c] = Interval::point(0.0); }
    if (s.nu >= 3) {
        std::vector<HP> Duu = deriv_net_u(Du, s.nu - 1, s.nv);
        eval_net(Duu, s.nu - 2, s.nv, U, V, huu);
    }
    if (s.nu >= 2 && s.nv >= 2) {
        std::vector<HP> Duv = deriv_net_v(Du, s.nu - 1, s.nv);
        eval_net(Duv, s.nu - 1, s.nv - 1, U, V, huv);
    }
    if (s.nv >= 3) {
        std::vector<HP> Dvv = deriv_net_v(Dv, s.nu, s.nv - 1);
        eval_net(Dvv, s.nu, s.nv - 2, U, V, hvv);
    }
    for (int c = 0; c < 3; ++c) {
        // Suu = (Huu - 2 Su Wu - S Wuu)/W
        Suu[c] = idiv(isub(isub(huu[c], imul(Interval::point(2.0), imul(Su[c], hu[3]))), imul(S[c], huu[3])), w);
        Suv[c] = idiv(isub(isub(isub(huv[c], imul(Su[c], hv[3])), imul(Sv[c], hu[3])), imul(S[c], huv[3])), w);
        Svv[c] = idiv(isub(isub(hvv[c], imul(Interval::point(2.0), imul(Sv[c], hv[3]))), imul(S[c], hvv[3])), w);
        if (Suu[c].isEmpty() || Suv[c].isEmpty() || Svv[c].isEmpty()) return false;
    }
    return true;
}

// 2x2 区间 Newton（Krawczyk）：F = (Hu, Hv)，Jacobian = Hessian
template <class Eval>
inline bool certify_crit2d(const Interval& U0, const Interval& V0, Eval& eval,
                           Interval* outU, Interval* outV, Interval* outH) {
    Interval U = U0, V = V0;
    for (int it = 0; it < 12; ++it) {
        Interval H, Hu, Hv, Huu, Huv, Hvv;
        if (!eval(U, V, H, Hu, Hv, Huu, Huv, Hvv)) return false;
        double mu = U.mid(), mv = V.mid();
        Interval Hm, Hum, Hvm, Huum, Huvm, Hvvm;
        if (!eval(Interval::point(mu), Interval::point(mv), Hm, Hum, Hvm, Huum, Huvm, Hvvm)) return false;
        double J00 = Huum.mid(), J01 = Huvm.mid(), J11 = Hvvm.mid();
        double det = J00 * J11 - J01 * J01;
        if (std::fabs(det) < 1e-13) return false;
        double i00 = J11 / det, i01 = -J01 / det, i10 = -J01 / det, i11 = J00 / det;
        // W = I - Jinv * Hess(X)
        Interval W00 = isub(Interval::point(1.0), iadd(imul(Interval::point(i00), Huu), imul(Interval::point(i01), Huv)));
        Interval W01 = isub(Interval::point(0.0), iadd(imul(Interval::point(i00), Huv), imul(Interval::point(i01), Hvv)));
        Interval W10 = isub(Interval::point(0.0), iadd(imul(Interval::point(i10), Huu), imul(Interval::point(i11), Huv)));
        Interval W11 = isub(Interval::point(1.0), iadd(imul(Interval::point(i10), Huv), imul(Interval::point(i11), Hvv)));
        // K = m - Jinv*F(m) + W*(X - m)
        Interval Ku = Interval(mu, mu), Kv = Interval(mv, mv);
        Ku = isub(Ku, iadd(imul(Interval::point(i00), Hum), imul(Interval::point(i01), Hvm)));
        Kv = isub(Kv, iadd(imul(Interval::point(i10), Hum), imul(Interval::point(i11), Hvm)));
        Interval du = isub(U, Interval::point(mu)), dv = isub(V, Interval::point(mv));
        Ku = iadd(Ku, iadd(imul(W00, du), imul(W01, dv)));
        Kv = iadd(Kv, iadd(imul(W10, du), imul(W11, dv)));
        const bool closed2 = (U.lo <= Ku.lo && Ku.hi <= U.hi && V.lo <= Kv.lo && Kv.hi <= V.hi);
        const bool shrink2 = (Ku.width() <= 0.9 * U.width() || U.width() < 1e-12) &&
                             (Kv.width() <= 0.9 * V.width() || V.width() < 1e-12);
        const bool inside = closed2 && shrink2;
        Interval nU(std::max(U.lo, Ku.lo), std::min(U.hi, Ku.hi));
        Interval nV(std::max(V.lo, Kv.lo), std::min(V.hi, Kv.hi));
        if (inside) {
            if (outU) *outU = nU;
            if (outV) *outV = nV;
            if (outH) {
                Interval Hc, Huc, Hvc, Huuc, Huvc, Hvvc;
                if (!eval(nU, nV, Hc, Huc, Hvc, Huuc, Huvc, Hvvc)) return false;
                *outH = Hc;
            }
            return true;
        }
        if (nU.isEmpty() || nV.isEmpty()) return false;
        if (!(nU.lo > U.lo || nU.hi < U.hi || nV.lo > V.lo || nV.hi < V.hi)) return false;
        U = nU; V = nV;
    }
    return false;
}

// 1D 区间 Newton：g(u) = Hu(u, v0)，g'(u) = Huu(u, v0)
template <class Eval>
inline bool certify_crit1d(const Interval& U0, const Interval& Vfixed, Eval& eval,
                           Interval* outU, Interval* outH, Interval* outHv) {
    Interval U = U0;
    for (int it = 0; it < 16; ++it) {
        Interval H, Hu, Hv, Huu, Huv, Hvv;
        if (!eval(U, Vfixed, H, Hu, Hv, Huu, Huv, Hvv)) return false;
        double mu = U.mid();
        Interval Hm, Hum, Hvm, Huum, Huvm, Hvvm;
        if (!eval(Interval::point(mu), Vfixed, Hm, Hum, Hvm, Huum, Huvm, Hvvm)) return false;
        if (!(Huum.lo > 0.0 || Huum.hi < 0.0)) return false;   // 导数符号不定
        Interval N = isub(Interval::point(mu), idiv(Hum, Huum));
        if (N.isEmpty()) return false;
        Interval nU(std::max(U.lo, N.lo), std::min(U.hi, N.hi));
        if (nU.isEmpty()) return false;
        const bool closed1 = (U.lo <= N.lo && N.hi <= U.hi);
        const bool shrink1 = (N.width() <= 0.9 * U.width() || U.width() < 1e-12);
        const bool inside = closed1 && shrink1;
        if (inside) {
            Interval Hc, Huc, Hvc, Huuc, Huvc, Hvvc;
            if (!eval(nU, Vfixed, Hc, Huc, Hvc, Huuc, Huvc, Hvvc)) return false;
            if (outU) *outU = nU;
            if (outH) *outH = Hc;
            if (outHv) *outHv = Hvc;
            return true;
        }
        if (!(nU.lo > U.lo || nU.hi < U.hi)) return false;
        U = nU;
    }
    return false;
}

inline void mid3(const Interval S[3], double out[3]) { for (int i = 0; i < 3; ++i) out[i] = S[i].mid(); }
inline void cross3_mid(const Interval a[3], const Interval b[3], double out[3]) {
    double x[3], y[3]; mid3(a, x); mid3(b, y);
    out[0] = x[1] * y[2] - x[2] * y[1];
    out[1] = x[2] * y[0] - x[0] * y[2];
    out[2] = x[0] * y[1] - x[1] * y[0];
}
inline double dot3_mid(const Interval a[3], const Interval b[3]) {
    double x[3], y[3]; mid3(a, x); mid3(b, y);
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

// 前向声明（M7 定义在文件末尾）：点包含认证
inline bool point_on_patch(const RBezierSurface& S, const double P[3], double s0, double t0,
                           Interval* outS, Interval* outT);
inline bool point_on_patch_flat(const RBezierSurface& S, const double P[3],
                                const Interval& Ss, const Interval& St,
                                Interval* outS, Interval* outT);
inline bool point_on_patch_anywhere(const RBezierSurface& S, const double P[3],
                                    Interval* outS, Interval* outT);

// 相切认证：返回 0=认证成功；1=法向不平行（非相切候选）；2=求值失败；4=高度不符
inline int certify_tangency_box_impl(const RBezierSurface& A, const RBezierSurface& B, const Box4& X, bool allow_inflate) {
    Interval S1[3], S1u[3], S1v[3], S1uu[3], S1uv[3], S1vv[3];
    Interval S2[3], S2u[3], S2v[3], S2uu[3], S2uv[3], S2vv[3];
    if (!eval_rational_d2(A, X.u, X.v, S1, S1u, S1v, S1uu, S1uv, S1vv)) return 2;
    if (!eval_rational_d2(B, X.s, X.t, S2, S2u, S2v, S2uu, S2uv, S2vv)) return 2;

    double n1[3], n2[3];
    cross3_mid(S1u, S1v, n1);
    cross3_mid(S2u, S2v, n2);
    double ln1 = std::sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
    double ln2 = std::sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
    if (ln1 < 1e-12 || ln2 < 1e-12) return 1;
    double cosang = std::fabs(n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (ln1 * ln2);
    if (cosang < 0.99) return 1;                       // 法向不平行 → 非相切候选

    double n[3] = { n1[0] / ln1, n1[1] / ln1, n1[2] / ln1 };
    if (n[0] * n2[0] + n[1] * n2[1] + n[2] * n2[2] < 0.0) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }
    double P0[3]; mid3(S1, P0);
    double c0 = n[0] * P0[0] + n[1] * P0[1] + n[2] * P0[2];

    auto height = [&](const Interval* s, const Interval* su, const Interval* sv,
                      const Interval* suu, const Interval* suv, const Interval* svv,
                      Interval& H, Interval& Hu, Interval& Hv,
                      Interval& Huu, Interval& Huv, Interval& Hvv) {
        auto d = [&](const Interval* v) -> Interval {
            return iadd(iadd(imul(Interval::point(n[0]), v[0]), imul(Interval::point(n[1]), v[1])),
                        imul(Interval::point(n[2]), v[2]));
        };
        H = isub(d(s), Interval::point(c0));
        Hu = d(su); Hv = d(sv); Huu = d(suu); Huv = d(suv); Hvv = d(svv);
    };
    auto evalH1 = [&](const Interval& U, const Interval& V, Interval& H, Interval& Hu, Interval& Hv,
                      Interval& Huu, Interval& Huv, Interval& Hvv) -> bool {
        Interval s[3], su[3], sv[3], suu[3], suv[3], svv[3];
        if (!eval_rational_d2(A, U, V, s, su, sv, suu, suv, svv)) return false;
        height(s, su, sv, suu, suv, svv, H, Hu, Hv, Huu, Huv, Hvv);
        return true;
    };
    auto evalH2 = [&](const Interval& U, const Interval& V, Interval& H, Interval& Hu, Interval& Hv,
                      Interval& Huu, Interval& Huv, Interval& Hvv) -> bool {
        Interval s[3], su[3], sv[3], suu[3], suv[3], svv[3];
        if (!eval_rational_d2(B, U, V, s, su, sv, suu, suv, svv)) return false;
        height(s, su, sv, suu, suv, svv, H, Hu, Hv, Huu, Huv, Hvv);
        return true;
    };

    // 伙伴面片是否"平面片"（高度 Hessian 恒为零 ⇒ 位于某平面内）
    bool flat2 = false;
    {
        Interval H, Hu, Hv, Huu, Huv, Hvv;
        if (evalH2(X.s, X.t, H, Hu, Hv, Huu, Huv, Hvv)) {
            double hh = std::max(std::max(std::fabs(Huu.mid()), std::fabs(Huv.mid())), std::fabs(Hvv.mid()));
            flat2 = (hh < 1e-9);
        }
    }
    if (flat2) {
        // 平面伙伴：用其自身法向作公共法向，高度函数更准确
        double q[3] = { n2[0] / ln2, n2[1] / ln2, n2[2] / ln2 };
        if (q[0] * n[0] + q[1] * n[1] + q[2] * n[2] < 0.0) { q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; }
        n[0] = q[0]; n[1] = q[1]; n[2] = q[2];
        c0 = n[0] * P0[0] + n[1] * P0[1] + n[2] * P0[2];
    }

    const double wmax = std::max(std::max(X.u.width(), X.v.width()), std::max(X.s.width(), X.t.width()));
    const double htol = std::max(1e-9, 2.0 * wmax);

    // ① 点相切：2x2 临界点认证
    Interval cU, cV, cH;
    if (certify_crit2d(X.u, X.v, evalH1, &cU, &cV, &cH)) {
        // 定理 2b 修正：认证临界点必须在 S2 的像上，否则"等高不重合"会误认证
        {
            Interval P1[3], d1[3], d2[3];
            if (!eval_rational_duv(A, cU, cV, P1, d1, d2)) return 4;
            double pv[3] = { P1[0].mid(), P1[1].mid(), P1[2].mid() };
            bool inc = flat2 ? point_on_patch_flat(B, pv, Interval(0.0, 1.0), Interval(0.0, 1.0), nullptr, nullptr)
                             : point_on_patch_anywhere(B, pv, nullptr, nullptr);
            if (!inc) return 4;
        }
        if (flat2) {
            if (cH.lo <= 0.0 && 0.0 <= cH.hi) return 0;              // 高度恰为零（在平面上）
            if (std::fabs(cH.mid()) <= htol) return 0;               // 容差内接触
            return 4;
        }
        Interval cS, cT, cH2;
        if (certify_crit2d(X.s, X.t, evalH2, &cS, &cT, &cH2)) {
            if (cH.lo <= cH2.hi && cH2.lo <= cH.hi) return 0;
            return 4;
        }
        return 4;
    }

    // ② 线/圆相切：1D 认证（两个方向都试：u 固定求 v，v 固定求 u）
    {
        // 先试 ∂h/∂v = 0（u 固定）——切触圆/水平切触线的典型情形
        auto evalSwapped = [&](const Interval& U, const Interval& V, Interval& H, Interval& Hu,
                               Interval& Hv, Interval& Huu, Interval& Huv, Interval& Hvv) -> bool {
            return evalH1(V, U, H, Hv, Hu, Hvv, Huv, Huu);   // 交换 u/v
        };
        Interval cV1, cH1v, cHu1;
        if (certify_crit1d(X.v, Interval::point(X.u.mid()), evalSwapped, &cV1, &cH1v, &cHu1)) {
            {
                Interval P1[3], d1[3], d2[3];
                if (eval_rational_duv(A, Interval::point(X.u.mid()), cV1, P1, d1, d2)) {
                    double pv[3] = { P1[0].mid(), P1[1].mid(), P1[2].mid() };
                    bool inc = flat2 ? point_on_patch_flat(B, pv, Interval(0.0, 1.0), Interval(0.0, 1.0), nullptr, nullptr)
                                     : point_on_patch_anywhere(B, pv, nullptr, nullptr);
                    if (inc) {
                        if (flat2) {
                            if (cH1v.lo <= 0.0 && 0.0 <= cH1v.hi) return 0;
                            if (std::fabs(cH1v.mid()) <= htol) return 0;
                        } else if (std::fabs(cH1v.mid()) <= htol) return 0;
                    }
                }
            }
        }
        Interval cU1, cH1, cHv1;
        if (certify_crit1d(X.u, Interval::point(X.v.mid()), evalH1, &cU1, &cH1, &cHv1)) {
            // 定理 2b 修正：切片认证点也必须在 S2 的像上
            {
                Interval P1[3], d1[3], d2[3];
                if (!eval_rational_duv(A, cU1, Interval::point(X.v.mid()), P1, d1, d2)) return 4;
                double pv[3] = { P1[0].mid(), P1[1].mid(), P1[2].mid() };
                bool inc = flat2 ? point_on_patch_flat(B, pv, Interval(0.0, 1.0), Interval(0.0, 1.0), nullptr, nullptr)
                             : point_on_patch_anywhere(B, pv, nullptr, nullptr);
            if (!inc) return 4;
            }
            if (flat2) {
                if (cH1.lo <= 0.0 && 0.0 <= cH1.hi) return 0;
                if (std::fabs(cH1.mid()) <= htol) return 0;
            } else if (std::fabs(cH1.mid()) <= htol) {
                return 0;
            }
        }
    }
    // ③ 边界情形：切触点落在盒边界上时，把盒向外膨胀 25% 再认证一次
    if (allow_inflate) {
        auto infl = [](const Interval& I, double frac) {
            double e = frac * I.width();
            double lo = std::max(0.0, I.lo - e), hi = std::min(1.0, I.hi + e);
            return Interval(lo, hi);
        };
        Box4 Y;
        Y.u = infl(X.u, 0.25); Y.v = infl(X.v, 0.25);
        Y.s = infl(X.s, 0.25); Y.t = infl(X.t, 0.25);
        if (certify_tangency_box_impl(A, B, Y, false) == 0) return 0;
    }
    return 4;
}

inline int certify_tangency_box(const RBezierSurface& A, const RBezierSurface& B, const Box4& X) {
    return certify_tangency_box_impl(A, B, X, true);
}

// ============================================================================
// M6：拓扑三分 —— 每个覆盖盒 ∈ {认证分支盒, 细化排除, undecided}
//      盒带邻接图 → 连通分量 / 端点计数
// ============================================================================

// 面邻接：恰好一维"仅接触"，其余三维有正重叠（避免角/边接触被当成连接）
inline bool boxes_face_adjacent(const Box4& a, const Box4& b, double tol) {
    const Interval* A[4] = { &a.u, &a.v, &a.s, &a.t };
    const Interval* B[4] = { &b.u, &b.v, &b.s, &b.t };
    int touch = 0, over = 0;
    for (int k = 0; k < 4; ++k) {
        double lo = std::max(A[k]->lo, B[k]->lo);
        double hi = std::min(A[k]->hi, B[k]->hi);
        if (lo > hi + tol) return false;              // 该维分离
        if (hi - lo > tol) ++over;
        else ++touch;
    }
    return touch == 1 && over == 3;
}

// 递归细分证明：整盒被排除（返回 true = 证明为空）
inline bool refine_excluded(const RBezierSurface& A, const RBezierSurface& B,
                            const Box4& box, int depth) {
    if (disjoint(hull_aabb(A), hull_aabb(B))) return true;
    if (depth == 0) return false;
    const double wu = box.u.width(), wv = box.v.width();
    const double ws = box.s.width(), wt = box.t.width();
    if (wu >= wv && wu >= ws && wu >= wt) {
        double mid = 0.5 * (box.u.lo + box.u.hi);
        RBezierSurface L, R; split_u(A, L, R);
        Box4 b1 = box, b2 = box;
        b1.u = Interval(box.u.lo, mid); b2.u = Interval(mid, box.u.hi);
        return refine_excluded(L, B, b1, depth - 1) && refine_excluded(R, B, b2, depth - 1);
    } else if (wv >= ws && wv >= wt) {
        double mid = 0.5 * (box.v.lo + box.v.hi);
        RBezierSurface L, R; split_v(A, L, R);
        Box4 b1 = box, b2 = box;
        b1.v = Interval(box.v.lo, mid); b2.v = Interval(mid, box.v.hi);
        return refine_excluded(L, B, b1, depth - 1) && refine_excluded(R, B, b2, depth - 1);
    } else if (ws >= wt) {
        double mid = 0.5 * (box.s.lo + box.s.hi);
        RBezierSurface L, R; split_u(B, L, R);
        Box4 b1 = box, b2 = box;
        b1.s = Interval(box.s.lo, mid); b2.s = Interval(mid, box.s.hi);
        return refine_excluded(A, L, b1, depth - 1) && refine_excluded(A, R, b2, depth - 1);
    } else {
        double mid = 0.5 * (box.t.lo + box.t.hi);
        RBezierSurface L, R; split_v(B, L, R);
        Box4 b1 = box, b2 = box;
        b1.t = Interval(box.t.lo, mid); b2.t = Interval(mid, box.t.hi);
        return refine_excluded(A, L, b1, depth - 1) && refine_excluded(A, R, b2, depth - 1);
    }
}

// 覆盖盒在 3D 空间中的像包络（区间求值两片子片）
inline AABB3 box_image_aabb(const RBezierSurface& A, const RBezierSurface& B, const Box4& X) {
    Interval s1[3], s2[3];
    Interval a1[3], a2[3], a3[3], b1[3], b2[3], b3[3];
    AABB3 out{ Interval::empty(), Interval::empty(), Interval::empty() };
    if (eval_rational_duv(A, X.u, X.v, s1, a1, a2) && eval_rational_duv(B, X.s, X.t, s2, b1, b2)) {
        out.x = Interval(std::min(s1[0].lo, s2[0].lo), std::max(s1[0].hi, s2[0].hi));
        out.y = Interval(std::min(s1[1].lo, s2[1].lo), std::max(s1[1].hi, s2[1].hi));
        out.z = Interval(std::min(s1[2].lo, s2[2].lo), std::max(s1[2].hi, s2[2].hi));
    }
    return out;
}

// 3D 像 AABB 面/体接触（用于沿曲线的连通性）
inline bool aabb_overlap3(const AABB3& a, const AABB3& b, double tol) {
    return (a.x.hi >= b.x.lo - tol && b.x.hi >= a.x.lo - tol) &&
           (a.y.hi >= b.y.lo - tol && b.y.hi >= a.y.lo - tol) &&
           (a.z.hi >= b.z.lo - tol && b.z.hi >= a.z.lo - tol);
}

// ============================================================================
// P2：安全细分条件 —— 交线可达半径（reach）的区间下界
//   Σ 的局部曲率 κ 上界：γ 为 Σ 的弧长参数化，γ'(s)=null(J)，
//     J γ'' = −D²F[γ',γ']  ⇒ |γ''| ≤ ‖J⁺‖·‖D²F‖
//     3D 像曲线曲率 κ ≤ ‖D²S₁‖ + ‖DS₁‖·|γ''|
//   reach r ≥ 1/κ_ub；安全细分要求 τ < r/2。
//   ‖J⁺‖ ≤ 1/√m，m 为 rank3_certified 的 Gershgorin 下界。
// ============================================================================
inline double reach_lower_bound(const RBezierSurface& A, const RBezierSurface& B, const Box4& X) {
    Interval S1[3], S1u[3], S1v[3], S1uu[3], S1uv[3], S1vv[3];
    Interval S2[3], S2u[3], S2v[3], S2uu[3], S2uv[3], S2vv[3];
    if (!eval_rational_d2(A, X.u, X.v, S1, S1u, S1v, S1uu, S1uv, S1vv)) return -1.0;
    if (!eval_rational_d2(B, X.s, X.t, S2, S2u, S2v, S2uu, S2uv, S2vv)) return -1.0;
    Interval J[3][4];
    for (int c = 0; c < 3; ++c) {
        J[c][0] = S1u[c]; J[c][1] = S1v[c];
        J[c][2] = Interval(-S2u[c].hi, -S2u[c].lo);
        J[c][3] = Interval(-S2v[c].hi, -S2v[c].lo);
    }
    double margin = 0.0;
    if (!rank3_certified(J, &margin) || margin <= 0.0) return -1.0;   // 非安全（秩亏）
    auto frob = [](const Interval* v, int n) {
        double acc = 0.0;
        for (int i = 0; i < n; ++i) {
            double m = std::max(std::fabs(v[i].lo), std::fabs(v[i].hi));
            acc += m * m;
        }
        return std::sqrt(acc);
    };
    // ‖DS₁‖, ‖D²S₁‖, ‖D²F‖
    Interval ds1[6] = { S1u[0], S1u[1], S1u[2], S1v[0], S1v[1], S1v[2] };
    Interval d2s1[9] = { S1uu[0], S1uu[1], S1uu[2], S1uv[0], S1uv[1], S1uv[2], S1vv[0], S1vv[1], S1vv[2] };
    Interval d2s2[9] = { S2uu[0], S2uu[1], S2uu[2], S2uv[0], S2uv[1], S2uv[2], S2vv[0], S2vv[1], S2vv[2] };
    double ds1_ub = frob(ds1, 6);
    double d2s1_ub = frob(d2s1, 9);
    double d2f_ub = d2s1_ub + frob(d2s2, 9);
    double sigma_min_lb = std::sqrt(margin);
    double gamma2_ub = d2f_ub / sigma_min_lb;          // |γ''| 上界
    double kappa_ub = d2s1_ub + ds1_ub * gamma2_ub;    // 曲线曲率上界
    if (kappa_ub <= 1e-12) return 1e300;
    return 1.0 / kappa_ub;
}

struct TopoStats {
    long cover = 0, certified = 0, empty = 0, undecided = 0;
    long band_components = 0;      // 全部覆盖盒的连通分量
    long cert_components = 0;      // 仅认证盒的连通分量
    long endpoints = 0;            // 认证盒中度数 = 1 的盒数
    long max_degree = 0;
    long safe_boxes = 0;           // 满足 P2 安全细分条件的认证盒
    double min_reach = 1e300;      // 认证盒上 reach 下界的最小值
    double tau = 0.0;
};

// 并查集
struct DSU {
    std::vector<int> p;
    explicit DSU(int n) : p(n) { for (int i = 0; i < n; ++i) p[i] = i; }
    int find(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) p[a] = b; }
};

inline TopoStats topology_analyze(const RBezierSurface& A, const RBezierSurface& B, double tau,
                                  int refine_depth = 3) {
    TopoStats st;
    st.tau = tau;
    auto r = cover(A, B, tau);
    st.cover = (long)r.boxes.size();
    std::vector<char> cert(r.boxes.size(), 0), empty(r.boxes.size(), 0);
    for (std::size_t i = 0; i < r.boxes.size(); ++i) {
        const Box4& X = r.boxes[i];
        if (certify_tangency_box(A, B, X) == 0 || krawczyk_slice(A, B, X) ||
            certify_box_deep(A, B, X)) {
            cert[i] = 1; ++st.certified;
            double rl = reach_lower_bound(A, B, X);
            if (rl > 0.0 && tau < 0.5 * rl) ++st.safe_boxes;
            if (rl > 0.0) st.min_reach = std::min(st.min_reach, rl);
        } else if (refine_excluded(A, B, X, refine_depth)) {
            empty[i] = 1; ++st.empty;
        } else {
            ++st.undecided;
        }
    }
    // 3D 像包络（沿曲线的连通性用 3D 重叠判断，避免 4D 对角盒只角接触）
    std::vector<AABB3> img(r.boxes.size());
    for (std::size_t i = 0; i < r.boxes.size(); ++i) img[i] = box_image_aabb(A, B, r.boxes[i]);
    const double img_tol = 1e-9;
    // 盒带连通分量（全部覆盖盒）
    {
        DSU d((int)r.boxes.size());
        for (std::size_t i = 0; i < r.boxes.size(); ++i)
            for (std::size_t j = i + 1; j < r.boxes.size(); ++j)
                if (aabb_overlap3(img[i], img[j], img_tol)) d.unite((int)i, (int)j);
        std::vector<char> seen(r.boxes.size(), 0);
        for (std::size_t i = 0; i < r.boxes.size(); ++i) {
            int k = d.find((int)i);
            if (!seen[k]) { seen[k] = 1; ++st.band_components; }
        }
    }
    // 仅认证盒：分量数 + 端点/度数
    {
        DSU d((int)r.boxes.size());
        std::vector<int> deg(r.boxes.size(), 0);
        for (std::size_t i = 0; i < r.boxes.size(); ++i) {
            if (!cert[i]) continue;
            for (std::size_t j = i + 1; j < r.boxes.size(); ++j) {
                if (!cert[j]) continue;
                if (aabb_overlap3(img[i], img[j], img_tol)) {
                    d.unite((int)i, (int)j);
                    ++deg[i]; ++deg[j];
                }
            }
        }
        std::vector<char> seen(r.boxes.size(), 0);
        for (std::size_t i = 0; i < r.boxes.size(); ++i) {
            if (!cert[i]) continue;
            int k = d.find((int)i);
            if (!seen[k]) { seen[k] = 1; ++st.cert_components; }
            // 端点代理：与参数域边界接触的认证盒
            { const Box4& X = r.boxes[i];
              if (X.u.lo <= 1e-12 || X.u.hi >= 1.0 - 1e-12 || X.v.lo <= 1e-12 || X.v.hi >= 1.0 - 1e-12 ||
                  X.s.lo <= 1e-12 || X.s.hi >= 1.0 - 1e-12 || X.t.lo <= 1e-12 || X.t.hi >= 1.0 - 1e-12)
                  ++st.endpoints; }
            st.max_degree = std::max(st.max_degree, (long)deg[i]);
        }
    }
    return st;
}

// ============================================================================
// M6.5：分支追踪 —— 认证盒 → 近似解种子 → 沿零空间方向延拓 → 有序分支
// ============================================================================

// 浮点 F、J（点区间求值取 mid）
inline bool eval_fj_float(const RBezierSurface& A, const RBezierSurface& B, const double p[4],
                          double F[3], double J[3][4]) {
    Interval s1[3], a1[3], a2[3], s2[3], b1[3], b2[3];
    if (!eval_rational_duv(A, Interval::point(p[0]), Interval::point(p[1]), s1, a1, a2)) return false;
    if (!eval_rational_duv(B, Interval::point(p[2]), Interval::point(p[3]), s2, b1, b2)) return false;
    for (int c = 0; c < 3; ++c) {
        F[c] = s1[c].mid() - s2[c].mid();
        J[c][0] = a1[c].mid(); J[c][1] = a2[c].mid();
        J[c][2] = -b1[c].mid(); J[c][3] = -b2[c].mid();
    }
    return true;
}

inline bool pinv_y(const double J[3][4], double Y[4][3]) {
    double Am[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double acc = 0.0;
            for (int k = 0; k < 4; ++k) acc += J[i][k] * J[j][k];
            Am[i][j] = acc;
        }
    double M[3][6];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) M[i][j] = Am[i][j];
        for (int j = 0; j < 3; ++j) M[i][3 + j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < 3; ++col) {
        int piv = col; double best = std::fabs(M[col][col]);
        for (int r = col + 1; r < 3; ++r) if (std::fabs(M[r][col]) > best) { best = std::fabs(M[r][col]); piv = r; }
        if (best < 1e-14) return false;
        if (piv != col) for (int j = 0; j < 6; ++j) std::swap(M[piv][j], M[col][j]);
        double d = M[col][col];
        for (int j = 0; j < 6; ++j) M[col][j] /= d;
        for (int r = 0; r < 3; ++r) if (r != col) {
            double f = M[r][col];
            for (int j = 0; j < 6; ++j) M[r][j] -= f * M[col][j];
        }
    }
    double Ainv[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Ainv[i][j] = M[i][3 + j];
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 3; ++i) {
            double acc = 0.0;
            for (int j = 0; j < 3; ++j) acc += J[j][k] * Ainv[j][i];
            Y[k][i] = acc;
        }
    return true;
}

// 零空间方向（I - YJ 的主特征向量，幂迭代）
inline bool null_dir(const double J[3][4], double t[4]) {
    double Y[4][3];
    if (!pinv_y(J, Y)) return false;
    double N[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            double acc = 0.0;
            for (int i = 0; i < 3; ++i) acc += Y[a][i] * J[i][b];
            N[a][b] = ((a == b) ? 1.0 : 0.0) - acc;
        }
    double v[4] = { 0.3, 0.5, 0.7, 0.9 };
    for (int it = 0; it < 80; ++it) {
        double w[4] = { 0, 0, 0, 0 };
        for (int a = 0; a < 4; ++a) for (int b = 0; b < 4; ++b) w[a] += N[a][b] * v[b];
        double nrm = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        if (nrm < 1e-14) return false;
        for (int a = 0; a < 4; ++a) v[a] = w[a] / nrm;
    }
    double nrm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
    for (int a = 0; a < 4; ++a) t[a] = v[a] / nrm;
    return true;
}

// Gauss-Newton 投影回解集
inline bool project_to_solution(const RBezierSurface& A, const RBezierSurface& B,
                                double p[4], int iters = 25) {
    for (int it = 0; it < iters; ++it) {
        double F[3], J[3][4], Y[4][3];
        if (!eval_fj_float(A, B, p, F, J)) return false;
        double nf = std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
        if (nf < 1e-12) return true;
        if (!pinv_y(J, Y)) return false;
        for (int a = 0; a < 4; ++a) {
            double d = 0.0;
            for (int i = 0; i < 3; ++i) d += Y[a][i] * F[i];
            p[a] -= d;
        }
    }
    double F[3], J[3][4];
    if (!eval_fj_float(A, B, p, F, J)) return false;
    double nf = std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
    return nf < 1e-9;
}

struct Branch {
    std::vector<std::array<double, 4>> pts;
    bool closed = false;
};

inline bool inside01(const double p[4], double tol = 1e-9) {
    for (int i = 0; i < 4; ++i) if (p[i] < -tol || p[i] > 1.0 + tol) return false;
    return true;
}

// 从种子双向追踪；返回是否成功
inline bool trace_from(const RBezierSurface& A, const RBezierSurface& B,
                       const double seed[4], double h, int max_steps, Branch& br) {
    br.pts.clear(); br.closed = false;
    // 先向前
    std::vector<std::array<double, 4>> fwd, bwd;
    for (int dir = 0; dir < 2; ++dir) {
        double p[4] = { seed[0], seed[1], seed[2], seed[3] };
        double t[4];
        if (!project_to_solution(A, B, p)) continue;
        double F[3], J[3][4];
        if (!eval_fj_float(A, B, p, F, J) || !null_dir(J, t)) continue;
        if (dir == 1) for (int a = 0; a < 4; ++a) t[a] = -t[a];
        std::vector<std::array<double, 4>> path;
        for (int step = 0; step < max_steps; ++step) {
            double q[4];
            for (int a = 0; a < 4; ++a) q[a] = p[a] + h * t[a];
            if (!inside01(q)) break;
            if (!project_to_solution(A, B, q)) break;
            // 更新切向并保持方向连续
            double F2[3], J2[3][4], t2[4];
            if (!eval_fj_float(A, B, q, F2, J2) || !null_dir(J2, t2)) break;
            double dot = 0.0;
            for (int a = 0; a < 4; ++a) dot += t2[a] * t[a];
            if (dot < 0.0) for (int a = 0; a < 4; ++a) t2[a] = -t2[a];
            for (int a = 0; a < 4; ++a) t[a] = t2[a];
            path.push_back({ q[0], q[1], q[2], q[3] });
            for (int a = 0; a < 4; ++a) p[a] = q[a];
        }
        if (dir == 0) fwd = path;
        else bwd = path;
    }
    // 合并：bwd 反向 + seed + fwd
    for (int i = (int)bwd.size() - 1; i >= 0; --i) br.pts.push_back(bwd[i]);
    br.pts.push_back({ seed[0], seed[1], seed[2], seed[3] });
    for (const auto& q : fwd) br.pts.push_back(q);
    if (br.pts.size() >= 3) {
        // 闭环判定：首尾接近
        const auto& a = br.pts.front(); const auto& b = br.pts.back();
        double d2 = 0.0;
        for (int i = 0; i < 4; ++i) d2 += (a[i] - b[i]) * (a[i] - b[i]);
        br.closed = (std::sqrt(d2) < 2.0 * h);
    }
    return br.pts.size() >= 2;
}

struct TraceStats {
    long seeds = 0, branches = 0, endpoints = 0, closed = 0;
    double total_len = 0.0;
};

inline TraceStats trace_branches(const RBezierSurface& A, const RBezierSurface& B, double tau,
                                 double h = -1.0, int max_steps = 4000) {
    TraceStats ts;
    if (h <= 0.0) h = std::max(1e-4, 2.0 * tau);
    auto r = cover(A, B, tau);
    // 认证盒 → 种子（盒中心投影）
    std::vector<std::array<double, 4>> seeds;
    for (const Box4& X : r.boxes) {
        if (certify_tangency_box(A, B, X) != 0 && !krawczyk_slice(A, B, X) && !certify_box_deep(A, B, X))
            continue;
        double p[4] = { X.u.mid(), X.v.mid(), X.s.mid(), X.t.mid() };
        if (!project_to_solution(A, B, p)) continue;
        bool dup = false;
        for (const auto& q : seeds) {
            double d2 = 0.0;
            for (int i = 0; i < 4; ++i) d2 += (p[i] - q[i]) * (p[i] - q[i]);
            if (std::sqrt(d2) < 2.5 * tau) { dup = true; break; }
        }
        if (!dup) seeds.push_back({ p[0], p[1], p[2], p[3] });
    }
    ts.seeds = (long)seeds.size();
    // 追踪：未覆盖的种子逐个追踪
    std::vector<char> used(seeds.size(), 0);
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        if (used[i]) continue;
        Branch br;
        if (!trace_from(A, B, seeds[i].data(), h, max_steps, br)) continue;
        ++ts.branches;
        if (br.closed) ++ts.closed; else ts.endpoints += 2;
        for (std::size_t k = 1; k < br.pts.size(); ++k) {
            double d2 = 0.0;
            for (int a = 0; a < 4; ++a) d2 += (br.pts[k][a] - br.pts[k - 1][a]) * (br.pts[k][a] - br.pts[k - 1][a]);
            ts.total_len += std::sqrt(d2);
        }
        // 标记被该分支覆盖的种子
        for (std::size_t j = 0; j < seeds.size(); ++j) {
            if (used[j] || j == i) continue;
            for (const auto& q : br.pts) {
                double d2 = 0.0;
                for (int a = 0; a < 4; ++a) d2 += (seeds[j][a] - q[a]) * (seeds[j][a] - q[a]);
                if (std::sqrt(d2) < 2.0 * h) { used[j] = 1; break; }
            }
        }
        used[i] = 1;
    }
    return ts;
}

// ============================================================================
// M7：2D 重合证书 —— 点包含（区间 Newton 2D）+ 切平面一致 + 网格覆盖
// ============================================================================

// 2x3 广义逆：Y = (JᵀJ)^{-1} Jᵀ，J 为 3x2
inline bool pinv32(const double J[3][2], double Y[2][3]) {
    double A00 = J[0][0]*J[0][0] + J[1][0]*J[1][0] + J[2][0]*J[2][0];
    double A01 = J[0][0]*J[0][1] + J[1][0]*J[1][1] + J[2][0]*J[2][1];
    double A11 = J[0][1]*J[0][1] + J[1][1]*J[1][1] + J[2][1]*J[2][1];
    double det = A00*A11 - A01*A01;
    if (std::fabs(det) < 1e-14) return false;
    double i00 = A11/det, i01 = -A01/det, i11 = A00/det;
    for (int j = 0; j < 3; ++j) {
        // Y = (J^T J)^{-1} J^T（2x3）；注意 J 的行是空间分量、列是方向
        Y[0][j] = i00*J[j][0] + i01*J[j][1];
        Y[1][j] = i01*J[j][0] + i11*J[j][1];
    }
    return true;
}

// 平面伙伴的点包含：
//   ① 选取 3x2 Jacobian 中最好的两行（2x2 非奇异）求解方系统（Krawczyk，存在+唯一）
//   ② 用平面距离检查第三坐标（平面片 ⇒ 法向恒定，距离=0 即包含）
inline bool point_on_patch_flat(const RBezierSurface& S, const double P[3],
                                const Interval& Ss, const Interval& St,
                                Interval* outS, Interval* outT) {
    // 选两行
    Interval Q[3], Su[3], Sv[3];
    if (!eval_rational_duv(S, Ss, St, Q, Su, Sv)) return false;
    int best_r1 = 0, best_r2 = 1;
    double best_det = -1.0;
    for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 3; ++b) {
            double det = Su[a].mid() * Sv[b].mid() - Su[b].mid() * Sv[a].mid();
            if (std::fabs(det) > best_det) { best_det = std::fabs(det); best_r1 = a; best_r2 = b; }
        }
    if (best_det < 1e-12) return false;
    // 2x2 Krawczyk（方系统）
    Interval U = Ss, V = St;
    for (int it = 0; it < 30; ++it) {
        Interval Qc[3], Suc[3], Svc[3];
        if (!eval_rational_duv(S, U, V, Qc, Suc, Svc)) return false;
        double mu = U.mid(), mv = V.mid();
        Interval Qm[3], Sum[3], Svm[3];
        if (!eval_rational_duv(S, Interval::point(mu), Interval::point(mv), Qm, Sum, Svm)) return false;
        double J00 = Sum[best_r1].mid(), J01 = Svm[best_r1].mid();
        double J10 = Sum[best_r2].mid(), J11 = Svm[best_r2].mid();
        double det = J00 * J11 - J01 * J10;
        if (std::fabs(det) < 1e-14) return false;
        double i00 = J11 / det, i01 = -J01 / det, i10 = -J10 / det, i11 = J00 / det;
        double F0 = Qm[best_r1].mid() - P[best_r1];
        double F1 = Qm[best_r2].mid() - P[best_r2];
        Interval W00 = isub(Interval::point(1.0), iadd(imul(Interval::point(i00), Suc[best_r1]), imul(Interval::point(i01), Svc[best_r1])));
        Interval W01 = isub(Interval::point(0.0), iadd(imul(Interval::point(i00), Suc[best_r2]), imul(Interval::point(i01), Svc[best_r2])));
        Interval W10 = isub(Interval::point(0.0), iadd(imul(Interval::point(i10), Suc[best_r1]), imul(Interval::point(i11), Svc[best_r1])));
        Interval W11 = isub(Interval::point(1.0), iadd(imul(Interval::point(i10), Suc[best_r2]), imul(Interval::point(i11), Svc[best_r2])));
        Interval Ku = isub(Interval(mu, mu), iadd(imul(Interval::point(i00), Interval::point(F0)), imul(Interval::point(i01), Interval::point(F1))));
        Interval Kv = isub(Interval(mv, mv), iadd(imul(Interval::point(i10), Interval::point(F0)), imul(Interval::point(i11), Interval::point(F1))));
        Interval du = isub(U, Interval::point(mu)), dv = isub(V, Interval::point(mv));
        Ku = iadd(Ku, iadd(imul(W00, du), imul(W01, dv)));
        Kv = iadd(Kv, iadd(imul(W10, du), imul(W11, dv)));
        bool inside = (U.lo <= Ku.lo && Ku.hi <= U.hi && V.lo <= Kv.lo && Kv.hi <= V.hi) &&
                      (Ku.width() <= 0.9 * U.width() || U.width() < 1e-12) &&
                      (Kv.width() <= 0.9 * V.width() || V.width() < 1e-12);
        Interval nU(std::max(U.lo, Ku.lo), std::min(U.hi, Ku.hi));
        Interval nV(std::max(V.lo, Kv.lo), std::min(V.hi, Kv.hi));
        if (inside) {
            // 第三坐标/平面距离检查：在认证小盒上评估 S 与 P 的距离（平面片法向恒定）
            Interval Qb[3], Qbu[3], Qbv[3];
            if (!eval_rational_duv(S, nU, nV, Qb, Qbu, Qbv)) return false;
            double d2 = 0.0;
            for (int c = 0; c < 3; ++c) {
                double lo = Qb[c].lo - P[c], hi = Qb[c].hi - P[c];
                double m = std::min(std::fabs(lo), std::fabs(hi));
                if (lo * hi > 0.0) m = std::fabs(lo) < std::fabs(hi) ? std::fabs(lo) : std::fabs(hi);
                d2 += m * m;
            }
            double scale = 1.0 + std::fabs(P[0]) + std::fabs(P[1]) + std::fabs(P[2]);
            if (std::sqrt(d2) > 1e-7 * scale) return false;
            if (outS) *outS = nU;
            if (outT) *outT = nV;
            return true;
        }
        if (nU.isEmpty() || nV.isEmpty()) return false;
        if (!(nU.lo > U.lo || nU.hi < U.hi || nV.lo > V.lo || nV.hi < V.hi)) return false;
        U = nU; V = nV;
    }
    return false;
}

// 点包含：P 是否在面片 S 的像上（2D 区间 Newton，闭包含 + 收缩）
inline bool point_on_patch(const RBezierSurface& S, const double P[3], double s0, double t0,
                           Interval* outS, Interval* outT) {
    Interval U(std::max(0.0, s0 - 0.06), std::min(1.0, s0 + 0.06));
    Interval V(std::max(0.0, t0 - 0.06), std::min(1.0, t0 + 0.06));
    for (int it = 0; it < 40; ++it) {
        Interval Q[3], Su[3], Sv[3];
        if (!eval_rational_duv(S, U, V, Q, Su, Sv)) return false;
        double mu = U.mid(), mv = V.mid();
        Interval Qm[3], Sum[3], Svm[3];
        if (!eval_rational_duv(S, Interval::point(mu), Interval::point(mv), Qm, Sum, Svm)) return false;
        double Fm[3], J[3][2];
        for (int c = 0; c < 3; ++c) { Fm[c] = Qm[c].mid() - P[c]; J[c][0] = Sum[c].mid(); J[c][1] = Svm[c].mid(); }
        double Y[2][3];
        if (!pinv32(J, Y)) return false;
        // 关键：F: R^2 -> R^3 为超定系统，YF=0 仅意味 F 在 J 的零空间中。
        // 必须要求"法向残差" r = F(m) - J d (d=YF) 足够小，否则会接受不在像上的点。
        {
            double d0 = 0.0, d1 = 0.0;
            for (int c = 0; c < 3; ++c) { d0 += Y[0][c] * Fm[c]; d1 += Y[1][c] * Fm[c]; }
            double r2 = 0.0;
            for (int c = 0; c < 3; ++c) {
                double rc = Fm[c] - (J[c][0] * d0 + J[c][1] * d1);
                r2 += rc * rc;
            }
            double scale = 1.0 + std::fabs(P[0]) + std::fabs(P[1]) + std::fabs(P[2]);
            if (std::sqrt(r2) > 1e-9 * scale) return false;
        }
        Interval W00 = Interval::point(1.0), W01 = Interval::point(0.0),
                 W10 = Interval::point(0.0), W11 = Interval::point(1.0);
        for (int c = 0; c < 3; ++c) {
            W00 = isub(W00, imul(Interval::point(Y[0][c]), Su[c]));
            W01 = isub(W01, imul(Interval::point(Y[0][c]), Sv[c]));
            W10 = isub(W10, imul(Interval::point(Y[1][c]), Su[c]));
            W11 = isub(W11, imul(Interval::point(Y[1][c]), Sv[c]));
        }
        Interval Ku(mu, mu), Kv(mv, mv);
        for (int c = 0; c < 3; ++c) {
            Ku = isub(Ku, imul(Interval::point(Y[0][c]), Interval::point(Fm[c])));
            Kv = isub(Kv, imul(Interval::point(Y[1][c]), Interval::point(Fm[c])));
        }
        Interval du = isub(U, Interval::point(mu)), dv = isub(V, Interval::point(mv));
        Ku = iadd(Ku, iadd(imul(W00, du), imul(W01, dv)));
        Kv = iadd(Kv, iadd(imul(W10, du), imul(W11, dv)));
        Interval nU(std::max(U.lo, Ku.lo), std::min(U.hi, Ku.hi));
        Interval nV(std::max(V.lo, Kv.lo), std::min(V.hi, Kv.hi));
        bool inside = (U.lo <= Ku.lo && Ku.hi <= U.hi && V.lo <= Kv.lo && Kv.hi <= V.hi) &&
                      (Ku.width() <= 0.9 * U.width() || U.width() < 1e-12) &&
                      (Kv.width() <= 0.9 * V.width() || V.width() < 1e-12);
        if (inside) {
            if (outS) *outS = nU;
            if (outT) *outT = nV;
            return true;
        }
        if (nU.isEmpty() || nV.isEmpty()) return false;
        if (!(nU.lo > U.lo || nU.hi < U.hi || nV.lo > V.lo || nV.hi < V.hi)) return false;
        U = nU; V = nV;
    }
    return false;
}

// 全参数域点包含：5x5 初值网格搜索（相切认证中预像可在对方任意处）
inline bool point_on_patch_anywhere(const RBezierSurface& S, const double P[3],
                                    Interval* outS, Interval* outT) {
    for (int i = 1; i <= 5; ++i)
        for (int j = 1; j <= 5; ++j) {
            double s0 = i / 6.0, t0 = j / 6.0;
            if (point_on_patch(S, P, s0, t0, outS, outT)) return true;
        }
    return false;
}

// ============================================================================
// P5：平面仿射片的精确 2D 重合判定
//   若两片均为平面仿射参数化（S(u,v)=O+u·U+v·V），则：
//   - 图像是平行四边形；
//   - "A 的像 ⊆ B 的像" ⇔ A 的 4 个角点在 B 的参数域内（凸集，顶点判据）；
//   - "部分 2D 重叠" ⇔ 两平行四边形交集面积 > 0（SAT + 裁剪求面积）。
//   返回：-1 非平面仿射（退回 ε 证书）；0 无 2D 重叠；1 包含（精确重合）；
//         2 部分 2D 重叠
// ============================================================================
struct AffinePlane {
    bool ok = false;
    double O[3], U[3], V[3], n[3];
};

inline AffinePlane affine_plane_of(const RBezierSurface& S, double scale) {
    AffinePlane ap;
    Interval Q[3], Qu[3], Qv[3], Quu[3], Quv[3], Qvv[3];
    if (!eval_rational_d2(S, Interval(0.0, 1.0), Interval(0.0, 1.0), Q, Qu, Qv, Quu, Quv, Qvv))
        return ap;
    double m2 = 0.0;
    for (int c = 0; c < 3; ++c) {
        m2 = std::max(m2, std::fabs(Quu[c].mid()));
        m2 = std::max(m2, std::fabs(Quv[c].mid()));
        m2 = std::max(m2, std::fabs(Qvv[c].mid()));
    }
    if (m2 > 1e-9 * scale) return ap;   // 非仿射
    // O 取参数原点值：S(0,0)=S(0.5,0.5)−0.5U−0.5V
    for (int c = 0; c < 3; ++c) {
        ap.U[c] = Qu[c].mid(); ap.V[c] = Qv[c].mid();
        ap.O[c] = Q[c].mid() - 0.5 * ap.U[c] - 0.5 * ap.V[c];
    }
    ap.n[0] = ap.U[1] * ap.V[2] - ap.U[2] * ap.V[1];
    ap.n[1] = ap.U[2] * ap.V[0] - ap.U[0] * ap.V[2];
    ap.n[2] = ap.U[0] * ap.V[1] - ap.U[1] * ap.V[0];
    double ln = std::sqrt(ap.n[0] * ap.n[0] + ap.n[1] * ap.n[1] + ap.n[2] * ap.n[2]);
    if (ln < 1e-14) return ap;
    for (int c = 0; c < 3; ++c) ap.n[c] /= ln;
    ap.ok = true;
    return ap;
}

// 平行四边形（B 的像）是否包含点 P；返回参数 (s,t)
inline bool affine_contains(const AffinePlane& B, const double P[3], double tol) {
    // 解 O + s U + t V = P（2x2 正规方程）
    double a = 0, b = 0, c = 0, r0 = 0, r1 = 0;
    double d[3] = { P[0] - B.O[0], P[1] - B.O[1], P[2] - B.O[2] };
    for (int k = 0; k < 3; ++k) {
        a += B.U[k] * B.U[k]; b += B.U[k] * B.V[k]; c += B.V[k] * B.V[k];
        r0 += B.U[k] * d[k]; r1 += B.V[k] * d[k];
    }
    double det = a * c - b * b;
    if (std::fabs(det) < 1e-30) return false;
    double s = (c * r0 - b * r1) / det;
    double t = (a * r1 - b * r0) / det;
    return s >= -tol && s <= 1.0 + tol && t >= -tol && t <= 1.0 + tol;
}

inline void affine_corner(const AffinePlane& A, double u, double v, double out[3]) {
    for (int k = 0; k < 3; ++k) out[k] = A.O[k] + u * A.U[k] + v * A.V[k];
}

inline int exact_coplanar_overlap(const RBezierSurface& A, const RBezierSurface& B, double scale) {
    AffinePlane pa = affine_plane_of(A, scale);
    AffinePlane pb = affine_plane_of(B, scale);
    if (!pa.ok || !pb.ok) return -1;
    // 法向平行 + 共面
    double cross = 0.0;
    cross += std::fabs(pa.n[1] * pb.n[2] - pa.n[2] * pb.n[1]);
    cross += std::fabs(pa.n[2] * pb.n[0] - pa.n[0] * pb.n[2]);
    cross += std::fabs(pa.n[0] * pb.n[1] - pa.n[1] * pb.n[0]);
    if (cross > 1e-6) return 0;
    double d = 0.0;
    for (int k = 0; k < 3; ++k) d += pa.n[k] * (pb.O[k] - pa.O[k]);
    if (std::fabs(d) > 1e-6 * scale) return 0;
    // 包含：A 的 4 角在 B 内（凸集顶点判据）
    bool a_in_b = true, b_in_a = true;
    double P[3];
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 1; ++j) {
            affine_corner(pa, (double)i, (double)j, P);
            if (!affine_contains(pb, P, 1e-7)) a_in_b = false;
            affine_corner(pb, (double)i, (double)j, P);
            if (!affine_contains(pa, P, 1e-7)) b_in_a = false;
        }
    if (a_in_b || b_in_a) return 1;
    // 部分重叠：SAT（分离轴）判断交集非空
    auto sat_disjoint = [](const AffinePlane& X, const AffinePlane& Y) {
        double axes[4][3] = { { X.n[0], X.n[1], X.n[2] }, { Y.n[0], Y.n[1], Y.n[2] }, {0,0,0}, {0,0,0} };
        // 平面内两条棱方向×法向得到面内轴
        double e1[3], e2[3];
        for (int k = 0; k < 3; ++k) { e1[k] = X.U[k]; e2[k] = X.V[k]; }
        axes[2][0] = e1[1] * X.n[2] - e1[2] * X.n[1];
        axes[2][1] = e1[2] * X.n[0] - e1[0] * X.n[2];
        axes[2][2] = e1[0] * X.n[1] - e1[1] * X.n[0];
        axes[3][0] = e2[1] * X.n[2] - e2[2] * X.n[1];
        axes[3][1] = e2[2] * X.n[0] - e2[0] * X.n[2];
        axes[3][2] = e2[0] * X.n[1] - e2[1] * X.n[0];
        auto proj = [](const AffinePlane& Z, const double ax[3], double& lo, double& hi) {
            lo = 1e300; hi = -1e300;
            for (int i = 0; i <= 1; ++i)
                for (int j = 0; j <= 1; ++j) {
                    double P2[3], s = 0.0;
                    affine_corner(Z, (double)i, (double)j, P2);
                    for (int k = 0; k < 3; ++k) s += P2[k] * ax[k];
                    lo = std::min(lo, s); hi = std::max(hi, s);
                }
        };
        for (const auto& ax : axes) {
            double nrm = std::sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
            if (nrm < 1e-14) continue;
            double axn[3] = { ax[0] / nrm, ax[1] / nrm, ax[2] / nrm };
            double la, ha, lb, hb;
            proj(X, axn, la, ha); proj(Y, axn, lb, hb);
            if (ha < lb - 1e-9 || hb < la - 1e-9) return true;
        }
        return false;
    };
    if (!sat_disjoint(pa, pb)) return 2;
    return 0;
}

// 采样网格上逐点做点包含 + 法向平行检查；返回命中比例
inline double coincidence_ratio(const RBezierSurface& A, const RBezierSurface& B,
                                double u0, double u1, double v0, double v1, int n) {
    int hit = 0, tot = 0;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            double u = u0 + (u1 - u0) * i / n, v = v0 + (v1 - v0) * j / n;
            Interval U = Interval::point(u), V = Interval::point(v);
            Interval Q[3], Qu[3], Qv[3];
            if (!eval_rational_duv(A, U, V, Q, Qu, Qv)) continue;
            double P[3] = { Q[0].mid(), Q[1].mid(), Q[2].mid() };
            // 粗搜 B 参数网格取初值
            double bs = 0.5, bt = 0.5, bd = 1e300;
            for (int a = 0; a <= 12; ++a)
                for (int b = 0; b <= 12; ++b) {
                    double s = a / 12.0, t = b / 12.0;
                    Interval Q2[3], Q2u[3], Q2v[3];
                    if (!eval_rational_duv(B, Interval::point(s), Interval::point(t), Q2, Q2u, Q2v)) continue;
                    double d = 0.0;
                    for (int c = 0; c < 3; ++c) { double e = Q2[c].mid() - P[c]; d += e * e; }
                    if (d < bd) { bd = d; bs = s; bt = t; }
                }
            ++tot;
            Interval osc, otc;
            if (!point_on_patch(B, P, bs, bt, &osc, &otc)) continue;
            // 法向平行（A 在 (u,v)，B 在包含处）
            Interval Qb[3], Qbu[3], Qbv[3];
            if (!eval_rational_duv(B, osc, otc, Qb, Qbu, Qbv)) continue;
            double na[3], nb[3];
            cross3_mid(Qu, Qv, na);
            cross3_mid(Qbu, Qbv, nb);
            double la = std::sqrt(na[0]*na[0] + na[1]*na[1] + na[2]*na[2]);
            double lb = std::sqrt(nb[0]*nb[0] + nb[1]*nb[1] + nb[2]*nb[2]);
            if (la < 1e-12 || lb < 1e-12) continue;
            double cosang = std::fabs(na[0]*nb[0] + na[1]*nb[1] + na[2]*nb[2]) / (la * lb);
            if (cosang > 0.999) ++hit;
        }
    return tot > 0 ? (double)hit / tot : 0.0;
}

// ============================================================================
// M8b：面片对分类器 —— disjoint / 1D 横截 / 2D 重合 / 共面接触(undecided)
//      把 M7 点包含证书接入主流程
// ============================================================================

struct PairClass {
    int type = 3;            // 0=disjoint, 1=1D certified, 2=2D coincidence, 3=undecided(contact), 4=coplanar overlap
    double coincidence = 0.0;
    long cover = 0, certified = 0, undecided = 0;
    long components = 0;
    int coplanar_overlap = 0;   // 1=contain, 2=partial (exact_coplanar_overlap)
};

// 共面接触检测：抽样若干覆盖盒中点，看法向是否平行（不要求盒子包含解）
inline bool coplanar_candidate(const RBezierSurface& A, const RBezierSurface& B,
                               const std::vector<Box4>& boxes, int samples = 6) {
    int checked = 0, parallel = 0;
    for (std::size_t i = 0; i < boxes.size() && checked < samples; i += std::max<std::size_t>(1, boxes.size() / samples)) {
        const Box4& X = boxes[i];
        Interval s1[3], a1[3], a2[3], s2[3], b1[3], b2[3];
        if (!eval_rational_duv(A, X.u, X.v, s1, a1, a2)) continue;
        if (!eval_rational_duv(B, X.s, X.t, s2, b1, b2)) continue;
        double n1[3], n2[3];
        cross3_mid(a1, a2, n1); cross3_mid(b1, b2, n2);
        double l1 = std::sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
        double l2 = std::sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
        if (l1 < 1e-12 || l2 < 1e-12) continue;
        double c = std::fabs(n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]) / (l1*l2);
        ++checked;
        if (c > 0.999) ++parallel;
    }
    return checked > 0 && parallel * 10 >= checked * 9;
}

// 面片笛卡尔尺寸（控制点归一化跨度），用于共面判定阈值
inline double patch_scale(const RBezierSurface& S) {
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    double wm = S.wmin > 0 ? S.wmin : 1e-9;
    double wM = S.wmax > 0 ? S.wmax : 1.0;
    for (const auto& p : S.P)
        for (int c = 0; c < 3; ++c) {
            double a = p[c].lo / wM, b = p[c].hi / wm;
            lo[c] = std::min(lo[c], a); hi[c] = std::max(hi[c], b);
        }
    return std::max(hi[0] - lo[0], std::max(hi[1] - lo[1], hi[2] - lo[2]));
}

inline PairClass classify_pair(const RBezierSurface& A, const RBezierSurface& B, double tau, long node_cap = 4000000) {
    PairClass pc;
    // ①' 共面重叠快速预检：cover 之前，避免共面细分失控
    {
        double sc = std::max(patch_scale(A), patch_scale(B));
        int cp = exact_coplanar_overlap(A, B, sc);
        if (cp > 0) { pc.type = 4; pc.coplanar_overlap = cp; return pc; }
    }
    auto r = cover(A, B, tau, 60, node_cap);
    pc.cover = (long)r.boxes.size();
    if (r.boxes.empty()) { pc.type = 0; return pc; }
    // 覆盖触顶或 box 超 certify 预算（接触/贴合面）→ 跳过 certify 标记未决，避免 certify 爆炸
    if (pc.cover >= node_cap / 2 || pc.cover > 2500) { pc.type = 3; return pc; }
    // ① 2D 重合证书（点包含 + 法向）
    pc.coincidence = coincidence_ratio(A, B, 0.15, 0.85, 0.15, 0.85, 2);
    if (pc.coincidence >= 0.9) { pc.type = 2; return pc; }
    // 部分重叠：分块检测（4 个象限），任一块全命中即 2D 重合
    {
        const double tiles[4][4] = { {0.05,0.5,0.05,0.5},{0.5,0.95,0.05,0.5},{0.05,0.5,0.5,0.95},{0.5,0.95,0.5,0.95} };
        for (int ti = 0; ti < 4; ++ti) {
            double r2 = coincidence_ratio(A, B, tiles[ti][0], tiles[ti][1], tiles[ti][2], tiles[ti][3], 2);
            if (r2 >= 0.9) { pc.type = 2; pc.coincidence = r2; return pc; }
        }
    }
    // ② 1D 认证
    auto st = certify_boxes(A, B, r.boxes);
    pc.certified = st.certified;
    pc.undecided = st.undecided;
    if (pc.certified > 0 || st.tangent > 0) {
        pc.type = 1;
        // 连通分量
        std::vector<AABB3> img(r.boxes.size());
        for (std::size_t i = 0; i < r.boxes.size(); ++i) img[i] = box_image_aabb(A, B, r.boxes[i]);
        DSU d((int)r.boxes.size());
        for (std::size_t i = 0; i < r.boxes.size(); ++i)
            if (st.box_certified[i])
                for (std::size_t j = i + 1; j < r.boxes.size(); ++j)
                    if (st.box_certified[j] && aabb_overlap3(img[i], img[j], 1e-9)) d.unite((int)i, (int)j);
        std::vector<char> seen(r.boxes.size(), 0);
        for (std::size_t i = 0; i < r.boxes.size(); ++i) {
            if (!st.box_certified[i]) continue;
            int k = d.find((int)i);
            if (!seen[k]) { seen[k] = 1; ++pc.components; }
        }
        return pc;
    }
    // ③ 共面/接触候选（2D 证书未过、1D 认证未过）
    pc.type = coplanar_candidate(A, B, r.boxes) ? 4 : 3;
    return pc;
}

} // namespace cert
} // namespace nurbs

#endif // NURBS_CERT_SSI_HPP
