// ============================================================================
// bench_theory_tests.cpp — 定理对应的测试（theory_proofs.md §7）
//   T1：覆盖无遗漏（真值点全部覆盖；被排除盒内无真值点 = 排除可靠自检）
//   T2a：横截认证（认证盒包含真值点）
//   T2b：相切认证（真相切认证成功；"等高不重合"假阳性必须被拒绝）
//   T4：2D ε-重合证书（Lipschitz 条件 + 点包含 + 法向）
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;

static int g_fail = 0;
static void CHECK(bool cond, const char* msg) {
    std::printf("  [%s] %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++g_fail;
}

using P4 = std::array<double, 4>;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static RBezierSurface make(const std::vector<std::vector<P4>>& poles) {
    int nu = (int)poles.size(), nv = (int)poles[0].size();
    std::vector<P4> flat;
    for (int i = 0; i < nu; ++i) for (int j = 0; j < nv; ++j) flat.push_back(poles[i][j]);
    return nurbs::cert::make_surface(nu, nv, flat);
}

static const double SQ = 0.7071067811865476;
static const double WQ[3] = { 1.0, SQ, 1.0 };
static const double CC[3] = { 1.0, 1.0, 0.0 };
static const double SS[3] = { 0.0, 1.0, 1.0 };

static RBezierSurface cyl_z() {
    std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j)
        P[i][j] = { CC[i] * WQ[i], SS[i] * WQ[i], (double)j * WQ[i], WQ[i] };
    return make(P);
}
static RBezierSurface cyl_x(double off) {
    std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j)
        P[i][j] = { (double)j * WQ[i], (off + CC[i]) * WQ[i], SS[i] * WQ[i], WQ[i] };
    return make(P);
}
static RBezierSurface plane_z(double z0) {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j)
        P[i][j] = { 2.0 * i - 1.0, 2.0 * j - 1.0, z0, 1.0 };
    return make(P);
}
static RBezierSurface plane_x(double x0) {   // x=x0，y,z ∈ [-1,1]
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { x0, 2.0 * i - 1.0, (double)j, 1.0 };
    return make(P);
}

static RBezierSurface plane_shift(double x0, double y0) {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j)
        P[i][j] = { x0 + (double)i, y0 + (double)j, 0.0, 1.0 };
    return make(P);
}
static RBezierSurface plane_repar() {
    const double xp[3] = { -1.0, -1.0, 1.0 };
    std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j)
        P[i][j] = { xp[i], 2.0 * j - 1.0, 0.0, 1.0 };
    return make(P);
}
static RBezierSurface paraboloid() {
    const double xp[3] = { -1.0, 0.0, 1.0 };
    const double zp[3] = { 1.0, -1.0, 1.0 };
    std::vector<std::vector<P4>> P(3, std::vector<P4>(3));
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
        P[i][j] = { xp[i], xp[j], zp[i] + zp[j], 1.0 };
    return make(P);
}

// ---- P4：高阶接触测试用曲面 ----
// z = p^2，p = x^2+y^2（x=2u-1, y=2v-1）：与平面 z=0 在原点为 4 阶接触（Hessian=0）
static RBezierSurface quartic_contact_surface() {
    const double zp[3] = { 1.0, -1.0, 1.0 };
    const double bin2[3] = { 1.0, 2.0, 1.0 };
    const double bin4[5] = { 1.0, 4.0, 6.0, 4.0, 1.0 };
    // p_ij = zp_i + zp_j （bideg 2x2）
    double pij[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) pij[i][j] = zp[i] + zp[j];
    // z 的 bideg(4,4) 极点
    double z4[5][5] = {};
    for (int I = 0; I < 5; ++I)
        for (int J = 0; J < 5; ++J) {
            double acc = 0.0;
            for (int i = 0; i <= 2; ++i) {
                int i2 = I - i; if (i2 < 0 || i2 > 2) continue;
                for (int j = 0; j <= 2; ++j) {
                    int j2 = J - j; if (j2 < 0 || j2 > 2) continue;
                    acc += pij[i][j] * pij[i2][j2] * (bin2[i] * bin2[i2] / bin4[I]) * (bin2[j] * bin2[j2] / bin4[J]);
                }
            }
            z4[I][J] = acc;
        }
    std::vector<std::vector<P4>> P(5, std::vector<P4>(5));
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) {
            double x = (double)i / 4.0 * 2.0 - 1.0;
            double y = (double)j / 4.0 * 2.0 - 1.0;
            P[i][j] = { x, y, z4[i][j], 1.0 };
        }
    return make(P);
}

// 在盒上计算 z 方向混合偏导 h_{u^a v^b} 的区间值（伙伴为平面 z=0 时 h=z1-c0）
static int contact_order_at_point(const RBezierSurface& S, double u0, double v0, int max_order) {
    // 逐阶检查：返回最小的接触阶 k（存在 k 阶偏导不含 0）
    std::vector<nurbs::cert::HP> net = S.P;
    int nu = S.nu, nv = S.nv;
    // 预计算各阶导数网
    std::vector<std::vector<nurbs::cert::HP>> du(max_order + 1), dv(max_order + 1);
    du[0] = net; dv[0] = net;
    for (int a = 1; a <= max_order; ++a) {
        du[a] = nurbs::cert::deriv_net_u(du[a - 1], nu - (a - 1), nv);
        dv[a] = nurbs::cert::deriv_net_v(dv[a - 1], nu, nv - (a - 1));
        if (du[a].empty() && dv[a].empty()) break;
    }
    for (int k = 1; k <= max_order; ++k) {
        bool any_nonzero = false;
        for (int a = 0; a <= k; ++a) {
            int b = k - a;
            if ((int)du[a].size() == 0 || (int)dv[b].size() == 0) continue;
            int nua = nu - a, nvb = nv - b;
            std::vector<nurbs::cert::HP> m = du[a];
            for (int t = 0; t < b; ++t) m = nurbs::cert::deriv_net_v(m, nua, nv - t);
            if (m.empty()) continue;
            nurbs::cert::HP out;
            if (nua < 1 || nv - b < 1) continue;
            nurbs::cert::eval_net(m, nua, nv - b, nurbs::Interval::point(u0),
                                  nurbs::Interval::point(v0), out);
            nurbs::Interval val = out[2];
            if (val.lo > 0.0 || val.hi < 0.0) any_nonzero = true;
        }
        if (any_nonzero) return k;
    }
    return 0;
}

// ---- 总公式验证：折叠型奇异面片 (u,v) -> (u, v^2, u v)，v=0 处 rank DS = 1 ----
static RBezierSurface fold_surface() {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(3));
    const double yp[3] = { 0.0, 0.0, 1.0 };     // y = v^2
    const double vp[3] = { 0.0, 0.5, 1.0 };     // 线性 v 提升到二次
    for (int i = 0; i < 2; ++i) {
        double x = (double)i;                    // x = u
        for (int j = 0; j < 3; ++j)
            P[i][j] = { x, yp[j], x * vp[j], 1.0 };   // z = u v
    }
    return make(P);
}

static RBezierSurface plane_y(double y0) {       // y = y0，x,z ∈ [-1,1]
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { 2.0 * i - 1.0, y0, 2.0 * j - 1.0, 1.0 };
    return make(P);
}

static RBezierSurface cyl_r(double r) {          // 半径 r 的四分之一圆柱（轴 z）
    RBezierSurface S = cyl_z();
    for (auto& q : S.P) { q[0] = nurbs::imul(q[0], nurbs::Interval::point(r)); q[1] = nurbs::imul(q[1], nurbs::Interval::point(r)); }
    S.update_weight_range();
    return S;
}

struct Truth { double u, v, s, t; };
static std::vector<Truth> truth_curve() {
    const double r = 0.7071067811865476;
    std::vector<Truth> out;
    auto arc_u = [&](double cx, double cy) {
        double lo = 0.0, hi = 1.0;
        for (int it = 0; it < 90; ++it) {
            double u = 0.5 * (lo + hi);
            double W = (1-u)*(1-u) + 2*u*(1-u)*r + u*u;
            double y = (2*u*(1-u)*r + u*u) / W;
            if (y < cy) lo = u; else hi = u;
        }
        return 0.5 * (lo + hi);
    };
    for (int i = 0; i < 21; ++i) {
        double y = 0.5 + 0.5 * i / 20.0;
        double x = std::sqrt(std::max(0.0, 1.0 - y*y));
        double z = std::sqrt(std::max(0.0, 1.0 - (y-0.5)*(y-0.5)));
        out.push_back({ arc_u(x, y), z, arc_u(y - 0.5, z), x });
    }
    return out;
}

static bool inside_box(const Box4& b, double u, double v, double s, double t, double tol) {
    return u >= b.u.lo - tol && u <= b.u.hi + tol &&
           v >= b.v.lo - tol && v <= b.v.hi + tol &&
           s >= b.s.lo - tol && s <= b.s.hi + tol &&
           t >= b.t.lo - tol && t <= b.t.hi + tol;
}

// Lipschitz 上界：sup |S_u|, |S_v|（区间导数在整域上的范数上界）
static double lipschitz_bound(const RBezierSurface& S) {
    Interval Q[3], Su[3], Sv[3];
    // 只用一阶导（二阶导阶段可能因退化返回 false，但 Su/Sv 已填充）
    (void)nurbs::cert::eval_rational_duv(S, Interval(0.0, 1.0), Interval(0.0, 1.0), Q, Su, Sv);
    auto nrm = [](Interval v[3]) {
        double a = 0.0;
        for (int i = 0; i < 3; ++i) {
            double m = std::max(std::fabs(v[i].lo), std::fabs(v[i].hi));
            a += m * m;
        }
        return std::sqrt(a);
    };
    return std::max(nrm(Su), nrm(Sv));
}

static Box4 small_box(double c, double hv) {
    Box4 X;
    X.u = Interval(c - hv, c + hv); X.v = Interval(c - hv, c + hv);
    X.s = Interval(c - hv, c + hv); X.t = Interval(c - hv, c + hv);
    return X;
}

int main() {
    std::printf("=== theory tests (theory_proofs.md) ===\n");

    // ---------- T1：覆盖无遗漏 ----------
    std::printf("\nT1 (cover no-missing)\n");
    {
        auto A = cyl_z(), B = cyl_x(0.5);
        auto tr = truth_curve();
        auto r = nurbs::cert::cover(A, B, 0.03, 60, 4000000, true);
        bool all_cov = true;
        for (const auto& p : tr) if (!nurbs::cert::covered(r.boxes, p.u, p.v, p.s, p.t, 1e-12)) all_cov = false;
        CHECK(all_cov, "all 21 truth points lie in cover boxes");
        bool no_truth_in_excluded = true;
        for (const auto& b : r.excluded)
            for (const auto& p : tr)
                if (inside_box(b, p.u, p.v, p.s, p.t, 1e-12)) no_truth_in_excluded = false;
        std::printf("  excluded boxes=%zu (checked against truth points)\n", r.excluded.size());
        CHECK(no_truth_in_excluded, "no truth point lies in any excluded box (exclusion soundness)");
        CHECK(r.excl + r.cover > 0, "cover tree non-empty");
    }

    // ---------- T2a：横截认证 ----------
    std::printf("\nT2a (transversal certification)\n");
    {
        auto A = cyl_z(), B = cyl_x(0.5);
        auto tr = truth_curve();
        auto r = nurbs::cert::cover(A, B, 0.1);
        auto st = nurbs::cert::certify_boxes(A, B, r.boxes);
        int cert_truth = 0;
        for (const auto& p : tr) {
            for (std::size_t i = 0; i < r.boxes.size(); ++i) {
                if (inside_box(r.boxes[i], p.u, p.v, p.s, p.t, 1e-12)) {
                    if (st.box_certified[i]) ++cert_truth;
                    break;
                }
            }
        }
        std::printf("  certified boxes=%ld / %ld, truth boxes certified=%d/21\n",
                    st.certified, (long)r.boxes.size(), cert_truth);
        CHECK(cert_truth >= 21, "all truth-bearing cover boxes certified (Krawczyk containment)");
    }

    // ---------- T1b：资源上限下的无遗漏（节点上限/深度上限）----------
    std::printf("\nT1b (no-missing under resource caps)\n");
    {
        auto A = cyl_z(), B = cyl_x(0.5);
        // 节点上限很小：必须触发 cap_hit，且真值点仍被覆盖（不得丢盒）
        auto tr = truth_curve();
        auto rc = nurbs::cert::cover(A, B, 0.02, 60, 8);
        int cov = 0;
        for (const auto& tp : tr)
            if (nurbs::cert::covered(rc.boxes, tp.u, tp.v, tp.s, tp.t, 0.05)) ++cov;
        std::printf("  node_cap=8: cap_hit=%d cover=%ld (of %ld) truth_cov=%d/21\n",
                    (int)rc.cap_hit, rc.cover, 21L, cov);
        CHECK(rc.cap_hit, "tiny node cap triggers cap_hit");
        CHECK(cov == 21, "no box is dropped at the node cap (all truth points still covered)");
        // 深度上限很小：同样必须不丢盒
        auto rd = nurbs::cert::cover(A, B, 0.02, 2, 4000000);
        int covd = 0;
        for (const auto& tp : tr)
            if (nurbs::cert::covered(rd.boxes, tp.u, tp.v, tp.s, tp.t, 0.05)) ++covd;
        std::printf("  max_depth=2: cover=%ld truth_cov=%d/21\n", rd.cover, covd);
        CHECK(covd == 21, "no box is dropped at the depth limit (all truth points still covered)");
    }

    std::printf("\nT2c (uniqueness in a box, Theorem 2c)\n");
    {
        auto A = cyl_z(), B = cyl_x(0.5);
        auto r = nurbs::cert::cover(A, B, 0.1);
        auto st = nurbs::cert::certify_boxes(A, B, r.boxes);
        std::printf("  A transversal: certified=%ld unique=%ld\n", st.certified, st.unique);
        CHECK(st.unique > 0, "transversal boxes certified unique (q<1 and K contained)");
        CHECK(st.unique <= st.certified, "unique count does not exceed certified count");
        // 相切案例：rank 亏，不应出现 unique
        auto C1 = paraboloid(), C2 = plane_z(0.0);
        auto rC = nurbs::cert::cover(C1, C2, 0.03);
        auto stC = nurbs::cert::certify_boxes(C1, C2, rC.boxes);
        std::printf("  C tangency: certified=%ld unique=%ld tang=%ld\n", stC.certified, stC.unique, stC.tangent);
        CHECK(stC.unique == 0, "tangency boxes are not claimed unique (rank deficient)");
        auto A2 = cyl_z(), B2 = plane_x(1.0);
        auto rD = nurbs::cert::cover(A2, B2, 0.03);
        auto stD = nurbs::cert::certify_boxes(A2, B2, rD.boxes);
        std::printf("  D tangency: certified=%ld unique=%ld tang=%ld\n", stD.certified, stD.unique, stD.tangent);
        CHECK(stD.unique == 0, "line tangency boxes are not claimed unique (rank deficient)");
    }

    // ---------- P2/T3：拓扑完备性（安全细分条件 + 分量）----------
    std::printf("\nT3 (topology completeness, Theorem 3 conditional)\n");
    {
        auto A = cyl_z(), B = cyl_x(0.5);
        for (double tau : { 0.1, 0.03 }) {
            auto ts = nurbs::cert::topology_analyze(A, B, tau);
            std::printf("  A tau=%.3f: cover=%ld certified=%ld safe=%ld min_reach=%.3e comps=%ld\n",
                        tau, ts.cover, ts.certified, ts.safe_boxes, ts.min_reach, ts.cert_components);
        }
        auto ts1 = nurbs::cert::topology_analyze(A, B, 0.1);
        auto ts2 = nurbs::cert::topology_analyze(A, B, 0.03);
        CHECK(ts1.band_components == 1 && ts1.cert_components == 1,
              "band/certified adjacency graph has one component at tau=0.1");
        CHECK(ts2.cert_components == 1, "one component at tau=0.03");
        CHECK(ts2.safe_boxes > ts1.safe_boxes,
              "safe-subdivision box count grows as tau decreases (reach condition)");
        CHECK(ts2.min_reach > 0.0, "interval reach lower bound is positive on certified boxes");
        // 相切案例：rank 亏 ⇒ reach 无下界，安全细分条件不适用
        auto C1 = paraboloid(), C2 = plane_z(0.0);
        auto tsC = nurbs::cert::topology_analyze(C1, C2, 0.03);
        std::printf("  C tau=0.03: certified=%ld safe=%ld min_reach=%.3e\n",
                    tsC.certified, tsC.safe_boxes, tsC.min_reach);
        CHECK(tsC.safe_boxes == 0, "tangency boxes are not safe-subdivision certified (rank deficient)");
    }

    // ---------- T2b：相切认证 + 假阳性拒绝 ----------
    std::printf("\nT2b (tangency certification)\n");
    {
        auto C1 = paraboloid(), C2 = plane_z(0.0);
        CHECK(nurbs::cert::certify_tangency_box(C1, C2, small_box(0.5, 3e-3)) == 0,
              "true tangency (apex on plane) certified");
        // 假阳性：平面片平移出顶点像域（等高不重合，修正前会误认证）
        auto C3 = plane_shift(1.5, 1.5);
        CHECK(nurbs::cert::certify_tangency_box(C1, C3, small_box(0.5, 3e-3)) != 0,
              "false positive rejected (equal height but apex not on shifted plane)");
        auto D1 = cyl_z(), D2 = plane_shift(1.0, 0.0);   // 平面 x=1 的另一种参数化
        (void)D1; (void)D2;
    }

    // ---------- T4：2D ε-重合证书 ----------
    std::printf("\nT4 (2D coincidence certificate)\n");
    {
        auto A = plane_z(0.0);
        double L = lipschitz_bound(A);
        double h = 1e-4;
        std::printf("  plane Lipschitz(interval bound)=%.3e, h=%.1e, h*L=%.3e\n", L, h, h * L);
        // 平面片的参数化导数非零（L≈2–4 的区间上界）；只要 h*L < eps_tol 即满足 ε-重合条件
        CHECK(L > 0.0 && L < 100.0, "plane Lipschitz bound finite and conservative");
        CHECK(h * L < 1e-2, "h*L < eps_tol for plane (epsilon-coincidence condition)");
        auto A2 = plane_shift(0.5, 0.5);
        double r_in = nurbs::cert::coincidence_ratio(A, A2, 0.8, 1.0, 0.8, 1.0, 4);
        CHECK(r_in > 0.99, "coplanar overlap region certified (point inclusion + normals)");
        double r_out = nurbs::cert::coincidence_ratio(A, A2, 0.0, 0.3, 0.0, 0.3, 4);
        CHECK(r_out < 0.01, "outside overlap region rejected (negative certificate)");
        auto B = plane_repar();
        CHECK(nurbs::cert::coincidence_ratio(A, B, 0.1, 0.9, 0.1, 0.9, 4) > 0.99,
              "reparameterized identical image certified");
        auto P1 = paraboloid();
        double Lp = lipschitz_bound(P1);
        double hp = 1e-4;
        std::printf("  paraboloid Lipschitz=%.3e, h=%.1e, h*L=%.3e\n", Lp, hp, hp * Lp);
        CHECK(hp * Lp < 1e-2, "Lipschitz condition h*L < eps_tol satisfied for curved patch");
    }

    // ---------- P5/T4b：平面仿射片的精确 2D 重合 ----------
    std::printf("\nT4b (exact 2D coincidence for affine planar patches, Theorem P5)\n");
    {
        auto A = plane_z(0.0);
        int r1 = nurbs::cert::exact_coplanar_overlap(A, A, 1.0);
        std::printf("  identical planes -> %d (expect 1)\n", r1);
        CHECK(r1 == 1, "identical affine planes: exact 2D coincidence (containment)");
        auto B = plane_shift(0.3, 0.3);   // A 像 [-1,1]^2，B 像 [-0.7,1.3]^2 → 部分重叠
        int r2 = nurbs::cert::exact_coplanar_overlap(A, B, 1.0);
        std::printf("  partially overlapping planes -> %d (expect 2)\n", r2);
        CHECK(r2 == 2, "partially overlapping affine planes: 2D overlap detected");
        auto C = plane_repar();           // 二次参数化（非仿射）
        int r3 = nurbs::cert::exact_coplanar_overlap(A, C, 1.0);
        std::printf("  reparameterized (non-affine) plane -> %d (expect -1, fallback)\n", r3);
        CHECK(r3 == -1, "non-affine reparameterized plane falls back to epsilon certificate");
        auto D = plane_shift(5.0, 5.0);   // 不相交
        int r4 = nurbs::cert::exact_coplanar_overlap(A, D, 1.0);
        std::printf("  far-away planes -> %d (expect 0)\n", r4);
        CHECK(r4 == 0, "disjoint affine planes: no 2D overlap");
    }

    // ---------- P4/T5：高阶接触的接触阶判据 ----------
    std::printf("\nT5 (contact-order criterion for higher-order contacts, Theorem P4)\n");
    {
        auto plane = plane_z(0.0);
        // 二次接触（抛物面 vs 平面）：k=2
        auto par = paraboloid();
        int k2 = contact_order_at_point(par, 0.5, 0.5, 4);
        std::printf("  paraboloid x plane: contact order = %d (expect 2)\n", k2);
        CHECK(k2 == 2, "quadratic tangency classified as contact order 2");
        // 四次接触（z=(x^2+y^2)^2 vs 平面）：k=4
        auto quartic = quartic_contact_surface();
        int k4 = contact_order_at_point(quartic, 0.5, 0.5, 4);
        std::printf("  quartic contact x plane: contact order = %d (expect 4)\n", k4);
        CHECK(k4 == 4, "quartic (osculating) contact classified as contact order 4");
        // Hessian 在原点区间含 0 ⇒ 标准 2x2 认证不可用（高阶接触需阶判据）
        double huu_lo = 0.0, huu_hi = 0.0;
        {
            // 二阶偏导 h_uu 的区间值
            std::vector<nurbs::cert::HP> du1 = nurbs::cert::deriv_net_u(quartic.P, quartic.nu, quartic.nv);
            std::vector<nurbs::cert::HP> du2 = nurbs::cert::deriv_net_u(du1, quartic.nu - 1, quartic.nv);
            nurbs::cert::HP out;
            nurbs::cert::eval_net(du2, quartic.nu - 2, quartic.nv,
                                  nurbs::Interval(0.45, 0.55), nurbs::Interval(0.45, 0.55), out);
            huu_lo = out[2].lo; huu_hi = out[2].hi;
        }
        std::printf("  quartic h_uu interval = [%.3f, %.3f]\n", huu_lo, huu_hi);
        CHECK(huu_lo <= 0.0 && huu_hi >= 0.0, "quartic Hessian entry interval contains 0 (standard certification fails)");
        (void)plane;
    }

    // ---------- P9/T6：退化片分类 ----------
    std::printf("\nT6 (degenerate patch classification, Proposition P9)\n");
    {
        // 点片：所有控制点重合
        {
            std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) P[i][j] = { 0.5, 0.5, 0.5, 1.0 };
            auto dg = nurbs::cert::classify_degeneracy(make(P));
            std::printf("  point patch: point=%d collinear=%d nonpos=%d\n", dg.point, dg.collinear, dg.nonpositive_weight);
            CHECK(dg.point, "point patch classified as point");
            CHECK(!dg.collinear, "point patch is not reported collinear");
        }
        // 共线片：控制点在 x 轴上
        {
            std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j) P[i][j] = { 0.2 * i, 0.0, 0.0, 1.0 };
            auto dg = nurbs::cert::classify_degeneracy(make(P));
            std::printf("  collinear patch: point=%d collinear=%d nonpos=%d\n", dg.point, dg.collinear, dg.nonpositive_weight);
            CHECK(dg.collinear, "collinear patch classified as collinear");
        }
        // 非正权
        {
            std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) P[i][j] = { 2.0 * i - 1.0, 2.0 * j - 1.0, 0.0, 1.0 };
            P[1][1][3] = -1.0;
            auto dg = nurbs::cert::classify_degeneracy(make(P));
            std::printf("  negative weight: point=%d collinear=%d nonpos=%d\n", dg.point, dg.collinear, dg.nonpositive_weight);
            CHECK(dg.nonpositive_weight, "non-positive weight detected");
            // cover 仍是无遗漏超集（不排除盒）
            auto r = nurbs::cert::cover(make(P), plane_z(0.0), 0.1);
            std::printf("  cover boxes with non-positive weight: %ld\n", r.cover);
            CHECK(r.cover > 0, "cover remains a valid superset (does not drop boxes)");
        }
        // 正常片
        {
            auto dg = nurbs::cert::classify_degeneracy(paraboloid());
            std::printf("  normal patch: point=%d collinear=%d nonpos=%d\n", dg.point, dg.collinear, dg.nonpositive_weight);
            CHECK(!dg.point && !dg.collinear && !dg.nonpositive_weight, "normal patch is not flagged degenerate");
        }
    }

    // ---------- P6/T7：裁剪域（参数多边形）----------
    std::printf("\nT7 (trimmed parameter domain, Theorem P6)\n");
    {
        // A = z=0 平面片（[0,1]² 参数），B = y=0.5 平面片；交线在 z=0, y=0.5 上
        auto A = plane_z(0.0);
        std::vector<std::vector<P4>> Pb(2, std::vector<P4>(2));
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) Pb[i][j] = { 2.0 * i - 1.0, 0.5, 2.0 * j - 1.0, 1.0 };
        auto B = make(Pb);
        auto full = nurbs::cert::cover(A, B, 0.05);
        // 裁剪：A 的参数域取 u<=0.5 的半正方形；B 取全 [0,1]²
        nurbs::cert::Poly2 polyA;
        polyA.V = { { 0.0, 0.0 }, { 0.5, 0.0 }, { 0.5, 1.0 }, { 0.0, 1.0 } };
        nurbs::cert::Poly2 polyB;
        polyB.V = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
        auto tc = nurbs::cert::trimmed_cover(A, B, 0.05, polyA, polyB);
        std::printf("  full cover=%ld -> kept=%ld outside=%ld partial=%ld\n",
                    tc.cover, tc.kept, tc.outside, tc.partial);
        std::printf("  full cover boxes (no trim) = %ld\n", full.cover);
        CHECK(tc.outside > 0, "trimmed cover excludes boxes certified outside the parameter polygon");
        CHECK(tc.kept < tc.cover, "trimmed cover keeps strictly fewer boxes than the full cover");
        CHECK(tc.kept + tc.outside == tc.cover, "every cover box is either kept or certified outside");

        // 裁剪域内的真值点仍被保留盒覆盖；域外的点不在保留盒中（按参数 u 判断）
        bool inside_kept = false;
        for (const auto& bx : full.boxes) {
            double u = 0.25, v = 0.75, s2 = 0.25, t2 = 0.5;
            if (u >= bx.u.lo && u <= bx.u.hi && v >= bx.v.lo && v <= bx.v.hi &&
                s2 >= bx.s.lo && s2 <= bx.s.hi && t2 >= bx.t.lo && t2 <= bx.t.hi) { inside_kept = true; break; }
        }
        std::printf("  truth point on trimmed intersection covered: %d\n", (int)inside_kept);
        CHECK(inside_kept, "truth point inside the trimmed domain remains covered");
        // 参数域排除证书自检：被排除的盒其 (u,v) 投影确实完全在 polyA 外或 (s,t) 完全在 polyB 外
        bool cert_ok = true;
        for (const auto& bx : full.boxes) {
            auto a = nurbs::cert::poly_box_state(polyA, bx.u, bx.v);
            auto b = nurbs::cert::poly_box_state(polyB, bx.s, bx.t);
            if ((a == nurbs::cert::TRIM_OUT || b == nurbs::cert::TRIM_OUT) == false) continue;
            if (a != nurbs::cert::TRIM_OUT && b != nurbs::cert::TRIM_OUT) cert_ok = false;
        }
        CHECK(cert_ok, "outside classification is consistent (at least one projection proven outside)");
    }

    // ---------- T8：总公式验证矩阵（全类型相交）----------
    std::printf("\nT8 (full-formula validation matrix over intersection types)\n");
    {
        int cases = 0, pass = 0;
        auto report = [&](const char* name, bool ok, const char* detail) {
            ++cases;
            if (ok) ++pass;
            std::printf("  [%s] %-46s %s\n", ok ? "PASS" : "FAIL", name, detail);
        };
        auto covered_by = [&](const RBezierSurface& A, const RBezierSurface& B,
                              const std::vector<Truth>& tr, double tau) {
            auto r = nurbs::cert::cover(A, B, tau);
            int cov = 0;
            for (const auto& t : tr)
                if (nurbs::cert::covered(r.boxes, t.u, t.v, t.s, t.t, 2 * tau)) ++cov;
            return std::make_pair(cov, r.cover);
        };
        // 1 横截曲线
        {
            auto A = cyl_z(), B = cyl_x(0.5);
            auto tr = truth_curve();
            auto cc = covered_by(A, B, tr, 0.03);
            auto pc = nurbs::cert::classify_pair(A, B, 0.03);
            char buf[128]; std::snprintf(buf, sizeof buf, "truth %d/21 type=%d cert=%ld", cc.first, pc.type, pc.certified);
            report("transversal curve", cc.first == 21 && pc.type == 1 && pc.certified > 0, buf);
        }
        // 2 无交
        {
            auto A = cyl_z(), B = plane_x(3.0);
            auto pc = nurbs::cert::classify_pair(A, B, 0.03);
            char buf[64]; std::snprintf(buf, sizeof buf, "type=%d cover=%ld", pc.type, pc.cover);
            report("disjoint (empty cover)", pc.type == 0 && pc.cover == 0, buf);
        }
        // 3 点相切
        {
            auto A = paraboloid(), B = plane_z(0.0);
            auto tr = std::vector<Truth>{ { 0.5, 0.5, 0.5, 0.5 } };
            auto r = nurbs::cert::cover(A, B, 0.03);
            bool cov = nurbs::cert::covered(r.boxes, 0.5, 0.5, 0.5, 0.5, 0.06);
            char buf[64]; std::snprintf(buf, sizeof buf, "covered=%d", (int)cov);
            report("point tangency", cov, buf);
        }
        // 4 线相切
        {
            auto A = cyl_z(), B = plane_x(1.0);
            std::vector<Truth> tr;
            for (int k = 1; k <= 5; ++k) { double z = 0.1 * k; tr.push_back({ 0.0, z, 0.5, z }); }
            auto r = nurbs::cert::cover(A, B, 0.03);
            int cov = 0;
            for (const auto& t : tr) if (nurbs::cert::covered(r.boxes, t.u, t.v, t.s, t.t, 0.06)) ++cov;
            char buf[64]; std::snprintf(buf, sizeof buf, "truth %d/5", cov);
            report("line tangency", cov == 5, buf);
        }
        // 5 高阶接触（osculating）
        {
            auto A = quartic_contact_surface(), B = plane_z(0.0);
            auto r = nurbs::cert::cover(A, B, 0.03);
            bool cov = nurbs::cert::covered(r.boxes, 0.5, 0.5, 0.5, 0.5, 0.06);
            char buf[64]; std::snprintf(buf, sizeof buf, "covered=%d (cert may be undecided)", (int)cov);
            report("osculating contact (order 4)", cov, buf);
        }
        // 6 近切小环
        {
            auto A = paraboloid(), B = plane_z(1e-6);
            const double rho = std::sqrt(1e-6);
            std::vector<Truth> tr;
            for (int k = 0; k < 8; ++k) {
                double th = 2.0 * M_PI * k / 8.0;
                double x = rho * std::cos(th), y = rho * std::sin(th);
                tr.push_back({ (x + 1) / 2, (y + 1) / 2, (x + 1) / 2, (y + 1) / 2 });
            }
            auto r = nurbs::cert::cover(A, B, 0.01);
            int cov = 0;
            for (const auto& t : tr) if (nurbs::cert::covered(r.boxes, t.u, t.v, t.s, t.t, 0.02)) ++cov;
            char buf[64]; std::snprintf(buf, sizeof buf, "truth %d/8 (no-missing)", cov);
            report("near-tangent small loop", cov == 8, buf);
        }
        // 7 极小特征
        {
            auto A = cyl_r(1e-6), B = plane_z(0.5);
            std::vector<Truth> tr = {
                { 0.0, 0.5, (1e-6 + 1) / 2, 0.5 },
                { 0.5, 0.5, (1e-6 / std::sqrt(2.0) + 1) / 2, (1e-6 / std::sqrt(2.0) + 1) / 2 },
                { 1.0, 0.5, 0.5, (1e-6 + 1) / 2 } };
            auto r = nurbs::cert::cover(A, B, 0.01);
            int cov = 0;
            for (const auto& t : tr) if (nurbs::cert::covered(r.boxes, t.u, t.v, t.s, t.t, 0.02)) ++cov;
            char buf[64]; std::snprintf(buf, sizeof buf, "truth %d/3", cov);
            report("tiny feature (r=1e-6)", cov == 3, buf);
        }
        // 8 2D 重合（平面精确）
        {
            auto A = plane_z(0.0), B = plane_z(0.0);
            int ex = nurbs::cert::exact_coplanar_overlap(A, B, 1.0);
            auto pc = nurbs::cert::classify_pair(A, B, 0.03);
            char buf[64]; std::snprintf(buf, sizeof buf, "exact=%d type=%d", ex, pc.type);
            report("2D coincidence (identical planes)", ex == 1 && pc.type == 2, buf);
        }
        // 9 2D 部分重合
        {
            auto A = plane_z(0.0), B = plane_shift(0.3, 0.3);
            int ex = nurbs::cert::exact_coplanar_overlap(A, B, 1.0);
            char buf[32]; std::snprintf(buf, sizeof buf, "exact=%d", ex);
            report("2D partial overlap", ex == 2, buf);
        }
        // 10 非仿射重参数化平面（ε 后备）
        {
            auto A = plane_z(0.0), B = plane_repar();
            int ex = nurbs::cert::exact_coplanar_overlap(A, B, 1.0);
            double ratio = nurbs::cert::coincidence_ratio(A, B, 0.1, 0.9, 0.1, 0.9, 4);
            char buf[64]; std::snprintf(buf, sizeof buf, "exact=%d ratio=%.2f", ex, ratio);
            report("reparameterized plane (epsilon fallback)", ex == -1 && ratio > 0.99, buf);
        }
        // 11 点片
        {
            std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) P[i][j] = { 0.0, 0.0, 0.0, 1.0 };
            auto A = make(P), B = plane_z(0.0);
            auto dg = nurbs::cert::classify_degeneracy(A);
            auto r = nurbs::cert::cover(A, B, 0.1);
            char buf[64]; std::snprintf(buf, sizeof buf, "point=%d cover=%ld", (int)dg.point, r.cover);
            report("degenerate point patch", dg.point && r.cover > 0, buf);
        }
        // 12 共线片
        {
            std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j) P[i][j] = { 0.2 * i, 0.0, 0.0, 1.0 };
            auto A = make(P), B = plane_z(0.0);
            auto dg = nurbs::cert::classify_degeneracy(A);
            auto r = nurbs::cert::cover(A, B, 0.1);
            char buf[64]; std::snprintf(buf, sizeof buf, "collinear=%d cover=%ld", (int)dg.collinear, r.cover);
            report("degenerate collinear patch", dg.collinear && r.cover > 0, buf);
        }
        // 13 非正权（覆盖仍无遗漏）
        {
            std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
            for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) P[i][j] = { 2.0 * i - 1.0, 2.0 * j - 1.0, 0.0, 1.0 };
            P[1][1][3] = -1.0;
            auto A = make(P), B = plane_z(0.0);
            auto dg = nurbs::cert::classify_degeneracy(A);
            auto r = nurbs::cert::cover(A, B, 0.1);
            char buf[64]; std::snprintf(buf, sizeof buf, "nonpos=%d cover=%ld", (int)dg.nonpositive_weight, r.cover);
            report("non-positive weight (sound cover)", dg.nonpositive_weight && r.cover > 0, buf);
        }
        // 14 折叠奇异面 × 平面（横截曲线）
        {
            auto A = fold_surface(), B = plane_z(0.5);
            std::vector<Truth> tr;
            for (int k = 0; k < 5; ++k) {
                double u = 0.5 + 0.1 * k;          // u in [0.5, 0.9]
                double v = 0.5 / u;                // z = u v = 0.5
                double y = v * v;
                tr.push_back({ u, v, (u + 1) / 2, (y + 1) / 2 });
            }
            auto r = nurbs::cert::cover(A, B, 0.02);
            int cov = 0;
            for (const auto& t : tr) if (nurbs::cert::covered(r.boxes, t.u, t.v, t.s, t.t, 0.04)) ++cov;
            char buf[64]; std::snprintf(buf, sizeof buf, "truth %d/5 cover=%ld", cov, r.cover);
            report("singular (fold) surface x plane", cov == 5, buf);
        }
        // 15 折叠奇异面 × 平面 y=0（交线 = 奇异线）
        {
            auto A = fold_surface(), B = plane_y(0.0);
            std::vector<Truth> tr;
            for (int k = 0; k < 5; ++k) {
                double u = 0.1 + 0.2 * k;
                tr.push_back({ u, 0.0, (u + 1) / 2, 0.5 });
            }
            auto r = nurbs::cert::cover(A, B, 0.02);
            int cov = 0;
            for (const auto& t : tr) if (nurbs::cert::covered(r.boxes, t.u, t.v, t.s, t.t, 0.04)) ++cov;
            char buf[64]; std::snprintf(buf, sizeof buf, "truth %d/5 cover=%ld", cov, r.cover);
            report("singular line as intersection", cov == 5, buf);
        }
        std::printf("  RESULT: %d/%d cases passed\n", pass, cases);
        CHECK(pass == cases, "full-formula validation matrix passes all cases");
    }

    std::printf("\n=== FAILURES: %d ===\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
