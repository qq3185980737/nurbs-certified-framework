// ============================================================================
// bench_cert_ssi.cpp — M3/M4：一般 NURBS（有理 Bezier 片）曲面求交的认证覆盖 + 横截认证
//   M3：4D 自适应细分 + 凸包包络排除 → 覆盖盒带 ⊇ Σ（定理 1 实现）
//   M4：覆盖盒的 3x4 广义逆 Krawczyk 认证 + Gershgorin 秩判定（定理 2 横截部分）
// 情形：横截交线 / 无交 / 0D 点相切 / 1D 线相切 / 近切小环
// ============================================================================
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <string>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;

using P4 = std::array<double, 4>;   // 齐次 (wx, wy, wz, w)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static RBezierSurface make(const std::vector<std::vector<P4>>& poles) {
    int nu = (int)poles.size(), nv = (int)poles[0].size();
    std::vector<P4> flat;
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) flat.push_back(poles[i][j]);
    return nurbs::cert::make_surface(nu, nv, flat);
}

static const double SQ = 0.7071067811865476;
static const double WQ[3] = { 1.0, SQ, 1.0 };
static const double CC[3] = { 1.0, 1.0, 0.0 };
static const double SS[3] = { 0.0, 1.0, 1.0 };

static RBezierSurface cyl_z() {
    std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) {
            double z = (double)j;
            P[i][j] = { CC[i] * WQ[i], SS[i] * WQ[i], z * WQ[i], WQ[i] };
        }
    return make(P);
}
static RBezierSurface cyl_x(double off) {
    std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) {
            double x = (double)j;
            P[i][j] = { x * WQ[i], (off + CC[i]) * WQ[i], SS[i] * WQ[i], WQ[i] };
        }
    return make(P);
}
static RBezierSurface plane_z(double z0) {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { 2.0 * i - 1.0, 2.0 * j - 1.0, z0, 1.0 };
    return make(P);
}
static RBezierSurface plane_x(double x0) {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { x0, 2.0 * i - 1.0, (double)j, 1.0 };
    return make(P);
}
static RBezierSurface paraboloid() {
    const double xp[3] = { -1.0, 0.0, 1.0 };
    const double zp[3] = { 1.0, -1.0, 1.0 };
    std::vector<std::vector<P4>> P(3, std::vector<P4>(3));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            P[i][j] = { xp[i], xp[j], zp[i] + zp[j], 1.0 };
    return make(P);
}

static RBezierSurface plane_shift(double x0, double y0) {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { x0 + (double)i, y0 + (double)j, 0.0, 1.0 };
    return make(P);
}
// 与 plane_z(0) 同像、不同参数化：x = 2s^2 - 1（极点 [-1,-1,1]）
static RBezierSurface plane_repar() {
    const double xp[3] = { -1.0, -1.0, 1.0 };
    std::vector<std::vector<P4>> P(3, std::vector<P4>(2));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { xp[i], 2.0 * j - 1.0, 0.0, 1.0 };
    return make(P);
}

struct Truth { double u, v, s, t; };

// 有理二次 Bezier 四分之一圆弧的参数反演：
//   (x(u), y(u)) = (Nx/W, Ny/W)，W=(1-u)^2+2u(1-u)r+u^2，r=sqrt(2)/2
//   注意：u 不是角度参数（atan2 得到的角度参数会带来 ~1e-2 的参数误差）
static double arc_u(double cx, double cy) {
    const double r = 0.7071067811865476;
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 90; ++it) {
        double u = 0.5 * (lo + hi);
        double W = (1 - u) * (1 - u) + 2 * u * (1 - u) * r + u * u;
        double y = (2 * u * (1 - u) * r + u * u) / W;
        if (y < cy) lo = u; else hi = u;
    }
    return 0.5 * (lo + hi);
}

static std::vector<Truth> truth_curve() {
    std::vector<Truth> out;
    for (int i = 0; i < 21; ++i) {
        double y = 0.5 + 0.5 * i / 20.0;
        double x = std::sqrt(std::max(0.0, 1.0 - y * y));
        double z = std::sqrt(std::max(0.0, 1.0 - (y - 0.5) * (y - 0.5)));
        double u1 = arc_u(x, y);            // S1 圆弧参数（cos 部分 x，sin 部分 y）
        double s1 = arc_u(y - 0.5, z);      // S2 圆弧参数（cos 部分 y-off，sin 部分 z）
        out.push_back({ u1, z, s1, x });
    }
    return out;
}
static std::vector<Truth> truth_tangent_line() {
    std::vector<Truth> out;
    for (double z : { 0.1, 0.3, 0.5, 0.7, 0.9 }) out.push_back({ 0.0, z, 0.5, z });
    return out;
}
static std::vector<Truth> truth_small_loop(double h) {
    std::vector<Truth> out;
    double r = std::sqrt(h);
    for (int k = 0; k < 12; ++k) {
        double th = 2.0 * M_PI * k / 12.0;
        double x = r * std::cos(th), y = r * std::sin(th);
        out.push_back({ (x + 1.0) / 2.0, (y + 1.0) / 2.0, (x + 1.0) / 2.0, (y + 1.0) / 2.0 });
    }
    return out;
}

static void run_case(const char* name, const RBezierSurface& A, const RBezierSurface& B,
                     const std::vector<Truth>& truth) {
    std::printf("\n=== %s ===\n", name);
    for (double tau : { 0.1, 0.03, 0.01 }) {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto r = nurbs::cert::cover(A, B, tau);
        auto st = nurbs::cert::certify_boxes(A, B, r.boxes);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        int cov = 0, cert = 0;
        std::string bad;
        int idx = 0;
        for (const Truth& p : truth) {
            for (std::size_t i = 0; i < r.boxes.size(); ++i) {
                const Box4& b = r.boxes[i];
                const double tol = 1e-12;
                if (p.u >= b.u.lo - tol && p.u <= b.u.hi + tol &&
                    p.v >= b.v.lo - tol && p.v <= b.v.hi + tol &&
                    p.s >= b.s.lo - tol && p.s <= b.s.hi + tol &&
                    p.t >= b.t.lo - tol && p.t <= b.t.hi + tol) {
                    ++cov;
                    if (st.box_certified[i]) ++cert;
                    else bad += std::to_string(idx) + " ";
                    break;
                }
            }
            ++idx;
        }
        std::printf("  tau=%.3f | cover=%5ld | rank3=%5ld cert=%5ld uniq=%5ld tang=%5ld | fail[zero=%ld inv=%ld rank=%ld K=%ld] | truth %2d/%2d cert %2d | %6.1f ms%s\n",
                    tau, r.cover, st.rank3, st.certified, st.unique, st.tangent,
                    st.fail[1], st.fail[2], st.fail[3], st.fail[4],
                    cov, (int)truth.size(), cert, ms, r.cap_hit ? " [CAP]" : "");
        if (!bad.empty()) std::printf("      uncertified truth idx: %s\n", bad.c_str());
    }
}

int main() {
    std::printf("M3/M4 certified cover + transversal certification\n");
    run_case("A 1D curve (two non-graph cylinders)", cyl_z(), cyl_x(0.5), truth_curve());
    run_case("B disjoint (offset 3)", cyl_z(), cyl_x(3.0), {});
    run_case("C 0D tangency (paraboloid vs plane z=0)", paraboloid(), plane_z(0.0),
             { {0.5, 0.5, 0.5, 0.5} });
    run_case("D 1D tangency (cylinder vs plane x=1)", cyl_z(), plane_x(1.0), truth_tangent_line());
    run_case("H near-tangent small loop h=1e-6", paraboloid(), plane_z(1e-6), truth_small_loop(1e-6));
    // ---- 单点诊断：真值点处的 F、J 与认证码 ----
    {
        std::printf("\n=== dump: F / J at truth points ===\n");
        auto A = cyl_z(); auto B = cyl_x(0.5);
        auto tr = truth_curve();
        double hv[4] = { 1e-2, 3e-3, 1e-3, 3e-4 };
        for (int idx = 1; idx <= 2; ++idx) {
            const Truth& p = tr[5 * idx];
            nurbs::Interval U = nurbs::Interval::point(p.u);
            nurbs::Interval V = nurbs::Interval::point(p.v);
            nurbs::Interval Sx = nurbs::Interval::point(p.s);
            nurbs::Interval T = nurbs::Interval::point(p.t);
            nurbs::Interval S1[3], S1u[3], S1v[3], S2[3], S2u[3], S2v[3];
            bool ok1 = nurbs::cert::eval_rational_duv(A, U, V, S1, S1u, S1v);
            bool ok2 = nurbs::cert::eval_rational_duv(B, Sx, T, S2, S2u, S2v);
            std::printf("  A truth[%2d] ok=%d/%d S1=(%.6f,%.6f,%.6f) S2=(%.6f,%.6f,%.6f)\n",
                        5 * idx, (int)ok1, (int)ok2,
                        S1[0].mid(), S1[1].mid(), S1[2].mid(),
                        S2[0].mid(), S2[1].mid(), S2[2].mid());
            std::printf("      F=(%.3e,%.3e,%.3e)\n",
                        (S1[0].mid() - S2[0].mid()), (S1[1].mid() - S2[1].mid()), (S1[2].mid() - S2[2].mid()));
            std::printf("      dS1u=(%.3f,%.3f,%.3f) dS1v=(%.3f,%.3f,%.3f)\n",
                        S1u[0].mid(), S1u[1].mid(), S1u[2].mid(), S1v[0].mid(), S1v[1].mid(), S1v[2].mid());
            std::printf("      dS2s=(%.3f,%.3f,%.3f) dS2t=(%.3f,%.3f,%.3f)\n",
                        S2u[0].mid(), S2u[1].mid(), S2u[2].mid(), S2v[0].mid(), S2v[1].mid(), S2v[2].mid());
            std::printf("      codes hw=1e-2:%d 3e-3:%d 1e-3:%d 3e-4:%d\n",
                        nurbs::cert::probe_box(A, B, p.u, p.v, p.s, p.t, hv[0]),
                        nurbs::cert::probe_box(A, B, p.u, p.v, p.s, p.t, hv[1]),
                        nurbs::cert::probe_box(A, B, p.u, p.v, p.s, p.t, hv[2]),
                        nurbs::cert::probe_box(A, B, p.u, p.v, p.s, p.t, hv[3]));
        }
    }

    // ---- M5 相切探针 ----
    {
        std::printf("\n=== probe: tangency certification at truth boxes ===\n");
        auto C1 = paraboloid(); auto C2 = plane_z(0.0);
        auto P = [](double c, double hv) { nurbs::cert::Box4 X;
            X.u = nurbs::Interval(c - hv, c + hv); X.v = nurbs::Interval(c - hv, c + hv);
            X.s = nurbs::Interval(c - hv, c + hv); X.t = nurbs::Interval(c - hv, c + hv); return X; };
        {
            nurbs::cert::Box4 X;
            X.u = nurbs::Interval(0.46875, 0.5); X.v = nurbs::Interval(0.46875, 0.5);
            X.s = nurbs::Interval(0.45, 0.55);   X.t = nurbs::Interval(0.45, 0.55);
            std::printf("  C dyadic-edge box (apex on upper faces) -> tangency code %d\n",
                        nurbs::cert::certify_tangency_box(C1, C2, X));
        }
        for (double hv : { 1e-2, 3e-3, 1e-3 })
            std::printf("  C point tangency hw=%.0e -> tangency code %d\n",
                        hv, nurbs::cert::certify_tangency_box(C1, C2, P(0.5, hv)));
        auto D1 = cyl_z(); auto D2 = plane_x(1.0);
        {
            nurbs::cert::Box4 X;
            X.u = nurbs::Interval(0.0, 4e-2); X.v = nurbs::Interval(0.45, 0.55);
            X.s = nurbs::Interval(0.45, 0.55); X.t = nurbs::Interval(0.45, 0.55);
            std::printf("  D line tangency (u in [0,0.04]) -> tangency code %d\n",
                        nurbs::cert::certify_tangency_box(D1, D2, X));
        }
    }

    // ---- M6 拓扑三分 ----
    std::printf("\n=== M6 topology (band / certified / empty / undecided / components / endpoints) ===\n");
    {
        for (double tau : { 0.1, 0.03, 0.01 }) {
            auto sA = nurbs::cert::topology_analyze(cyl_z(), cyl_x(0.5), tau);
            std::printf("  A tau=%.3f | cover=%4ld cert=%3ld empty=%3ld undec=%3ld | bandComp=%2ld certComp=%2ld endpoints=%2ld maxdeg=%ld\n",
                        tau, sA.cover, sA.certified, sA.empty, sA.undecided,
                        sA.band_components, sA.cert_components, sA.endpoints, sA.max_degree);
        }
        for (double tau : { 0.1, 0.03 }) {
            auto sD = nurbs::cert::topology_analyze(cyl_z(), plane_x(1.0), tau);
            std::printf("  D tau=%.3f | cover=%4ld cert=%3ld empty=%3ld undec=%3ld | bandComp=%2ld certComp=%2ld endpoints=%2ld maxdeg=%ld\n",
                        tau, sD.cover, sD.certified, sD.empty, sD.undecided,
                        sD.band_components, sD.cert_components, sD.endpoints, sD.max_degree);
        }
    }

    // ---- M6.5 分支追踪 ----
    std::printf("\n=== M6.5 branch tracing (seeds / branches / endpoints / closed / length) ===\n");
    {
        for (double tau : { 0.1, 0.03 }) {
            auto tr = nurbs::cert::trace_branches(cyl_z(), cyl_x(0.5), tau);
            std::printf("  A tau=%.3f | seeds=%3ld branches=%ld endpoints=%ld closed=%ld len=%.3f\n",
                        tau, tr.seeds, tr.branches, tr.endpoints, tr.closed, tr.total_len);
        }
        for (double tau : { 0.1, 0.03 }) {
            auto tr = nurbs::cert::trace_branches(cyl_z(), plane_x(1.0), tau);
            std::printf("  D tau=%.3f | seeds=%3ld branches=%ld endpoints=%ld closed=%ld len=%.3f\n",
                        tau, tr.seeds, tr.branches, tr.endpoints, tr.closed, tr.total_len);
        }
    }

    // ---- M7 2D 重合证书 ----
    std::printf("\n=== M7 2D coincidence certificate (grid point inclusion + normal) ===\n");
    {
        // E1：共面重叠（B 平移 0.5）；A 坐标 [0,0.5]^2 落在 B 的覆盖区内
        auto A = plane_z(0.0);
        auto B = plane_shift(0.5, 0.5);
        double r_in = nurbs::cert::coincidence_ratio(A, B, 0.78, 1.0, 0.78, 1.0, 5);
        double r_out = nurbs::cert::coincidence_ratio(A, B, 0.0, 0.35, 0.0, 0.35, 4);
        std::printf("  E1 coplanar overlap | inside ratio=%.2f  outside ratio=%.2f\n", r_in, r_out);
        // E2：非共面重参数化同像（x = 2s^2-1）
        auto A2 = plane_z(0.0);
        auto B2 = plane_repar();
        double r2 = nurbs::cert::coincidence_ratio(A2, B2, 0.05, 0.95, 0.05, 0.95, 4);
        std::printf("  E2 reparameterized identical image | ratio=%.2f\n", r2);
    }

    std::printf("\ndone.\n");
    return 0;
}
