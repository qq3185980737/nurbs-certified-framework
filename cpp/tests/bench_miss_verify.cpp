// ============================================================================
// bench_miss_verify.cpp — 独立校验 OCCT 漏检案例的本方法结果
//   对每个案例：取解析真值 3D 交点 P → 用点包含认证求 P 在两张面上的参数
//   → 检查 cover(τ) 是否覆盖该 4D 参数点 → 打印验证结论
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;

static RBezierSurface make_patch(int nu, int nv, const std::vector<std::array<double, 4>>& cart_w) {
    RBezierSurface s;
    s.nu = nu; s.nv = nv;
    s.P.resize((std::size_t)nu * nv);
    for (std::size_t k = 0; k < s.P.size(); ++k) {
        double w = cart_w[k][3];
        s.P[k][0] = Interval::point(cart_w[k][0] * w);
        s.P[k][1] = Interval::point(cart_w[k][1] * w);
        s.P[k][2] = Interval::point(cart_w[k][2] * w);
        s.P[k][3] = Interval::point(w);
    }
    s.update_weight_range();
    return s;
}

static bool inside_box(const Box4& b, double u, double v, double s, double t, double tol) {
    return u >= b.u.lo - tol && u <= b.u.hi + tol &&
           v >= b.v.lo - tol && v <= b.v.hi + tol &&
           s >= b.s.lo - tol && s <= b.s.hi + tol &&
           t >= b.t.lo - tol && t <= b.t.hi + tol;
}

// 无初值点包含：粗网格最近点 + 点包含（边界 0/1 也在网格内）
static bool find_params(const RBezierSurface& S, const double P[3], double* ou, double* ov) {
    double best_s = 0.5, best_t = 0.5, best_d = 1e300;
    const int N = 12;
    for (int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j) {
            double s = (double)i / N, t = (double)j / N;
            Interval Q[3], Qu[3], Qv[3];
            if (!nurbs::cert::eval_rational_duv(S, Interval::point(s), Interval::point(t), Q, Qu, Qv)) continue;
            double d = 0.0;
            for (int c = 0; c < 3; ++c) { double e = Q[c].mid() - P[c]; d += e * e; }
            if (d < best_d) { best_d = d; best_s = s; best_t = t; }
        }
    Interval os, ot;
    if (nurbs::cert::point_on_patch(S, P, best_s, best_t, &os, &ot)) {
        if (ou) *ou = os.mid();
        if (ov) *ov = ot.mid();
        return true;
    }
    // 数值回退：Gauss-Newton（3x2 广义逆）
    double u = best_s, v = best_t;
    for (int it = 0; it < 60; ++it) {
        Interval Q[3], Qu[3], Qv[3];
        if (!nurbs::cert::eval_rational_duv(S, Interval::point(u), Interval::point(v), Q, Qu, Qv)) return false;
        double F[3], Ju[3], Jv[3];
        for (int c = 0; c < 3; ++c) { F[c] = Q[c].mid() - P[c]; Ju[c] = Qu[c].mid(); Jv[c] = Qv[c].mid(); }
        double a = 0, b = 0, cc = 0, g0 = 0, g1 = 0;
        for (int k = 0; k < 3; ++k) { a += Ju[k] * Ju[k]; b += Ju[k] * Jv[k]; cc += Jv[k] * Jv[k]; g0 += Ju[k] * F[k]; g1 += Jv[k] * F[k]; }
        double det = a * cc - b * b;
        if (std::fabs(det) < 1e-30) break;
        double du = (cc * g0 - b * g1) / det;
        double dv = (-b * g0 + a * g1) / det;
        u -= du; v -= dv;
        if (u < -1e-9) u = 0.0; if (u > 1.0 + 1e-9) u = 1.0;
        if (v < -1e-9) v = 0.0; if (v > 1.0 + 1e-9) v = 1.0;
    }
    // 验证残差
    Interval Q[3], Qu[3], Qv[3];
    if (nurbs::cert::eval_rational_duv(S, Interval::point(u), Interval::point(v), Q, Qu, Qv)) {
        double d = 0.0;
        for (int c = 0; c < 3; ++c) { double e = Q[c].mid() - P[c]; d += e * e; }
        if (std::sqrt(d) < 1e-9) { if (ou) *ou = u; if (ov) *ov = v; return true; }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: %s <data.txt> [tau]\n", argv[0]); return 1; }
    const double tau = (argc >= 3) ? std::atof(argv[2]) : 0.02;
    std::ifstream f(argv[1]);
    if (!f) { std::printf("cannot open\n"); return 1; }

    std::vector<std::vector<RBezierSurface>> surfs;
    std::string tok;
    while (f >> tok) {
        if (tok == "NSURF") { int n; f >> n; surfs.resize(n); }
        else if (tok == "SURF") {
            int fid, pi, du, dv, nu, nv;
            f >> fid >> pi >> du >> dv >> nu >> nv;
            std::vector<std::array<double, 4>> pts((std::size_t)nu * nv);
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j) {
                    double x, y, z, w; f >> x >> y >> z >> w;
                    pts[(std::size_t)i * nv + j] = { x, y, z, w };
                }
            std::string e; f >> e;
            surfs[fid].push_back(make_patch(nu, nv, pts));
        } else { std::string line; std::getline(f, line); }
    }
    std::printf("parsed %d faces, tau=%.4f\n", (int)surfs.size(), tau);

    struct Case { int a, b; double P[3]; const char* name; };
    const double eps = 1e-6;
    std::vector<Case> cases = {
        { 0, 1, { 1.0, 0.0, 0.0 }, "sphere r=1 x plane z=0 (equator point)" },
        { 4, 5, { 1e-6, 0.0, 0.0 }, "tiny sphere r=1e-6 x plane" },
        { 6, 9, { 1e-5, 0.0, 0.0 }, "tiny cylinder r=1e-5 x plane" },
        { 12, 13, { 0.35, -0.2, 0.0 }, "identical planes (coincident point)" },
        { 10, 11, { 1.0, 0.0, 0.0 }, "identical spheres (coincident point)" },
    };
    // 环面-平面切触圆：8 个角度上的真值点
    static const char* torus_names[8] = {
        "torus circle pt 0", "torus circle pt 45", "torus circle pt 90", "torus circle pt 135",
        "torus circle pt 180", "torus circle pt 225", "torus circle pt 270", "torus circle pt 315" };
    for (int k = 0; k < 8; ++k) {
        double th = 2.0 * M_PI * k / 8.0;
        cases.push_back({ 2, 3, { 2.0 * std::cos(th), 2.0 * std::sin(th), 0.0 }, torus_names[k] });
    }

    int verified = 0, covered = 0, total = 0;
    for (const Case& c : cases) {
        std::printf("\n[%s] pair %d-%d, truth P=(%.3e,%.3e,%.3e)\n", c.name, c.a, c.b, c.P[0], c.P[1], c.P[2]);
        if (c.a >= (int)surfs.size() || c.b >= (int)surfs.size()) { std::printf("  pair out of range\n"); continue; }
        bool ok = false;
        for (const auto& A : surfs[c.a]) {
            for (const auto& B : surfs[c.b]) {
                double u, v, s, t;
                if (!find_params(A, c.P, &u, &v)) continue;
                if (!find_params(B, c.P, &s, &t)) continue;
                ++total;
                std::printf("  point inclusion: S1 params (%.4f,%.4f), S2 params (%.4f,%.4f)\n", u, v, s, t);
                auto r = nurbs::cert::cover(A, B, tau);
                bool cov = false;
                for (const Box4& bx : r.boxes) if (inside_box(bx, u, v, s, t, 2 * tau)) { cov = true; break; }
                std::printf("  cover boxes=%ld -> truth point %s\n", r.cover, cov ? "COVERED" : "NOT COVERED");
                ++verified;
                if (cov) ++covered;
                ok = true;
                break;
            }
            if (ok) break;
        }
        if (!ok) std::printf("  point inclusion FAILED on one of the patches\n");
    }
    std::printf("\nRESULT verified=%d covered=%d (total param combos=%d)\n", verified, covered, total);
    return 0;
}
