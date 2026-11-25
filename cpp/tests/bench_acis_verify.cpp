// ============================================================================
// bench_acis_verify.cpp — ACIS 曲面类型对照：真值点是否被覆盖（无遗漏）
//   解析 SURF 块；对给定 (fidA,fidB,真值3D点)：
//     在 A/B 的任意片对上求参数（粗网格+点包含），再用该片对的 cover(τ) 检查覆盖
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::Interval;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static RBezierSurface make_patch(int nu, int nv, const std::vector<std::array<double, 4>>& cw) {
    RBezierSurface s; s.nu = nu; s.nv = nv; s.P.resize((std::size_t)nu * nv);
    for (std::size_t k = 0; k < s.P.size(); ++k) {
        double w = cw[k][3];
        s.P[k][0] = Interval::point(cw[k][0] * w);
        s.P[k][1] = Interval::point(cw[k][1] * w);
        s.P[k][2] = Interval::point(cw[k][2] * w);
        s.P[k][3] = Interval::point(w);
    }
    s.update_weight_range();
    return s;
}

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
        if (ou) *ou = os.mid(); if (ov) *ov = ot.mid();
        return true;
    }
    // 数值回退：Gauss-Newton
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
        double du = (cc * g0 - b * g1) / det, dv = (-b * g0 + a * g1) / det;
        u -= du; v -= dv;
        if (u < 0.0) u = 0.0; if (u > 1.0) u = 1.0;
        if (v < 0.0) v = 0.0; if (v > 1.0) v = 1.0;
    }
    Interval Q[3], Qu[3], Qv[3];
    if (nurbs::cert::eval_rational_duv(S, Interval::point(u), Interval::point(v), Q, Qu, Qv)) {
        double d = 0.0;
        for (int c = 0; c < 3; ++c) { double e = Q[c].mid() - P[c]; d += e * e; }
        if (std::sqrt(d) < 1e-9) { if (ou) *ou = u; if (ov) *ov = v; return true; }
    }
    return false;
}

struct Case { int a, b; const char* name; std::vector<std::array<double, 3>> pts; };

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: %s <acis_types.txt> [tau]\n", argv[0]); return 1; }
    const double tau = (argc >= 3) ? std::atof(argv[2]) : 0.02;
    std::ifstream f(argv[1]);
    std::vector<std::vector<RBezierSurface>> surfs;
    std::string tok;
    while (f >> tok) {
        if (tok == "NSURF") { int n; f >> n; surfs.resize(n); }
        else if (tok == "SURF") {
            int fid, pi, du, dv, nu, nv;
            f >> fid >> pi >> du >> dv >> nu >> nv;
            std::vector<std::array<double, 4>> pts((std::size_t)nu * nv);
            for (int i = 0; i < nu; ++i) for (int j = 0; j < nv; ++j) {
                double x, y, z, w; f >> x >> y >> z >> w;
                pts[(std::size_t)i * nv + j] = { x, y, z, w };
            }
            std::string e; f >> e;
            surfs[fid].push_back(make_patch(nu, nv, pts));
        } else { std::string line; std::getline(f, line); }
    }
    std::printf("parsed %d faces, tau=%.3f\n", (int)surfs.size(), tau);

    auto circle = [](double z, double r, int n) {
        std::vector<std::array<double, 3>> v;
        for (int k = 0; k < n; ++k) {
            double th = 2.0 * M_PI * k / n;
            v.push_back({ r * std::cos(th), r * std::sin(th), z });
        }
        return v;
    };
    std::vector<Case> cases = {
        { 0, 2, "plane z=0   x cylinder r=1", circle(0.0, 1.0, 8) },
        { 1, 2, "plane z=0.5 x cylinder r=1", circle(0.5, 1.0, 8) },
        { 0, 3, "plane z=0   x sphere r=1",   circle(0.0, 1.0, 8) },
        { 1, 3, "plane z=0.5 x sphere r=1",   circle(0.5, std::sqrt(1.0 - 0.25), 8) },
        { 1, 4, "plane z=0.5 x torus",        circle(0.5, 2.5, 8) },
        { 0, 4, "plane z=0   x torus (tangency circle r=2)", circle(0.0, 2.0, 8) },
        { 1, 5, "plane z=0.5 x cone",         circle(0.5, 0.5 * std::tan(30.0 * M_PI / 180.0), 8) },
    };

    int total_pts = 0, covered = 0, failed_find = 0;
    for (const Case& c : cases) {
        std::printf("\n[%s] face %d-%d\n", c.name, c.a, c.b);
        int cov = 0, nf = 0;
        for (const auto& P : c.pts) {
            ++total_pts;
            bool found = false;
            for (const auto& A : surfs[c.a]) {
                for (const auto& B : surfs[c.b]) {
                    double u, v, s2, t2;
                    if (!find_params(A, P.data(), &u, &v)) continue;
                    if (!find_params(B, P.data(), &s2, &t2)) continue;
                    auto r = nurbs::cert::cover(A, B, tau);
                    if (nurbs::cert::covered(r.boxes, u, v, s2, t2, 2 * tau)) {
                        ++cov; ++covered; found = true;
                    }
                    break;
                }
                if (found) break;
            }
            if (!found) { ++nf; ++failed_find; }
        }
        std::printf("  truth covered %d/%d\n", cov, (int)c.pts.size());
    }
    std::printf("\nRESULT total_truth_pts=%d covered=%d not_covered=%d\n", total_pts, covered, failed_find);
    return (covered == total_pts) ? 0 : 1;
}
