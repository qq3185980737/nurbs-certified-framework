// ============================================================================
// bench_probe_pair.cpp — 诊断：某面片对的覆盖盒认证失败原因直方图
// 用法: bench_probe_pair <data.txt> <i> <j> <tau>
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <map>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;

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

int main(int argc, char** argv) {
    if (argc < 5) { std::printf("usage: %s data i j tau\n", argv[0]); return 1; }
    const char* path = argv[1];
    int ti = std::atoi(argv[2]), tj = std::atoi(argv[3]);
    double tau = std::atof(argv[4]);
    std::ifstream f(path);
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
    if (ti >= (int)surfs.size() || tj >= (int)surfs.size()) { std::printf("bad pair\n"); return 1; }
    std::printf("pair %d-%d: patches %zu x %zu, tau=%.3f\n", ti, tj, surfs[ti].size(), surfs[tj].size(), tau);
    for (const auto& A : surfs[ti]) {
        for (const auto& B : surfs[tj]) {
            auto r = nurbs::cert::cover(A, B, tau);
            std::map<int, long> hist;
            long rank3 = 0, eval_fail = 0;
            for (const Box4& X : r.boxes) {
                Interval S1[3], S1u[3], S1v[3], S2[3], S2u[3], S2v[3];
                bool ok = nurbs::cert::eval_rational_duv(A, X.u, X.v, S1, S1u, S1v) &&
                          nurbs::cert::eval_rational_duv(B, X.s, X.t, S2, S2u, S2v);
                if (!ok) { ++eval_fail; hist[-1]++; continue; }
                Interval J[3][4];
                for (int c = 0; c < 3; ++c) {
                    J[c][0] = S1u[c]; J[c][1] = S1v[c];
                    J[c][2] = Interval(-S2u[c].hi, -S2u[c].lo);
                    J[c][3] = Interval(-S2v[c].hi, -S2v[c].lo);
                }
                if (nurbs::cert::rank3_certified(J)) ++rank3;
                int code = nurbs::cert::krawczyk_slice_code(A, B, X);
                hist[code]++;
            }
            std::printf("  cover=%ld rank3=%ld eval_fail=%ld\n", r.cover, rank3, eval_fail);
            for (auto& kv : hist) std::printf("    code %2d : %ld\n", kv.first, kv.second);
        }
    }
    return 0;
}
