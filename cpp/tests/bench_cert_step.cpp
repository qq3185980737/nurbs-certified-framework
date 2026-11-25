// ============================================================================
// bench_cert_step.cpp — M8：真实 STEP 面片（bracket.step 导出的有理 Bezier 片）
//   逐对运行：4D 覆盖 + 横截/相切认证 + 拓扑三分；统计相交/不交、分支、时间
// 用法: bench_cert_step <data.txt> [tau]
// 数据格式: NSURF n; SURF fid pi du dv nu nv; <nu*nv 行 "x y z w">; ENDSURF
// ============================================================================
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <set>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::Interval;
using nurbs::cert::Box4;

static RBezierSurface make_patch(int nu, int nv, const std::vector<std::array<double, 4>>& cart_w) {
    RBezierSurface s;
    s.nu = nu; s.nv = nv;
    s.P.resize((std::size_t)nu * nv);
    for (std::size_t k = 0; k < s.P.size(); ++k) {
        double w = cart_w[k][3];
        s.P[k][0] = nurbs::Interval::point(cart_w[k][0] * w);
        s.P[k][1] = nurbs::Interval::point(cart_w[k][1] * w);
        s.P[k][2] = nurbs::Interval::point(cart_w[k][2] * w);
        s.P[k][3] = nurbs::Interval::point(w);
    }
    return s;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: %s <data.txt> [tau]\n", argv[0]);
        return 1;
    }
    const char* path = argv[1];
    double tau = (argc >= 3) ? std::atof(argv[2]) : 0.02;
    std::ifstream f(path);
    if (!f) { std::printf("cannot open %s\n", path); return 1; }

    std::vector<std::vector<RBezierSurface>> surfs;
    std::string tok;
    std::set<std::pair<int,int>> adjacent;
    while (f >> tok) {
        if (tok == "NSURF") { int n; f >> n; surfs.resize(n); }
        else if (tok == "NADJ") {
            int n; f >> n;
            for (int k = 0; k < n; ++k) {
                int i, j; f >> i >> j;
                if (i > j) { int tmp = i; i = j; j = tmp; }
                adjacent.insert({i, j});
            }
        }
        else if (tok == "NURB") {
            // NURB fid conv du dv nu nv / U knots / V knots / poles(w 含权)
            int fid, conv, du, dv, nu, nv;
            f >> fid >> conv >> du >> dv >> nu >> nv;
            (void)conv;
            if (fid >= (int)surfs.size()) surfs.resize(fid + 1);
            nurbs::cert::NurbsPatch N;
            N.du = du; N.dv = dv; N.nu = nu; N.nv = nv;
            N.U.resize((std::size_t)(nu + du + 1));
            for (int i = 0; i < nu + du + 1; ++i) f >> N.U[i];
            N.V.resize((std::size_t)(nv + dv + 1));
            for (int i = 0; i < nv + dv + 1; ++i) f >> N.V[i];
            N.P.resize((std::size_t)nu * nv);
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j) {
                    double x, y, z, w; f >> x >> y >> z >> w;
                    nurbs::cert::HP hp;
                    hp[0] = nurbs::Interval::point(x * w);
                    hp[1] = nurbs::Interval::point(y * w);
                    hp[2] = nurbs::Interval::point(z * w);
                    hp[3] = nurbs::Interval::point(w);
                    N.P[(std::size_t)i * nv + j] = hp;
                }
            std::string e; f >> e;
            for (auto& patch : nurbs::cert::to_bezier_patches(N)) surfs[fid].push_back(patch);
        }
        else if (tok == "SURF") {
            int fid, pi, du, dv, nu, nv;
            f >> fid >> pi >> du >> dv >> nu >> nv;
            std::vector<std::array<double, 4>> pts((std::size_t)nu * nv);
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j) {
                    double x, y, z, w;
                    f >> x >> y >> z >> w;
                    pts[(std::size_t)i * nv + j] = { x, y, z, w };
                }
            std::string e; f >> e;
            surfs[fid].push_back(make_patch(nu, nv, pts));
        } else if (tok == "NCURVE" || tok == "NPAIR") {
            int n; f >> n;
        } else {
            std::string line;
            std::getline(f, line);
        }
    }

    int only_i = -1, only_j = -1;
    if (argc >= 5) { only_i = std::atoi(argv[3]); only_j = std::atoi(argv[4]); }
    int ns = (int)surfs.size();
    std::printf("parsed %d faces, tau=%.4f\n", ns, tau);
    long pairs = 0, intersecting = 0, disjoint = 0, branches = 0, cert_ok = 0, undec = 0;
    long adjacent_skipped = 0;
    long t_1d = 0, t_2d = 0, t_copl = 0, t_und = 0, t_cover = 0, t_tang = 0, t_bound = 0;
    double ms_total = 0.0;
    for (int a = 0; a < ns; ++a) {
        for (int b = a + 1; b < ns; ++b) {
            if (only_i >= 0 && (a != only_i || b != only_j)) continue;
            if (only_i < 0 && adjacent.count({a, b})) { ++adjacent_skipped; continue; }
            for (const auto& A : surfs[a]) {
                for (const auto& B : surfs[b]) {
                    ++pairs;
                    auto t0 = std::chrono::high_resolution_clock::now();
                    auto pc = nurbs::cert::classify_pair(A, B, tau);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    ms_total += std::chrono::duration<double, std::milli>(t1 - t0).count();
                    if (pc.type == 0) ++disjoint;
                    else {
                        ++intersecting;
                        if (pc.type == 1) { ++t_1d; branches += pc.components; }
                        else if (pc.type == 2) ++t_2d;
                        else if (pc.type == 4) ++t_copl;
                        else if (pc.type == 5) ++t_cover;
                        else if (pc.type == 6) ++t_tang;
                        else if (pc.type == 7) ++t_bound;
                        else ++t_und;
                        cert_ok += pc.certified;
                        undec += pc.undecided;
                        std::printf("  pair face %2d-%2d: type=%d cover=%5ld cert=%4ld comps=%ld coinc=%.2f\n",
                                    a, b, pc.type, pc.cover, pc.certified, pc.components, pc.coincidence);
                    }
                }
            }
        }
    }
    std::printf("\nRESULT face-pairs=%ld | disjoint=%ld | 1D_certified=%ld (comps=%ld) | 2D_coincidence=%ld | coplanar_contact=%ld | tangent_candidate=%ld | cover_limit=%ld | patch_boundary=%ld | undecided=%ld | adjacent_skipped=%ld\n",
                pairs, disjoint, t_1d, branches, t_2d, t_copl, t_tang, t_cover, t_bound, t_und, adjacent_skipped);
    std::printf("       certified_boxes=%ld undecided_boxes=%ld time=%.1f ms (%.2f ms/pair)\n",
                cert_ok, undec, ms_total, ms_total / std::max(1L, pairs));
    return 0;
}
