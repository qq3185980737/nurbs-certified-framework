// ============================================================================
// bench_tangency_chain.cpp — P3 延拓链验证：切触曲线的认证接触盒链
//   用法: bench_tangency_chain <data.txt> <i> <j> <tau>
//   步骤：
//     1) cover(τ)；
//     2) 每个覆盖盒做相切认证（certify_tangency_box==0 → 接触盒）；
//     3) 接触盒的 3D 像 AABB 邻接图 → 连通分量；
//     4) 对切触曲线上的 16 个角度真值点，检查是否被某个接触盒的 3D 像 AABB 覆盖；
//   输出：接触盒数 / 分量数 / 角度覆盖数
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <functional>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;
using nurbs::cert::AABB3;

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

static bool aabb_contains(const AABB3& a, double x, double y, double z, double tol) {
    return x >= a.x.lo - tol && x <= a.x.hi + tol &&
           y >= a.y.lo - tol && y <= a.y.hi + tol &&
           z >= a.z.lo - tol && z <= a.z.hi + tol;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::printf("usage: %s data i j tau\n", argv[0]); return 1; }
    std::ifstream f(argv[1]);
    int ti = std::atoi(argv[2]), tj = std::atoi(argv[3]);
    double tau = std::atof(argv[4]);
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

    long total_contact = 0, total_boxes = 0;
    bool any_covered[16] = { false };
    for (const auto& A : surfs[ti]) {
        for (const auto& B : surfs[tj]) {
            auto r = nurbs::cert::cover(A, B, tau);
            if (r.cover == 0) continue;
            total_boxes += r.cover;
            std::vector<AABB3> contact;
            for (const Box4& X : r.boxes) {
                if (nurbs::cert::certify_tangency_box(A, B, X) != 0) continue;
                contact.push_back(nurbs::cert::box_image_aabb(A, B, X));
            }
            total_contact += (long)contact.size();
            // 连通分量（3D 像 AABB 重叠）
            std::vector<int> parent(contact.size());
            for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = (int)i;
            std::function<int(int)> find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
            for (std::size_t i = 0; i < contact.size(); ++i)
                for (std::size_t j = i + 1; j < contact.size(); ++j)
                    if (nurbs::cert::aabb_overlap3(contact[i], contact[j], 1e-9))
                        parent[find((int)i)] = find((int)j);
            int comps = 0;
            for (std::size_t i = 0; i < contact.size(); ++i) if (find((int)i) == (int)i) ++comps;
            // 16 角度覆盖
            for (int k = 0; k < 16; ++k) {
                double th = 2.0 * M_PI * k / 16.0;
                double px = 2.0 * std::cos(th), py = 2.0 * std::sin(th), pz = 0.0;
                for (const AABB3& a : contact)
                    if (aabb_contains(a, px, py, pz, 2.0 * tau)) { any_covered[k] = true; break; }
            }
            std::printf("  patch pair: cover=%ld contact=%zu comps=%d\n", r.cover, contact.size(), comps);
        }
    }
    int cov = 0;
    for (int k = 0; k < 16; ++k) if (any_covered[k]) ++cov;
    std::printf("RESULT total_boxes=%ld contact_boxes=%ld angular_coverage=%d/16\n",
                total_boxes, total_contact, cov);
    return 0;
}
