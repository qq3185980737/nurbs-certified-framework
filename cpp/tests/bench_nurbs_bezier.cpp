// ============================================================================
// bench_nurbs_bezier.cpp — P7 验证：NURBS→Bézier 精确转换
//   1) 由 Bézier 圆柱片构造 NURBS（插入内部节点 0.5，重数 2）
//   2) to_bezier_patches 应还原出 2 个 Bézier 片
//   3) 在随机参数处比较 NURBS de Boor 求值与 Bézier 分片求值（3D 点）
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <random>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::NurbsPatch;
using nurbs::cert::HP;
using nurbs::Interval;

static const double SQ = 0.7071067811865476;
static const double WQ[3] = { 1.0, SQ, 1.0 };
static const double CC[3] = { 1.0, 1.0, 0.0 };
static const double SS[3] = { 0.0, 1.0, 1.0 };

static RBezierSurface cyl_z() {
    std::vector<std::vector<std::array<double, 4>>> P(3, std::vector<std::array<double, 4>>(2));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { CC[i] * WQ[i], SS[i] * WQ[i], (double)j * WQ[i], WQ[i] };
    int nu = 3, nv = 2;
    std::vector<std::array<double, 4>> flat;
    for (int i = 0; i < nu; ++i) for (int j = 0; j < nv; ++j) flat.push_back(P[i][j]);
    return nurbs::cert::make_surface(nu, nv, flat);
}

static void eval_patch_pt(const RBezierSurface& S, double u, double v, double out[3]) {
    Interval S3[3], Su[3], Sv[3];
    nurbs::cert::eval_rational_duv(S, Interval::point(u), Interval::point(v), S3, Su, Sv);
    for (int c = 0; c < 3; ++c) out[c] = S3[c].mid();
}

int main() {
    auto bz = cyl_z();
    // 构造 NURBS：从 Bézier 出发，插入内部节点 0.5（重数 2）
    NurbsPatch N;
    N.du = 2; N.dv = 1; N.nu = 3; N.nv = 2;
    N.P = bz.P;
    N.U = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
    N.V = { 0.0, 0.0, 1.0, 1.0 };
    // 插两次 u=0.5
    for (int r = 0; r < 2; ++r) {
        std::vector<HP> Pn; std::vector<double> Un;
        nurbs::cert::insert_knot_u(N.P, N.nu, N.nv, N.U, N.du, 0.5, Pn, Un);
        N.P = Pn; N.U = Un; ++N.nu;
    }
    std::printf("NURBS: du=%d dv=%d nu=%d nv=%d, U size=%zu\n", N.du, N.dv, N.nu, N.nv, N.U.size());
    auto patches = nurbs::cert::to_bezier_patches(N);
    std::printf("Bézier patches: %zu (expect 2)\n", patches.size());

    // 对比 NURBS 求值与分片求值
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    double max_err = 0.0;
    for (int trial = 0; trial < 200; ++trial) {
        double u = uni(rng), v = uni(rng);
        HP h;
        nurbs::cert::nurbs_eval_h(N, u, v, h);
        double nu_pt[3] = { h[0].mid() / h[3].mid(), h[1].mid() / h[3].mid(), h[2].mid() / h[3].mid() };
        // 找所在参数片
        int su = (u < 0.5) ? 0 : 1;
        double u0 = (su == 0) ? 0.0 : 0.5, u1 = (su == 0) ? 0.5 : 1.0;
        int patch_idx = su * 1 + 0;   // segV = 1
        double local_u = (u - u0) / (u1 - u0);
        double p[3];
        eval_patch_pt(patches[patch_idx], local_u, v, p);
        double d = 0.0;
        for (int c = 0; c < 3; ++c) d += (p[c] - nu_pt[c]) * (p[c] - nu_pt[c]);
        max_err = std::max(max_err, std::sqrt(d));
    }
    std::printf("max |NURBS - Bezier patch| = %.3e\n", max_err);
    bool ok = (patches.size() == 2) && (max_err < 1e-12);
    std::printf("%s\n", ok ? "PASS: exact NURBS->Bezier conversion" : "FAIL");
    return ok ? 0 : 1;
}
