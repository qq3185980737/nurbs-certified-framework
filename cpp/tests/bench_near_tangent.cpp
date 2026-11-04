// ============================================================================
// bench_near_tangent.cpp — P8：近切小环的 τ 阈值验证
//   案例：抛物面 z=(2u-1)^2+(2v-1)^2 与平面 z=h 相交于半径 sqrt(h) 的小圆
//   对 h=1e-6：true radius = 1e-3；测试 τ = 1e-3, 5e-4, 2.5e-4, 1.25e-4 下
//   覆盖/秩判定/认证/真值覆盖 是否随 τ 减小而恢复
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;

using P4 = std::array<double, 4>;

static RBezierSurface make(const std::vector<std::vector<P4>>& poles) {
    int nu = (int)poles.size(), nv = (int)poles[0].size();
    std::vector<P4> flat;
    for (int i = 0; i < nu; ++i) for (int j = 0; j < nv; ++j) flat.push_back(poles[i][j]);
    return nurbs::cert::make_surface(nu, nv, flat);
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

static RBezierSurface plane_z(double z0) {
    std::vector<std::vector<P4>> P(2, std::vector<P4>(2));
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            P[i][j] = { 2.0 * i - 1.0, 2.0 * j - 1.0, z0, 1.0 };
    return make(P);
}

int main() {
    const double h = 1e-6;
    const double radius = std::sqrt(h);          // 1e-3
    auto A = paraboloid();
    auto B = plane_z(h);
    // 真值点：圆上 12 点 (x,y,h)，参数 u=(x+1)/2, v=(y+1)/2（两片相同）
    std::vector<std::array<double, 4>> truth;
    for (int k = 0; k < 12; ++k) {
        double th = 2.0 * M_PI * k / 12.0;
        double x = radius * std::cos(th), y = radius * std::sin(th);
        truth.push_back({ (x + 1.0) / 2.0, (y + 1.0) / 2.0, (x + 1.0) / 2.0, (y + 1.0) / 2.0 });
    }
    std::printf("h=%.1e, true radius=%.3e; feature/4 = %.3e\n", h, radius, radius / 4.0);
    std::printf("tau        cover    rank3   cert   tang   truth_cov   time_ms\n");
    for (double tau : { 1e-3, 5e-4, 2.5e-4, 1.25e-4 }) {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto r = nurbs::cert::cover(A, B, tau);
        auto st = nurbs::cert::certify_boxes(A, B, r.boxes);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        int cov = 0;
        for (const auto& p : truth)
            if (nurbs::cert::covered(r.boxes, p[0], p[1], p[2], p[3], 2 * tau)) ++cov;
        std::printf("%-10.2e %-8ld %-7ld %-6ld %-6ld %2d/12       %.1f\n",
                    tau, r.cover, st.rank3, st.certified, st.tangent, cov, ms);
    }
    return 0;
}
