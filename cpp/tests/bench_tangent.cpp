// ============================================================================
// bench_tangent.cpp — ② 相切交点精化（对标 JCAD 2025 有理插值处理相切）
//   相切/近奇异交点：LM 线性收敛、Krawczyk 不适用 → 阻尼 Gauss-Newton 精化
//   把切点定位到机器精度。对比：统一求解器(相切距离) vs refineTangentPoint 精化后。
// 确定性构造多个相切场景。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace nurbs;

namespace {
NURBSCurve fromPoints(const std::vector<Point2>& pts, int degree) {
    NURBSCurve c;
    c.degree = degree;
    c.controlPoints.resize(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        c.controlPoints[i].x = pts[i].x; c.controlPoints[i].y = pts[i].y; c.controlPoints[i].w = 1.0;
    }
    c.knots = generateOpenUniformKnotVector((int)pts.size(), degree);
    c.uMin = c.knots[degree];
    c.uMax = c.knots[c.knots.size() - 1 - degree];
    validateCurve(c);
    return c;
}
} // namespace

int main() {
    SolverOptions opt; opt.recordHistory = false;
    printf("==== ② 相切交点精化（阻尼 Gauss-Newton）====\n");

    // 场景 1：抛物线 y=u(1-u) 顶点与水平线相切
    {
        NURBSCurve c1 = fromPoints({ {0, 0}, {0.5, 0.5}, {1, 0} }, 2);  // y=u(1-u)
        NURBSCurve c2 = fromPoints({ {0, 0.25}, {1, 0.25} }, 1);        // y=0.25
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        double du = 0;
        for (auto& r : sol) if (r.state != GeometryState::SEPARATED) du = r.distance;
        double u = 0.5, v = 0.5;
        double refined = refineTangentPoint(c1, c2, u, v);
        printf("[抛物线顶点相切] 求解器距离 %.2e → 精化后 %.2e（u=%.6f,v=%.6f）改善 %.0f×\n",
               du, refined, u, v, du / (refined > 0 ? refined : 1e-300));
    }

    // 场景 2：圆角相切（构造两曲线在某点切线相同）
    {
        // c1: (u, u^2) 顶点，c2: 反抛物线在下方相切（倒置共享顶点法向）
        NURBSCurve c1 = fromPoints({ {0, 0}, {0.5, 0.5}, {1, 0} }, 2);   // y=u(1-u), 顶点(0.5,0.25)
        NURBSCurve c2 = fromPoints({ {0, 0.5}, {0.5, 0.0}, {1, 0.5} }, 2); // 倒抛物线 y=0.5-2(u-0.5)^2 顶点(0.5,0.5)
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        double du = 1e9;
        for (auto& r : sol) if (r.state != GeometryState::SEPARATED) du = std::min(du, r.distance);
        double u = 0.5, v = 0.5;
        double refined = refineTangentPoint(c1, c2, u, v);
        printf("[两抛物线相切] 求解器距离 %.2e → 精化后 %.2e（u=%.6f,v=%.6f）改善 %.0f×\n",
               du, refined, u, v, du / (refined > 0 ? refined : 1e-300));
    }

    // 场景 3：随机曲线对中检出的近奇异解精化（统计）
    {
        struct RNG { uint32_t s; explicit RNG(uint32_t x) : s(x) {} double nextF() { s += 0x6D2B79F5u; uint64_t t = s; t = (t ^ (t >> 15)) * (t | 1u); t ^= (t + ((t ^ (t >> 7)) * (t | 61u)) ^ ((t ^ (t >> 14)) << 32)); uint32_t r = (uint32_t)(t ^ (t >> 14)); return (double)r / 4294967296.0; } };
        RNG rng(321);
        double sumBefore = 0, sumAfter = 0; int n = 0;
        for (int i = 0; i < 200; ++i) {
            std::vector<Point2> p1(4 + (i % 3)), p2(4 + ((i + 1) % 3));
            for (size_t k = 0; k < p1.size(); ++k) p1[k] = { 0.1 * (i % 5) + 0.7 * k / (p1.size() - 1), 0.6 * rng.nextF() - 0.3 };
            for (size_t k = 0; k < p2.size(); ++k) p2[k] = { 0.1 * (i % 5) + 0.05 + 0.7 * k / (p2.size() - 1), 0.6 * rng.nextF() - 0.3 };
            NURBSCurve c1 = fromPoints(p1, 3), c2 = fromPoints(p2, 3);
            BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
            auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
            for (auto& r : sol) if (r.state == GeometryState::INTERSECTION) {
                double u = r.u, v = r.v;
                double after = refineTangentPoint(c1, c2, u, v);
                sumBefore += r.distance; sumAfter += after; n++;
            }
        }
        printf("[随机 200 对 %d 交点] 平均距离 %.2e → 精化后 %.2e（改善 %.1f×）\n",
               n, n ? sumBefore / n : 0, n ? sumAfter / n : 0,
               n ? sumBefore / (sumAfter > 0 ? sumAfter : 1e-300) : 0);
    }
    return 0;
}
