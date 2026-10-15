// ============================================================================
// bench_prove_gs.cpp — ① GS-Newton 区间证明 vs Krawczyk 对比
//   目标：证明 GS 比 Krawczyk 更紧（同一盒宽下更高覆盖率）或更快，
//   从而允许全局枚举用更大的叶盒（更少盒、更快）。
//   δ 扫描：对每个横截解构造 [u±δ]×[v±δ] 盒，比较 Krawczyk / GS 证明覆盖率与耗时。
// 随机曲线（mulberry32，确定性）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <chrono>

using namespace nurbs;

namespace {
struct RNG {
    uint32_t s;
    explicit RNG(uint32_t seed) : s(seed) {}
    double nextF() {
        s += 0x6D2B79F5u;
        uint64_t t = (uint64_t)s;
        t = (t ^ (t >> 15)) * (t | 1u);
        t ^= (t + ((t ^ (t >> 7)) * (t | 61u)) ^ ((t ^ (t >> 14)) << 32));
        uint32_t r = (uint32_t)(t ^ (t >> 14));
        return (double)r / 4294967296.0;
    }
};
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
    RNG rng(777);
    // 生成确定性解集（随机曲线对求交）
    struct Case { NURBSCurve c1, c2; std::vector<SolveResult> sol; };
    std::vector<Case> cases;
    int totalInter = 0;
    for (int i = 0; i < 200; ++i) {
        std::vector<Point2> p1(4 + (i % 3)), p2(4 + ((i + 1) % 3));
        for (size_t k = 0; k < p1.size(); ++k) p1[k] = { 0.1 * (i % 5) + 0.7 * k / (p1.size() - 1), 0.6 * rng.nextF() - 0.3 };
        for (size_t k = 0; k < p2.size(); ++k) p2[k] = { 0.1 * (i % 5) + 0.05 + 0.7 * k / (p2.size() - 1), 0.6 * rng.nextF() - 0.3 };
        NURBSCurve c1 = fromPoints(p1, 3), c2 = fromPoints(p2, 3);
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        Case cs{ c1, c2, {} };
        for (auto& r : sol) if (r.state == GeometryState::INTERSECTION) cs.sol.push_back(r);
        totalInter += (int)cs.sol.size();
        if (!cs.sol.empty()) cases.push_back(cs);
    }
    printf("==== ① GS-Newton vs Krawczyk 区间证明对比 ====\n");
    printf("横截解总数 %d（%d 对非空）\n", totalInter, (int)cases.size());

    const double deltas[] = { 1e-2, 5e-3, 2e-3, 1e-3, 1e-4, 1e-6 };
    printf("\nδ   | Krawczyk 覆盖率 | GS 覆盖率 | Krawczyk 耗时(µs) | GS 耗时(µs)\n");
    printf("----+----------------+-----------+-------------------+------------\n");
    for (double d : deltas) {
        int k = 0, g = 0, gExcl = 0;
        long long tK = 0, tG = 0;
        for (auto& cs : cases) for (auto& r : cs.sol) {
            Interval ub(r.u - d, r.u + d), vb(r.v - d, r.v + d);
            auto t0 = std::chrono::high_resolution_clock::now();
            bool pk = proveUniqueInBox(cs.c1, cs.c2, ub, vb);
            auto t1 = std::chrono::high_resolution_clock::now();
            Interval gub(r.u - d, r.u + d), gvb(r.v - d, r.v + d);
            int pg = proveUniqueGS(cs.c1, cs.c2, gub, gvb);
            auto t2 = std::chrono::high_resolution_clock::now();
            if (pk) k++;
            if (pg == 1) g++; else if (pg == -1) gExcl++;
            tK += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
            tG += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        }
        printf("%.0e | %14.1f%% | %12.1f%% | %17.2f | %11.2f\n", d,
               100.0 * k / totalInter, 100.0 * g / totalInter,
               tK / (double)totalInter / 1000.0, tG / (double)totalInter / 1000.0);
        if (d == 1e-4)
            printf("    (GS 排除误判 %d 个，应=0；%s)\n", gExcl,
                   gExcl == 0 ? "排除测试未误伤横截解 ✓" : "警告：GS 误排除横截解");
    }
    return 0;
}
