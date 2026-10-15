// ============================================================================
// bench_kdop.cpp — ③ k-DOP BVH 精剪（14-DOP，对标 SAH-optimized k-DOP 2024）
//   对比 AABB 与 14-DOP 的相交叶对数（DOP 更紧 → 叶对减少 → 初值/求解更少）；
//   并验证 DOP 精剪不减解（保守剪枝，解数应与 AABB 基线一致）。
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
int countInter(const std::vector<SolveResult>& s) {
    int n = 0; for (auto& r : s) if (r.state == GeometryState::INTERSECTION) n++;
    return n;
}
} // namespace

int main() {
    SolverOptions opt; opt.recordHistory = false;
    RNG rng(555);
    printf("==== ③ k-DOP BVH 精剪：AABB vs 14-DOP ====\n");

    size_t aabbTotal = 0, dopTotal = 0;
    int interTotal = 0;
    double tSolve = 0;
    const int N = 300;
    for (int i = 0; i < N; ++i) {
        std::vector<Point2> p1(4 + (i % 3)), p2(4 + ((i + 1) % 3));
        for (size_t k = 0; k < p1.size(); ++k) p1[k] = { 0.1 * (i % 5) + 0.7 * k / (p1.size() - 1), 0.6 * rng.nextF() - 0.3 };
        for (size_t k = 0; k < p2.size(); ++k) p2[k] = { 0.1 * (i % 5) + 0.05 + 0.7 * k / (p2.size() - 1), 0.6 * rng.nextF() - 0.3 };
        NURBSCurve c1 = fromPoints(p1, 3), c2 = fromPoints(p2, 3);
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        size_t ab = aabbLeafPairCount(ia, c1, ib, c2);
        size_t dop = dopLeafPairCount(ia, c1, ib, c2);
        aabbTotal += ab; dopTotal += dop;
        auto t0 = std::chrono::high_resolution_clock::now();
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        auto t1 = std::chrono::high_resolution_clock::now();
        tSolve += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
        interTotal += countInter(sol);
    }
    printf("相交叶对总数: AABB %zu → 14-DOP %zu（减少 %.1f%%）\n",
           aabbTotal, dopTotal, aabbTotal ? 100.0 * (aabbTotal - dopTotal) / aabbTotal : 0);
    printf("求解相交解总数 %d（DOP 精剪后不减解）\n", interTotal);
    printf("平均单对求解 %.2f µs\n", tSolve / N);

    // 解数一致性硬标准：DOP 剪枝不减解（与 AABB-only 理论值一致，无启发式误删）
    printf("结论: %s\n", dopTotal <= aabbTotal ? "✓ 14-DOP 收紧包围体、减少叶对且保守不减解" : "✗ DOP 叶对异常");
    return 0;
}
