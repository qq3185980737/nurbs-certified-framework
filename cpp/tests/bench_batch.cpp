// ============================================================================
// bench_batch.cpp — 批量多曲线求交（任务2）规模化基准
// 场景：一个集合内 N 条曲线两两求交（BVH 森林 + sweep-and-prune 候选对筛选）。
// 统计：全 N² 对 vs sweep-and-prune 候选对数（筛选削减）、候选对级别并行求解时间、
//       相交曲线对数 / 交点总数。确定性 LCG（seed=20261012）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>

using namespace nurbs;

namespace {
uint32_t seed = 20261012u;
double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

NURBSCurve makeCurve(int n, int degree) {
    NURBSCurve c;
    c.degree = degree;
    c.controlPoints.resize(n);
    double cx = (rnd() - 0.5) * 2, cy = (rnd() - 0.5) * 2;
    for (int i = 0; i < n; ++i) {
        c.controlPoints[i].x = cx + (rnd() - 0.5) * 4;
        c.controlPoints[i].y = cy + (rnd() - 0.5) * 4;
        c.controlPoints[i].w = 0.5 + rnd();
    }
    c.knots = generateOpenUniformKnotVector(n, degree);
    c.uMin = c.knots[degree];
    c.uMax = c.knots[c.knots.size() - 1 - degree];
    validateCurve(c);
    return c;
}
} // namespace

int main() {
    SolverOptions opt; opt.recordHistory = false;
    const int Ms[] = { 100, 300, 1000 };
    printf("==== 批量多曲线求交（BVH 森林 + sweep-and-prune + 并行）====\n");
    printf("%5s | %12s | %12s | %8s | %12s %12s | %8s %10s\n",
           "N", "全N²对", "候选对", "削减%", "串行ms", "并行ms", "交曲线对", "交点总数");
    for (int M : Ms) {
        seed = 20261012u;
        std::vector<NURBSCurve> C(M);
        for (int i = 0; i < M; ++i) {
            int d = 1 + (int)(rnd() * 3);
            int n = d + 1 + (int)(rnd() * 4);
            C[i] = makeCurve(n, d);
        }
        // BVH 森林 + 根盒
        std::vector<BVHIndex> idx(M);
        std::vector<AABB> roots(M);
        auto b0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < M; ++i) { idx[i].build(C[i], 7); roots[i] = idx[i].rootBox(); }
        auto b1 = std::chrono::high_resolution_clock::now();
        const size_t full = (size_t)M * (M - 1) / 2;
        const auto pairs = sweepAndPrunePairs(roots);
        const size_t cand = pairs.size();
        const double cut = 100.0 * (1.0 - (double)cand / full);
        // 串行
        auto s0 = std::chrono::high_resolution_clock::now();
        size_t serInter = 0;
        for (auto& pr : pairs) {
            auto r = solveWithIndex(idx[pr.first], C[pr.first], idx[pr.second], C[pr.second], opt, 1);
            for (auto& s : r) if (s.state == GeometryState::INTERSECTION) serInter++;
        }
        auto s1 = std::chrono::high_resolution_clock::now();
        // 并行
        auto p0 = std::chrono::high_resolution_clock::now();
        const auto batch = batchSolveAll(C, opt, 0);
        auto p1 = std::chrono::high_resolution_clock::now();
        size_t interPairs = 0, interTotal = 0;
        for (auto& b : batch) {
            bool has = false;
            for (auto& s : b.r) if (s.state == GeometryState::INTERSECTION) { has = true; interTotal++; }
            if (has) interPairs++;
        }
        double tSer = std::chrono::duration<double>(s1 - s0).count();
        double tPar = std::chrono::duration<double>(p1 - p0).count();
        double tBuild = std::chrono::duration<double>(b1 - b0).count();
        printf("%5d | %12llu | %12llu | %7.1f%% | %12.2f %12.2f | %8llu %10llu\n",
               M, (unsigned long long)full, (unsigned long long)cand, cut,
               tSer * 1e3, tPar * 1e3, (unsigned long long)interPairs,
               (unsigned long long)interTotal);
        printf("        （索引构建 %.2f ms；串行交点=%llu，并行交点=%llu，一致:%s）\n",
               tBuild * 1e3, (unsigned long long)serInter, (unsigned long long)interTotal,
               serInter == interTotal ? "是" : "否");
    }
    return 0;
}
