// ============================================================================
// bench_prove.cpp — A: 保证性求交（区间/Krawczyk 解存在性证明）覆盖率基准
// 对每个 BVH 加速求解得到的 INTERSECTION 解 (u,v)，构造参数邻域盒，
// 用严格舍入区间扩展 + Krawczyk 算子证明该盒内存在唯一解。
// 指标：横截解中获得存在性证明的比例（预期近 100%）；
//       相切/退化解 Krawczyk 失败是正确行为（Jacobian 奇异，衔接解分类）。
// 确定性 LCG seed=20261012。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include "nurbs/interval_arith.hpp"
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
    const int N = 300;
    std::vector<NURBSCurve> C(2 * N);
    for (int i = 0; i < 2 * N; ++i) {
        int d = 1 + (int)(rnd() * 3);
        int n = d + 1 + (int)(rnd() * 4);
        C[i] = makeCurve(n, d);
    }
    printf("==== A: 保证性求交 — Krawczyk 解存在性证明覆盖率（300 对）====\n");
    // 预求解一次，收集所有横截解
    struct Sol { NURBSCurve c1, c2; double u, v; };
    std::vector<Sol> sols;
    for (int i = 0; i < N; ++i) {
        BVHIndex ia, ib; ia.build(C[2 * i], 7); ib.build(C[2 * i + 1], 7);
        auto sol = solveWithIndex(ia, C[2 * i], ib, C[2 * i + 1], opt, 1);
        for (auto& r : sol) if (r.state == GeometryState::INTERSECTION)
            sols.push_back({ C[2 * i], C[2 * i + 1], r.u, r.v });
    }
    printf("横截交点总数: %d\n\n", (int)sols.size());
    // 参数化盒半宽 δ 扫描
    const double d2[] = { 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };
    printf("%12s | %8s | %8s\n", "盒半宽 δ", "证明数", "覆盖率");
    for (double delta : d2) {
        int proven = 0;
        for (auto& s : sols) {
            if (proveUniqueInBox(s.c1, s.c2, Interval(s.u - delta, s.u + delta),
                                 Interval(s.v - delta, s.v + delta))) proven++;
        }
        printf("%12g | %8d | %7.1f%%\n", delta, proven,
               100.0 * proven / (sols.empty() ? 1 : (int)sols.size()));
    }
    return 0;
}
