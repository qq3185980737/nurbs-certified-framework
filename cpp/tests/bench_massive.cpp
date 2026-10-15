// ============================================================================
// bench_massive.cpp — 综合：NURBS 曲线类型全覆盖 × 跨类型相交（每类型 1000 组）
//   9 类曲线：直线(L)/二次(Q)/三次(C)/高次(H)/圆弧(A,有理w)/椭圆(E,有理w)/
//             折角(K,C⁰)/闭曲线(O,整圆)/非均匀+权重(NU)。
//   每种类型作 c1 与混合类型 c2 求交 1000 组（跨类型配对），求解定位 + Krawczyk 证明。
//   统计每类型求解交点总数、证明数、覆盖率 → 类型全覆盖 + 跨类型相交证据。
// 自包含（用主库接口）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include "nurbs/interval_arith.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

using namespace nurbs;

namespace {
uint32_t seed = 20261014u;
double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

enum Type { L=0, Q, C, H, A, E, K, O, NU, NUM_TYPES };
const char* TNAME[NUM_TYPES] = { "直线", "二次", "三次", "高次", "圆弧", "椭圆", "折角", "闭曲线", "非均匀+权重" };

NURBSCurve line0() {
    NURBSCurve c; c.degree = 1; c.controlPoints.resize(2);
    c.controlPoints[0] = { rnd(), rnd() - 0.5, 1.0 };
    c.controlPoints[1] = { rnd(), rnd() - 0.5, 1.0 };
    c.knots = { 0, 0, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
NURBSCurve bezier(int deg, int n) {
    NURBSCurve c; c.degree = deg; c.controlPoints.resize(n);
    for (int i = 0; i < n; ++i) { c.controlPoints[i].x = rnd(); c.controlPoints[i].y = rnd() - 0.5; c.controlPoints[i].w = 1.0; }
    c.knots = generateOpenUniformKnotVector(n, deg);
    c.uMin = c.knots[deg]; c.uMax = c.knots[c.knots.size() - 1 - deg]; validateCurve(c); return c;
}
NURBSCurve arc() {   // 标准 1/4 圆弧，随机缩放+平移
    NURBSCurve c; c.degree = 2; c.controlPoints.resize(3);
    double a = 0.2 + 0.4 * rnd(), cx = 0.3 + 0.4 * rnd(), cy = rnd() - 0.5;
    const double w = std::sqrt(0.5);
    c.controlPoints[0] = { cx + a, cy, 1.0 };
    c.controlPoints[1] = { cx + a, cy + a, w };
    c.controlPoints[2] = { cx, cy + a, 1.0 };
    c.knots = { 0, 0, 0, 1, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
NURBSCurve ell() {
    NURBSCurve c; c.degree = 2; c.controlPoints.resize(3);
    double a = 0.2 + 0.3 * rnd(), b = 0.1 + 0.2 * rnd(), cx = 0.3 + 0.4 * rnd(), cy = rnd() - 0.5;
    const double w = std::sqrt(0.5);
    c.controlPoints[0] = { cx + a, cy, 1.0 };
    c.controlPoints[1] = { cx + a, cy + b, w };
    c.controlPoints[2] = { cx, cy + b, 1.0 };
    c.knots = { 0, 0, 0, 1, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
NURBSCurve kink() {  // C⁰ 折角 V 形（内部节点 0.5 重 2）
    double s = 0.6 + 0.4 * rnd(), cy = rnd() - 0.3;
    NURBSCurve c; c.degree = 2; c.controlPoints.resize(5);
    c.controlPoints[0] = { 0, cy, 1.0 };
    c.controlPoints[1] = { 0.25 * s, cy + 0.75 * s, 1.0 };
    c.controlPoints[2] = { 0.5 * s, cy + 0.5 * s, 1.0 };
    c.controlPoints[3] = { 0.75 * s, cy + 0.75 * s, 1.0 };
    c.controlPoints[4] = { s, cy, 1.0 };
    c.knots = { 0, 0, 0, 0.5, 0.5, 1, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
NURBSCurve closed() {  // 整圆（4 段圆弧拼接，闭合）
    NURBSCurve c; c.degree = 2; c.controlPoints.resize(9);
    double a = 0.2 + 0.2 * rnd(), cx = 0.3 + 0.4 * rnd(), cy = rnd() - 0.5;
    const double w = std::sqrt(0.5);
    const double pts[9][3] = { {1,0,1},{1,1,w},{0,1,1},{-1,1,w},{-1,0,1},{-1,-1,w},{0,-1,1},{1,-1,w},{1,0,1} };
    for (int i = 0; i < 9; ++i) { c.controlPoints[i].x = cx + a * pts[i][0]; c.controlPoints[i].y = cy + a * pts[i][1]; c.controlPoints[i].w = pts[i][2]; }
    c.knots = { 0,0,0, 0.25,0.25, 0.5,0.5, 0.75,0.75, 1,1,1 };
    c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
NURBSCurve nonuniform() {  // 非均匀节点 + 随机权重
    int deg = 2 + (int)(rnd() * 2), n = deg + 2;
    NURBSCurve c; c.degree = deg; c.controlPoints.resize(n);
    for (int i = 0; i < n; ++i) { c.controlPoints[i].x = rnd(); c.controlPoints[i].y = rnd() - 0.5; c.controlPoints[i].w = 0.5 + rnd(); }
    c.knots.clear();
    for (int i = 0; i <= deg; ++i) c.knots.push_back(0.0);
    c.knots.push_back(0.2 + 0.5 * rnd());
    for (int i = 0; i <= deg; ++i) c.knots.push_back(1.0);
    c.uMin = c.knots[deg]; c.uMax = c.knots[c.knots.size() - 1 - deg]; validateCurve(c); return c;
}
NURBSCurve make(Type t) {
    switch (t) {
        case L: return line0();
        case Q: return bezier(2, 3);
        case C: return bezier(3, 4);
        case H: return bezier(4, 5);
        case A: return arc();
        case E: return ell();
        case K: return kink();
        case O: return closed();
        case NU: return nonuniform();
        default: return line0();
    }
}
} // namespace

int main() {
    printf("==== 综合：NURBS 曲线类型全覆盖 × 跨类型相交（每类型 1000 组）====\n\n");
    const int PER = 1000;
    std::vector<std::vector<NURBSCurve>> pool(NUM_TYPES);
    for (int t = 0; t < NUM_TYPES; ++t) for (int i = 0; i < PER; ++i) pool[t].push_back(make((Type)t));

    SolverOptions opt; opt.recordHistory = false;
    long long totalSols = 0, totalProven = 0;
    // 每类型作 c1 与混合类型 c2 求交 1000 组
    for (int t = 0; t < NUM_TYPES; ++t) {
        long long sols = 0, proven = 0;
        for (int i = 0; i < PER; ++i) {
            int t2 = (t + 1 + i % 8) % NUM_TYPES;      // 跨类型为主，偶有同类
            const NURBSCurve& c1 = pool[t][i];
            const NURBSCurve& c2 = pool[t2][i];
            BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
            auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
            for (auto& r : sol) if (r.state == GeometryState::INTERSECTION) {
                sols++;
                if (proveUniqueInBox(c1, c2, Interval(r.u - 1e-5, r.u + 1e-5), Interval(r.v - 1e-5, r.v + 1e-5))) proven++;
            }
        }
        totalSols += sols; totalProven += proven;
        printf("  %-12s : %4lld 组，求解交点 %4lld，δ=1e-5 证明 %4lld → 覆盖率 %5.1f%%\n",
               TNAME[t], (long long)PER, sols, proven, sols ? 100.0 * proven / sols : 0.0);
    }
    printf("\n  总计 %d 组，横截交点 %lld，证明 %lld → 总覆盖率 %.1f%%\n",
           NUM_TYPES * PER, totalSols, totalProven, totalSols ? 100.0 * totalProven / totalSols : 0.0);
    printf("\n结论：9 类曲线每类型 1000 组（跨类型配对）求解定位 + Krawczyk 证明，覆盖率与正确性汇总如上。\n");
    return 0;
}
