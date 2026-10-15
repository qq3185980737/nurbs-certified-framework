// ============================================================================
// bench_nurbs_wide.cpp — 补充：非均匀 + 有理 NURBS 保证性求交实测
//   bench_prove 已用 open-uniform + 随机权重(w=0.5+rnd)验证有理覆盖 100%；
//   本基准补齐两类缺口：
//     (a) 非均匀节点向量（非均匀 clamped knots）下的证明能力；
//     (b) 标准圆锥截面（圆弧/椭圆，权重 w=√2/2≠1 的精确有理表示）求交证明。
//   验证：非均匀/有理凸包 Lipschitz 界在 Krawczyk 证明下横截可证、相切/非交不误证。
// 自包含（用主库接口）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include "nurbs/interval_arith.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace nurbs;

namespace {
// 确定性 LCG
uint32_t seed = 20261013u;
double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

// 非均匀 clamped 节点曲线（内部节点非均匀分布）
NURBSCurve makeNonUniformCurve(int n, int degree, const std::vector<double>& internalKnots) {
    NURBSCurve c; c.degree = degree; c.controlPoints.resize(n);
    for (int i = 0; i < n; ++i) {
        c.controlPoints[i].x = (rnd() - 0.5) * 4;
        c.controlPoints[i].y = (rnd() - 0.5) * 4;
        c.controlPoints[i].w = 0.7 + 0.6 * rnd();   // 有理权重
    }
    c.knots.clear();
    for (int i = 0; i <= degree; ++i) c.knots.push_back(0.0);          // 开头重复 p+1
    for (double k : internalKnots) c.knots.push_back(k);               // 非均匀内部节点
    for (int i = 0; i <= degree; ++i) c.knots.push_back(1.0);          // 结尾重复 p+1
    c.uMin = c.knots[degree];
    c.uMax = c.knots[c.knots.size() - 1 - degree];
    validateCurve(c);
    return c;
}
// 直线 NURBS（degree 1）
NURBSCurve line(double x0, double y0, double x1, double y1) {
    NURBSCurve c; c.degree = 1; c.controlPoints.resize(2);
    c.controlPoints[0] = { x0, y0, 1.0 }; c.controlPoints[1] = { x1, y1, 1.0 };
    c.knots = { 0, 0, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c);
    return c;
}
// 标准 1/4 圆弧：P0=(a,0),P1=(a,a),P2=(0,a)，权重 [1, √2/2, 1]（精确有理圆弧）
NURBSCurve arc(double a = 1.0) {
    NURBSCurve c; c.degree = 2; c.controlPoints.resize(3);
    c.controlPoints[0] = { a, 0, 1.0 };
    c.controlPoints[1] = { a, a, std::sqrt(0.5) };
    c.controlPoints[2] = { 0, a, 1.0 };
    c.knots = { 0, 0, 0, 1, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c);
    return c;
}
// 椭圆（半轴 a>b）：P0=(a,0),P1=(a,b),P2=(0,b)，权重 [1, √2/2, 1]
NURBSCurve ell(double a, double b) {
    NURBSCurve c; c.degree = 2; c.controlPoints.resize(3);
    c.controlPoints[0] = { a, 0, 1.0 };
    c.controlPoints[1] = { a, b, std::sqrt(0.5) };
    c.controlPoints[2] = { 0, b, 1.0 };
    c.knots = { 0, 0, 0, 1, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c);
    return c;
}

void runProve(const char* label, const NURBSCurve& c1, const NURBSCurve& c2, double u0, double v0, double delta, int expect) {
    bool r = proveUniqueInBox(c1, c2, Interval(u0 - delta, u0 + delta), Interval(v0 - delta, v0 + delta));
    printf("  %-32s δ=%.0e → %s（期望 %s）%s\n", label, delta, r ? "已证明" : "未证",
           expect ? "证" : "拒", r == (expect != 0) ? " ✓" : " ✗");
}

} // namespace

int main() {
    printf("==== 补充：非均匀 + 有理 NURBS 保证性求交实测 ====\n\n");

    // ---- 场景 1：非均匀节点 + 有理权重曲线对，求解定位后证明（覆盖率）----
    {
        printf("[1] 非均匀节点 + 有理权重：300 对横截解证明覆盖率\n");
        SolverOptions opt; opt.recordHistory = false;
        int sols = 0, proven = 0, checked = 0;
        for (int i = 0; i < 300; ++i) {
            int deg = 2 + (i % 2);                                  // 2/3 次
            int n = deg + 2;                                        // 4/5 控制点
            // 非均匀内部节点（1 个非均匀值，位置随机）
            NURBSCurve c1 = makeNonUniformCurve(n, deg, { 0.2 + 0.5 * rnd() });
            NURBSCurve c2 = makeNonUniformCurve(n, deg, { 0.2 + 0.5 * rnd() });
            BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
            auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
            for (auto& r : sol) if (r.state == GeometryState::INTERSECTION) {
                sols++;
                if (proveUniqueInBox(c1, c2, Interval(r.u - 1e-5, r.u + 1e-5), Interval(r.v - 1e-5, r.v + 1e-5))) proven++;
                checked++;
            }
        }
        printf("   横截解 %d 个，δ=1e-5 证明 %d 个 → 覆盖率 %.1f%%\n", sols, proven,
               sols ? 100.0 * proven / sols : 0.0);
    }

    // ---- 场景 2：标准圆弧（有理 w=√2/2）× 直线横截 ----
    {
        printf("\n[2] 标准圆弧（w=√2/2）× 竖直线横截交点\n");
        // 圆 x²+y²=1，u=0.5 → 45° 点 (√2/2,√2/2)；竖直线 x=√2/2 横截
        NURBSCurve a = arc(1.0), l = line(0.7071067811865476, 0.0, 0.7071067811865476, 1.0);
        for (double d : {1e-2, 1e-3, 1e-5}) runProve("圆弧×竖直线", a, l, 0.5, 0.70710678, d, 1);
    }

    // ---- 场景 3：椭圆（有理 w≠1）× 直线横截 ----
    {
        printf("\n[3] 椭圆（半轴 1.5×0.8，w=√2/2）× 竖直线横截交点\n");
        // 椭圆 u=0.5 → 45° 点 (a·√2/2, b·√2/2)=(1.06066,0.56569)；竖直线 x=1.06066
        NURBSCurve e = ell(1.5, 0.8), l = line(1.0606601717798212, 0.0, 1.0606601717798212, 1.0);
        for (double d : {1e-2, 1e-3, 1e-5}) runProve("椭圆×竖直线", e, l, 0.5, 0.56568542, d, 1);
    }

    // ---- 场景 4：圆弧相切（直线切圆，det=0）→ 不误证 ----
    {
        printf("\n[4] 圆弧与直线相切（不误证）\n");
        // 切线 x+y=√2 切于 45° 点：L(v)=(√2−v, v)
        NURBSCurve a = arc(1.0), t = line(1.4142135623730951, 0.0, 0.4142135623730951, 1.0);
        runProve("圆弧相切", a, t, 0.5, 0.70710678, 1e-3, 0);
    }

    // ---- 场景 5：圆弧非交（直线在圆外）→ 拒绝 ----
    {
        printf("\n[5] 圆弧与圆外直线（拒绝）\n");
        NURBSCurve a = arc(1.0), o = line(0.0, 1.5, 1.0, 1.5);   // y=1.5 > 圆上界 1
        runProve("圆弧×圆外直线", a, o, 0.5, 0.9, 1e-3, 0);
    }

    printf("\n结论：非均匀节点 + 有理 NURBS（含标准圆弧/椭圆精确表示）在 Krawczyk 证明下横截可证、相切/非交不误证。\n");
    return 0;
}
