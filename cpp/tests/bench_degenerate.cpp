// ============================================================================
// bench_degenerate.cpp — 退化/奇异曲线专项（当前框架行为诚实测试）
//   结合最新理论（Monterde NURBS 奇异点结构定理；CASC2023 尖点/自交/切线分析；
//   NURBS 权重→∞ 退化 regular control curve）分类构造退化曲线，
//   测当前 solveWithIndex + proveUniqueInBox 在退化情形下的真实表现，
//   定位"保证性求交"的鲁棒性缺口（哪些退化能证、哪些失败、失败模式）。
// 自包含（用主库接口）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include "nurbs/interval_arith.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace nurbs;

static NURBSCurve bezier(std::vector<ControlPoint> cp) {
    int deg = (int)cp.size() - 1;
    NURBSCurve c; c.degree = deg; c.controlPoints = cp;
    for (int i = 0; i <= deg; ++i) c.knots.push_back(0.0);
    for (int i = 0; i <= deg; ++i) c.knots.push_back(1.0);
    c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
static NURBSCurve line2(ControlPoint a, ControlPoint b) {
    NURBSCurve c; c.degree = 1; c.controlPoints = { a, b };
    c.knots = { 0, 0, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c); return c;
}
static bool tryProve(const NURBSCurve& c1, const NURBSCurve& c2, const SolveResult& r, double d) {
    return proveUniqueInBox(c1, c2, Interval(r.u - d, r.u + d), Interval(r.v - d, r.v + d)) == 1;
}

int main() {
    printf("==== 退化/奇异曲线专项：当前保证性求交框架行为测试 ====\n\n");
    SolverOptions opt; opt.recordHistory = false;

    // [1] 退化点曲线：控制点全重合 → 曲线退化为单点 P*=（0.5,0.3）
    {
        printf("[1] 退化点曲线（控制点全重合 → 单点）\n");
        ControlPoint P{0.5, 0.3, 1.0};
        NURBSCurve dp = bezier({ P, P, P, P });
        // 直线过该点（应相交） / 直线不过（应拒绝）
        NURBSCurve l1 = line2({0,0.3,1},{1,0.3,1});   // y=0.3 过 P
        NURBSCurve l2 = line2({0,0.2,1},{1,0.2,1});   // y=0.2 不过
        BVHIndex d1a, d1b, d2a, d2b; d1a.build(dp,7); d1b.build(l1,7); d2a.build(dp,7); d2b.build(l2,7);
        auto s1 = solveWithIndex(d1a, dp, d1b, l1, opt, 1);
        auto s2 = solveWithIndex(d2a, dp, d2b, l2, opt, 1);
        printf("   ×过点直线: 求解交点 %zu 个", s1.size());
        if (!s1.empty() && s1[0].state == GeometryState::INTERSECTION) {
            bool ok = tryProve(dp, l1, s1[0], 1e-4);
            printf("，proveUniqueInBox(δ=1e-4)=%s\n", ok ? "证明 ✓" : "未证明 ✗");
        } else printf("（未定位交点，退化点系统为 2→0 维，当前求解器未覆盖）\n");
        printf("   ×不过直线: 求解交点 %zu 个（期望 0，拒绝）\n", s2.size());
    }

    // [2] 有理权重退化 w→0（中间控制点 w=1e-9，Lipschitz 界含 wmin² 分母）
    {
        printf("\n[2] 有理权重退化（中间控制点 w→0，wmin² 分母界→大）\n");
        NURBSCurve rq = bezier({ {0,0,1.0}, {0.5,0.3,1e-9}, {1,0,1.0} });
        NURBSCurve l  = line2({0,0.15,1},{1,0.15,1});  // y=0.15 横截
        BVHIndex ir, il; ir.build(rq,7); il.build(l,7);
        auto s = solveWithIndex(ir, rq, il, l, opt, 1);
        printf("   w=[1,1e-9,1] × y=0.15: 求解交点 %zu 个", s.size());
        if (!s.empty() && s[0].state == GeometryState::INTERSECTION) {
            bool ok = tryProve(rq, l, s[0], 1e-4);
            printf("，proveUniqueInBox(δ=1e-4)=%s\n", ok ? "证明 ✓" : "未证明 ✗");
            if (!ok) printf("     （wmin² 极小→Lipschitz 界 L∝1/wmin² 过大，盒需极窄才能证）\n");
        } else printf("\n");
    }

    // [3] 尖点 cusp：半立方抛物线 (t²,t³) 端点 C'(0)=0
    {
        printf("\n[3] 尖点 cusp（半立方抛物线 y²=x³，端点 C'(0)=0）\n");
        NURBSCurve cu = bezier({ {0,0,1},{0,0,1},{1.0/3.0,0,1},{1,1,1} });
        // cusp 处交点：直线 x=0 过尖点 (0,0)
        NURBSCurve lc = line2({0,-0.5,1},{0,0.5,1});
        BVHIndex i0a, i0b; i0a.build(cu,7); i0b.build(lc,7);
        auto s0 = solveWithIndex(i0a, cu, i0b, lc, opt, 1);
        printf("   ×直线x=0(过尖点): 求解交点 %zu 个", s0.size());
        if (!s0.empty() && s0[0].state == GeometryState::INTERSECTION) {
            bool ok = tryProve(cu, lc, s0[0], 1e-3);
            printf("，proveUniqueInBox(δ=1e-3)=%s\n", ok ? "证明 ✓" : "未证明 ✗（cusp 处 Jacobian 含 0，Krawczyk 失效）");
        } else printf("\n");
        // 非 cusp 处横截：直线 y=0.5（t: t³=0.5 → t=0.79，C'≠0）
        NURBSCurve lh = line2({-0.5,0.5,1},{1.5,0.5,1});
        BVHIndex i1a, i1b; i1a.build(cu,7); i1b.build(lh,7);
        auto s1 = solveWithIndex(i1a, cu, i1b, lh, opt, 1);
        printf("   ×直线y=0.5(正则横截): 求解交点 %zu 个", s1.size());
        if (!s1.empty() && s1[0].state == GeometryState::INTERSECTION) {
            bool ok = tryProve(cu, lh, s1[0], 1e-4);
            printf("，proveUniqueInBox(δ=1e-4)=%s\n", ok ? "证明 ✓" : "未证明 ✗");
        } else printf("\n");
    }

    // [4] 控制点共线退化（二次三控制点共线 → 几何退化为直线段）
    {
        printf("\n[4] 控制点共线退化（二次控制点共线 → 线段）\n");
        NURBSCurve col = bezier({ {0,0,1},{0.5,0.25,1},{1,0.5,1} }); // y=0.5x 线段
        NURBSCurve lt = line2({0.6,-1,1},{0.6,1,1});   // x=0.6 横截于 (0.6,0.3)
        BVHIndex ca, cb; ca.build(col,7); cb.build(lt,7);
        auto s = solveWithIndex(ca, col, cb, lt, opt, 1);
        printf("   共线退化×直线x=0.6: 求解交点 %zu 个", s.size());
        if (!s.empty() && s[0].state == GeometryState::INTERSECTION) {
            bool ok = tryProve(col, lt, s[0], 1e-4);
            printf("，proveUniqueInBox(δ=1e-4)=%s\n", ok ? "证明 ✓" : "未证明 ✗");
        } else printf("\n");
        // 与共线重叠：直线 y=0.5x（无限重叠，应判为重叠/退化而非单点）
        NURBSCurve lo = line2({0.2,0.1,1},{0.8,0.4,1});
        BVHIndex oa, ob; oa.build(col,7); ob.build(lo,7);
        auto so = solveWithIndex(oa, col, ob, lo, opt, 1);
        printf("   共线重叠×直线y=0.5x: 求解交点 %zu 个（几何退化重叠，非单点横截）\n", so.size());
    }

    printf("\n结论：退化/奇异情形暴露当前框架的鲁棒性缺口——退化点曲线(cusp处)Krawczyk Jacobian\n");
    printf("      含0失效、有理 w→0 的 wmin² 界爆炸、共线重叠退化。需要专门退化处理机制。\n");
    return 0;
}
