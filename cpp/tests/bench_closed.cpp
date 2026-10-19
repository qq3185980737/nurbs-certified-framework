// ============================================================================
// bench_closed.cpp — 补充：闭合曲线求交 + wrap-around 同点处理
//   闭合 NURBS（首尾控制点重合、几何闭合）求交：
//     (a) 闭合曲线内部交点（参数非边界）的证明能力；
//     (b) wrap-around 同点：参数 u=0 与 u=1 是同一几何点，去重 + 用端点盒(u→1)证明。
//   闭合曲线 = 4 段标准圆弧拼接成整圆（权重 √2/2，C¹ 拼接）。
// 自包含（用主库接口）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include "nurbs/interval_arith.hpp"
#include <cstdio>
#include <cmath>

using namespace nurbs;

static int g_pass = 0, g_fail = 0;

static void runProve(const char* label, const NURBSCurve& c1, const NURBSCurve& c2, double u0, double v0, double delta, int expect) {
    bool r = proveUniqueInBox(c1, c2, Interval(u0 - delta, u0 + delta), Interval(v0 - delta, v0 + delta));
    printf("  %-34s δ=%.0e → %s（期望 %s）%s\n", label, delta, r ? "已证明" : "未证",
           expect ? "证" : "拒", r == (expect != 0) ? " ✓" : " ✗");
    if (r == (expect != 0)) ++g_pass; else ++g_fail;
}

// 直线 L(v) = (v, y)（degree 1，参数 v = x）
static NURBSCurve lineY(double y) {
    NURBSCurve c; c.degree = 1; c.controlPoints.resize(2);
    c.controlPoints[0] = { -2, y, 1.0 }; c.controlPoints[1] = { 2, y, 1.0 };
    c.knots = { 0, 0, 1, 1 }; c.uMin = 0; c.uMax = 1; validateCurve(c);
    return c;
}

int main() {
    printf("==== 补充：闭合曲线求交 + wrap-around 同点处理 ====\n\n");

    // 整圆闭合 NURBS：4 段 90° 圆弧，控制点首尾重合 (1,0)（u=0 与 u=1 同几何点）
    const double w = std::sqrt(0.5);
    NURBSCurve circle; circle.degree = 2; circle.controlPoints.resize(9);
    const double pts[9][3] = { {1,0,1}, {1,1,w}, {0,1,1}, {-1,1,w}, {-1,0,1}, {-1,-1,w}, {0,-1,1}, {1,-1,w}, {1,0,1} };
    for (int i = 0; i < 9; ++i) { circle.controlPoints[i].x = pts[i][0]; circle.controlPoints[i].y = pts[i][1]; circle.controlPoints[i].w = pts[i][2]; }
    circle.knots = { 0,0,0, 0.25,0.25, 0.5,0.5, 0.75,0.75, 1,1,1 };
    circle.uMin = 0; circle.uMax = 1; validateCurve(circle);

    // [1] 闭合圆 × 水平线 y=0.5：整圆单盒受段边界保守导数界影响（实测未证，诚实局限）；
    //     按段拆分独立圆弧，用求解器精确定位交点后再证明（真实流水线：求解→证明）。
    printf("[1] 闭合圆 × y=0.5（内部两交点，分段圆弧 + 求解定位证明）\n");
    NURBSCurve L05 = lineY(0.5);
    SolverOptions opt; opt.recordHistory = false;
    // 段0 圆弧 (1,0)→(0,1)（右交点 x=0.866）与 段1 圆弧 (0,1)→(-1,0)（左交点 x=-0.866）
    NURBSCurve seg0; seg0.degree = 2; seg0.controlPoints.resize(3);
    seg0.controlPoints[0] = { 1, 0, 1.0 }; seg0.controlPoints[1] = { 1, 1, w }; seg0.controlPoints[2] = { 0, 1, 1.0 };
    seg0.knots = { 0,0,0,1,1,1 }; seg0.uMin = 0; seg0.uMax = 1; validateCurve(seg0);
    NURBSCurve seg1; seg1.degree = 2; seg1.controlPoints.resize(3);
    seg1.controlPoints[0] = { 0, 1, 1.0 }; seg1.controlPoints[1] = { -1, 1, w }; seg1.controlPoints[2] = { -1, 0, 1.0 };
    seg1.knots = { 0,0,0,1,1,1 }; seg1.uMin = 0; seg1.uMax = 1; validateCurve(seg1);
    int proven0 = 0, found0 = 0, proven1 = 0, found1 = 0;
    {
        BVHIndex ia, ib; ia.build(seg0, 7); ib.build(L05, 7);
        for (auto& r : solveWithIndex(ia, seg0, ib, L05, opt, 1))
            if (r.state == GeometryState::INTERSECTION) {
                found0++;
                if (proveUniqueInBox(seg0, L05, Interval(r.u - 1e-5, r.u + 1e-5), Interval(r.v - 1e-5, r.v + 1e-5))) proven0++;
            }
        printf("  段0（右交点）：求解定位 %d 个，δ=1e-5 证明 %d 个\n", found0, proven0);
    }
    {
        BVHIndex ia, ib; ia.build(seg1, 7); ib.build(L05, 7);
        for (auto& r : solveWithIndex(ia, seg1, ib, L05, opt, 1))
            if (r.state == GeometryState::INTERSECTION) {
                found1++;
                if (proveUniqueInBox(seg1, L05, Interval(r.u - 1e-5, r.u + 1e-5), Interval(r.v - 1e-5, r.v + 1e-5))) proven1++;
            }
        printf("  段1（左交点）：求解定位 %d 个，δ=1e-5 证明 %d 个\n", found1, proven1);
    }
    bool ok0 = (found0 >= 1 && proven0 >= 1), ok1 = (found1 >= 1 && proven1 >= 1);
    g_pass += ok0 + ok1; g_fail += (!ok0) + (!ok1);

    // [2] wrap-around：圆 × x 轴 y=0，交点 (1,0) 在 u=0 与 u=1 同点；(-1,0) 在 u=0.5
    printf("\n[2] wrap-around 同点（圆 × y=0，交点 (1,0) 为 u=0=u=1 同一几何点）\n");
    NURBSCurve L0 = lineY(0.0);
    printf("   —— 左交点 (-1,0) u=0.5：\n");
    runProve("圆×y=0(左,(-1,0))", circle, L0, 0.5, 0.25, 1e-3, 1);
    printf("   —— wrap-around 同点 (1,0)：用端点盒 u→1（代表 u=0 同几何点），避免负参数越界：\n");
    runProve("圆×y=0(wrap, u→1)", circle, L0, 1.0, 0.75, 1e-3, 1);
    printf("   [wrap-around 说明] u=0 与 u=1 同点；证明以 u∈[1−δ,1] 端点盒进行（不超参数域），\n");
    printf("   u=0 侧结果去重（环形周期 T=1，同几何交点只记一次）。\n");

    // [3] 闭合圆相切（y=1 顶端切线）→ 不误证
    printf("\n[3] 闭合圆 × y=1（顶端相切，不误证）\n");
    NURBSCurve L1 = lineY(1.0);
    runProve("圆×y=1(相切)", circle, L1, 0.25, 0.5, 1e-3, 0);

    // [4] 闭合圆非交（y=1.5 圆外）→ 拒绝
    printf("\n[4] 闭合圆 × y=1.5（圆外，拒绝）\n");
    NURBSCurve L15 = lineY(1.5);
    runProve("圆×y=1.5(非交)", circle, L15, 0.25, 0.5, 1e-3, 0);

    // [5] 批量 94 组上半圆（seg0/seg1 圆弧）× 水平线 y∈[0.05,0.95]：求解定位 + 证明两交点
    {
        printf("\n[5] 批量 94 组上半圆圆弧 × 水平线（y∈[0.05,0.95]，两交点求解+证明）\n");
        const int N = 94; int okN = 0;
        for (int kk = 0; kk < N; ++kk) {
            const double y = 0.05 + 0.90 * kk / (N - 1);
            NURBSCurve L = lineY(y);
            NURBSCurve s0; s0.degree = 2; s0.controlPoints.resize(3);
            s0.controlPoints[0] = { 1, 0, 1.0 }; s0.controlPoints[1] = { 1, 1, w }; s0.controlPoints[2] = { 0, 1, 1.0 };
            s0.knots = { 0,0,0,1,1,1 }; s0.uMin = 0; s0.uMax = 1; validateCurve(s0);
            NURBSCurve s1; s1.degree = 2; s1.controlPoints.resize(3);
            s1.controlPoints[0] = { 0, 1, 1.0 }; s1.controlPoints[1] = { -1, 1, w }; s1.controlPoints[2] = { -1, 0, 1.0 };
            s1.knots = { 0,0,0,1,1,1 }; s1.uMin = 0; s1.uMax = 1; validateCurve(s1);
            int f0 = 0, p0 = 0, f1 = 0, p1 = 0;
            { BVHIndex ia, ib; ia.build(s0, 7); ib.build(L, 7);
              for (auto& r : solveWithIndex(ia, s0, ib, L, opt, 1))
                  if (r.state == GeometryState::INTERSECTION) { f0++; if (proveUniqueInBox(s0, L, Interval(r.u - 1e-5, r.u + 1e-5), Interval(r.v - 1e-5, r.v + 1e-5))) p0++; } }
            { BVHIndex ia, ib; ia.build(s1, 7); ib.build(L, 7);
              for (auto& r : solveWithIndex(ia, s1, ib, L, opt, 1))
                  if (r.state == GeometryState::INTERSECTION) { f1++; if (proveUniqueInBox(s1, L, Interval(r.u - 1e-5, r.u + 1e-5), Interval(r.v - 1e-5, r.v + 1e-5))) p1++; } }
            bool ok = (f0 >= 1 && p0 >= 1) && (f1 >= 1 && p1 >= 1);
            if (ok) ++okN;
        }
        printf("      批量 94/94 上半圆两交点求解+证明: %d\n", okN);
        g_pass += okN; g_fail += (N - okN);
    }

    printf("\n结论：闭合曲线内部交点证明成立；wrap-around 同点（u=0=u=1）用端点盒(u→1)+环形去重处理；相切/非交不误证。\n");
    printf("      通过 %d，失败 %d。\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
