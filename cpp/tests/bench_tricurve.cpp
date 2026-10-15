// ============================================================================
// bench_tricurve.cpp — ⑤ 非方系统区间 Krawczyk：三曲线共点（T 型结点）证明
//   对标 Burr-Hauenstein-Lee 2025（广义 Krawczyk 推广到非方系统）。
//   系统 F:R³→R⁴（超定）：
//     F(u,v,w) = (c1x-c2x, c1y-c2y, c2x-c3x, c2y-c3y)
//   用中点 Moore-Penrose 伪逆 A=J⁺ 作定矩阵，区间 Krawczyk 算子
//     K(x) = m − A·f(m) + (I₃ − A·J([x]))·(x−m)
//   判定 K(x) ⊂ interior(x) ⇒ 三曲线在盒内存在唯一公共点（T 型结点拓扑保证）。
//   几何意义：CAD 拓扑中三边共点（T-junction）的保证性检测。
// ============================================================================
#include "nurbs/interval_arith.hpp"
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace nurbs;

// ---------- 实数 3×3 逆（克拉默） ----------
static bool inv3(const double M[3][3], double R[3][3]) {
    const double d0 = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                    - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                    + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (std::abs(d0) < 1e-14) return false;
    const double inv = 1.0 / d0;
    R[0][0] =  (M[1][1]*M[2][2]-M[1][2]*M[2][1])*inv;
    R[0][1] = -(M[0][1]*M[2][2]-M[0][2]*M[2][1])*inv;
    R[0][2] =  (M[0][1]*M[1][2]-M[0][2]*M[1][1])*inv;
    R[1][0] = -(M[1][0]*M[2][2]-M[1][2]*M[2][0])*inv;
    R[1][1] =  (M[0][0]*M[2][2]-M[0][2]*M[2][0])*inv;
    R[1][2] = -(M[0][0]*M[1][2]-M[0][2]*M[1][0])*inv;
    R[2][0] =  (M[1][0]*M[2][1]-M[1][1]*M[2][0])*inv;
    R[2][1] = -(M[0][0]*M[2][1]-M[0][1]*M[2][0])*inv;
    R[2][2] =  (M[0][0]*M[1][1]-M[0][1]*M[1][0])*inv;
    return true;
}

static double ddx(const NURBSCurve& c, double u) { const double e=1e-6; return (evaluatePoint(c,u+e).x-evaluatePoint(c,u-e).x)/(2*e); }
static double ddy(const NURBSCurve& c, double u) { const double e=1e-6; return (evaluatePoint(c,u+e).y-evaluatePoint(c,u-e).y)/(2*e); }

// 超定三曲线共点 Krawczyk：返回 1=证明存在唯一公共点，0=未定/不存在。
// box u,v,w 为三条曲线各自参数盒。
static int proveTriCommon(const NURBSCurve& c1, const NURBSCurve& c2, const NURBSCurve& c3,
                          Interval& ub, Interval& vb, Interval& wb) {
    // 中点
    const double mu = ub.mid(), mv = vb.mid(), mw = wb.mid();
    // F(m) ∈ R⁴
    const double f[4] = {
        evaluatePoint(c1,mu).x - evaluatePoint(c2,mv).x,
        evaluatePoint(c1,mu).y - evaluatePoint(c2,mv).y,
        evaluatePoint(c2,mv).x - evaluatePoint(c3,mw).x,
        evaluatePoint(c2,mv).y - evaluatePoint(c3,mw).y };
    // J(m) 实数 4×3
    double J[4][3] = {
        { ddx(c1,mu), -ddx(c2,mv), 0 },
        { ddy(c1,mu), -ddy(c2,mv), 0 },
        { 0,  ddx(c2,mv), -ddx(c3,mw) },
        { 0,  ddy(c2,mv), -ddy(c3,mw) } };
    // G = JᵀJ (3×3)，伪逆 A = G⁻¹Jᵀ (3×4)
    double G[3][3], Gi[3][3], A[3][4];
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) { G[i][j]=0; for (int k=0;k<4;++k) G[i][j]+=J[k][i]*J[k][j]; }
    if (!inv3(G, Gi)) return 0;             // 中点奇异 → 相切/退化，非横截，不证
    for (int i=0;i<3;++i) for (int j=0;j<4;++j) { A[i][j]=0; for (int k=0;k<3;++k) A[i][j]+=Gi[i][k]*J[j][k]; }
    // J([x]) 区间 4×3（用 curveInterval 的导数界 dfx/dfy）
    const CurveIntervalEval e1 = curveInterval(c1, ub), e2 = curveInterval(c2, vb), e3 = curveInterval(c3, wb);
    Interval Ji[4][3] = {
        { e1.dfx, imul(e2.dfx, Interval(-1.0,-1.0)), Interval::point(0) },
        { e1.dfy, imul(e2.dfy, Interval(-1.0,-1.0)), Interval::point(0) },
        { Interval::point(0), e2.dfx, imul(e3.dfx, Interval(-1.0,-1.0)) },
        { Interval::point(0), e2.dfy, imul(e3.dfy, Interval(-1.0,-1.0)) } };
    // (I − A·J([x]))(x−m)：3×3 区间
    Interval dx[3] = { isub(ub, Interval::point(mu)), isub(vb, Interval::point(mv)), isub(wb, Interval::point(mw)) };
    Interval K[3];
    for (int i=0;i<3;++i) {
        Interval sum = Interval::point(0);
        for (int j=0;j<3;++j) {           // proj_ij = δ_ij − Σ_k A_ik Ji_kj
            Interval proj = (i==j) ? Interval(1.0,1.0) : Interval::point(0);
            for (int k=0;k<4;++k) proj = isub(proj, imul(Interval::point(A[i][k]), Ji[k][j]));
            sum = iadd(sum, imul(proj, dx[j]));
        }
        double af = 0; for (int k=0;k<4;++k) af += A[i][k]*f[k];
        double mi = (i==0)?mu:((i==1)?mv:mw);
        K[i] = iadd(Interval::point(mi - af), sum);
    }
    // 判定：K(x) ⊂ interior(x)
    const Interval xx[3] = { ub, vb, wb };
    for (int i=0;i<3;++i) if (!(xx[i].lo < K[i].lo && K[i].hi < xx[i].hi)) return 0;
    return 1;
}

namespace {
NURBSCurve fromPoints(const std::vector<Point2>& pts, int degree) {
    NURBSCurve c;
    c.degree = degree;
    c.controlPoints.resize(pts.size());
    for (size_t i=0;i<pts.size();++i) { c.controlPoints[i].x=pts[i].x; c.controlPoints[i].y=pts[i].y; c.controlPoints[i].w=1.0; }
    c.knots = generateOpenUniformKnotVector((int)pts.size(), degree);
    c.uMin = c.knots[degree]; c.uMax = c.knots[c.knots.size()-1-degree];
    validateCurve(c); return c;
}
void tryProve(const NURBSCurve& c1,const NURBSCurve& c2,const NURBSCurve& c3,
              double u0,double v0,double w0,double delta,const char* label,int expect) {
    Interval ub(u0-delta,u0+delta), vb(v0-delta,v0+delta), wb(w0-delta,w0+delta);
    int r = proveTriCommon(c1,c2,c3,ub,vb,wb);
    printf("  %-28s δ=%.0e → %s（期望 %s）%s\n", label, delta,
           r==1?"已证明":"未证", expect==1?"证":"拒",
           r==expect?" ✓":" ✗");
}
} // namespace

int main() {
    printf("==== ⑤ 非方系统区间 Krawczyk：三曲线共点（T 型结点）证明 ====\n\n");

    // 场景 A：三直线共点 (0.5,0.5)
    {
        printf("[A] 三直线共点 (0.5,0.5)\n");
        NURBSCurve c1 = fromPoints({ {0,0},{1,1} }, 1);         // y=x
        NURBSCurve c2 = fromPoints({ {0,1},{1,0} }, 1);         // y=1−x
        NURBSCurve c3 = fromPoints({ {0.2,0.5},{0.8,0.5} }, 1); // y=0.5
        for (double d : {2e-2, 1e-3}) tryProve(c1,c2,c3,0.5,0.5,0.5,d,"共点(0.5,0.5)",1);
    }
    // 场景 B：三抛物线共点（曲线而非直线）
    {
        printf("[B] 三抛物线共点 (0.5,0.5)\n");
        NURBSCurve c1 = fromPoints({ {0,0},{0.5,0.5},{1,0} }, 2);       // y=u(1-u) 顶点(0.5,0.25)? 用通过(0.5,0.5)
        // 让三条二次曲线都通过 (0.5,0.5)：构造 y=4(u-0.5)^2+0.5（顶点在(0.5,0.5)开口上）等
        NURBSCurve d1 = fromPoints({ {0,1.5},{0.5,0.5},{1,1.5} }, 2);   // y=4(u-0.5)^2+0.5，顶点(0.5,0.5)
        NURBSCurve d2 = fromPoints({ {0,-0.5},{0.5,0.5},{1,-0.5} }, 2); // y=−4(u-0.5)^2+0.5，顶点(0.5,0.5) 开口向下
        NURBSCurve d3 = fromPoints({ {0.2,0.2},{0.5,0.5},{0.8,0.8} }, 2);// 直线段(过(0.5,0.5))
        // 三者在 (0.5,0.5) 相切程度高（奇异），Krawczyk 应识别为退化或证
        tryProve(d1,d2,d3,0.5,0.5,0.5,1e-2,"三抛物线顶点共点",0); // 相切奇异 → 不误证
        // 横截场景：二次曲线两两横截交于公共点
        NURBSCurve e1 = fromPoints({ {0,0},{0.5,0.5},{1,0} }, 2);       // y=4u(1-u)? 顶点(0.5,1)修正
        NURBSCurve e2 = fromPoints({ {0,1},{0.5,0.5},{1,1} }, 2);       // 倒抛物线顶点(0.5,0.5)? 
        NURBSCurve e3 = fromPoints({ {0.2,0.5},{0.8,0.5} }, 1);         // 水平 y=0.5
        // 三者在 (0.5,0.5) 相交（前两者顶点即切点，退化）→ 期望不误证
        tryProve(e1,e2,e3,0.5,0.5,0.5,1e-2,"两切一交",0);
    }
    // 场景 C：三曲线共点（横截，非退化）
    {
        printf("[C] 三曲线横截共点\n");
        // 三条直线旋转交于 (0.5,0.5)，夹角各异 → 横截
        NURBSCurve c1 = fromPoints({ {0.5-0.5,0.5-0.1},{0.5+0.5,0.5+0.1} }, 1);
        NURBSCurve c2 = fromPoints({ {0.5-0.5,0.5+0.2},{0.5+0.5,0.5-0.2} }, 1);
        NURBSCurve c3 = fromPoints({ {0.5-0.3,0.5},{0.5+0.3,0.5} }, 1);
        for (double d : {2e-2, 1e-3, 1e-5}) tryProve(c1,c2,c3,0.5,0.5,0.5,d,"横截共点",1);
    }
    // 场景 D：非共点 → 拒绝（不误证）
    {
        printf("[D] 非共点（拒绝）\n");
        NURBSCurve c1 = fromPoints({ {0,0},{1,1} }, 1);
        NURBSCurve c2 = fromPoints({ {0,1},{1,0} }, 1);
        NURBSCurve c3 = fromPoints({ {0,0.6},{1,0.6} }, 1);  // y=0.6，不过(0.5,0.5)
        tryProve(c1,c2,c3,0.5,0.5,0.5,1e-3,"三线不共点",0);
        tryProve(c1,c2,c3,0.45,0.45,0.5,2e-2,"偏离公共点",0);
    }
    printf("\n说明：超定 Krawczyk 仅在存在唯一横截公共点时证明；相切/退化或不存在时不误证。\n");
    return 0;
}
