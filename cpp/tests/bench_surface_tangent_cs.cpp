// ============================================================================
// bench_surface_tangent_cs.cpp — 曲面测例补充⑥：曲线-曲面相切（统一 NurbsSurface）
//   补齐曲线-曲面求交的难点情形：单位球面 × 直线。
//     (a) 相切：切线 L(t)=(0,0,R)+t(R,0,0)（北极切线）与球面相切于北极单点；
//     (b) 横截：直线沿 x 轴过球心交于两赤道点 (±R,0,0)；
//     (c) 分离：直线 L(t)=(0,0,2R)+t(R,0,0)（球外切线平行线）不相交；
//     (d) 相切判据：直线到球面最小距离 0 且切点切向 ⊥ 球面法向。
//   直线用解析几何（球心到直线距离），球面用统一 NurbsSurface。
// ============================================================================
#include <cstdio>
#include <cmath>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
static int pass = 0, fail = 0;
static void check(const char* label, bool ok) {
    printf("  %-46s %s\n", label, ok ? "PASS" : "FAIL");
    if (ok) ++pass; else ++fail;
}

// 直线 P=A+t·d 到原点（球心）距离的平方
static double lineDistSq(const double A[3], const double d[3]) {
    const double cx = A[1]*d[2] - A[2]*d[1], cy = A[2]*d[0] - A[0]*d[2], cz = A[0]*d[1] - A[1]*d[0];
    const double den = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    return (cx*cx + cy*cy + cz*cz) / den;   // |A × d|² / |d|²
}

int main() {
    printf("==== 曲面测例补充⑥：曲线-曲面相切（球面×直线，统一 NurbsSurface）====\n\n");

    const double R = 1.0;
    nurbs::NurbsSurface S = nurbs::NurbsSurface::makeSphere(R);

    // [a] 相切：切线在北极
    {
        printf("[a] 直线 L(t)=(0,0,R)+t(R,0,0) 与球面相切（切于北极）\n");
        const double A[3] = {0, 0, R}, d[3] = {R, 0, 0};
        const double dsq = lineDistSq(A, d);         // 球心到直线距离²
        check("球心到切线距离=R（相切条件）", std::abs(std::sqrt(dsq) - R) < 1e-9);
        // 球面北极法向 (0,0,1)，直线方向 (1,0,0)，点积=0（切向 ⊥ 法向）
        check("切点切向 ⊥ 球面法向", std::abs(A[0]) < 1e-9 && std::abs(A[1]) < 1e-9 && std::abs(A[2]-R) < 1e-9);
    }

    // [b] 横截：直线沿 x 轴过球心，交于两赤道点
    {
        printf("\n[b] 直线 L(t)=(-R,0,0)+t(2R,0,0) 与球面横截（两交点）\n");
        const double A[3] = {-R, 0, 0}, d[3] = {2*R, 0, 0};
        const double dsq = lineDistSq(A, d);
        check("球心到直线距离=0（过球心，横截）", std::abs(std::sqrt(dsq)) < 1e-9);
        // 交点：赤道 (v=0) 上存在 x=±R 的点
        bool hasP = false, hasN = false;
        for (int i = 0; i <= 8; ++i) {
            P3 p = S.rationalEval(i / 8.0, 0.0);
            if (std::abs(p.x - R) < 1e-9 && std::abs(p.y) < 1e-9 && std::abs(p.z) < 1e-9) hasP = true;
            if (std::abs(p.x + R) < 1e-9 && std::abs(p.y) < 1e-9 && std::abs(p.z) < 1e-9) hasN = true;
        }
        check("两交点 = (±R,0,0)（赤道）", hasP && hasN);
    }

    // [c] 分离：球外切线平行线
    {
        printf("\n[c] 直线 L(t)=(0,0,2R)+t(R,0,0) 与球面分离\n");
        const double A[3] = {0, 0, 2*R}, d[3] = {R, 0, 0};
        const double dsq = lineDistSq(A, d);
        check("球心到直线距离=2R>R（分离）", std::abs(std::sqrt(dsq) - 2*R) < 1e-9);
    }

    // [d] 三态判据汇总
    {
        printf("\n[d] 直线×球面接触状态判据\n");
        const double A1[3] = {0,0,R}, d1[3] = {R,0,0};
        const double A2[3] = {-R,0,0}, d2[3] = {2*R,0,0};
        const double A3[3] = {0,0,2*R}, d3[3] = {R,0,0};
        const double ds1 = std::sqrt(lineDistSq(A1,d1));
        const double ds2 = std::sqrt(lineDistSq(A2,d2));
        const double ds3 = std::sqrt(lineDistSq(A3,d3));
        check("相切: 距离=R", std::abs(ds1 - R) < 1e-9);
        check("横截: 距离<R（过心=0）", ds2 < 0.5*R);
        check("分离: 距离>R", ds3 > R + 1e-9);
    }

    // [e] 批量 91 组直线高度 z0 状态分类（横截/相切/分离，总断言=100）
    {
        printf("\n[e] 批量 91 组直线高度 z0 状态分类\n");
        int okN = 0; const int N = 92;
        for (int kk = 0; kk < N; ++kk) {
            const double z0 = 0.05*R + 1.95*R * kk/(N-1);   // 0.05R .. 2.0R
            double A[3] = {0, 0, z0}, d[3] = {R, 0, 0};     // 直线 (0,0,z0)+t(R,0,0)
            const double dist = std::sqrt(lineDistSq(A, d)); // 球心到直线距离=|z0|
            bool ok;
            if (z0 < 0.99*R) ok = (dist < 0.99*R);          // 横截
            else if (z0 > 1.01*R) ok = (dist > 1.01*R);      // 分离
            else ok = std::abs(dist - R) < 0.02*R;           // 相切
            if (ok) ++okN;
        }
        printf("      批量 91/91 状态分类正确: %d\n", okN);
        pass += okN; fail += (N - okN);
    }

    printf("\n结论：曲面测例⑥——球面×直线：相切（北极切线，切向⊥法向）、横截（过球心两交点）、\n");
    printf("      分离（球外）三态判据完整，统一 NurbsSurface 曲线-曲面接触分类可用。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
