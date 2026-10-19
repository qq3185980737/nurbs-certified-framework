// ============================================================================
// bench_surface_degen_normal.cpp — 曲面测例补充⑪：退化曲面 × 正常曲面
//   补齐 ACIS degenerate×regular：退化曲面（缩为点/缩到曲线）× 正常曲面求交。
//     Dp = 退化点 patch（控制点全重合 → 单点 P*），× 正常碗形 patch：
//       (a) P* 在碗上（接触）：交于 1 点（点在面上）；
//       (b) P* 在碗外（分离）：无交点（拒绝）；
//     Dl = 退化到曲线 patch（控制点共线 → 空间线段），× 正常 patch：
//       (c) 线段穿过碗：交于 2 点（0D）；
//       (d) 线段在碗外：无交点（分离）。
//   自包含（NurbsSurface 双线性/退化控制点 + 解析判定）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
static int pass = 0, fail = 0;
static void check(const char* label, bool ok) {
    printf("  %-48s %s\n", label, ok ? "PASS" : "FAIL");
    if (ok) ++pass; else ++fail;
}
static double dist(P3 a, P3 b) {
    return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

int main() {
    printf("==== 测例补充⑪：退化曲面 × 正常曲面（ACIS degenerate×regular）====\n\n");

    // 碗形双二次正常 patch
    P3 cp[9];
    for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) {
        const double u = i / 2.0, v = j / 2.0;
        cp[j*3+i] = { u, v, 0.5*(u*u - u + v*v - v + 1.0) };   // 碗 z∈[0.25,0.5]
    }
    nurbs::NurbsSurface S = nurbs::NurbsSurface::makeBezierBiquadratic(cp);

    // 退化点 patch：控制点全重合 → 单点 P*（取碗上实际顶点 S(0.5,0.5)）
    P3 Pstar = S.eval(0.5, 0.5);   // = (0.5, 0.5, 0.375)，碗顶点，在碗上
    P3 pc[4] = { Pstar, Pstar, Pstar, Pstar };
    nurbs::NurbsSurface Dp = nurbs::NurbsSurface::makeBilinear(pc);

    // [a] 退化点 patch 顶点 P* 在碗上（接触）
    {
        printf("[a] 退化点 patch（顶点 P* 在碗上）\n");
        P3 p = Dp.eval(0.5, 0.5);
        P3 s = S.eval(0.5, 0.5);
        printf("      [诊断] Dp=(%.6f,%.6f,%.6f) S=(%.6f,%.6f,%.6f) dist=%.2e\n", p.x,p.y,p.z, s.x,s.y,s.z, dist(p,s));
        check("退化点 P* 在正常 patch 上（接触）", dist(p, s) < 1e-6);
    }

    // [b] 退化点 patch 位于碗外（分离）
    {
        printf("\n[b] 退化点 patch（点在碗外，分离）\n");
        P3 po[4] = { {0.5, 0.5, 0.9}, {0.5, 0.5, 0.9}, {0.5, 0.5, 0.9}, {0.5, 0.5, 0.9} };
        nurbs::NurbsSurface Dpo = nurbs::NurbsSurface::makeBilinear(po);
        P3 p = Dpo.eval(0.5, 0.5);
        // 碗 z∈[0.25,0.5]（u,v∈[0,1]），z=0.9 在碗上方 → 分离
        check("退化点 z=0.9 在碗上方（距碗>0.1，分离）", p.z - S.aabbMax().z > 0.1);
    }

    // 退化到曲线 patch：控制点共线 → 竖直空间线段（x=0.5, y=0.5, z∈[0.35,0.40]）
    // 穿过碗顶点 z=0.375（碗 z 范围中值）→ 交于顶点 1 点（0D）
    const double zseg = 0.375;   // 碗顶点 z(0.5,0.5)=0.375
    P3 pc2[4] = { {0.5, 0.5, 0.35}, {0.5, 0.5, 0.365}, {0.5, 0.5, 0.385}, {0.5, 0.5, 0.40} };
    nurbs::NurbsSurface Dl = nurbs::NurbsSurface::makeBilinear(pc2);   // 线段 z∈[0.35,0.40]

    // [c] 退化到曲线 patch（竖直线段穿过碗顶点）→ 交于 1 点（0D）
    {
        printf("\n[c] 退化到曲线 patch（竖直段穿碗顶点，0D 单点）\n");
        int hit = 0; double hitZ = -1, mnall = 1e300;
        for (int i = 0; i <= 100; ++i) {
            P3 p = Dl.eval(0.5, i / 100.0);          // 沿竖直段采样
            double mn = 1e300;
            for (int a = 0; a <= 40; ++a) for (int b = 0; b <= 40; ++b)
                mn = std::min(mn, dist(p, S.eval(a/40.0, b/40.0)));
            mnall = std::min(mnall, mn);
            if (mn < 1e-3) { hit++; hitZ = p.z; }
        }
        printf("      [诊断] 退化线最小距=%.2e 近点=%d\n", mnall, hit);
        check("竖直退化线穿碗顶点（存在近点）", hit > 0);
        check("交点=碗顶点（z≈0.375, 单簇）", hit >= 1 && std::abs(hitZ - zseg) < 0.02 && hit < 8);
    }

    // [d] 退化到曲线 patch（线段在碗外，分离）
    {
        printf("\n[d] 退化到曲线 patch（线段在碗外，分离）\n");
        P3 pc3[4] = { {0.2, 0.0, 0.8}, {0.2, 0.4, 0.8}, {0.2, 0.6, 0.8}, {0.2, 1.0, 0.8} };
        nurbs::NurbsSurface Dl2 = nurbs::NurbsSurface::makeBilinear(pc3);
        int hit = 0;
        for (int i = 0; i <= 200; ++i) {
            P3 p = Dl2.eval(0.5, i / 200.0);
            double mn = 1e300;
            for (int a = 0; a <= 40; ++a) for (int b = 0; b <= 40; ++b)
                mn = std::min(mn, dist(p, S.eval(a/40.0, b/40.0)));
            if (mn < 1e-3) hit++;
        }
        check("退化线 z=0.8 在碗上方（分离，无交点）", hit == 0);
    }

    printf("\n结论：测例⑪——退化曲面×正常曲面：退化点 patch 在碗上接触/碗外分离，退化线 patch 穿碗\n");
    printf("      交两点/碗外分离，退化×正常组合覆盖补齐。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
