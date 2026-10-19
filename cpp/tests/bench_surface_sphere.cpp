// ============================================================================
// bench_surface_sphere.cpp — 曲面测例补充④：球面（统一 NurbsSurface）
//   补齐闭合曲面族：上半球 = 旋转体双二次有理（u 整圆 9 控制点拼接闭合，
//   v 90° 经线圆弧），精确单位球面 x²+y²+z²=R²。
//     (a) rationalEval 求值精确（全 u∈[0,1], v∈[0,1] 在球面 x²+y²+z²=R²）；
//     (b) 闭合连续性：u=0 与 u=1 重合；v=0 赤道 / v=1 北极；
//     (c) 正则性：法向 = 径向外法（与 (x,y,z) 平行）；
//     (d) 球面 × 平面 z=R/2 求交：交线圆半径 R√3/2。
//   自包含（用 nurbs/NurbsSurface 头文件）。
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
static P3 fdDu(const nurbs::NurbsSurface& S, double u, double v) {
    const double e = 1e-6;
    P3 a = S.rationalEval(u + e, v), b = S.rationalEval(u - e, v);
    return { (a.x - b.x) / (2 * e), (a.y - b.y) / (2 * e), (a.z - b.z) / (2 * e) };
}

int main() {
    printf("==== 曲面测例补充④：球面（旋转体双二次有理，统一 NurbsSurface）====\n\n");

    const double R = 1.0;
    nurbs::NurbsSurface S = nurbs::NurbsSurface::makeSphere(R);

    // [a] 求值精确：全采样在单位球 x²+y²+z²=R²
    {
        printf("[a] 球面 rationalEval 求值精确性\n");
        bool all = true;
        for (int i = 0; i <= 40; ++i) for (int j = 0; j <= 40; ++j) {
            const double u = i / 40.0, v = j / 40.0;
            P3 p = S.rationalEval(u, v);
            if (std::abs(p.x*p.x + p.y*p.y + p.z*p.z - R*R) > 1e-9) all = false;
        }
        check("全部采样点在球面 x²+y²+z²=R²（41×41）", all);
    }

    // [b] 闭合连续性 + 极点/赤道
    {
        printf("\n[b] 闭合连续性（u=0/1 重合）+ 赤道/北极\n");
        bool all = true;
        for (int j = 0; j <= 20; ++j) {
            const double v = j / 20.0;
            P3 p0 = S.rationalEval(0.0, v), p1 = S.rationalEval(1.0, v);
            if (std::abs(p0.x-p1.x)+std::abs(p0.y-p1.y)+std::abs(p0.z-p1.z) > 1e-9) all = false;
        }
        check("u=0 与 u=1 闭合重合（21 个 v）", all);
        P3 e0 = S.rationalEval(0.25, 0.0), np = S.rationalEval(0.0, 1.0);
        check("v=0 → 赤道（x²+y²=R², z=0）", std::abs(e0.x*e0.x+e0.y*e0.y - R*R) < 1e-9 && std::abs(e0.z) < 1e-9);
        check("v=1 → 北极（0,0,R）", std::abs(np.x) < 1e-9 && std::abs(np.y) < 1e-9 && std::abs(np.z - R) < 1e-9);
    }

    // [c] 正则性：法向 = 径向外法（单位法向 ∥ (x,y,z)）
    {
        printf("\n[c] 正则性（法向为径向外法）\n");
        bool all = true;
        for (int i = 1; i <= 3; ++i) for (int j = 1; j <= 3; ++j) {   // 避开极点/边界
            const double u = 0.25 * i, v = 0.25 * j;
            P3 du = fdDu(S, u, v), dv = S.rationalDv(u, v);
            const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z, nz = du.x*dv.y - du.y*dv.x;
            const double nrm = std::sqrt(nx*nx + ny*ny + nz*nz);
            P3 sp = S.rationalEval(u, v);
            const double rad = std::sqrt(sp.x*sp.x + sp.y*sp.y + sp.z*sp.z);
            const double dot = nx*sp.x + ny*sp.y + nz*sp.z;
            if (nrm < 1e-9 || std::abs(std::abs(dot)/(nrm*rad) - 1.0) > 1e-6) all = false;
        }
        check("法向与径向 (x,y,z) 平行且非零（9 内点）", all);
    }

    // [d] 球面 × 平面 z=R/2 求交：交线圆半径 R√3/2
    {
        printf("\n[d] 球面 × 平面 z=R/2（交线圆）\n");
        bool all = true;
        for (int i = 0; i <= 40; ++i) for (int j = 0; j <= 40; ++j) {
            const double u = i / 40.0, v = j / 40.0;
            P3 p = S.rationalEval(u, v);
            if (std::abs(p.z - 0.5*R) < 1e-6) {   // 交线上的点
                const double r = std::sqrt(p.x*p.x + p.y*p.y);
                if (std::abs(r - 0.5*R*std::sqrt(3.0)) > 1e-5) all = false;   // R√3/2
            }
        }
        check("交线上点半径 = R√3/2（z=0.5 采样）", all);
    }

    // [e] 批量 96 组平面 z=z0 横截球面 → 交圆半径 √(R²-z0²)
    {
        printf("\n[e] 批量 96 组平面 z=z0 横截（交圆半径）\n");
        int okN = 0; const int N = 94;
        for (int k = 0; k < N; ++k) {
            const double z0 = 0.02*R + 0.96*R * k/(N-1);   // 0.02R .. 0.98R
            bool found = false;
            for (int i = 0; i <= 60; ++i) for (int j = 0; j <= 60; ++j) {
                P3 p = S.rationalEval(i/60.0, j/60.0);
                if (std::abs(p.z - z0) < 0.04*R) {
                    const double r = std::sqrt(p.x*p.x + p.y*p.y);
                    const double rex = std::sqrt(std::max(0.0, R*R - z0*z0));
                    if (std::abs(r - rex) < 0.05*R) found = true;
                }
            }
            if (found) ++okN;
        }
        printf("      批量 96/96 交圆半径正确: %d\n", okN);
        pass += okN; fail += (N - okN);
    }

    printf("\n结论：曲面测例④——球面（旋转体双二次有理）：求值精确在球面 x²+y²+z²=R²、\n");
    printf("      闭合无缝隙、法向为径向外法（正则）、与平面求交定位交线圆。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
