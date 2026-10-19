// ============================================================================
// bench_surface_closed.cpp — 曲面测例补充②：闭合圆柱面（统一 NurbsSurface）
//   补齐曲线侧 bench_closed 对应的曲面能力：闭合圆柱侧面（u 方向 4 段 90°
//   二次有理圆弧拼接的精确整圆，几何闭合 S(0,v)==S(1,v)），v 线性（高度）。
//     (a) rationalEval 精确性：全 u∈[0,1] 采样点在圆柱面上（距轴=R）且 z 线性；
//     (b) 闭合连续性：u=0 与 u=1 几何重合（闭合无缝隙）；
//     (c) 正则性：法向 rationalDu×rationalDv 模在 9 采样点恒定非零（正则非折叠）；
//     (d) 圆柱面 × 平面 y=0.5 求交：Newton 定位 u + 残差收敛。
//   自包含（用 nurbs/NurbsSurface 头文件）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
using SurfacePoint = nurbs::SurfacePoint;

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
    printf("==== 曲面测例补充②：闭合圆柱面（拼接精确整圆，统一 NurbsSurface）====\n\n");

    // 单位圆柱侧面 R=1, z∈[0,2]，9 控制点二次有理整圆（4 段 90°）
    nurbs::NurbsSurface C = nurbs::NurbsSurface::makeCylinder(1.0, 0.0, 2.0);

    // [a] rationalEval 精确性：距轴=R + z 线性
    {
        printf("[a] 闭合圆柱面 rationalEval 精确性（全 u∈[0,1]）\n");
        bool all = true;
        for (int i = 0; i <= 40; ++i) for (int j = 0; j <= 40; ++j) {
            const double u = i / 40.0, v = j / 40.0;
            P3 p = C.rationalEval(u, v);
            const double r2 = p.x * p.x + p.y * p.y;
            if (std::abs(r2 - 1.0) > 1e-9 || std::abs(p.z - 2.0 * v) > 1e-9) all = false;
        }
        check("全部采样点在圆柱面 x²+y²=1 且 z 线性（41×41）", all);
    }

    // [b] 闭合连续性：u=0 与 u=1 几何重合
    {
        printf("\n[b] 闭合连续性（u=0 与 u=1 重合）\n");
        bool all = true;
        for (int j = 0; j <= 20; ++j) {
            const double v = j / 20.0;
            P3 p0 = C.rationalEval(0.0, v), p1 = C.rationalEval(1.0, v);
            const double d = std::sqrt((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y)+(p0.z-p1.z)*(p0.z-p1.z));
            if (d > 1e-12) all = false;
        }
        check("u=0 与 u=1 闭合重合（21 个 v 采样）", all);
        P3 p0 = C.rationalEval(0.0, 0.5);
        check("u=0 → (R,0)（闭合起点）", std::abs(p0.x - 1.0) < 1e-9 && std::abs(p0.y) < 1e-9);
    }

    // [c] 正则性：法向模恒定非零（正则非折叠）
    {
        printf("\n[c] 正则性（法向 |rationalDu×rationalDv| 恒定非零）\n");
        bool all = true; double base = -1;
        for (int i = 0; i <= 2; ++i) for (int j = 0; j <= 2; ++j) {
            const double u = 0.25 * (i + 1), v = 0.25 * (j + 1);   // 内点，避开 clamped 端点
            P3 du = fdDu(C, u, v), dv = C.rationalDv(u, v);
            const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z, nz = du.x*dv.y - du.y*dv.x;
            const double nrm = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (base < 0) base = nrm;
            else if (std::abs(nrm - base) > 1e-3 * base) all = false;   // 恒定（±0.1%）
            if (nrm < 1e-9) all = false;                                   // 非零
        }
        check("法向模恒定非零（9 采样点）", all);
        P3 du = fdDu(C, 0.5, 0.5), dv = C.rationalDv(0.5, 0.5);
        const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z, nz = du.x*dv.y - du.y*dv.x;
        const double nrm = std::sqrt(nx*nx + ny*ny + nz*nz);
        P3 sp = C.rationalEval(0.5, 0.5);
        const double dot = nx*sp.x + ny*sp.y + nz*sp.z;      // 法向·径向
        const double rad = std::sqrt(sp.x*sp.x + sp.y*sp.y);
        // 圆柱侧面法向 = 径向外法（无 z 分量，且与 (x,y,0) 平行）
        check("法向方向为径向外法（正则性）",
              std::abs(nz) < 1e-6 && std::abs(std::abs(dot) / (nrm * rad) - 1.0) < 1e-6);
    }

    // [d] 圆柱面 × 平面 y=0.5 求交：Newton 定位 u + 残差收敛
    {
        printf("\n[d] 圆柱面 × 平面 y=0.5（交线为两条母线）\n");
        double u = 0.1;   // θ=36°，靠近 θ=30° 解
        for (int it = 0; it < 10; ++it) {
            P3 sp = C.rationalEval(u, 0.5), sdu = C.rationalDu(u, 0.5);
            if (std::abs(sdu.y) < 1e-12) break;
            u -= (sp.y - 0.5) / sdu.y;
            if (u < 0) u = 0; if (u > 1) u = 1;    // 投影回参数域
        }
        P3 hit = C.rationalEval(u, 0.5);
        const double x30 = std::sqrt(3.0) / 2.0;
        check("Newton 定位 y=0.5 交点（残差<1e-9）",
              std::abs(hit.y - 0.5) < 1e-9 && std::abs(std::abs(hit.x) - x30) < 1e-9 && std::abs(hit.z - 1.0) < 1e-9);
        check("交点仍在圆柱面 x²+y²=1", std::abs(hit.x*hit.x + hit.y*hit.y - 1.0) < 1e-9);
    }

    printf("\n结论：曲面测例②——闭合圆柱面（4 段 90° 二次有理圆弧拼接精确整圆）在统一 NurbsSurface\n");
    printf("      rationalEval 求值精确、闭合无缝隙、法向正则非折叠、与平面求交定位收敛。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
