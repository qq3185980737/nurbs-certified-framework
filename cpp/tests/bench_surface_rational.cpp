// ============================================================================
// bench_surface_rational.cpp — 曲面测例补充①：有理权重 + 非均匀节点（统一 NurbsSurface）
//   补齐曲线侧 bench_nurbs_wide 对应的曲面能力：
//     (a) 有理四分之一圆柱侧面（u 方向二次圆弧 w=[1,√2/2,1]，v 线性）——
//         rationalEval 精确求值、rationalDu/rationalDv 与有限差分对照；
//     (b) 权重一致性：w≡1 时 rationalEval 与 eval 完全一致；
//     (c) 非均匀 v 节点（vKnots=[0,0,0.35,1,1]）下 rationalEval 精确性（z 仍线性）；
//     (d) 有理曲面 × 竖直线求交：Newton 定位 + 残差收敛（保证性框架的数值入口）。
//   有理求值走齐次 4D de Boor；权重字段 w 由 nurbs::NurbsSurface 承载。
// 自包含（用 nurbs/NurbsSurface 头文件）。
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

// 解析圆柱侧面点：S(u,v)=(R cos(pi u/2), R sin(pi u/2), z0 + H v)
static P3 exactQuarterCyl(double R, double z0, double H, double u, double v) {
    const double a = 3.14159265358979323846 * u * 0.5;
    return { R * std::cos(a), R * std::sin(a), z0 + H * v };
}

// 有限差分偏导
static P3 fdDu(const nurbs::NurbsSurface& S, double u, double v) {
    const double e = 1e-6;
    P3 a = S.rationalEval(u + e, v), b = S.rationalEval(u - e, v);
    return { (a.x - b.x) / (2 * e), (a.y - b.y) / (2 * e), (a.z - b.z) / (2 * e) };
}
static P3 fdDv(const nurbs::NurbsSurface& S, double u, double v) {
    const double e = 1e-6;
    P3 a = S.rationalEval(u, v + e), b = S.rationalEval(u, v - e);
    return { (a.x - b.x) / (2 * e), (a.y - b.y) / (2 * e), (a.z - b.z) / (2 * e) };
}

int main() {
    printf("==== 曲面测例补充①：有理权重 + 非均匀节点（统一 NurbsSurface）====\n\n");

    // [a] 有理四分之一圆柱侧面：rationalEval 精确性（圆弧 + 线性高度）
    {
        printf("[a] 有理四分之一圆柱侧面（u 二次圆弧 w=[1,√2/2,1]，v 线性）求值精确性\n");
        nurbs::NurbsSurface S = nurbs::NurbsSurface::makeRationalQuarterCylinder(1.0, 0.0, 2.0);
        bool all = true;
        for (int i = 0; i <= 20; ++i) for (int j = 0; j <= 20; ++j) {
            const double u = i / 20.0, v = j / 20.0;
            P3 p = S.rationalEval(u, v), e = exactQuarterCyl(1.0, 0.0, 2.0, u, v);
            const double r2 = p.x * p.x + p.y * p.y;
            if (std::abs(r2 - 1.0) > 1e-9 || std::abs(p.z - e.z) > 1e-9) all = false;
        }
        check("在圆弧上 x²+y²=1 且 z 线性（21×21 采样）", all);
        // 二次有理圆弧（w=√2/2）参数化：u=0.5 是中点 45°，端点 u=0/1 为 (R,0)/(0,R)
        P3 pm = S.rationalEval(0.5, 0.5), p0 = S.rationalEval(0.0, 0.5), p1 = S.rationalEval(1.0, 0.5);
        check("u=0.5 → 45° 中点（√2/2,√2/2），u=0/1 → 端点",
              std::abs(pm.x - std::sqrt(0.5)) < 1e-9 && std::abs(pm.y - std::sqrt(0.5)) < 1e-9
              && std::abs(p0.x - 1.0) < 1e-9 && std::abs(p0.y) < 1e-9
              && std::abs(p1.x) < 1e-9 && std::abs(p1.y - 1.0) < 1e-9 && std::abs(p0.z - 1.0) < 1e-9);
    }

    // [b] rationalDu / rationalDv 与有限差分对照
    {
        printf("\n[b] 有理偏导 rationalDu/rationalDv 与有限差分对照（9 采样点）\n");
        nurbs::NurbsSurface S = nurbs::NurbsSurface::makeRationalQuarterCylinder(1.0, 0.0, 2.0);
        bool all = true;
        for (int i = 1; i <= 3; ++i) for (int j = 1; j <= 3; ++j) {
            const double u = 0.25 * i, v = 0.25 * j;   // 内点，避开 clamped 端点中心差分
            P3 du = S.rationalDu(u, v), fdu = fdDu(S, u, v);
            P3 dv = S.rationalDv(u, v), fdv = fdDv(S, u, v);
            const double e1 = std::abs(du.x - fdu.x) + std::abs(du.y - fdu.y) + std::abs(du.z - fdu.z);
            const double e2 = std::abs(dv.x - fdv.x) + std::abs(dv.y - fdv.y) + std::abs(dv.z - fdv.z);
            if (e1 > 1e-5 || e2 > 1e-5) all = false;
        }
        check("rationalDu/Dv 与有限差分一致（9 内点）", all);
        // 圆柱侧面法向应为径向：normal = du × dv ∝ (cos a, sin a, 0)
        P3 du = S.rationalDu(0.5, 0.5), dv = S.rationalDv(0.5, 0.5);
        const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z;
        const double ang = std::atan2(ny, nx);
        check("法向 du×dv 指向径向（45°）", std::abs(std::abs(ang) - 0.7853981633974483) < 1e-6);
    }

    // [c] 权重一致性：w≡1 时 rationalEval 与 eval 完全一致
    {
        printf("\n[c] 权重一致性（w≡1 → rationalEval == eval）\n");
        // 用 makeBilinear 平面（w=1）
        SurfacePoint cp[4] = { {0,0,0}, {1,0,0}, {0,1,0.3}, {1,1,0.7} };
        nurbs::NurbsSurface S = nurbs::NurbsSurface::makeBilinear(cp);
        bool all = true;
        for (int i = 0; i <= 10; ++i) for (int j = 0; j <= 10; ++j) {
            const double u = i / 10.0, v = j / 10.0;
            P3 a = S.rationalEval(u, v), b = S.eval(u, v);
            if (std::abs(a.x - b.x) + std::abs(a.y - b.y) + std::abs(a.z - b.z) > 1e-12) all = false;
        }
        check("w≡1 平面 rationalEval == eval（11×11）", all);
    }

    // [d] 非均匀 v 节点：rationalEval 精确性（z 仍线性 z0+H·v）
    {
        printf("\n[d] 非均匀 v 节点（vKnots=[0,0,0.35,1,1]）求值精确性\n");
        nurbs::NurbsSurface S; S.nu = 3; S.nv = 3; S.pu = 2; S.pv = 1;
        S.uKnots = {0,0,0,1,1,1}; S.vKnots = {0,0,0.35,1,1};
        const double R = 1.0, z0 = 0.0, H = 2.0, r2 = std::sqrt(0.5);
        const double vz[3] = { 0.0, 0.35 * H, H };
        const double ax[3] = { R, R, 0.0 }, ay[3] = { 0.0, R, R }, wu[3] = {1.0, r2, 1.0};
        S.C.resize(9); S.w.assign(9, 1.0);
        for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) {
            S.C[j*3+i] = { ax[i], ay[i], z0 + vz[j] };
            S.w[j*3+i] = wu[i];
        }
        bool all = true;
        for (int i = 0; i <= 20; ++i) for (int j = 0; j <= 20; ++j) {
            const double u = i / 20.0, v = j / 20.0;
            P3 p = S.rationalEval(u, v);
            const double r2 = p.x * p.x + p.y * p.y;
            const double ez = z0 + H * v;             // z 应仍线性（非均匀节点下）
            if (std::abs(r2 - R * R) > 1e-9 || std::abs(p.z - ez) > 1e-9) all = false;
        }
        check("非均匀 v 节点求值精确（21×21，在圆上 + z 线性）", all);
        P3 q = S.rationalEval(0.5, 0.5);
        check("u=0.5 → 45° 点（x=y=√2/2）",
              std::abs(q.x - std::sqrt(0.5)) < 1e-9 && std::abs(q.y - std::sqrt(0.5)) < 1e-9);
    }

    // [e] 有理曲面 × 竖直线求交：Newton 定位 + 残差收敛
    {
        printf("\n[e] 有理四分之一圆柱侧面 × 竖直线 x=y=√2/2（交点 u=0.5, z 可变）\n");
        nurbs::NurbsSurface S = nurbs::NurbsSurface::makeRationalQuarterCylinder(1.0, 0.0, 2.0);
        // 竖直线 L(w)=(x0,y0,w)，求 F(u,v)=S(u,v)-(x0,y0,0) 的 u 使 x=y=x0
        const double x0 = std::sqrt(0.5);      // = cos45°
        // 残差 g(u) = S.x(u,v)-x0，在 u=0.5 应为 0
        P3 p = S.rationalEval(0.5, 0.3);
        check("竖直线 x=y=√2/2 与圆柱交于 u=0.5（残差）",
              std::abs(p.x - x0) < 1e-12 && std::abs(p.y - x0) < 1e-12);
        // Newton 在 u 上求 g(u)=0（固定 z 由直线给出，v 解出）
        // 简化：沿 u 找使 x=x0 的解（u=0.5）
        double u = 0.4;
        for (int it = 0; it < 8; ++it) {
            P3 sp = S.rationalEval(u, 0.5), sdu = S.rationalDu(u, 0.5);
            u -= (sp.x - x0) / sdu.x;
        }
        check("Newton 定位 u→0.5（残差<1e-12）", std::abs(u - 0.5) < 1e-12);
    }

    printf("\n结论：曲面测例①——有理权重（w=[1,√2/2,1]）+ 非均匀节点在统一 NurbsSurface 下\n");
    printf("      rationalEval/du/dv 求值与偏导精确、w≡1 与无理一致、非均匀节点 z 线性、求交定位收敛。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
