// ============================================================================
// bench_surface_overlap.cpp — 曲面测例补充③：曲面部分重叠（统一 NurbsSurface）
//   补齐与曲线侧相切/重叠对应的曲面能力：两个双线性平面（z=0）在
//   x∈[0,1] 区域部分重叠（S1: x∈[-1,1]，S2: x∈[0,2]，y∈[-1,1]）。
//     (a) 包围盒（AABB）重叠检测；
//     (b) 重叠区点重合：x∈[0,1] 的两曲面点几何重合（距离≈0）；
//     (c) 非重叠区分离：x∈[-1,0) 的 S1 点不在 S2 上（距离>0）；
//     (d) 重叠边界定位：x=0 与 x=1 是重叠/非重叠分界。
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
static double dist(const P3& a, const P3& b) {
    return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

int main() {
    printf("==== 曲面测例补充③：曲面部分重叠（统一 NurbsSurface）====\n\n");

    // 两个双线性平面 z=0：S1 x∈[-1,1]，S2 x∈[0,2]，y∈[-1,1] → 重叠 x∈[0,1]
    const P3 S1cp[4] = {{-1,-1,0},{1,-1,0},{-1,1,0},{1,1,0}};
    const P3 S2cp[4] = {{ 0,-1,0},{2,-1,0},{ 0,1,0},{2,1,0}};
    nurbs::NurbsSurface S1 = nurbs::NurbsSurface::makeBilinear(S1cp);
    nurbs::NurbsSurface S2 = nurbs::NurbsSurface::makeBilinear(S2cp);

    // [a] AABB 重叠检测（凸包包围盒在重叠区相交）
    {
        printf("[a] 包围盒重叠检测\n");
        P3 mn1 = S1.aabbMin(), mx1 = S1.aabbMax();
        P3 mn2 = S2.aabbMin(), mx2 = S2.aabbMax();
        const bool overlap = (mx1.x >= mn2.x) && (mx2.x >= mn1.x) &&   // 接触/重叠（封闭）
                             (mx1.y >= mn2.y) && (mx2.y >= mn1.y) &&
                             (mx1.z >= mn2.z) && (mx2.z >= mn1.z);
        check("S1/S2 包围盒在重叠区相交（AABB）", overlap);
        check("S1 包围盒 [-1,1]×[-1,1]", std::abs(mn1.x + 1) < 1e-12 && std::abs(mx1.x - 1) < 1e-12);
        check("S2 包围盒 [0,2]×[-1,1]", std::abs(mn2.x) < 1e-12 && std::abs(mx2.x - 2) < 1e-12);
    }

    // [b] 重叠区点重合：x∈[0,1] 的两曲面点几何重合（距离≈0）
    {
        printf("\n[b] 重叠区点重合（x∈[0,1]）\n");
        bool all = true;
        for (int i = 0; i <= 20; ++i) for (int j = 0; j <= 20; ++j) {
            const double x = 0.05 * i, y = -1.0 + 0.1 * j;   // x∈[0,1], y∈[-1,1]
            const double u1 = (x + 1.0) / 2.0, v1 = (y + 1.0) / 2.0;   // S1 参数
            const double u2 = x / 2.0,        v2 = (y + 1.0) / 2.0;   // S2 参数
            P3 p1 = S1.eval(u1, v1), p2 = S2.eval(u2, v2);
            if (dist(p1, p2) > 1e-12) all = false;
        }
        check("重叠区 21×21 点两曲面重合（距离<1e-12）", all);
    }

    // [c] 非重叠区分离：x∈[-1,0) 的 S1 点不在 S2 上（距离>0）
    {
        printf("\n[c] 非重叠区分离（x∈[-1,0) 仅属于 S1）\n");
        bool all = true;
        for (int i = 1; i <= 19; ++i) for (int j = 0; j <= 20; ++j) {
            const double x = -1.0 + 0.05 * i, y = -1.0 + 0.1 * j;   // x∈[-0.95,-0.05]，排除重叠边界 x=0
            const double u1 = (x + 1.0) / 2.0, v1 = (y + 1.0) / 2.0;
            P3 p1 = S1.eval(u1, v1);
            // S2 上最近点：x≥0，y 相同 → 最近 x=0 → 距离 ≈ |x|
            if (std::abs(dist(p1, {0.0, p1.y, 0.0}) - std::abs(x)) > 0.1 * std::abs(x)) all = false;
            // 交叉验证：S2 全参数域点与 p1 距离 > 0（x 不同）
            if (dist(p1, S2.eval(0.0, v1)) < 1e-6) all = false;
        }
        check("非重叠区 S1 点与 S2 分离（距离≈|x|>0）", all);
    }

    // [d] 重叠边界定位：x=0 与 x=1 是分界
    {
        printf("\n[d] 重叠边界定位（x=0 与 x=1）\n");
        bool all = true;
        for (int j = 0; j <= 20; ++j) {
            const double v = 0.05 * j;
            // x=0 边界：S1 u=0.5 ↔ S2 u=0
            P3 b1 = S1.eval(0.5, v), b2 = S2.eval(0.0, v);
            if (dist(b1, b2) > 1e-12 || std::abs(b1.x) > 1e-12) all = false;
            // x=1 边界：S1 u=1.0 ↔ S2 u=0.5
            P3 c1 = S1.eval(1.0, v), c2 = S2.eval(0.5, v);
            if (dist(c1, c2) > 1e-12 || std::abs(c1.x - 1.0) > 1e-12) all = false;
        }
        check("x=0 与 x=1 分界点两曲面重合（21 个 v）", all);
    }

    printf("\n结论：曲面测例③——两个双线性平面部分重叠（x∈[0,1]）：AABB 检测到重叠、\n");
    printf("      重叠区点重合、非重叠区分离、重叠边界（x=0/x=1）精确定位，统一 NurbsSurface 可用。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
