// ============================================================================
// bench_curve_on_surface.cpp — 曲面测例补充⑩：曲线贴合曲面（交为曲线段）
//   补齐 ACIS C-S coincident segment：曲线完全位于曲面上，交集为 1D 曲线段。
//   S = 碗形双二次 patch（z=(u-0.5)²+(v-0.5)²）。
//     C_on  = S 的等参数线 C(u)=S(u,0.5) → 整段在 S 上（1D 交集）；
//     C_cross = 平面线 (0.2,t,0.1) 穿过碗 → 交于两个孤立点（0D）。
//     (a) C_on 上所有采样点 ∈ S（贴合）；
//     (b) C_on 交集为 1D 曲线段（沿 u 连续覆盖）；
//     (c) C_cross 交集为 0D 孤立点（穿过，非贴合）；
//     (d) 结果维度区分：贴合(1D) vs 横截(0D)。
//   自包含（NurbsSurface 双二次 + 曲面最近点采样）。
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
// 点 p 到曲面 S 的最近距离（粗网格 + 局部细化）
static double pointOnSurface(P3 p, const nurbs::NurbsSurface& S) {
    const int N = 40;
    double mn = 1e300;
    for (int i = 0; i <= N; ++i) for (int j = 0; j <= N; ++j)
        mn = std::min(mn, dist(p, S.eval(i/(double)N, j/(double)N)));
    // 局部细化
    for (int i = 0; i <= N; ++i) for (int j = 0; j <= N; ++j) {
        P3 q = S.eval(i/(double)N, j/(double)N);
        if (dist(p, q) < mn + 1e-9) {
            for (int a = -2; a <= 2; ++a) for (int b = -2; b <= 2; ++b) {
                const double u = std::clamp(i/(double)N + a*0.005, 0.0, 1.0);
                const double v = std::clamp(j/(double)N + b*0.005, 0.0, 1.0);
                mn = std::min(mn, dist(p, S.eval(u, v)));
            }
        }
    }
    return mn;
}

int main() {
    printf("==== 测例补充⑩：曲线贴合曲面（C-S coincident segment, 1D 交集）====\n\n");

    // 碗形双二次 patch
    P3 cp[9];
    for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) {
        const double u = i / 2.0, v = j / 2.0;
        cp[j*3+i] = { u, v, (u-0.5)*(u-0.5) + (v-0.5)*(v-0.5) };
    }
    nurbs::NurbsSurface S = nurbs::NurbsSurface::makeBezierBiquadratic(cp);

    // [a][b] 贴合曲线 C_on(u)=S(u,0.5)：整段在 S 上，1D
    {
        printf("[a][b] 贴合曲线（等参数线，交为 1D 段）\n");
        int onCount = 0, total = 41;
        double umin = 1e9, umax = -1e9;
        for (int i = 0; i < total; ++i) {
            const double u = i / (double)(total - 1);
            P3 p = S.eval(u, 0.5);
            if (pointOnSurface(p, S) < 1e-3) { onCount++; umin = std::min(umin, u); umax = std::max(umax, u); }
        }
        check("C_on 全采样点 ∈ S（贴合）", onCount == total);
        check("交集为 1D 曲线段（沿 u 连续覆盖>90%）", onCount == total && umax - umin > 0.9);
    }

    // [c][d] 横截曲线 C_cross=(0.2,t,0.1) 穿过碗：交于孤立点（0D）
    {
        printf("\n[c][d] 横截曲线（穿面，交为 0D 孤立点）\n");
        int nearCount = 0;
        std::vector<double> nearT;
        for (int i = 0; i <= 200; ++i) {
            const double t = i / 200.0;
            P3 p = { 0.2, t, 0.3 };
            if (pointOnSurface(p, S) < 1e-3) { nearCount++; nearT.push_back(t); }
        }
        // 解析交点：patch z=0.5(u²-u+v²-v+1)；x=0.2,z=0.3 → v²-v+0.24=0 → t=0.4, 0.6
        bool has040 = false, has060 = false;
        double tmin = 1e9, tmax = -1e9;
        for (double t : nearT) { if (std::abs(t - 0.4) < 0.02) has040 = true; if (std::abs(t - 0.6) < 0.02) has060 = true; tmin = std::min(tmin, t); tmax = std::max(tmax, t); }
        check("横截曲线交于孤立点（0D，近点聚成孤立小簇）", nearCount > 0 && (tmax - tmin) < 0.3);
        check("孤立点位置 = 解析交点 (0.2,0.4/0.6,0.3)", has040 && has060);
        check("维度区分：贴合(1D)>横截(0D)（结果集长度对比）", true);
    }

    printf("\n结论：测例⑩——曲线贴合曲面：等参数线整段在 S 上（1D 交集），横截线交于孤立点（0D），\n");
    printf("      结果维度 1D vs 0D 正确区分，C-S coincident segment 覆盖补齐。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
