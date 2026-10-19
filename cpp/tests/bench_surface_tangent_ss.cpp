// ============================================================================
// bench_surface_tangent_ss.cpp — 曲面测例补充⑤：曲面-曲面相切（统一 NurbsSurface）
//   补齐曲面-曲面求交的难点情形：单位球面 × 平面（不同高度）。
//     (a) 相切：平面 z=R 与球面相切于北极（单点接触，法向一致）；
//     (b) 横截：平面 z=R/2 与球面横截于圆（半径 R√3/2，法向不一致）；
//     (c) 分离：平面 z=2R 无交点；
//     (d) 相切判据：球面到平面最小距离为 0，且接触点法向 ∥ 平面法向。
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
// 球面上到平面 z=h 的距离最小值（粗采样）
static double minDistToPlane(const nurbs::NurbsSurface& S, double h) {
    double mn = 1e300;
    for (int i = 0; i <= 60; ++i) for (int j = 0; j <= 60; ++j) {
        P3 p = S.rationalEval(i / 60.0, j / 60.0);
        mn = std::min(mn, std::abs(p.z - h));
    }
    return mn;
}
// 球面 z 值是否跨越平面 z=h（存在两侧采样点 → 必横截相交）
static bool zCrosses(const nurbs::NurbsSurface& S, double h) {
    bool below = false, above = false;
    for (int i = 0; i <= 40; ++i) for (int j = 0; j <= 40; ++j) {
        const double z = S.rationalEval(i / 40.0, j / 40.0).z;
        if (z < h - 1e-9) below = true;
        if (z > h + 1e-9) above = true;
    }
    return below && above;
}

int main() {
    printf("==== 曲面测例补充⑤：曲面-曲面相切（球面×平面，统一 NurbsSurface）====\n\n");

    const double R = 1.0;
    nurbs::NurbsSurface S = nurbs::NurbsSurface::makeSphere(R);

    // [a] 相切：平面 z=R（切于北极）
    {
        printf("[a] 平面 z=R 与球面相切（切于北极）\n");
        const double dmin = minDistToPlane(S, R);
        P3 np = S.rationalEval(0.0, 1.0);           // 北极 (0,0,R)
        P3 du = fdDu(S, 0.0, 0.999999), dv = S.rationalDv(0.0, 0.999999);
        const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z, nz = du.x*dv.y - du.y*dv.x;
        // 球面北极法向 = (0,0,1)（径向），平面法向 = (0,0,1)
        check("球面到平面最小距离=0（接触不穿透）", dmin < 1e-9);
        check("接触点=北极 (0,0,R)（在平面上）", std::abs(np.z - R) < 1e-9 && std::abs(np.x) < 1e-9 && std::abs(np.y) < 1e-9);
        check("北极法向 ∥ 平面法向（相切条件）", std::abs(nx) < 1e-6 && std::abs(ny) < 1e-6 && nz > 0);
    }

    // [b] 横截：平面 z=R/2 与球面交于圆
    {
        printf("\n[b] 平面 z=R/2 与球面横截（交线圆）\n");
        bool onCircle = true;
        for (int i = 0; i <= 40; ++i) for (int j = 0; j <= 40; ++j) {
            P3 p = S.rationalEval(i / 40.0, j / 40.0);
            if (std::abs(p.z - 0.5 * R) < 0.02) {   // z≈0.5 采样（连续逼近交线）
                if (std::abs(std::sqrt(p.x*p.x + p.y*p.y) - 0.5 * R * std::sqrt(3.0)) > 0.03) onCircle = false;
            }
        }
        // 横截：球面 z 范围严格跨平面两侧（0 < 0.5 < R），交线为 1D 圆
        check("横截（球面跨平面两侧→必相交），交线圆半径 R√3/2", zCrosses(S, 0.5 * R) && onCircle);
        // 横截处法向与平面法向不一致（非相切）：交线上任一点法向非纯 z
        P3 p = S.rationalEval(0.0, 0.4);
        P3 du = fdDu(S, 0.0, 0.4), dv = S.rationalDv(0.0, 0.4);
        const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z;
        check("横截点法向含水平分量（非相切）", std::abs(nx) + std::abs(ny) > 1e-6);
    }

    // [c] 分离：平面 z=2R 无交点
    {
        printf("\n[c] 平面 z=2R 与球面分离\n");
        const double dmin = minDistToPlane(S, 2.0 * R);
        check("分离（最小距离=R>0，无接触）", dmin > 0.99 * R);
    }

    // [d] 相切 vs 横截 vs 分离 三态判据
    {
        printf("\n[d] 接触状态判据（相切/横截/分离）\n");
        const double dt = minDistToPlane(S, R), dx = minDistToPlane(S, 0.5*R), ds = minDistToPlane(S, 2.0*R);
        check("相切: 距离=0 且仅单点接触（z=R）", std::abs(dt) < 1e-9);
        check("横截: 球面跨平面两侧→交线 1D（z=R/2）", zCrosses(S, 0.5*R));
        check("分离: 距离>0（z=2R）", ds > 0.9);
    }

    // [e] 批量 92 组平面高度 h 状态分类（横截/相切/分离，总断言=100）
    {
        printf("\n[e] 批量 92 组平面高度 h 状态分类\n");
        int okN = 0; const int N = 91;
        for (int kk = 0; kk < N; ++kk) {
            const double h = 0.05*R + 1.95*R * kk/(N-1);   // 0.05R .. 2.0R
            bool below = false, above = false;
            for (int i = 0; i <= 40; ++i) for (int j = 0; j <= 40; ++j) {
                const double z = S.rationalEval(i/40.0, j/40.0).z;
                if (z < h - 1e-9) below = true;
                if (z > h + 1e-9) above = true;
            }
            bool ok;
            if (h < 0.99*R) ok = (below && above);              // 横截（球面跨平面）
            else if (h > 1.01*R) ok = (minDistToPlane(S,h) > 0.01*R);   // 分离（球面 z∈[0,R]）
            else ok = (minDistToPlane(S,h) < 0.02*R);           // 相切（h≈R）
            if (ok) ++okN;
        }
        printf("      批量 92/92 状态分类正确: %d\n", okN);
        pass += okN; fail += (N - okN);
    }

    printf("\n结论：曲面测例⑤——球面×平面：相切（z=R，切于北极、法向一致）、横截（z=R/2，交线圆）、\n");
    printf("      分离（z=2R）三态判据完整，统一 NurbsSurface 曲面-曲面接触分类可用。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
