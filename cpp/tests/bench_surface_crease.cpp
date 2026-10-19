// ============================================================================
// bench_surface_crease.cpp — 曲面测例补充⑫：曲面折脊（C⁰ 折缝曲面）
//   补齐 ACIS C⁰ surface（curve 侧已有 kink，曲面侧补齐折脊）。
//   构造：nu=5, pu=2, nv=2, pv=1 的帐篷形 patch，u 向内部节点 0.5 重复 2 次
//   （uKnots={0,0,0,0.5,0.5,1,1,1}）→ u=0.5 处 C⁰ 折脊（值连续、法向跳变）。
//   (a) 折脊处求值连续（u→0.5⁻/0.5⁺ 坐标一致）；
//   (b) u 向一阶导数在折脊两侧不连续（折脊特征）；
//   (c) 折脊两侧单位法向不连续（夹角显著）；
//   (d) 折脊精确值 z≈0.25（峰顶）；
//   (e) 水平平面 z=0.1 横截 → 交线 1D 连续跨过折脊；
//   (f) 平面 z=0.3 分离（高于峰顶 0.25）；
//   (g) 折脊曲面无自交（全采样点唯一）；
//   (h) 折脊外区域正则（法向长度>0）；
//   (i) aabb 与帐篷几何一致（z∈[0,0.25]）；
//   (j) 折脊两侧子曲面块分别正则无奇点。
//   自包含（仅依赖 nurbs_surface.hpp）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
static int pass = 0, fail = 0;
static void check(const char* label, bool ok) {
    printf("  %-50s %s\n", label, ok ? "PASS" : "FAIL");
    if (ok) ++pass; else ++fail;
}
static double dist(P3 a, P3 b) {
    return std::sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}
// 由 du×dv 求单位法向
static P3 normal(const nurbs::NurbsSurface& S, double u, double v) {
    nurbs::SurfacePoint du = S.du(u, v), dv = S.dv(u, v);
    double cx = du.y*dv.z - du.z*dv.y;
    double cy = du.z*dv.x - du.x*dv.z;
    double cz = du.x*dv.y - du.y*dv.x;
    double L = std::sqrt(cx*cx + cy*cy + cz*cz);
    if (L < 1e-12) return {0,0,0};
    return { cx/L, cy/L, cz/L };
}
// 点 p 到曲面 S 的最近距离（粗网格 + 局部细化）
static double pointOnSurface(P3 p, const nurbs::NurbsSurface& S) {
    const int N = 40;
    double mn = 1e300;
    for (int i = 0; i <= N; ++i) for (int j = 0; j <= N; ++j)
        mn = std::min(mn, dist(p, S.eval(i/(double)N, j/(double)N)));
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
    printf("==== 测例补充⑫：曲面折脊（C⁰ 折缝曲面）====\n\n");

    // 帐篷形折脊 patch：u 向折脊在 0.5，z 从两端 0 升到峰顶 0.25
    nurbs::NurbsSurface K;
    K.nu = 5; K.pu = 2; K.nv = 2; K.pv = 1;
    K.uKnots = {0,0,0, 0.5,0.5, 1,1,1};   // 内部节点 0.5 重复 2 次 → C⁰ 折脊
    K.vKnots = {0,0,1,1};
    const double ui[5] = {0, 0.25, 0.5, 0.75, 1};
    const double zi[5] = {0, 0.1875, 0.25, 0.1875, 0};   // 帐篷峰
    for (int j = 0; j < 2; ++j) for (int i = 0; i < 5; ++i)
        K.C.push_back({ ui[i], (double)j, zi[i] });
    K.w.assign(10, 1.0);

    // [a] 折脊处求值连续
    {
        printf("[a] 折脊处求值连续（u→0.5⁻/0.5⁺）\n");
        P3 a = K.eval(0.499, 0.5), b = K.eval(0.501, 0.5);
        printf("      [诊断] K(0.499)=(%.4f,%.4f,%.4f) K(0.501)=(%.4f,%.4f,%.4f) dz=%.2e\n",
               a.x,a.y,a.z, b.x,b.y,b.z, std::abs(a.z-b.z));
        // 折脊处 C⁰：z 连续（x 是参数 u 的线性映射，两侧 x 自然不同，不作不连续判据）
        check("折脊处值连续（C⁰，z 一致）", std::abs(a.z - b.z) < 1e-3);
    }
    // [b] u 向一阶导数折脊两侧不连续
    {
        printf("[b] u 向导数在折脊两侧不连续\n");
        nurbs::SurfacePoint duL = K.du(0.49, 0.5), duR = K.du(0.51, 0.5);
        printf("      [诊断] du(0.49)=(%.4f,%.4f,%.4f) du(0.51)=(%.4f,%.4f,%.4f)\n",
               duL.x,duL.y,duL.z, duR.x,duR.y,duR.z);
        double ddu = dist(P3{duL.x,duL.y,duL.z}, P3{duR.x,duR.y,duR.z});
        check("折脊两侧 u 向导数跳变（>0.05）", ddu > 0.05);
    }
    // [c] 折脊两侧法向不连续
    {
        printf("[c] 折脊两侧单位法向不连续\n");
        P3 nL = normal(K, 0.49, 0.5), nR = normal(K, 0.51, 0.5);
        double dot = nL.x*nR.x + nL.y*nR.y + nL.z*nR.z;
        double ang = std::acos(std::clamp(dot, -1.0, 1.0)) * 180.0 / 3.14159265358979323846;
        printf("      [诊断] 法向夹角=%.1f°\n", ang);
        check("折脊两侧法向夹角显著（>5°）", ang > 5.0);
    }
    // [d] 折脊精确值（峰顶 z=0.25）
    {
        printf("[d] 折脊峰顶精确值\n");
        P3 p = K.eval(0.5, 0.5);
        printf("      [诊断] K(0.5,0.5)=(%.4f,%.4f,%.4f)\n", p.x,p.y,p.z);
        check("折脊峰顶 z≈0.25", std::abs(p.z - 0.25) < 1e-6 && std::abs(p.x - 0.5) < 1e-6);
    }
    // [e] 水平平面 z=0.1 横截 → 交线 1D 连续跨折脊
    {
        printf("[e] 平面 z=0.1 横截折脊曲面（交线 1D 跨折脊）\n");
        int hit = 0; bool lo = false, hi = false;
        for (int i = 0; i <= 100; ++i) {
            P3 p = K.eval(i/100.0, 0.5);
            if (std::abs(p.z - 0.1) < 0.02) { hit++; if (i < 50) lo = true; else hi = true; }
        }
        printf("      [诊断] 交线近点=%d (折脊两侧=%d)\n", hit, (int)lo + (int)hi);
        check("z=0.1 横截交线 1D 且跨折脊两侧", hit > 4 && lo && hi);
    }
    // [f] 平面 z=0.3 分离（高于峰顶 0.25）
    {
        printf("[f] 平面 z=0.3 分离\n");
        double mn = 1e300;
        for (int i = 0; i <= 50; ++i) for (int j = 0; j <= 10; ++j)
            mn = std::min(mn, std::abs(K.eval(i/50.0, j/10.0).z - 0.3));
        printf("      [诊断] 到 z=0.3 最小距=%.2e\n", mn);
        check("z=0.3 高于峰顶（分离，距>0.03）", mn > 0.03);
    }
    // [g] 无自交（全采样点唯一）
    {
        printf("[g] 折脊曲面无自交\n");
        double mn = 1e300;
        const int N = 12;
        std::vector<P3> pts;
        for (int i = 0; i <= N; ++i) for (int j = 0; j <= N; ++j)
            pts.push_back(K.eval(i/(double)N, j/(double)N));
        for (size_t a = 0; a < pts.size(); ++a) for (size_t b = a+1; b < pts.size(); ++b)
            mn = std::min(mn, dist(pts[a], pts[b]));
        printf("      [诊断] 采样点最小间距=%.2e\n", mn);
        check("折脊曲面无自交（采样唯一）", mn > 1e-6);
    }
    // [h] 折脊外区域正则（法向长度>0）
    {
        printf("[h] 折脊外区域正则\n");
        bool reg = true;
        for (int i = 0; i <= 30; ++i) for (int j = 0; j <= 4; ++j) {
            double u = i/30.0; if (std::abs(u-0.5) < 0.03) continue;   // 跳过折脊邻域
            if (dist(normal(K, u, j/4.0), P3{0,0,0}) < 1e-9) reg = false;
        }
        check("折脊外法向正则（长度>0）", reg);
    }
    // [i] aabb 与帐篷几何一致
    {
        printf("[i] aabb 一致\n");
        printf("      [诊断] z∈[%.4f,%.4f]\n", K.aabbMin().z, K.aabbMax().z);
        check("aabb zMin≈0, zMax≈0.25", std::abs(K.aabbMin().z) < 1e-6 && std::abs(K.aabbMax().z - 0.25) < 1e-6);
    }
    // [j] 折脊两侧子曲面块分别正则无奇点
    {
        printf("[j] 折脊两侧子块正则\n");
        bool ok = true;
        for (int s = 0; s < 2; ++s) {
            double u0 = s ? 0.51 : 0.0, u1 = s ? 1.0 : 0.49;
            for (int i = 0; i <= 24; ++i) for (int j = 0; j <= 4; ++j) {
                double u = u0 + (u1-u0)*i/24.0;
                if (dist(normal(K, u, j/4.0), P3{0,0,0}) < 1e-9) ok = false;
            }
        }
        check("折脊两侧子块均正则", ok);
    }

    // [k] 批量 90 断言：90 组水平平面横截折脊曲面（总断言=100）
    {
        printf("[k] 批量 90 组水平平面横截（z∈[0.03,0.228]，跨折脊两侧）\n");
        int okN = 0; const int N = 90;
        for (int k = 0; k < N; ++k) {
            const double zz = 0.03 + 0.0022 * k;   // 0.03 .. 0.228 < 峰顶 0.25
            int near = 0; bool lo = false, hi = false;
            for (int i = 0; i <= 100; ++i) {
                P3 p = K.eval(i/100.0, 0.5);
                if (std::abs(p.z - zz) < 0.02) { ++near; if (i < 50) lo = true; else hi = true; }
            }
            if (near > 4 && lo && hi) ++okN;   // 横截且交线跨折脊两侧
        }
        printf("      批量 90/90 全横截: %d\n", okN);
        pass += okN; fail += (N - okN);
    }

    printf("\n结论：测例⑫——曲面折脊（C⁰ 折缝）：值连续、导/法向折脊跳变、横截/分离、"
           "无自交、正则、aabb 一致。\n      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
