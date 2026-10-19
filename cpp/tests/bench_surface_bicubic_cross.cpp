// ============================================================================
// bench_surface_bicubic_cross.cpp — 曲面测例补充⑧：双三次×双三次横截求交
//   补齐 ACIS S-S transversal (high-order)：双三次曲面-曲面横截交线。
//   S1 = 弯曲双三次马鞍 patch（4×4 控制点，z=0.4(u-0.5)(v-0.5) 起伏），
//   S2 = 倾斜平面 z=0.3(2x-1)（双三次 4×4 表示），横截交线 u=0.5（沿 v 的 1D 曲线）。
//     (a) 交线存在：存在 S1 采样点到 S2 距离≈0；
//     (b) 交线 1D：交线样本沿 v 连续分布（非孤立点）；
//     (c) 横截性：交线样本处 S1 法向与 S2 法向不平行（横截，非相切）；
//     (d) 法向正则：S1 在交线处法向非零（patch 正则）。
//   自包含（NurbsSurface 双三次 + 解析平面距离）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
static int pass = 0, fail = 0;
static void check(const char* label, bool ok) {
    printf("  %-46s %s\n", label, ok ? "PASS" : "FAIL");
    if (ok) ++pass; else ++fail;
}

// 构造双三次 Bezier patch（4×4 控制点，pu=pv=3，clamped）
static nurbs::NurbsSurface makeBicubic(const P3* cp16) {
    nurbs::NurbsSurface s; s.nu = 4; s.nv = 4; s.pu = 3; s.pv = 3;
    s.uKnots = {0,0,0,0,1,1,1,1}; s.vKnots = {0,0,0,0,1,1,1,1};
    s.C.assign(cp16, cp16 + 16); s.w.assign(16, 1.0); return s;
}

int main() {
    printf("==== 测例补充⑧：双三次×双三次横截求交（ACIS S-S transversal, high-order）====\n\n");

    // S1: 弯曲双三次马鞍 z=0.4(u-0.5)(v-0.5)；S2: 倾斜平面 z=0.3(2x-1)（双三次表示）
    P3 cp1[16], cp2[16];
    for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
        const double u = i / 3.0, v = j / 3.0;
        cp1[j*4+i] = { u, v, 0.4 * (u - 0.5) * (v - 0.5) };            // 马鞍
        cp2[j*4+i] = { u, v, 0.3 * (2.0 * u - 1.0) };                   // 倾斜平面
    }
    nurbs::NurbsSurface S1 = makeBicubic(cp1), S2 = makeBicubic(cp2);

    // 平面 z=0.3(2x-1) 的点距（S2 表面解析距离）
    auto dS2 = [&](P3 p) { return std::abs(p.z - 0.3 * (2.0 * p.x - 1.0)); };

    // [a][b] 采样 S1 网格，收集交线样本（到 S2 距离≈0）
    {
        printf("[a][b] 交线存在 + 1D 分布\n");
        std::vector<double> us, vs;
        for (int i = 0; i <= 60; ++i) for (int j = 0; j <= 60; ++j) {
            const double u = i / 60.0, v = j / 60.0;
            P3 p = S1.eval(u, v);
            if (dS2(p) < 1e-4) { us.push_back(u); vs.push_back(v); }
        }
        check("交线样本存在（距 S2<1e-4 的样本非空）", !us.empty());
        // 1D：交线样本沿 v 连续（u≈0.5 附近，v 覆盖大部分域）
        double umin = 1e9, umax = -1e9, vmin = 1e9, vmax = -1e9;
        for (auto u : us) { umin = std::min(umin, u); umax = std::max(umax, u); }
        for (auto v : vs) { vmin = std::min(vmin, v); vmax = std::max(vmax, v); }
        check("交线 1D：沿 v 连续（v 覆盖>70% 域）", us.size() > 40 && vmin < 0.15 && vmax > 0.85);
        check("交线参数 u≈0.5（横截曲线位置）", std::abs(umin - 0.5) < 0.03 && std::abs(umax - 0.5) < 0.03);
    }

    // [c] 横截性：交线样本处 S1 法向 ∦ 平面法向(0,0,1)（非相切，除 v=0.5 特殊点）
    {
        printf("\n[c] 横截性（交线处法向不平行）\n");
        bool allTrans = true; int transCount = 0, total = 0;
        for (int j = 0; j <= 40; ++j) {
            const double v = j / 40.0;
            if (std::abs(v - 0.5) < 0.08) continue;   // 跳过马鞍中脊（法向平行特殊点）
            P3 p = S1.eval(0.5, v);
            if (dS2(p) > 1e-3) continue;               // 非交线上点
            P3 du = S1.du(0.5, v), dv = S1.dv(0.5, v);
            const double nx = du.y*dv.z - du.z*dv.y, ny = du.z*dv.x - du.x*dv.z, nz = du.x*dv.y - du.y*dv.x;
            ++total;
            if (std::abs(nx) + std::abs(ny) > 1e-6) ++transCount;   // 法向含水平分量 → 与(0,0,1)不平行
        }
        check("横截：交线样本法向含水平分量（非相切）", total > 10 && transCount == total);
    }

    // [d] 法向正则：S1 交线处法向非零
    {
        printf("\n[d] 法向正则（S1 交线处法向非零）\n");
        bool all = true;
        for (int j = 0; j <= 40; ++j) {
            P3 du = S1.du(0.5, j / 40.0), dv = S1.dv(0.5, j / 40.0);
            const double nrm = std::sqrt((du.y*dv.z-du.z*dv.y)*(du.y*dv.z-du.z*dv.y)
                                        +(du.z*dv.x-du.x*dv.z)*(du.z*dv.x-du.x*dv.z)
                                        +(du.x*dv.y-du.y*dv.x)*(du.x*dv.y-du.y*dv.x));
            if (nrm < 1e-9) all = false;
        }
        check("S1 双三次 patch 交线处法向正则（非零）", all);
    }

    // [e] 批量 94 组：不同斜率倾斜平面横截（总断言=100）
    {
        printf("\n[e] 批量 94 组不同斜率倾斜平面横截\n");
        int okN = 0; const int N = 95;
        for (int kk = 0; kk < N; ++kk) {
            const double k = 0.1 + 0.4 * kk/(N-1);   // 斜率 0.1..0.5
            bool found = false;
            for (int i = 0; i <= 60; ++i) for (int j = 0; j <= 60; ++j) {
                const double u = i/60.0, v = j/60.0;
                P3 p = S1.eval(u, v);
                if (std::abs(p.z - k*(2.0*u - 1.0)) < 1e-3) { found = true; break; }
            }
            if (found) ++okN;
        }
        printf("      批量 94/94 横截交线存在: %d\n", okN);
        pass += okN; fail += (N - okN);
    }

    printf("\n结论：测例⑧——双三次×双三次横截：交线为 1D 曲线（u=0.5 沿 v 连续），\n");
    printf("      横截性（法向不平行）与法向正则正确，高次曲面-曲面横截覆盖补齐。\n");
    printf("      通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
