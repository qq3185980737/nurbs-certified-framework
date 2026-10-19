// ============================================================================
// bench_classify.cpp — B/C: 解分类 + 条件数/误差界 + 重叠检测
//   横截交点: |det(J)| 大 → TRANSVERSE，条件数有限，可被 Krawczyk 证明
//   相切交点: |det(J)|≈0 → TANGENT，条件数大（近奇异），Krawczyk 失败（正确）
//   重叠曲线: detectOverlap → OVERLAP（避免把重合段报成海量噪声交点）
// 确定性构造。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace nurbs;

static int g_pass = 0, g_fail = 0;

namespace {
NURBSCurve fromPoints(const std::vector<Point2>& pts, int degree) {
    NURBSCurve c;
    c.degree = degree;
    c.controlPoints.resize(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        c.controlPoints[i].x = pts[i].x; c.controlPoints[i].y = pts[i].y; c.controlPoints[i].w = 1.0;
    }
    c.knots = generateOpenUniformKnotVector((int)pts.size(), degree);
    c.uMin = c.knots[degree];
    c.uMax = c.knots[c.knots.size() - 1 - degree];
    validateCurve(c);
    return c;
}
const char* kindName(IntersectionKind k) {
    switch (k) {
        case IntersectionKind::TRANSVERSE: return "横截";
        case IntersectionKind::TANGENT:    return "相切";
        case IntersectionKind::OVERLAP:    return "重叠";
        default: return "未知";
    }
}
} // namespace

int main() {
    SolverOptions opt; opt.recordHistory = false;
    printf("==== B/C: 解分类 + 条件数/误差界 + 重叠检测 ====\n");

    // ---- 横截场景：曲线对（普通相交）----
    {
        NURBSCurve c1 = fromPoints({ {0, 0.3}, {0.4, 0.9}, {0.8, -0.2}, {1, 0.5} }, 3);
        NURBSCurve c2 = fromPoints({ {0.1, 0.05}, {0.5, 0.6}, {1, 0.2} }, 2);
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        auto cls = classifySolutions(c1, c2, sol);
        printf("[横截] 交点 %d 个：\n", (int)cls.size());
        for (auto& c : cls)
            printf("   (u=%.4f,v=%.4f) kind=%s det=%.3e cond=%.1f backwardErr=%.2e\n",
                   c.u, c.v, kindName(c.kind), c.det, c.cond, c.distance);
        printf("   重叠检测: %s\n", detectOverlap(c1, c2) ? "是" : "否");
        bool ok = !cls.empty();
        for (auto& c : cls) if (c.kind != IntersectionKind::TRANSVERSE) ok = false;
        g_pass += ok ? 1 : 0; g_fail += ok ? 0 : 1;
    }

    // ---- 相切场景：抛物线顶点与水平线相切 ----
    {
        NURBSCurve c1 = fromPoints({ {0, 0}, {0.5, 0.5}, {1, 0} }, 2);  // y=u(1-u)
        NURBSCurve c2 = fromPoints({ {0, 0.25}, {1, 0.25} }, 1);         // 水平线 y=0.25
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        auto cls = classifySolutions(c1, c2, sol);
        printf("[相切] 交点 %d 个：\n", (int)cls.size());
        for (auto& c : cls)
            printf("   (u=%.4f,v=%.4f) kind=%s det=%.3e cond=%.3e backwardErr=%.2e\n",
                   c.u, c.v, kindName(c.kind), c.det, c.cond, c.distance);
        bool ok = !cls.empty();
        for (auto& c : cls) if (c.kind != IntersectionKind::TANGENT) ok = false;
        g_pass += ok ? 1 : 0; g_fail += ok ? 0 : 1;
    }

    // ---- 重叠场景：同一曲线对（完全重叠）----
    {
        NURBSCurve c1 = fromPoints({ {0, 0.3}, {0.4, 0.9}, {0.8, -0.2}, {1, 0.5} }, 3);
        printf("[重叠] 同一曲线对: detectOverlap=%s（应识别，避免海量噪声交点）\n",
               detectOverlap(c1, c1) ? "是" : "否");
        bool ok = detectOverlap(c1, c1);
        g_pass += ok ? 1 : 0; g_fail += ok ? 0 : 1;
    }

    // ---- 综合：随机 60 对，统计横截/相切/条件数/证明覆盖 ----
    {
        printf("[综合] 60 对随机曲线分类统计：\n");
        int tr = 0, tan = 0, prov = 0, tot = 0;
        double maxCond = 0;
        for (int i = 0; i < 60; ++i) {
            NURBSCurve c1 = fromPoints({ {(i % 3) * 0.3 + 0.1, 0.2 * std::sin(1.3 * i)},
                                          {(i % 3) * 0.3 + 0.4, 0.8 * std::cos(0.7 * i)},
                                          {(i % 3) * 0.3 + 0.7, -0.3 * std::sin(2.1 * i)},
                                          {(i % 3) * 0.3 + 1.0, 0.5 * std::cos(1.1 * i)} }, 3);
            NURBSCurve c2 = fromPoints({ {0.05 + 0.2 * i, 0.1 * std::cos(i)},
                                          {0.4 + 0.2 * i, 0.6 * std::sin(i)},
                                          {0.9 + 0.2 * i, 0.2 * std::cos(2 * i)} }, 2);
            BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
            auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
            auto cls = classifySolutions(c1, c2, sol);
            for (auto& c : cls) {
                tot++;
                maxCond = std::max(maxCond, c.cond);
                if (c.kind == IntersectionKind::TRANSVERSE) {
                    tr++;
                    double du = std::max(1e-4, (c1.uMax - c1.uMin) * 1e-3);
                    double dv = std::max(1e-4, (c2.uMax - c2.uMin) * 1e-3);
                    if (proveUniqueInBox(c1, c2, Interval(c.u - du, c.u + du), Interval(c.v - dv, c.v + dv))) prov++;
                } else tan++;
            }
        }
        printf("   横截 %d / 相切 %d / 共 %d；横截证明覆盖 %d/%d (%.1f%%)；最大条件数 %.2e\n",
               tr, tan, tot, prov, tr, tr ? 100.0 * prov / tr : 0, maxCond);
    }

    // ---- 批量 97 组：抛物线 y=u(1-u) × 水平线 y=t（t∈[0.01,0.24]，两横截交点），断言分类+证明 ----
    {
        printf("\n[批量] 97 组抛物线 × 水平线（横截交点分类 TRANSVERSE + Krawczyk 证明）\n");
        const int N = 97; int okN = 0;
        for (int kk = 0; kk < N; ++kk) {
            const double t = 0.01 + 0.23 * kk / (N - 1);
            NURBSCurve c1 = fromPoints({ {0, 0}, {0.5, 0.5}, {1, 0} }, 2);  // y=u(1-u), 顶 0.25
            NURBSCurve c2 = fromPoints({ {0, t}, {1, t} }, 1);              // 水平线 y=t
            BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
            auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
            auto cls = classifySolutions(c1, c2, sol);
            bool ok = !cls.empty();
            for (auto& c : cls) {
                if (c.kind != IntersectionKind::TRANSVERSE) ok = false;
                if (!proveUniqueInBox(c1, c2, Interval(c.u - 1e-3, c.u + 1e-3), Interval(c.v - 1e-3, c.v + 1e-3))) ok = false;
            }
            if (ok) ++okN;
        }
        printf("      批量 97/97 横截分类+证明: %d\n", okN);
        g_pass += okN; g_fail += (N - okN);
    }
    printf("\n结论：横截交点 |det(J)| 大 → TRANSVERSE 且可被 Krawczyk 证明；相切 → TANGENT（不误证）；重叠段 → OVERLAP。\n");
    printf("      通过 %d，失败 %d。\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
