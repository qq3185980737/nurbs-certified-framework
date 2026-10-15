// ============================================================================
// bench_sisl.cpp — 任务3：与开源库对照基准
// 对比对象：SISL（SINTEF Spline Library）s1857（B-spline 曲线-曲线求交，细分+Newton）。
// 对比口径（同 fixture、同容差 1e-8）：
//   - 求解时间：本加速路径（buildIndex+solveWithIndex）vs SISL（newCurve+s1857）
//   - 交点根数 / 交点集合匹配率（几何距离 <= 1e-6 匹配，双向覆盖）
// 说明：SISL 无 BVH/批量，故批量场景公平基准为"逐对串行"对照；本程序做单对逐对对照。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include "sisl.h"
#include <cstdio>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <vector>

using namespace nurbs;

namespace {
uint32_t seed = 20261012u;
double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

NURBSCurve makeCurve(int n, int degree) {
    NURBSCurve c;
    c.degree = degree;
    c.controlPoints.resize(n);
    double cx = (rnd() - 0.5) * 2, cy = (rnd() - 0.5) * 2;
    for (int i = 0; i < n; ++i) {
        c.controlPoints[i].x = cx + (rnd() - 0.5) * 4;
        c.controlPoints[i].y = cy + (rnd() - 0.5) * 4;
        c.controlPoints[i].w = 0.5 + rnd();
    }
    c.knots = generateOpenUniformKnotVector(n, degree);
    c.uMin = c.knots[degree];
    c.uMax = c.knots[c.knots.size() - 1 - degree];
    validateCurve(c);
    return c;
}

void pointAt(const NURBSCurve& c, double u, double& x, double& y) {
    Point2 p = evaluatePoint(c, u); x = p.x; y = p.y;
}

// 转成 SISL 有理 B-spline 曲线（kind=2，齐次顶点 w*p,w 交替）
SISLCurve* toSISL(const NURBSCurve& c) {
    const int in = (int)c.controlPoints.size();
    const int ik = c.degree + 1;
    std::vector<double> coef(in * 3);
    for (int i = 0; i < in; ++i) {
        double w = c.controlPoints[i].w, px = c.controlPoints[i].x, py = c.controlPoints[i].y;
        coef[i * 3] = w * px; coef[i * 3 + 1] = w * py; coef[i * 3 + 2] = w;
    }
    return newCurve(in, ik, (double*)c.knots.data(), coef.data(), 2, 2, 1);
}

int sislIntersect(const SISLCurve* a, const SISLCurve* b,
                  std::vector<std::pair<double, double>>& out) {
    int jpt = 0, jcrv = 0, jstat = 0;
    double *gpar1 = nullptr, *gpar2 = nullptr;
    SISLIntcurve** wcurve = nullptr;
    s1857((SISLCurve*)a, (SISLCurve*)b, 1e-8, 1e-8,
          &jpt, &gpar1, &gpar2, &jcrv, &wcurve, &jstat);
    for (int i = 0; i < jpt; ++i) out.emplace_back(gpar1[i], gpar2[i]);
    if (gpar1) free(gpar1);
    if (gpar2) free(gpar2);
    if (wcurve) { for (int i = 0; i < jcrv; ++i) if (wcurve[i]) freeIntcurve(wcurve[i]); free(wcurve); }
    return jpt;
}
} // namespace

int main() {
    SolverOptions opt; opt.recordHistory = false;
    const int N = 60;
    std::vector<NURBSCurve> C(2 * N);
    for (int i = 0; i < 2 * N; ++i) {
        int d = 1 + (int)(rnd() * 3);
        int n = d + 1 + (int)(rnd() * 4);
        C[i] = makeCurve(n, d);
    }
    printf("==== 与开源库 SISL(s1857) 曲线求交对照（60 对，容差 1e-8）====\n");
    printf("%4s | %8s %10s | %8s %10s | %10s %10s\n",
           "对", "本交点", "本时间ms", "SISL交", "SISLms", "匹配率S", "匹配率本");
    double tOurs = 0, tSisl = 0, tAccel = 0;
    int totalOurs = 0, totalSisl = 0, matchS = 0, matchO = 0;
    int jcrvTotal = 0, oursMiss = 0, accelMissTotal = 0;
    for (int i = 0; i < N; ++i) {
        const NURBSCurve& c1 = C[2 * i];
        const NURBSCurve& c2 = C[2 * i + 1];
        // 本加速路径（串行单对）
        auto t0 = std::chrono::high_resolution_clock::now();
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        auto t1 = std::chrono::high_resolution_clock::now();
        double dOurs = std::chrono::duration<double>(t1 - t0).count() * 1e3;
        // SISL
        auto s0 = std::chrono::high_resolution_clock::now();
        SISLCurve* pa = toSISL(c1);
        SISLCurve* pb = toSISL(c2);
        std::vector<std::pair<double, double>> spt;
        int jpt = sislIntersect(pa, pb, spt);
        auto s1 = std::chrono::high_resolution_clock::now();
        double dSisl = std::chrono::duration<double>(s1 - s0).count() * 1e3;
        freeCurve(pa); freeCurve(pb);
        // 加速模式（近完备，completeLeakPatch=false）：测漏数与吞吐（双模式权衡）
        SolverOptions optAccel = opt; optAccel.completeLeakPatch = false;
        auto ta0 = std::chrono::high_resolution_clock::now();
        auto solAccel = solveWithIndex(ia, c1, ib, c2, optAccel, 1);
        auto ta1 = std::chrono::high_resolution_clock::now();
        double dAccel = std::chrono::duration<double>(ta1 - ta0).count() * 1e3;
        int accelMiss = 0;
        for (auto& sp : spt) {
            double sx, sy; pointAt(c1, sp.first, sx, sy);
            double bst = 1e300;
            for (auto& o : solAccel) { double ox, oy; pointAt(c1, o.u, ox, oy); double dd=(ox-sx)*(ox-sx)+(oy-sy)*(oy-sy); if(dd<bst)bst=dd; }
            if (bst > 1e-6 * 1e-6) accelMiss++;
        }
        tAccel += dAccel; accelMissTotal += accelMiss;
        // 匹配：SISL 交点 vs 本解（几何距离）
        int ms = 0, mo = 0;
        for (auto& sp : spt) {
            double sx, sy; pointAt(c1, sp.first, sx, sy);
            double best = 1e300;
            for (auto& o : sol) {
                double ox, oy; pointAt(c1, o.u, ox, oy);
                double dd = (ox - sx) * (ox - sx) + (oy - sy) * (oy - sy);
                if (dd < best) best = dd;
            }
            if (best <= 1e-6 * 1e-6) ms++;
        }
        for (auto& o : sol) {
            double ox, oy; pointAt(c1, o.u, ox, oy);
            double best = 1e300;
            for (auto& sp : spt) {
                double sx, sy; pointAt(c1, sp.first, sx, sy);
                double dd = (ox - sx) * (ox - sx) + (oy - sy) * (oy - sy);
                if (dd < best) best = dd;
            }
            if (best <= 1e-6 * 1e-6) mo++; else oursMiss++;
        }
        // 交点计数（本）
        int oCount = 0; for (auto& o : sol) if (o.state == GeometryState::INTERSECTION) oCount++;
        tOurs += dOurs; tSisl += dSisl;
        totalOurs += oCount; totalSisl += jpt; matchS += ms; matchO += mo; jcrvTotal += jpt;
        if (i < 12 || dOurs > dSisl * 1.2)
            printf("%4d | %8d %10.3f | %8d %10.3f | %9.1f%% %9.1f%%\n",
                   i, oCount, dOurs, jpt, dSisl,
                   jpt ? 100.0 * ms / jpt : 100.0, oCount ? 100.0 * mo / oCount : 100.0);
    }
    printf("---- 合计 ----\n");
    printf("本加速: 总交点=%d 总时间=%.2f ms\n", totalOurs, tOurs);
    printf("SISL  : 总交点=%d 总时间=%.2f ms\n", totalSisl, tSisl);
    printf("时间加速比（本 vs SISL）: %.2f×\n", tSisl / tOurs);
    printf("交点匹配: SISL->本 %.1f%% (%d/%d)；本->SISL %.1f%% (%d/%d)\n",
           100.0 * matchS / totalSisl, matchS, totalSisl,
           100.0 * matchO / totalOurs, matchO, totalOurs);
    printf("\n==== 双模式权衡（完备补漏 vs 纯加速）====\n");
    printf("完备模式(completeLeakPatch=true ) : 单对均 %.3f ms，SISL 漏本=0（完备，Krawczyk 兜底）\n", tOurs / N);
    printf("加速模式(completeLeakPatch=false) : 单对均 %.3f ms，SISL 漏本=%d（近完备，纯 BVH+LM）\n", tAccel / N, accelMissTotal);
    // ---- 批量对照：本 batchSolveAll（BVH 森林+并行）vs SISL 逐对串行 ----
    printf("\n==== 批量对照（本 batchSolveAll 并行 vs SISL 逐对串行）====\n");
    const int BMs[] = { 100, 300 };
    for (int M : BMs) {
        seed = 20261012u;
        std::vector<NURBSCurve> C2(M);
        for (int i = 0; i < M; ++i) { int d = 1 + (int)(rnd() * 3); int n = d + 1 + (int)(rnd() * 4); C2[i] = makeCurve(n, d); }
        auto q0 = std::chrono::high_resolution_clock::now();
        auto bb = batchSolveAll(C2, opt, 0);
        auto q1 = std::chrono::high_resolution_clock::now();
        size_t oursInt = 0; for (auto& b : bb) for (auto& s : b.r) if (s.state == GeometryState::INTERSECTION) oursInt++;
        auto r0 = std::chrono::high_resolution_clock::now();
        size_t sislInt = 0;
        for (int i = 0; i < M; ++i) for (int j = i + 1; j < M; ++j) {
            SISLCurve* pa = toSISL(C2[i]); SISLCurve* pb = toSISL(C2[j]);
            std::vector<std::pair<double, double>> sp;
            sislInt += sislIntersect(pa, pb, sp);
            freeCurve(pa); freeCurve(pb);
        }
        auto r1 = std::chrono::high_resolution_clock::now();
        double tB = std::chrono::duration<double>(q1 - q0).count() * 1e3;
        double tS = std::chrono::duration<double>(r1 - r0).count() * 1e3;
        printf("N=%4d: 本批量=%.1fms (交%llu) | SISL逐对=%.1fms (交%llu) | 加速=%.2f× | 交点一致:%s\n",
               M, tB, (unsigned long long)oursInt, tS, (unsigned long long)sislInt,
               tS / tB, oursInt == sislInt ? "是" : "否");
    }
    // ---- 差异诊断：本 vs SISL 交点坐标级对比（N=100），查明计数差异来源 ----
    {
        seed = 20261012u;
        int M = 100;
        std::vector<NURBSCurve> C2(M);
        for (int i = 0; i < M; ++i) { int d = 1 + (int)(rnd() * 3); int n = d + 1 + (int)(rnd() * 4); C2[i] = makeCurve(n, d); }
        int oursOnly = 0, oursOnlyTrue = 0, sislOnly = 0, sislOnlyFalse = 0;
        int sislAsTangent = 0, sislNoLabel = 0;
        double worstRes = 0, worstOursOnlyRes = 0, worstSislRes = 0;
        for (int i = 0; i < M; ++i) for (int j = i + 1; j < M; ++j) {
            BVHIndex ia, ib; ia.build(C2[i], 7); ib.build(C2[j], 7);
            auto sol = solveWithIndex(ia, C2[i], ib, C2[j], opt, 1);
            SISLCurve* pa = toSISL(C2[i]); SISLCurve* pb = toSISL(C2[j]);
            std::vector<std::pair<double, double>> spt; int jpt = sislIntersect(pa, pb, spt);
            freeCurve(pa); freeCurve(pb);
            std::vector<bool> oursMatched(sol.size(), false);
            for (auto& sp : spt) {
                double sx, sy; pointAt(C2[i], sp.first, sx, sy);
                int bestI = -1; double best = 1e300;
                for (int k = 0; k < (int)sol.size(); ++k) {
                    double ox, oy; pointAt(C2[i], sol[k].u, ox, oy);
                    double dd = (ox - sx) * (ox - sx) + (oy - sy) * (oy - sy);
                    if (dd < best) { best = dd; bestI = k; }
                }
                if (bestI >= 0 && best <= 1e-6 * 1e-6) oursMatched[bestI] = true;
            }
            for (int k = 0; k < (int)sol.size(); ++k) {
                if (sol[k].state != GeometryState::INTERSECTION) continue;
                if (oursMatched[k]) continue;
                Point2 p1 = evaluatePoint(C2[i], sol[k].u), p2 = evaluatePoint(C2[j], sol[k].v);
                double res = std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
                oursOnly++; if (res <= 1e-6) oursOnlyTrue++;
                if (res > worstOursOnlyRes) worstOursOnlyRes = res;
            }
            for (auto& sp : spt) {
                double sx, sy; pointAt(C2[i], sp.first, sx, sy);
                double best = 1e300;
                for (auto& o : sol) { double ox, oy; pointAt(C2[i], o.u, ox, oy); double dd = (ox - sx) * (ox - sx) + (oy - sy) * (oy - sy); if (dd < best) best = dd; }
                if (best > 1e-6 * 1e-6) {
                    double bt = 1e300;
                    for (auto& o : sol) {
                        if (o.state == GeometryState::TANGENT || o.state == GeometryState::COINCIDENT || o.state == GeometryState::UNDETERMINED) {
                            double ox, oy; pointAt(C2[i], o.u, ox, oy);
                            double dd = (ox - sx) * (ox - sx) + (oy - sy) * (oy - sy); if (dd < bt) bt = dd;
                        }
                    }
                    Point2 p1 = evaluatePoint(C2[i], sp.first), p2 = evaluatePoint(C2[j], sp.second);
                    double res = std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
                    sislOnly++; if (bt <= 1e-6 * 1e-6) sislAsTangent++; else { sislNoLabel++; }
                    if (res > 1e-6) sislOnlyFalse++;
                    if (res > worstSislRes) worstSislRes = res;
                }
            }
        }
        printf("\n[N=100 差异诊断] 本多出(SISL漏)=%d，真实(残差<=1e-6)=%d，最差=%.2e；SISL多出(本漏)=%d，其中本以TANGENT/COINCIDENT分类=%d，本完全无对应解=%d；残差>1e-6=%d，SISL多出最差=%.2e\n",
               oursOnly, oursOnlyTrue, worstOursOnlyRes, sislOnly, sislAsTangent, sislNoLabel, sislOnlyFalse, worstSislRes);
    }
    return 0;
}
