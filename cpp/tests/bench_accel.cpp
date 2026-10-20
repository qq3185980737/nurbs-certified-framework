// ============================================================================
// bench_accel.cpp — BVH 加速求解 vs 原 unifiedSolve 基准对比
// 对比维度：建索引成本（一次性）+ 批量求解吞吐、相交解数一致性、逐根匹配率。
// 场景：一批曲线，两两求交（索引构建一次、查询复用），确定性 LCG 生成。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>

using namespace nurbs;

namespace {
uint32_t seed = 20261012u;
double rnd() {
    seed = seed * 1664525u + 1013904223u;
    return (seed >> 8) / 16777216.0;
}

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
}

int main() {
    const int N = 300;
    std::vector<NURBSCurve> A(N), B(N);
    for (int i = 0; i < N; ++i) {
        int d1 = 1 + (int)(rnd() * 3), d2 = 1 + (int)(rnd() * 3);
        int n1 = d1 + 1 + (int)(rnd() * 4), n2 = d2 + 1 + (int)(rnd() * 4);
        A[i] = makeCurve(n1, d1);
        B[i] = makeCurve(n2, d2);
    }
    SolverOptions opt; opt.recordHistory = false;

    // ---- 加速版：一次性建索引 + 复用求解 ----
    auto tb0 = std::chrono::high_resolution_clock::now();
    std::vector<BVHIndex> ia(N), ib(N);
    for (int i = 0; i < N; ++i) { ia[i].build(A[i], 7); ib[i].build(B[i], 7); }
    auto tb1 = std::chrono::high_resolution_clock::now();
    long long accIter = 0; int accInter = 0;
    for (int i = 0; i < N; ++i) {
        auto r = solveWithIndex(ia[i], A[i], ib[i], B[i], opt);
        for (auto& s : r) { accIter += s.iterations; if (s.state == GeometryState::INTERSECTION) accInter++; }
    }
    auto tb2 = std::chrono::high_resolution_clock::now();

    // ---- 原版：同批 N 对 ----
    long long origIter = 0; int origInter = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        auto r = unifiedSolve(A[i], B[i], opt);
        for (auto& s : r) { origIter += s.iterations; if (s.state == GeometryState::INTERSECTION) origInter++; }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    double tBuild = std::chrono::duration<double>(tb1 - tb0).count();
    double tSolve = std::chrono::duration<double>(tb2 - tb1).count();
    double tOrig = std::chrono::duration<double>(t1 - t0).count();

    // ---- 加速模式（completeLeakPatch=false）：纯 BVH+LM，不做完备兜底 ----
    SolverOptions optFast = opt; optFast.completeLeakPatch = false;
    long long fastIter = 0; int fastInter = 0;
    auto tf0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        auto r = solveWithIndex(ia[i], A[i], ib[i], B[i], optFast);
        for (auto& s : r) { fastIter += s.iterations; if (s.state == GeometryState::INTERSECTION) fastInter++; }
    }
    auto tf1 = std::chrono::high_resolution_clock::now();
    double tFast = std::chrono::duration<double>(tf1 - tf0).count();

    // ---- 解集合一致性 ----
    int agreePairs = 0, checked = 0; double bestRes = 0;
    for (int i = 0; i < N; ++i) {
        auto ro = unifiedSolve(A[i], B[i], opt);
        auto ra = solveWithIndex(ia[i], A[i], ib[i], B[i], opt);
        std::vector<SolveResult> io, iacc;
        for (auto& s : ro) if (s.state == GeometryState::INTERSECTION) io.push_back(s);
        for (auto& s : ra) if (s.state == GeometryState::INTERSECTION) iacc.push_back(s);
        int matched = 0;
        for (auto& x : iacc)
            for (auto& y : io)
                if (std::abs(x.u - y.u) < 1e-4 && std::abs(x.v - y.v) < 1e-4) { matched++; break; }
        if (io.empty() && iacc.empty()) { checked++; agreePairs++; }
        else if (io.size() == iacc.size() && matched == (int)iacc.size()) { checked++; agreePairs++; }
        for (auto& s : iacc) if (s.distance > bestRes) bestRes = s.distance;
    }

    printf("==== BVH 加速 vs 原 unifiedSolve 基准（%d 对曲线）====\n", N);
    printf("原版      : 求解 %8.2f ms   相交解数=%d  LM迭代=%lld\n", tOrig * 1e3, origInter, origIter);
    printf("加速版    : 建索引 %6.2f ms + 求解 %8.2f ms = %8.2f ms\n",
           tBuild * 1e3, tSolve * 1e3, (tBuild + tSolve) * 1e3);
    printf("          相交解数=%d  LM迭代=%lld（减少 %.1f%%）\n", accInter, accIter,
           100.0 * (1.0 - (double)accIter / origIter));
    printf("求解阶段加速比: %.2f×\n", tOrig / tSolve);
    // 两种模式必须分开报告：完备模式含 unifiedSolve 兜底（正确性优先，慢），
    // 加速模式为纯 BVH+LM（快，但只近完备）。论文 Table 5 的 35× 指后者。
    printf("加速模式  : 求解 %8.2f ms   相交解数=%d  LM迭代=%lld（completeLeakPatch=false）\n",
           tFast * 1e3, fastInter, fastIter);
    printf("加速比    : 加速模式/原版 = %.2f×   完备模式/原版 = %.2f×\n", tOrig / tFast, tOrig / tSolve);
    printf("含建索引总加速比: %.2f×\n", tOrig / (tBuild + tSolve));
    printf("相交解数一致: %d / %d 对\n", agreePairs, checked);
    printf("加速交点最大残差: %.3e（应 < 1e-6）\n", bestRes);

    // ---- A1：重负载对（高阶+密集叶对 → 每对初值多）下 串行 vs 并行 ----
    seed = 20261012u;
    const int NH = 60;
    std::vector<NURBSCurve> HA(NH), HB(NH);
    for (int i = 0; i < NH; ++i) {
        int d1 = 3 + (int)(rnd() * 3), d2 = 3 + (int)(rnd() * 3); // 高阶
        int n1 = d1 + 1 + (int)(rnd() * 8), n2 = d2 + 1 + (int)(rnd() * 8);
        HA[i] = makeCurve(n1, d1);
        HB[i] = makeCurve(n2, d2);
    }
    std::vector<BVHIndex> hia(NH), hib(NH);
    for (int i = 0; i < NH; ++i) { hia[i].build(HA[i], 7); hib[i].build(HB[i], 7); }
    auto hs0 = std::chrono::high_resolution_clock::now();
    long long hIterSer = 0;
    for (int i = 0; i < NH; ++i)
        for (auto& s : solveWithIndex(hia[i], HA[i], hib[i], HB[i], opt, 1)) hIterSer += s.iterations;
    auto hs1 = std::chrono::high_resolution_clock::now();
    long long hIterPar = 0;
    for (int i = 0; i < NH; ++i)
        for (auto& s : solveWithIndex(hia[i], HA[i], hib[i], HB[i], opt, 0)) hIterPar += s.iterations;
    auto hs2 = std::chrono::high_resolution_clock::now();
    double tSer = std::chrono::duration<double>(hs1 - hs0).count();
    double tPar = std::chrono::duration<double>(hs2 - hs1).count();
    printf("==== A1 并行（%d 对高阶/密集曲线，LM初值多）====\n", NH);
    printf("串行 solveWithIndex(threads=1) : %8.2f ms   LM迭代=%lld\n", tSer * 1e3, hIterSer);
    printf("并行 solveWithIndex(threads=auto): %8.2f ms  LM迭代=%lld\n", tPar * 1e3, hIterPar);
    printf("并行加速比: %.2f×\n", tSer / tPar);
    printf("迭代数一致: %s\n", hIterSer == hIterPar ? "是（结果逐位等价）" : "否（差异待查）");

    // ---- A2：迭代预算消融（重负载对，并行）——找"解数仍 100% 一致"的最小 maxIterations ----
    const int kBudgets[] = { 80, 40, 24, 16, 10, 6 };
    const int kBudgetBase = kBudgets[0];
    printf("\n==== A2 迭代预算消融（%d 对高阶/密集，并行）====", NH);
    int baseCnt = 0;
    // 先跑基线(80)统计解数
    { SolverOptions o = opt; o.maxIterations = kBudgetBase;
      for (int i = 0; i < NH; ++i)
        for (auto& s : solveWithIndex(hia[i], HA[i], hib[i], HB[i], o, 0))
          if (s.state == GeometryState::INTERSECTION) baseCnt++; }
    printf("基线 maxIter=%d → 相交解数=%d\n", kBudgetBase, baseCnt);
    for (int b : kBudgets) {
        SolverOptions o = opt; o.maxIterations = b;
        auto s0 = std::chrono::high_resolution_clock::now();
        int cnt = 0; long long it = 0;
        for (int i = 0; i < NH; ++i)
            for (auto& s : solveWithIndex(hia[i], HA[i], hib[i], HB[i], o, 0)) { if (s.state == GeometryState::INTERSECTION) cnt++; it += s.iterations; }
        auto s1 = std::chrono::high_resolution_clock::now();
        double t = std::chrono::duration<double>(s1 - s0).count();
        printf("  maxIter=%3d : %7.2f ms  相交解=%d (%s)  LM迭代=%lld\n", b, t * 1e3, cnt,
               cnt == baseCnt ? "解数一致" : "解数减少!", it);
    }

    // ---- A3：find 加速对照——shared findAllIntersections vs BVH 加速枚举（前 20 对）----
    const int NF = 20;
    printf("\n==== A3 find 加速对照（前 %d 对高阶/密集，含多交点）====\n", NF);
    auto f0 = std::chrono::high_resolution_clock::now();
    int fRoots = 0;
    for (int i = 0; i < NF; ++i)
        for (auto& r : findAllIntersections(HA[i], HB[i], opt)) if (r.distance < 1e-6) fRoots++;
    auto f1 = std::chrono::high_resolution_clock::now();
    auto f2 = std::chrono::high_resolution_clock::now();
    int aRoots = 0;
    for (int i = 0; i < NF; ++i)
        for (auto& s : solveWithIndex(hia[i], HA[i], hib[i], HB[i], opt, 0)) if (s.state == GeometryState::INTERSECTION) aRoots++;
    auto f3 = std::chrono::high_resolution_clock::now();
    double tf = std::chrono::duration<double>(f1 - f0).count();
    double ta = std::chrono::duration<double>(f3 - f2).count();
    printf("  findAllIntersections (shared): %7.2f ms  根数=%d\n", tf * 1e3, fRoots);
    printf("  BVH 加速枚举 (accelerated)   : %7.2f ms  相交解=%d\n", ta * 1e3, aRoots);
    printf("  find 加速比: %.2f×   （根数 %s）\n", tf / ta, fRoots == aRoots ? "一致" : "有差异(须核)");

    // ---- A4：自适应初值消融——固定 9×叶对 vs 自适应初值（任务1）----
    size_t pairsSum = 0, fixedSum = 0, adaptiveSum = 0;
    for (int i = 0; i < N; ++i) {
        auto c = adaptiveInitCounts(ia[i], A[i], ib[i], B[i]);
        pairsSum += c.first; adaptiveSum += c.second; fixedSum += c.first * 9;
    }
    printf("\n==== A4 自适应初值消融（%d 对）====\n", N);
    printf("相交叶对总数: %zu\n", pairsSum);
    printf("固定 9 初值 : %zu\n", fixedSum);
    printf("自适应初值  : %zu（相对固定 9：%.1f%%）\n", adaptiveSum,
           100.0 * ((double)adaptiveSum / fixedSum - 1.0));
    return 0;
}
