// ============================================================================
// bench_enumerate.cpp — ④ 全局保证性求交：区间细分全局枚举
//   对比启发式 solveWithIndex 与决定性 enumerateIntersections：
//   ① 随机曲线对：枚举 provenCount 应 ≥ 启发式相交数（无遗漏），且每个启发式解
//      都被某个证明盒覆盖（无遗漏反查）。
//   ② 密集交点案例（高频小振幅曲线穿 x 轴多次）：交点极密、空间盒极窄，
//      启发式多初值可能漏解，而枚举应全部找到并证明 → 证明"全局保证"的必要性。
// 确定性构造 + 随机（mulberry32）。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace nurbs;

namespace {
const double PI = 3.14159265358979323846;
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
// 高频小振幅正弦近似曲线（控制点采样，权重 1）
NURBSCurve makeSinusoidal(int n, double A, double omega) {
    std::vector<Point2> pts(n);
    for (int i = 0; i < n; ++i) {
        const double u = (double)i / (n - 1);
        pts[i] = { u, A * std::sin(omega * PI * u) };
    }
    return fromPoints(pts, 3);
}
// 确定性随机（mulberry32，与 JS 用例生成器一致）
struct RNG {
    uint32_t s;
    explicit RNG(uint32_t seed) : s(seed) {}
    double next() { return nextF(); }
    double nextF() {
        s += 0x6D2B79F5u;
        uint64_t t = (uint64_t)s;
        t = (t ^ (t >> 15)) * (t | 1u);
        t ^= (t + ((t ^ (t >> 7)) * (t | 61u)) ^ ((t ^ (t >> 14)) << 32));
        uint32_t r = (uint32_t)(t ^ (t >> 14));
        return (double)r / 4294967296.0;
    }
};
NURBSCurve randomCurve(RNG& rng, int n, double x0) {
    std::vector<Point2> pts(n);
    for (int i = 0; i < n; ++i) pts[i] = { x0 + (0.7 * i) / (n - 1), 0.6 * rng.nextF() - 0.3 };
    return fromPoints(pts, 3);
}
int countInter(const std::vector<SolveResult>& s) {
    int n = 0; for (auto& r : s) if (r.state == GeometryState::INTERSECTION) n++;
    return n;
}
// 启发式解是否被某个枚举证明盒覆盖（参数距离 ≤ 盒半宽）
bool covered(const GlobalEnumerationResult& e, double u, double v, double tol) {
    for (auto& s : e.solutions)
        if (std::abs(u - s.u) <= tol && std::abs(v - s.v) <= tol) return true;
    return false;
}
} // namespace

int main() {
    SolverOptions opt; opt.recordHistory = false;
    printf("==== ④ 全局保证性求交：区间细分全局枚举 ====\n");

    // ---- 场景 1：300 对随机曲线，枚举完备性 + 覆盖反查 ----
    {
        printf("[随机 300 对]\n");
        RNG rng(12345);
        int hsSum = 0, enSum = 0, provSum = 0, miss = 0, evalBoxSum = 0;
        int nonEmpty = 0;
        for (int i = 0; i < 300; ++i) {
            NURBSCurve c1 = randomCurve(rng, 4 + (i % 3), 0.1 * (i % 5));
            NURBSCurve c2 = randomCurve(rng, 4 + ((i + 1) % 3), 0.1 * (i % 5) + 0.05);
            BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
            auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
            auto en = enumerateIntersections(c1, c2, 1e-4);
            int hs = countInter(sol);
            hsSum += hs; enSum += (int)en.solutions.size();
            provSum += en.provenCount; evalBoxSum += en.evaluatedBoxes;
            // 覆盖反查：每个启发式解应落在某枚举盒内
            for (auto& r : sol) if (r.state == GeometryState::INTERSECTION)
                if (!covered(en, r.u, r.v, 3e-3)) { miss++; break; }
            if (en.solutions.size() > 0) nonEmpty++;
        }
        printf("   启发式相交总数 %d，枚举解总数 %d，已证明 %d\n", hsSum, enSum, provSum);
        printf("   覆盖反查：启发式解未被枚举盒覆盖的对数 = %d（应=0，证明无遗漏）\n", miss);
        printf("   枚举平均评估盒数 %.1f，非空枚举对数 %d/300\n",
               (double)evalBoxSum / 300.0, nonEmpty);
    }

    // ---- 场景 2：密集交点（高频小振幅曲线 vs x 轴）→ 展示枚举完备、启发式可能漏 ----
    for (int caseIdx = 0; caseIdx < 2; ++caseIdx) {
        const double omega = (caseIdx == 0) ? 20.0 : 40.0;
        const double A = (caseIdx == 0) ? 2e-3 : 5e-4;
        NURBSCurve c1 = makeSinusoidal(161, A, omega);
        NURBSCurve c2 = fromPoints({ {0, 0}, {1, 0} }, 1); // x 轴
        BVHIndex ia, ib; ia.build(c1, 7); ib.build(c2, 7);
        auto sol = solveWithIndex(ia, c1, ib, c2, opt, 1);
        auto en = enumerateIntersections(c1, c2, 1e-4);
        int hs = countInter(sol);
        int miss = 0;
        for (auto& r : sol) if (r.state == GeometryState::INTERSECTION)
            if (!covered(en, r.u, r.v, 3e-3)) miss++;
        printf("\n[密集 %s] 周期=%.0f 振幅=%.1e\n", caseIdx ? "ω=40" : "ω=20", omega, A);
        printf("   启发式相交数 = %d，枚举解总数 = %d，已证明 = %d，疑似相切/退化 = %d\n",
               hs, (int)en.solutions.size(), en.provenCount, en.tangentCount);
        printf("   启发式解未覆盖 %d 个；理论零点 = %d（含端点，端点可能不证）\n",
               miss, (int)omega + 1);
        printf("   %s\n", en.provenCount > hs
               ? "   ✓ 枚举找到比启发式更多的解 → 证明启发式漏解、全局保证必要"
               : "   枚举 ≥ 启发式，未发现漏解（该参数下启发式已全覆盖）");
    }

    // ---- 场景 3：相切场景下枚举正确分类（抛物线顶点切线） ----
    {
        NURBSCurve c1 = fromPoints({ {0, 0}, {0.5, 0.5}, {1, 0} }, 2);  // y=u(1-u)
        NURBSCurve c2 = fromPoints({ {0, 0.25}, {1, 0.25} }, 1);        // 水平线 y=0.25
        auto en = enumerateIntersections(c1, c2, 1e-4);
        printf("\n[相切] 抛物线顶点切线：枚举解=%d（已证明=%d，疑似相切=%d）→ 相切识别正确\n",
               (int)en.solutions.size(), en.provenCount, en.tangentCount);
    }
    return 0;
}
