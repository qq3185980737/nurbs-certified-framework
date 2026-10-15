// ============================================================================
// test_robust.cpp — 鲁棒性 / 防御性测试（C++ 独立，不依赖 JS 夹具）
// ----------------------------------------------------------------------------
// 目标：算法库面对 CAD 真实退化/极端输入时不崩溃、不越界、不死循环，
//       且对结构非法输入在入口（validateCurve）明确抛异常。
// 覆盖：退化点 / 共线重叠 / 端点相接 / 折返自交 / 极端尺度 / 高重节点，
//       以及非法结构（knot 长度错、负权重、NaN、参数域倒置、空控制点）。
// ============================================================================
#include "nurbs/nurbs_core.hpp"
#include "nurbs/unified_solver.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <stdexcept>

using namespace nurbs;

static int g_pass = 0, g_fail = 0;
#define CHECK(name, cond) do { \
    if (cond) { g_pass++; } \
    else { g_fail++; std::printf("[robust] FAIL: %s\n", name); } \
} while (0)

// 标准 clamped 均匀 knot：长度 n+p+1
static std::vector<double> clampedKnots(int n, int p) {
    std::vector<double> k;
    for (int i = 0; i <= p; ++i) k.push_back(0.0);
    const int mid = n - p - 1;
    for (int i = 1; i <= mid; ++i) k.push_back(i / (double)(mid + 1));
    for (int i = 0; i <= p; ++i) k.push_back(1.0);
    return k; // length = n+p+1
}
static NURBSCurve mkCurve(const std::vector<ControlPoint>& pts, const std::vector<double>& knots, int deg) {
    NURBSCurve c;
    c.degree = deg;
    c.controlPoints = pts;
    c.knots = knots;
    c.uMin = knots[deg];
    c.uMax = knots[knots.size() - 1 - deg];
    return c;
}
static ControlPoint C(double x, double y, double w = 1.0) { return ControlPoint{ x, y, w }; }

static bool isFinite(double x) { return std::isfinite(x); }

int main() {
    // ---- 1) 退化几何：evaluatePoint 必须有限且不崩溃 ----
    {
        const NURBSCurve pointCurve = mkCurve({ C(2,3), C(2,3), C(2,3), C(2,3) }, clampedKnots(4, 3), 3);
        const NURBSCurve collinear = mkCurve({ C(0,0), C(1,0), C(2,0), C(3,0) }, clampedKnots(4, 3), 3);
        const NURBSCurve fold = mkCurve({ C(0,0), C(2,2), C(4,0), C(2,-2), C(0,0) }, clampedKnots(5, 4), 4);
        const NURBSCurve big = mkCurve({ C(0,0), C(1e6,0), C(2e6,1e6), C(3e6,0) }, clampedKnots(4, 3), 3);
        const NURBSCurve tiny = mkCurve({ C(0,0), C(1e-6,0), C(2e-6,1e-6), C(3e-6,0) }, clampedKnots(4, 3), 3);
        const std::vector<NURBSCurve> curves = { pointCurve, collinear, fold, big, tiny };
        bool allFinite = true;
        for (const auto& c : curves) {
            const double us[3] = { c.uMin, (c.uMin + c.uMax) * 0.5, c.uMax };
            for (double u : us) {
                Point2 p;
                try { p = evaluatePoint(c, u); }
                catch (...) { allFinite = false; break; }
                if (!isFinite(p.x) || !isFinite(p.y)) { allFinite = false; break; }
            }
        }
        CHECK("退化曲线 evaluatePoint 有限", allFinite);
    }

    // ---- 2) 退化几何：一阶/二阶导不崩溃（允许点退化时切向量为 0/NaN）----
    {
        const NURBSCurve pointCurve = mkCurve({ C(2,3), C(2,3), C(2,3), C(2,3) }, clampedKnots(4, 3), 3);
        bool noCrash = true;
        try { evaluateFirstDerivative(pointCurve, 0.5); }
        catch (...) { noCrash = false; }
        try { evaluateSecondDerivative(pointCurve, 0.5); }
        catch (...) { noCrash = false; }
        CHECK("退化点曲线导数不崩溃", noCrash);
    }

    // ---- 3) 退化对：unifiedSolve / findAllIntersections 不崩溃、结果有限、迭代有界 ----
    {
        const NURBSCurve pointCurve = mkCurve({ C(2,3), C(2,3), C(2,3), C(2,3) }, clampedKnots(4, 3), 3);
        const NURBSCurve collinear = mkCurve({ C(0,0), C(1,0), C(2,0), C(3,0) }, clampedKnots(4, 3), 3);
        const NURBSCurve endA = mkCurve({ C(0,0), C(2,1), C(4,0) }, clampedKnots(3, 2), 2);
        const NURBSCurve endB = mkCurve({ C(4,0), C(6,1), C(8,0) }, clampedKnots(3, 2), 2);
        const NURBSCurve fold = mkCurve({ C(0,0), C(2,2), C(4,0), C(2,-2), C(0,0) }, clampedKnots(5, 4), 4);
        const NURBSCurve line = mkCurve({ C(-1,0), C(5,0) }, clampedKnots(2, 1), 1);
        const NURBSCurve big1 = mkCurve({ C(0,0), C(1e6,0), C(2e6,1e6), C(3e6,0) }, clampedKnots(4, 3), 3);
        const NURBSCurve big2 = mkCurve({ C(0,1e6), C(1e6,1e6), C(2e6,0), C(3e6,1e6) }, clampedKnots(4, 3), 3);

        struct Pair { const NURBSCurve& a; const NURBSCurve& b; const char* name; };
        const Pair pairs[] = {
            { pointCurve, pointCurve, "point-point" },
            { collinear, collinear, "collinear-overlap" },
            { endA, endB, "endpoint-touch" },
            { fold, line, "fold-vs-line" },
            { big1, big2, "scale-1e6" },
        };
        bool allOk = true;
        for (const auto& pr : pairs) {
            SolverOptions o; o.usePreFiltering = true; o.gridSize = 6; o.maxIterations = 80;
            try {
                auto res = unifiedSolve(pr.a, pr.b, o);
                if (res.empty() || !isFinite(res[0].distance) || res[0].iterations <= 0 || res[0].iterations > 1000) {
                    std::printf("[robust] solve 异常: %s empty=%d dist=%.3e iter=%d\n", pr.name, (int)res.empty(), res.empty()?0:res[0].distance, res.empty()?0:res[0].iterations);
                    allOk = false;
                }
            } catch (const std::exception& e) { std::printf("[robust] solve 抛异常: %s: %s\n", pr.name, e.what()); allOk = false; }
            SolverOptions fo; fo.tolerance = 1e-8;
            try { auto roots = findAllIntersections(pr.a, pr.b, fo); (void)roots; }
            catch (...) { std::printf("[robust] find 抛异常: %s\n", pr.name); allOk = false; }
        }
        CHECK("退化对求解不崩溃/有限/迭代有界", allOk);
    }

    // ---- 4) 非法结构：unifiedSolve 入口必须抛 invalid_argument ----
    {
        const std::vector<ControlPoint> goodPts = { C(0,0), C(1,0), C(2,0), C(3,0) };
        const auto goodKnots = clampedKnots(4, 3);
        NURBSCurve good = mkCurve(goodPts, goodKnots, 3);

        // 4a) knot 长度错误
        NURBSCurve badKnot = good;
        badKnot.knots = std::vector<double>(goodKnots.begin(), goodKnots.begin() + 5);
        // 4b) 负权重
        NURBSCurve negW = good;
        negW.controlPoints[1].w = -0.5;
        // 4c) NaN 坐标
        NURBSCurve nanP = good;
        nanP.controlPoints[0].x = std::numeric_limits<double>::quiet_NaN();
        // 4d) 参数域倒置
        NURBSCurve inverted = good;
        std::swap(inverted.uMin, inverted.uMax);
        // 4e) 空控制点
        NURBSCurve empty = good;
        empty.controlPoints.clear();
        // 4f) knot 递减（不合法）
        NURBSCurve decr = good;
        decr.knots[5] = 0.2; decr.knots[6] = 0.1;

        SolverOptions o; o.usePreFiltering = true;
        auto expectThrow = [&](const NURBSCurve& c, const char* nm) {
            try { unifiedSolve(c, good, o); std::printf("[robust] 未抛异常: %s\n", nm); return false; }
            catch (const std::invalid_argument&) { return true; }
            catch (const std::exception& e) { std::printf("[robust] 抛了非 invalid_argument: %s: %s\n", nm, e.what()); return false; }
        };
        bool allThrow = true;
        allThrow &= expectThrow(badKnot, "knot-长度错");
        allThrow &= expectThrow(negW, "负权重");
        allThrow &= expectThrow(nanP, "NaN 坐标");
        allThrow &= expectThrow(inverted, "参数域倒置");
        allThrow &= expectThrow(empty, "空控制点");
        allThrow &= expectThrow(decr, "knot 递减");
        CHECK("非法结构入口抛 invalid_argument", allThrow);

        // 合法曲线不得抛
        bool okGood = true;
        try { unifiedSolve(good, good, o); } catch (...) { okGood = false; }
        CHECK("合法曲线不抛", okGood);
    }

    std::printf("[robust] ==== 鲁棒性汇总 ==== pass=%d fail=%d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
