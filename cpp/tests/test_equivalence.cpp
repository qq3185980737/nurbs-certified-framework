// ============================================================================
// test_equivalence.cpp — JS↔C++ 移植等价性验证
// ----------------------------------------------------------------------------
// 读取 gen_fixture.mjs 生成的 fixture.txt（JS 用同一种子 20261012 与同批用例生成），
// 用本 C++ 库重算 evaluatePoint / unifiedSolve / findAllIntersections，逐位对比。
// 同时校验 mulberry32 双实现一致性。全部精确一致则退出码 0，否则 1。
//
// 运行：test_equivalence <path-to-fixture.txt>
// ============================================================================
#include "nurbs/nurbs_core.hpp"
#include "nurbs/unified_solver.hpp"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <cmath>
#include <array>
#include <algorithm>
#include <chrono>

using Clock = std::chrono::steady_clock;
static double msSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

using namespace nurbs;

// ---- mulberry32（C++ 版，与 JS 逐位一致）----
struct Mulberry32 {
    std::uint32_t a;
    explicit Mulberry32(std::uint32_t seed) : a(seed) {}
    double next() {
        a += 0x6D2B79F5u;
        std::uint32_t t = (a ^ (a >> 15)) * (1u | a);
        t = (t + (t ^ (t >> 7)) * (61u | t)) ^ t;
        const std::uint32_t x = t ^ (t >> 14);
        return static_cast<double>(x) / 4294967296.0;
    }
};

static int stateToInt(GeometryState s) {
    switch (s) {
        case GeometryState::INTERSECTION: return 0;
        case GeometryState::TANGENT: return 1;
        case GeometryState::SEPARATED: return 2;
        case GeometryState::COINCIDENT: return 3;
        case GeometryState::UNDETERMINED: return 4;
    }
    return -1;
}

struct FixtureCurve { int ncp; int deg; double uMin, uMax; std::vector<ControlPoint> pts; std::vector<double> knots; };
struct EvalRec { int cid; double u, x, y; };
struct SolveRec { int c1, c2; int st; double u, v, dist; int iter; int conv; int nres; };
struct FindRec { int c1, c2; std::vector<std::array<double,3>> roots; };

static NURBSCurve buildCurve(const FixtureCurve& fc) {
    NURBSCurve c;
    c.degree = fc.deg;
    c.uMin = fc.uMin; c.uMax = fc.uMax;
    c.knots = fc.knots;
    c.controlPoints = fc.pts;
    return c;
}

int main(int argc, char** argv) {
    const char* path = (argc > 1) ? argv[1] : "tests/fixture.txt";
    std::ifstream in(path);
    if (!in) { std::fprintf(stderr, "[eq] 无法打开夹具: %s\n", path); return 2; }

    std::string tag;
    in >> tag; // NURBS_FIXTURE

    // RND
    in >> tag; // RND
    int nRnd = 0; in >> nRnd;
    std::vector<double> rndExp(nRnd);
    for (int i = 0; i < nRnd; ++i) in >> rndExp[i];

    // CURVES
    in >> tag; // NC
    int nCurves = 0; in >> nCurves;
    std::vector<FixtureCurve> curves(nCurves);
    for (int i = 0; i < nCurves; ++i) {
        in >> tag; // CURVE
        int id; in >> id;
        FixtureCurve& fc = curves[id];
        in >> fc.ncp >> fc.deg >> fc.uMin >> fc.uMax;
        int nk; in >> nk;
        fc.pts.resize(fc.ncp);
        for (int k = 0; k < fc.ncp; ++k) in >> fc.pts[k].x >> fc.pts[k].y >> fc.pts[k].w;
        in >> tag; // K
        fc.knots.resize(nk);
        for (int k = 0; k < nk; ++k) in >> fc.knots[k];
    }

    // EVALS
    in >> tag; // NE
    int nEval = 0; in >> nEval;
    std::vector<EvalRec> evals(nEval);
    for (int i = 0; i < nEval; ++i) {
        in >> tag; // EVAL
        in >> evals[i].cid >> evals[i].u >> evals[i].x >> evals[i].y;
    }

    // SOLVES
    in >> tag; // NS
    int nSolve = 0; in >> nSolve;
    std::vector<SolveRec> solves(nSolve);
    for (int i = 0; i < nSolve; ++i) {
        in >> tag; // SOLVE
        in >> solves[i].c1 >> solves[i].c2 >> solves[i].st >> solves[i].u >> solves[i].v
           >> solves[i].dist >> solves[i].iter >> solves[i].conv >> solves[i].nres;
    }

    // FINDS
    in >> tag; // NF
    int nFind = 0; in >> nFind;
    std::vector<FindRec> finds(nFind);
    for (int i = 0; i < nFind; ++i) {
        in >> tag; // FIND
        in >> finds[i].c1 >> finds[i].c2;
        int nr = 0; in >> nr;
        finds[i].roots.resize(nr);
        for (int k = 0; k < nr; ++k) in >> finds[i].roots[k][0] >> finds[i].roots[k][1] >> finds[i].roots[k][2];
    }

    int pass = 0, fail = 0;

    // 1) mulberry32 一致性
    Mulberry32 rng(20261012u);
    int rndBad = 0;
    for (int i = 0; i < nRnd; ++i) if (rng.next() != rndExp[i]) rndBad++;
    if (rndBad == 0) { pass++; } else { fail++; std::printf("[eq] RND 不一致: %d/%d\n", rndBad, nRnd); }

    // 2) evaluatePoint 逐位
    int evalExact = 0; double evalMaxDiff = 0;
    auto tEval = Clock::now();
    for (const auto& e : evals) {
        NURBSCurve c = buildCurve(curves[e.cid]);
        Point2 p;
        try { p = evaluatePoint(c, e.u); } catch (...) { continue; }
        const double dx = std::abs(p.x - e.x), dy = std::abs(p.y - e.y);
        evalMaxDiff = std::max(evalMaxDiff, std::max(dx, dy));
        if (p.x == e.x && p.y == e.y) evalExact++;
    }
    if (evalExact == nEval) { pass++; } else { fail++; std::printf("[eq] evaluatePoint 精确 %d/%d, maxDiff=%.3e\n", evalExact, nEval, evalMaxDiff); }
    const double tEvalMs = msSince(tEval);

    // 3) unifiedSolve 逐位
    int solveExact = 0, solveFieldFail = 0;
    auto tSolve = Clock::now();
    for (const auto& s : solves) {
        NURBSCurve c1 = buildCurve(curves[s.c1]);
        NURBSCurve c2 = buildCurve(curves[s.c2]);
        SolverOptions o; o.usePreFiltering = true; o.gridSize = 6; o.maxIterations = 80;
        std::vector<SolveResult> res = unifiedSolve(c1, c2, o);
        if (res.empty()) { solveFieldFail++; continue; }
        SolveResult best = res[0];
        for (const auto& r : res) if (r.distance < best.distance) best = r;
        const bool ok = (stateToInt(best.state) == s.st) && (best.u == s.u) && (best.v == s.v) &&
                        (best.distance == s.dist) && (best.iterations == s.iter) &&
                        (static_cast<int>(best.converged) == s.conv) &&
                        (static_cast<int>(res.size()) == s.nres);
        if (ok) solveExact++; else { solveFieldFail++; if (solveFieldFail <= 5) {
            std::printf("[eq] SOLVE %d: JS(st=%d u=%.17g v=%.17g d=%.17g it=%d conv=%d n=%d) C++(st=%d u=%.17g v=%.17g d=%.17g it=%d conv=%d n=%d)\n",
                s.c1, s.st, s.u, s.v, s.dist, s.iter, s.conv, s.nres,
                stateToInt(best.state), best.u, best.v, best.distance, best.iterations, (int)best.converged, (int)res.size());
        } }
    }
    if (solveFieldFail == 0) { pass++; } else { fail++; std::printf("[eq] unifiedSolve 精确 %d/%d\n", solveExact, nSolve); }
    const double tSolveMs = msSince(tSolve);

    // 4) findAllIntersections 逐位
    int findExact = 0, findFieldFail = 0;
    auto tFind = Clock::now();
    for (const auto& f : finds) {
        NURBSCurve c1 = buildCurve(curves[f.c1]);
        NURBSCurve c2 = buildCurve(curves[f.c2]);
        SolverOptions o; o.tolerance = 1e-8;
        std::vector<Root> roots = findAllIntersections(c1, c2, o);
        bool ok = (roots.size() == f.roots.size());
        if (ok) for (size_t k = 0; k < roots.size(); ++k)
            if (roots[k].u != f.roots[k][0] || roots[k].v != f.roots[k][1] || roots[k].distance != f.roots[k][2]) { ok = false; break; }
        if (ok) findExact++; else { findFieldFail++; if (findFieldFail <= 5) {
            std::printf("[eq] FIND %d/%d: JS n=%d C++ n=%d\n", f.c1, f.c2, (int)f.roots.size(), (int)roots.size());
        } }
    }
    if (findFieldFail == 0) { pass++; } else { fail++; std::printf("[eq] findAllIntersections 精确 %d/%d\n", findExact, nFind); }
    const double tFindMs = msSince(tFind);

    std::printf("[eq] ==== 等价性汇总 ====\n");
    std::printf("[eq] RND=%d  eval=%d/%d  solve=%d/%d  find=%d/%d\n",
                (rndBad==0), evalExact, nEval, solveExact, nSolve, findExact, nFind);
    std::printf("[eq] 耗时: eval=%d 次=%.1fms  solve=%d 组=%.1fms  find=%d 组=%.1fms\n",
                nEval, tEvalMs, nSolve, tSolveMs, nFind, tFindMs);
    std::printf("[eq] pass=%d fail=%d\n", pass, fail);
    return fail == 0 ? 0 : 1;
}
