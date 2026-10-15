// ============================================================================
// nurbs_c_api.cpp — 纯 C 接口实现（extern "C" 包装 C++ 核心）
// ============================================================================
#include "nurbs/nurbs_c_api.h"
#include "nurbs/nurbs_core.hpp"
#include "nurbs/unified_solver.hpp"
#include "nurbs/accelerated_solver.hpp"
#include <new>
#include <cstdlib>
#include <stdexcept>

using namespace nurbs;

// 在实现侧补全不透明句柄：持有 C++ 核心曲线 / BVH 索引
struct NurbsCurve {
    NURBSCurve impl;
};
struct NurbsBvhIndex {
    BVHIndex impl;
};

namespace {

int mapError() {
    try { throw; }
    catch (const std::invalid_argument&) { return NURBS_ERR_INVALID; }
    catch (const std::exception&) { return NURBS_ERR_RUNTIME; }
    catch (...) { return NURBS_ERR_RUNTIME; }
}

int stateToInt(GeometryState s) {
    switch (s) {
        case GeometryState::INTERSECTION: return NURBS_STATE_INTERSECTION;
        case GeometryState::TANGENT: return NURBS_STATE_TANGENT;
        case GeometryState::SEPARATED: return NURBS_STATE_SEPARATED;
        case GeometryState::COINCIDENT: return NURBS_STATE_COINCIDENT;
        default: return NURBS_STATE_UNDETERMINED;
    }
}

// 把 vector<SolveResult> 拷贝为 malloc 数组（供统一求解 / BVH 求解复用）
int copyResults(const std::vector<SolveResult>& res, NurbsSolveResult** out_results) {
    NurbsSolveResult* arr = static_cast<NurbsSolveResult*>(
        std::malloc(res.empty() ? 1 : res.size() * sizeof(NurbsSolveResult)));
    if (!arr) return NURBS_ERR_RUNTIME;
    for (size_t i = 0; i < res.size(); ++i) {
        arr[i].state = stateToInt(res[i].state);
        arr[i].u = res[i].u; arr[i].v = res[i].v;
        arr[i].point1_x = res[i].point1.x; arr[i].point1_y = res[i].point1.y;
        arr[i].point2_x = res[i].point2.x; arr[i].point2_y = res[i].point2.y;
        arr[i].distance = res[i].distance;
        arr[i].iterations = res[i].iterations;
        arr[i].converged = res[i].converged ? 1 : 0;
    }
    *out_results = arr;
    return NURBS_OK;
}

NURBSCurve toCurve(const NurbsControlPoint* ctrl, int n, const double* knots, int degree) {
    NURBSCurve c;
    c.degree = degree;
    c.controlPoints.resize(n);
    for (int i = 0; i < n; ++i) {
        c.controlPoints[i].x = ctrl[i].x;
        c.controlPoints[i].y = ctrl[i].y;
        c.controlPoints[i].w = ctrl[i].w;
    }
    c.knots.assign(knots, knots + degree + n + 1);
    c.uMin = c.knots[degree];
    c.uMax = c.knots[c.knots.size() - 1 - degree];
    return c;
}

} // anonymous namespace

extern "C" {

int nurbs_curve_create(const NurbsControlPoint* ctrl, int nctrl,
                       const double* knots, int nknots, int degree,
                       NurbsCurve** out) {
    if (!ctrl || nctrl <= 0 || !knots || nknots != nctrl + degree + 1 || degree < 0 || !out)
        return NURBS_ERR_INVALID;
    try {
        NURBSCurve c = toCurve(ctrl, nctrl, knots, degree);
        validateCurve(c); // 创建即校验
        NurbsCurve* h = new NurbsCurve;
        h->impl = std::move(c);
        *out = h;
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

void nurbs_curve_destroy(NurbsCurve* c) {
    delete c;
}

int nurbs_evaluate_point(const NurbsCurve* c, double u, double* ox, double* oy) {
    if (!c || !ox || !oy) return NURBS_ERR_INVALID;
    try {
        Point2 p = evaluatePoint(c->impl, u);
        *ox = p.x; *oy = p.y;
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

int nurbs_evaluate_first_derivative(const NurbsCurve* c, double u, double* ox, double* oy) {
    if (!c || !ox || !oy) return NURBS_ERR_INVALID;
    try {
        Point2 p = evaluateFirstDerivative(c->impl, u);
        *ox = p.x; *oy = p.y;
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

int nurbs_evaluate_second_derivative(const NurbsCurve* c, double u, double* ox, double* oy) {
    if (!c || !ox || !oy) return NURBS_ERR_INVALID;
    try {
        Point2 p = evaluateSecondDerivative(c->impl, u);
        *ox = p.x; *oy = p.y;
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

static int doUnifiedSolve(const NurbsCurve* a, const NurbsCurve* b,
                          const SolverOptions& opts,
                          NurbsSolveResult** out_results, int* out_n) {
    if (!a || !b || !out_results || !out_n) return NURBS_ERR_INVALID;
    try {
        std::vector<SolveResult> res = unifiedSolve(a->impl, b->impl, opts);
        int rc = copyResults(res, out_results);
        if (rc != NURBS_OK) return rc;
        *out_n = static_cast<int>(res.size());
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

int nurbs_unified_solve(const NurbsCurve* a, const NurbsCurve* b,
                        NurbsSolveResult** out_results, int* out_n) {
    SolverOptions opts;
    opts.usePreFiltering = true; opts.gridSize = 6; opts.maxIterations = 80; opts.tolerance = 0;
    return doUnifiedSolve(a, b, opts, out_results, out_n);
}

int nurbs_unified_solve_ex(const NurbsCurve* a, const NurbsCurve* b,
                           int use_prefiltering, int grid_size, int max_iterations, double tolerance,
                           NurbsSolveResult** out_results, int* out_n) {
    SolverOptions opts;
    opts.usePreFiltering = use_prefiltering != 0;
    opts.gridSize = grid_size > 0 ? grid_size : 0;
    opts.maxIterations = max_iterations > 0 ? max_iterations : 0;
    opts.tolerance = tolerance > 0 ? tolerance : 0;
    return doUnifiedSolve(a, b, opts, out_results, out_n);
}

void nurbs_solve_results_free(NurbsSolveResult* r) {
    std::free(r);
}

int nurbs_find_all_intersections(const NurbsCurve* a, const NurbsCurve* b,
                                 NurbsRoot** out_roots, int* out_n) {
    if (!a || !b || !out_roots || !out_n) return NURBS_ERR_INVALID;
    try {
        SolverOptions opts; opts.tolerance = 1e-8;
        std::vector<Root> roots = findAllIntersections(a->impl, b->impl, opts);
        NurbsRoot* arr = static_cast<NurbsRoot*>(
            std::malloc(roots.empty() ? 1 : roots.size() * sizeof(NurbsRoot)));
        if (!arr) return NURBS_ERR_RUNTIME;
        for (size_t i = 0; i < roots.size(); ++i) {
            arr[i].u = roots[i].u; arr[i].v = roots[i].v; arr[i].distance = roots[i].distance;
            arr[i].point1_x = roots[i].point1.x; arr[i].point1_y = roots[i].point1.y;
            arr[i].point2_x = roots[i].point2.x; arr[i].point2_y = roots[i].point2.y;
            arr[i].iterations = roots[i].iterations;
        }
        *out_roots = arr;
        *out_n = static_cast<int>(roots.size());
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

void nurbs_roots_free(NurbsRoot* r) {
    std::free(r);
}

void nurbs_set_geometry_tolerance(double v)      { setThresholds(v, std::nullopt, std::nullopt, std::nullopt, std::nullopt); }
void nurbs_set_angular_tolerance_deg(double deg) { setThresholds(std::nullopt, deg * 3.141592653589793 / 180.0, std::nullopt, std::nullopt, std::nullopt); }
void nurbs_set_boundary_tolerance(double v)      { setThresholds(std::nullopt, std::nullopt, v, std::nullopt, std::nullopt); }
void nurbs_set_param_boundary_tol(double v)      { setThresholds(std::nullopt, std::nullopt, std::nullopt, v, std::nullopt); }
void nurbs_set_convergence_threshold(double v)   { setThresholds(std::nullopt, std::nullopt, std::nullopt, std::nullopt, v); }

// ---------------- BVH 加速求解 ----------------
int nurbs_bvh_index_create(const NurbsCurve* curve, NurbsBvhIndex** out) {
    if (!curve || !out) return NURBS_ERR_INVALID;
    try {
        NurbsBvhIndex* h = new NurbsBvhIndex;
        h->impl.build(curve->impl, 7);
        *out = h;
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

void nurbs_bvh_index_destroy(NurbsBvhIndex* idx) {
    delete idx;
}

int nurbs_bvh_solve(const NurbsBvhIndex* ia, const NurbsCurve* ca,
                    const NurbsBvhIndex* ib, const NurbsCurve* cb,
                    NurbsSolveResult** out_results, int* out_n) {
    if (!ia || !ca || !ib || !cb || !out_results || !out_n) return NURBS_ERR_INVALID;
    try {
        SolverOptions opts;
        opts.usePreFiltering = true; opts.gridSize = 6; opts.maxIterations = 80; opts.tolerance = 0;
        std::vector<SolveResult> res = solveWithIndex(ia->impl, ca->impl, ib->impl, cb->impl, opts);
        int rc = copyResults(res, out_results);
        if (rc != NURBS_OK) return rc;
        *out_n = static_cast<int>(res.size());
        return NURBS_OK;
    } catch (...) { return mapError(); }
}

} // extern "C"
