// ============================================================================
// unified_solver.cpp — 统一数值迭代求解引擎实现（C++）
// 逐函数对照 src/solver/unified_solver.js 移植，运算顺序与类型保持一致。
// ============================================================================
#include "nurbs/unified_solver.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <limits>

namespace nurbs {

namespace {
// ---- 阈值体系（全局可变，对应 JS 模块级 let）----
double GEOMETRY_TOLERANCE = 1e-8;
const double kPi = 3.141592653589793;
double ANGULAR_TOLERANCE = 12 * kPi / 180;
double BOUNDARY_INTERSECTION_TOLERANCE = 1e-4;
double PARAM_BOUNDARY_TOL = 1e-3;
double CONVERGENCE_THRESHOLD = 1e-12;
const double DAMPING_INITIAL = 0.01;
const double DAMPING_FACTOR = 10.0;
const double MAX_DAMPING = 1e10;
const double MIN_DAMPING = 1e-14;

double vectorNorm(double x, double y) { return std::sqrt(x * x + y * y); }

double distanceP(const Point2& p1, const Point2& p2) {
    const double dx = p1.x - p2.x, dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

double dotProduct(const Point2& v1, const Point2& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

// 2x2 矩阵求逆；奇异返回 false（out 不写）
bool invert2x2(double m00, double m01, double m10, double m11,
               double& i00, double& i01, double& i10, double& i11) {
    const double det = m00 * m11 - m01 * m10;
    if (std::abs(det) < kEpsilon) return false;
    const double invDet = 1.0 / det;
    i00 = m11 * invDet; i01 = -m01 * invDet;
    i10 = -m10 * invDet; i11 = m00 * invDet;
    return true;
}

double estimateConditionNumber(const std::vector<std::vector<double>>& matrix) {
    const double normF = std::sqrt(matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1] +
                                   matrix[1][0] * matrix[1][0] + matrix[1][1] * matrix[1][1]);
    double i00, i01, i10, i11;
    if (!invert2x2(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1], i00, i01, i10, i11))
        return std::numeric_limits<double>::infinity();
    const double normInvF = std::sqrt(i00 * i00 + i01 * i01 + i10 * i10 + i11 * i11);
    return normF * normInvF;
}

// 与 JS Number.prototype.toFixed(6) 对齐（round half away from zero），用于去重键
std::string toFixed6(double x) {
    if (!std::isfinite(x)) return x > 0 ? "Infinity" : (x < 0 ? "-Infinity" : "NaN");
    bool neg = x < 0; if (neg) x = -x;
    const double scaled = x * 1e6;
    const double r = std::round(scaled);
    const long long intPart = static_cast<long long>(r / 1e6);
    const long long frac = static_cast<long long>(r) - intPart * 1000000;
    char buf[96];
    std::snprintf(buf, sizeof buf, "%lld.%06lld", intPart, frac);
    std::string s = buf;
    if (neg) s = "-" + s;
    return s;
}
} // anonymous namespace

// ============================================================================
// 阈值管理
// ============================================================================
ThresholdDefaults thresholdDefaults() {
    ThresholdDefaults t;
    t.GEOMETRY_TOLERANCE = GEOMETRY_TOLERANCE;
    t.ANGULAR_TOLERANCE = ANGULAR_TOLERANCE;
    t.BOUNDARY_INTERSECTION_TOLERANCE = BOUNDARY_INTERSECTION_TOLERANCE;
    t.PARAM_BOUNDARY_TOL = PARAM_BOUNDARY_TOL;
    t.CONVERGENCE_THRESHOLD = CONVERGENCE_THRESHOLD;
    return t;
}

void setThresholds(std::optional<double> g, std::optional<double> ang,
                   std::optional<double> b, std::optional<double> p,
                   std::optional<double> c) {
    if (g) GEOMETRY_TOLERANCE = *g;
    if (ang) ANGULAR_TOLERANCE = *ang;
    if (b) BOUNDARY_INTERSECTION_TOLERANCE = *b;
    if (p) PARAM_BOUNDARY_TOL = *p;
    if (c) CONVERGENCE_THRESHOLD = *c;
}

// ============================================================================
// 目标函数
// ============================================================================
ObjectiveResult objectiveFunction(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                  double u, double v) {
    // 合并求值：一次扫描同时出 点+一阶导。evaluatePointAndDerivative 与
    // evaluatePoint/evaluateFirstDerivative 分别调用逐位一致（公式相同，
    // 仅省共享的 de Boor 扫描），故本函数输出与原版完全相同（等价保持），
    // 但 LM 每步求值从 8 次 de Boor 扫描降为 4 次。
    const auto p1 = evaluatePointAndDerivative(curve1, u);
    const auto p2 = evaluatePointAndDerivative(curve2, v);
    const Point2& point1 = p1.first;  const Point2& deriv1 = p1.second;
    const Point2& point2 = p2.first;  const Point2& deriv2 = p2.second;

    const double rx = point1.x - point2.x;
    const double ry = point1.y - point2.y;
    const double f = rx * rx + ry * ry;

    const double j11 = deriv1.x, j12 = -deriv2.x;
    const double j21 = deriv1.y, j22 = -deriv2.y;

    const double h00 = j11 * j11 + j21 * j21;
    const double h01 = j11 * j12 + j21 * j22;
    const double h11 = j12 * j12 + j22 * j22;

    const double g0 = j11 * rx + j21 * ry;
    const double g1 = j12 * rx + j22 * ry;

    ObjectiveResult o;
    o.f = f;
    o.gradient = { 2 * g0, 2 * g1 };
    o.hessian = { {h00, h01}, {h01, h11} };
    o.point1 = point1; o.point2 = point2;
    o.deriv1 = deriv1; o.deriv2 = deriv2;
    o.distance = std::sqrt(f);
    return o;
}

// ============================================================================
// Levenberg-Marquardt 迭代求解器
// ============================================================================
SolveResult levenbergMarquardt(const NURBSCurve& curve1, const NURBSCurve& curve2,
                               double u0, double v0, const SolverOptions& options) {
    const int maxIter = options.maxIterations > 0 ? options.maxIterations : kMaxIterations;
    const double tol = options.tolerance > 0 ? options.tolerance : CONVERGENCE_THRESHOLD;
    const bool recordHistory = options.recordHistory;

    double u = u0, v = v0;
    double damping = DAMPING_INITIAL;

    std::vector<SolveResult::IterRecord> history;
    int iterCount = 0;
    bool converged = false;
    double conditionNum = 0;
    ObjectiveResult obj = objectiveFunction(curve1, curve2, u, v);
    double prevF = std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < maxIter; ++iter) {
        if (recordHistory) {
            history.push_back({ iter, u, v, obj.distance, obj.f,
                                vectorNorm(obj.gradient[0], obj.gradient[1]) });
        }
        ++iterCount;

        const auto& hessian = obj.hessian;
        const double grad0 = obj.gradient[0], grad1 = obj.gradient[1];

        const double dh00 = hessian[0][0] + damping;
        const double dh01 = hessian[0][1];
        const double dh10 = hessian[1][0];
        const double dh11 = hessian[1][1] + damping;

        double i00, i01, i10, i11;
        if (!invert2x2(dh00, dh01, dh10, dh11, i00, i01, i10, i11)) {
            damping = std::min(MAX_DAMPING, damping * DAMPING_FACTOR);
            obj = objectiveFunction(curve1, curve2, u, v);
            continue;
        }

        const double du = -(i00 * grad0 + i01 * grad1) * 0.5;
        const double dv = -(i10 * grad0 + i11 * grad1) * 0.5;

        const double un = std::max(curve1.uMin, std::min(curve1.uMax, u + du));
        const double vn = std::max(curve2.uMin, std::min(curve2.uMax, v + dv));
        ObjectiveResult candidate = objectiveFunction(curve1, curve2, un, vn);

        if (candidate.f <= obj.f + 1e-30) {
            u = un; v = vn; obj = std::move(candidate);
            damping = std::max(MIN_DAMPING, damping / DAMPING_FACTOR);
        } else {
            damping = std::min(MAX_DAMPING, damping * DAMPING_FACTOR);
            continue;
        }

        const double gradNorm = vectorNorm(grad0, grad1);
        if (obj.distance < tol) { converged = true; break; }
        if (gradNorm < 1e-8) { converged = true; break; }
        if (std::abs(du) < 1e-12 && std::abs(dv) < 1e-12) { converged = true; break; }
        if (prevF != std::numeric_limits<double>::infinity() && prevF > 0 &&
            obj.f >= prevF * (1 - 1e-10)) { converged = true; break; }
        prevF = obj.f;
    }

    if (converged) conditionNum = estimateConditionNumber(obj.hessian);

    const double dist = obj.distance;

    GeometryState state = GeometryState::UNDETERMINED;
    std::string message;

    const auto atParamBoundary = [&](double uu, const NURBSCurve& c) {
        return std::abs(uu - c.uMin) < PARAM_BOUNDARY_TOL ||
               std::abs(uu - c.uMax) < PARAM_BOUNDARY_TOL;
    };

    if (!converged) {
        state = GeometryState::UNDETERMINED;
        message = "算法未在最大迭代次数内收敛";
    } else if (dist < GEOMETRY_TOLERANCE ||
               (dist < BOUNDARY_INTERSECTION_TOLERANCE &&
                (atParamBoundary(u, curve1) || atParamBoundary(v, curve2)))) {
        const auto boundaryTan = [&](const NURBSCurve& c, double param)
                                 -> std::optional<Point2> {
            const auto& pts = c.controlPoints;
            if (pts.size() < 2) return std::nullopt;
            if (std::abs(param - c.uMin) < PARAM_BOUNDARY_TOL)
                return Point2{ pts[1].x - pts[0].x, pts[1].y - pts[0].y };
            if (std::abs(param - c.uMax) < PARAM_BOUNDARY_TOL) {
                const size_t n = pts.size();
                return Point2{ pts[n - 1].x - pts[n - 2].x, pts[n - 1].y - pts[n - 2].y };
            }
            return std::nullopt;
        };
        std::optional<Point2> tan1, tan2;
        const auto b1 = atParamBoundary(u, curve1) ? boundaryTan(curve1, u) : std::nullopt;
        const auto b2 = atParamBoundary(v, curve2) ? boundaryTan(curve2, v) : std::nullopt;
        const Point2 t1 = b1 ? *b1 : obj.deriv1;
        const Point2 t2 = b2 ? *b2 : obj.deriv2;

        const double tan1Norm = vectorNorm(t1.x, t1.y);
        const double tan2Norm = vectorNorm(t2.x, t2.y);

        if (tan1Norm > kEpsilon && tan2Norm > kEpsilon) {
            const double cosAngle = (t1.x * t2.x + t1.y * t2.y) / (tan1Norm * tan2Norm);
            const double angle = std::acos(std::max(-1.0, std::min(1.0, cosAngle)));
            if (angle < ANGULAR_TOLERANCE || std::abs(angle - kPi) < ANGULAR_TOLERANCE) {
                state = GeometryState::TANGENT;
                message = "检测到相切接触";
            } else {
                state = GeometryState::INTERSECTION;
                message = dist < GEOMETRY_TOLERANCE ? "检测到交点" : "检测到边界交点";
            }
        } else {
            state = GeometryState::INTERSECTION;
            message = "检测到交点 (奇异点)";
        }
    } else {
        state = GeometryState::SEPARATED;
        message = "曲线分离，最短距离已求得";
    }

    SolveResult r;
    r.state = state;
    r.u = u; r.v = v;
    r.point1 = obj.point1; r.point2 = obj.point2;
    r.distance = dist;
    r.iterations = iterCount;
    r.residual = obj.f;
    r.conditionNumber = conditionNum;
    r.converged = converged;
    r.message = message;
    r.iterationHistory = std::move(history);
    return r;
}

// ============================================================================
// 多初值策略
// ============================================================================
std::vector<Point2> generateGridSamples(int n1, int n2) {
    std::vector<Point2> samples;
    samples.reserve(n1 * n2);
    for (int i = 0; i < n1; ++i) {
        const double u = static_cast<double>(i) / (n1 - 1);
        for (int j = 0; j < n2; ++j) {
            const double v = static_cast<double>(j) / (n2 - 1);
            samples.push_back({ u, v });
        }
    }
    return samples;
}

static double estimateSegmentRadius(const NURBSCurve& curve, double uStart, double uEnd) {
    const auto& controlPoints = curve.controlPoints;
    const auto& knots = curve.knots;
    const int degree = curve.degree;
    const int spanStart = findSpanIndex(degree, uStart, knots);
    const int spanEnd = findSpanIndex(degree, uEnd, knots);

    std::vector<const ControlPoint*> relevant;
    const int lo = std::max(0, spanStart - degree);
    const int hi = std::min(static_cast<int>(controlPoints.size()) - 1, spanEnd + degree);
    for (int i = lo; i <= hi; ++i) relevant.push_back(&controlPoints[i]);
    if (relevant.empty()) return 1.0;

    double cx = 0, cy = 0;
    for (const auto* pt : relevant) { cx += pt->x; cy += pt->y; }
    cx /= relevant.size(); cy /= relevant.size();

    double maxXDist = 0, maxYDist = 0;
    for (const auto* pt : relevant) {
        maxXDist = std::max(maxXDist, std::abs(pt->x - cx));
        maxYDist = std::max(maxYDist, std::abs(pt->y - cy));
    }
    return std::max(maxXDist, maxYDist) * 1.5;
}

std::vector<Region> findCandidateRegions(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                         int subdivisions) {
    struct Seg { double start, end; int index; };
    std::vector<Seg> uSegments, vSegments;
    for (int i = 0; i < subdivisions; ++i) {
        uSegments.push_back({ static_cast<double>(i) / subdivisions,
                              static_cast<double>(i + 1) / subdivisions, i });
        vSegments.push_back({ static_cast<double>(i) / subdivisions,
                              static_cast<double>(i + 1) / subdivisions, i });
    }
    struct BoxSeg { double minX, minY, maxX, maxY; Seg segment; };
    std::vector<BoxSeg> uBoxes, vBoxes;
    for (const auto& seg : uSegments) {
        const double midU = (seg.start + seg.end) / 2;
        const Point2 pt = evaluatePoint(curve1, midU);
        const double radius = estimateSegmentRadius(curve1, seg.start, seg.end);
        uBoxes.push_back({ pt.x - radius, pt.y - radius, pt.x + radius, pt.y + radius, seg });
    }
    for (const auto& seg : vSegments) {
        const double midV = (seg.start + seg.end) / 2;
        const Point2 pt = evaluatePoint(curve2, midV);
        const double radius = estimateSegmentRadius(curve2, seg.start, seg.end);
        vBoxes.push_back({ pt.x - radius, pt.y - radius, pt.x + radius, pt.y + radius, seg });
    }
    std::vector<Region> regions;
    for (size_t i = 0; i < uBoxes.size(); ++i)
        for (size_t j = 0; j < vBoxes.size(); ++j) {
            const AABB b1{ uBoxes[i].minX, uBoxes[i].minY, uBoxes[i].maxX, uBoxes[i].maxY };
            const AABB b2{ vBoxes[j].minX, vBoxes[j].minY, vBoxes[j].maxX, vBoxes[j].maxY };
            if (aabbIntersect(b1, b2))
                regions.push_back({ uSegments[i].start, uSegments[i].end,
                                    vSegments[j].start, vSegments[j].end });
        }
    return regions;
}

// ============================================================================
// 统一求解器
// ============================================================================
SolveResult globalOptimization(const NURBSCurve& curve1, const NURBSCurve& curve2,
                               const SolverOptions& options) {
    const std::vector<Point2> denseGrid = generateGridSamples(10, 10);
    SolveResult* bestConverged = nullptr;
    double minConvDistance = std::numeric_limits<double>::infinity();
    SolveResult* bestAny = nullptr;
    double minAnyDistance = std::numeric_limits<double>::infinity();

    for (const auto& sample : denseGrid) {
        SolverOptions o = options;
        o.recordHistory = false;
        SolveResult result = levenbergMarquardt(curve1, curve2, sample.x, sample.y, o);
        if (result.distance < minAnyDistance) {
            minAnyDistance = result.distance;
            bestAny = new SolveResult(std::move(result));
        }
        if (result.converged && result.distance < minConvDistance) {
            minConvDistance = result.distance;
            bestConverged = new SolveResult(std::move(result));
        }
    }
    if (bestConverged || bestAny) {
        SolveResult pick = bestConverged ? *bestConverged : *bestAny;
        delete bestConverged; delete bestAny;
        return pick;
    }
    delete bestConverged; delete bestAny;
    SolveResult fallback;
    fallback.state = GeometryState::UNDETERMINED;
    fallback.u = 0.5; fallback.v = 0.5;
    fallback.point1 = evaluatePoint(curve1, 0.5);
    fallback.point2 = evaluatePoint(curve2, 0.5);
    fallback.distance = std::numeric_limits<double>::infinity();
    fallback.iterations = 0;
    fallback.residual = std::numeric_limits<double>::infinity();
    fallback.converged = false;
    fallback.message = "全局优化未能找到有效解";
    return fallback;
}

std::vector<SolveResult> unifiedSolve(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                      const SolverOptions& options) {
    const bool usePreFiltering = options.usePreFiltering;
    const int gridSize = options.gridSize > 0 ? options.gridSize : 5;
    const int maxIterations = options.maxIterations > 0 ? options.maxIterations : kMaxIterations;
    const double tolerance = options.tolerance > 0 ? options.tolerance : GEOMETRY_TOLERANCE;

    validateCurve(curve1);
    validateCurve(curve2);

    if (usePreFiltering) {
        const AABB box1 = computeAABB(curve1);
        const AABB box2 = computeAABB(curve2);
        if (!aabbIntersect(box1, box2)) {
            SolverOptions o;
            o.maxIterations = maxIterations;
            o.tolerance = tolerance;
            std::vector<SolveResult> out;
            out.push_back(globalOptimization(curve1, curve2, o));
            return out;
        }
    }

    struct Guess { double u, v; };
    std::vector<Guess> initialGuesses;

    if (usePreFiltering) {
        const auto regions = findCandidateRegions(curve1, curve2, gridSize);
        for (const auto& region : regions)
            initialGuesses.push_back({ (region.uStart + region.uEnd) / 2,
                                       (region.vStart + region.vEnd) / 2 });
    } else {
        const auto samples = generateGridSamples(gridSize, gridSize);
        for (const auto& s : samples) initialGuesses.push_back({ s.x, s.y });
    }

    const double cornerU[4] = { curve1.uMin, curve1.uMin, curve1.uMax, curve1.uMax };
    const double cornerV[4] = { curve2.uMin, curve2.uMax, curve2.uMin, curve2.uMax };
    for (int k = 0; k < 4; ++k) {
        const double cu = cornerU[k], cv = cornerV[k];
        bool dup = false;
        for (const auto& g : initialGuesses)
            if (std::abs(g.u - cu) < 1e-9 && std::abs(g.v - cv) < 1e-9) { dup = true; break; }
        if (!dup) initialGuesses.push_back({ cu, cv });
    }

    std::vector<SolveResult> results;
    std::vector<std::string> foundKeys;

    for (const auto& guess : initialGuesses) {
        SolverOptions o;
        o.maxIterations = maxIterations / 3;   // floor 整数除法，与 Math.floor 一致
        o.tolerance = tolerance;
        o.recordHistory = false;
        SolveResult result = levenbergMarquardt(curve1, curve2, guess.u, guess.v, o);
        if (result.converged) {
            const std::string key = toFixed6(result.u) + "," + toFixed6(result.v);
            bool seen = false;
            for (const auto& k : foundKeys) if (k == key) { seen = true; break; }
            if (!seen) {
                foundKeys.push_back(key);
                results.push_back(std::move(result));
            }
        }
    }

    // JS Array.prototype.sort 是稳定排序；等值项顺序影响贪心合并，故用 stable_sort
    std::stable_sort(results.begin(), results.end(),
              [](const SolveResult& a, const SolveResult& b) { return a.distance < b.distance; });

    if (results.empty()) {
        SolverOptions o;
        o.maxIterations = maxIterations;
        o.tolerance = tolerance;
        results.push_back(globalOptimization(curve1, curve2, o));
    }
    return results;
}

// ============================================================================
// 全交点求解：递归细分隔离 + 逐区间 LM 抛光 + 容差聚类
// ============================================================================
AABB pieceBox(const NURBSCurve& curve, double a, double b) {
    const int K = 16;
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    for (int i = 0; i <= K; ++i) {
        const double u = a + (b - a) * i / K;
        const Point2 p = evaluatePoint(curve, u);
        minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
    }
    return { minX, minY, maxX, maxY };
}

std::vector<Root> findAllIntersections(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                       const SolverOptions& options) {
    const double tolerance = options.tolerance > 0 ? options.tolerance : GEOMETRY_TOLERANCE;
    const int maxDepth = options.maxDepth > 0 ? options.maxDepth : 22;
    const double domainW = std::max(curve1.uMax - curve1.uMin, curve2.uMax - curve2.uMin);
    const double domainW0 = domainW != 0 ? domainW : 1;
    const double paramTol = options.paramTol > 0 ? options.paramTol : domainW0 * 1e-3;

    struct Cell { double u0, u1, v0, v1; int depth; };
    std::vector<Cell> work;
    work.push_back({ curve1.uMin, curve1.uMax, curve2.uMin, curve2.uMax, 0 });
    std::vector<Cell> seeds;
    int guard = 0;

    while (!work.empty()) {
        if (++guard > 1000000) break;
        Cell cell = work.back();
        work.pop_back();
        const AABB b1 = pieceBox(curve1, cell.u0, cell.u1);
        const AABB b2 = pieceBox(curve2, cell.v0, cell.v1);
        if (!aabbIntersect(b1, b2)) continue;

        const double wU = cell.u1 - cell.u0, wV = cell.v1 - cell.v0;
        if (cell.depth >= maxDepth || (wU < paramTol && wV < paramTol)) {
            seeds.push_back(cell);
            continue;
        }
        if (wU >= wV) {
            const double um = (cell.u0 + cell.u1) / 2;
            work.push_back({ cell.u0, um, cell.v0, cell.v1, cell.depth + 1 });
            work.push_back({ um, cell.u1, cell.v0, cell.v1, cell.depth + 1 });
        } else {
            const double vm = (cell.v0 + cell.v1) / 2;
            work.push_back({ cell.u0, cell.u1, cell.v0, vm, cell.depth + 1 });
            work.push_back({ cell.u0, cell.u1, vm, cell.v1, cell.depth + 1 });
        }
    }

    std::vector<Root> roots;
    for (const auto& cell : seeds) {
        const double u = (cell.u0 + cell.u1) / 2, v = (cell.v0 + cell.v1) / 2;
        SolverOptions o;
        o.maxIterations = 100;
        o.tolerance = tolerance;
        o.recordHistory = false;
        SolveResult r = levenbergMarquardt(curve1, curve2, u, v, o);
        if (r.converged && r.distance < tolerance)
            roots.push_back({ r.u, r.v, r.distance, r.point1, r.point2, r.iterations });
    }

    std::stable_sort(roots.begin(), roots.end(),
              [](const Root& a, const Root& b) { return a.u < b.u; });
    std::vector<Root> out;
    for (const auto& rr : roots) {
        bool merged = false;
        for (auto& o : out) {
            if (std::abs(o.u - rr.u) < paramTol * 6 && std::abs(o.v - rr.v) < paramTol * 6) {
                if (rr.distance < o.distance) {
                    o.u = rr.u; o.v = rr.v; o.distance = rr.distance;
                    o.point1 = rr.point1; o.point2 = rr.point2;
                }
                merged = true;
                break;
            }
        }
        if (!merged) out.push_back(rr);
    }
    std::stable_sort(out.begin(), out.end(),
              [](const Root& a, const Root& b) { return a.u < b.u; });
    return out;
}

} // namespace nurbs
