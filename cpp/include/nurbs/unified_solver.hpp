// ============================================================================
// unified_solver.hpp — 统一数值迭代求解引擎（C++ 移植）
// 从 src/solver/unified_solver.js 逐函数移植，实现"一套框架自动判别"：
// 相交 / 相切 / 分离 / 重合 / 无法判定，基于 Newton-Raphson + LM 阻尼。
// ============================================================================
#ifndef NURBS_UNIFIED_SOLVER_HPP
#define NURBS_UNIFIED_SOLVER_HPP

#include "nurbs_core.hpp"
#include <string>
#include <vector>
#include <optional>

namespace nurbs {

// ---------- 几何状态枚举 ----------
enum class GeometryState {
    INTERSECTION,   // 相交（距离 ≈ 0）
    TANGENT,        // 相切 / 退化接触
    SEPARATED,      // 分离（存在最短距离）
    COINCIDENT,     // 重合（无穷多交点）
    UNDETERMINED    // 无法判定（算法未收敛）
};

// ---------- 阈值体系（全局可配，对应 JS 模块级 let） ----------
// 与 JS 默认值完全一致
struct ThresholdDefaults {
    double GEOMETRY_TOLERANCE = 1e-8;              // 几何判定容差（距离）
    double ANGULAR_TOLERANCE = 12 * 3.141592653589793 / 180;  // 相切角度容差（弧度）12°
    double BOUNDARY_INTERSECTION_TOLERANCE = 1e-4; // 边界接触容差
    double PARAM_BOUNDARY_TOL = 1e-3;              // 参数域"边界"判定宽度
    double CONVERGENCE_THRESHOLD = 1e-12;          // 收敛判据阈值
};

// 返回当前阈值快照（供恢复）
ThresholdDefaults thresholdDefaults();

// 覆盖阈值：仅覆盖 std::nullopt 之外的非空项（对应 JS setThresholds 对象按字段可选）
void setThresholds(std::optional<double> GEOMETRY_TOLERANCE = std::nullopt,
                   std::optional<double> ANGULAR_TOLERANCE = std::nullopt,
                   std::optional<double> BOUNDARY_INTERSECTION_TOLERANCE = std::nullopt,
                   std::optional<double> PARAM_BOUNDARY_TOL = std::nullopt,
                   std::optional<double> CONVERGENCE_THRESHOLD = std::nullopt);

// ---------- 求解结果 ----------
struct SolveResult {
    GeometryState state = GeometryState::UNDETERMINED;
    double u = 0, v = 0;
    Point2 point1{}, point2{};
    double distance = 0;
    int iterations = 0;
    double residual = 0;
    double conditionNumber = 0;
    bool converged = false;
    std::string message;

    struct IterRecord { int iteration; double u, v, distance, residual, gradientNorm; };
    std::vector<IterRecord> iterationHistory;
};

// ---------- 求解器选项 ----------
// 语义与 JS options 对象一致；用 0 / false 表示"未指定"，按各函数默认值处理。
struct SolverOptions {
    int maxIterations = 0;        // <=0 => 使用默认 MAX_ITERATIONS
    double tolerance = 0;         // <=0 => 由各函数决定默认容差
    bool recordHistory = true;
    bool usePreFiltering = true;
    int gridSize = 0;             // <=0 => 5
    // findAllIntersections 专属
    int maxDepth = 0;             // <=0 => 22
    double paramTol = 0;          // <=0 => domainW*1e-3
    // 加速路径 solveWithIndex 的完备性补漏开关（路线 A）：
    //   true （默认）= 完备模式：对可能相交的曲线对用主路径 unifiedSolve
    //                    （Krawczyk 完备证明）兜底，漏→0，正确性优先（较慢）。
    //   false          = 加速模式：纯 BVH+LM 近完备，漏相切/重根/端点难例，
    //                    吞吐高。供"完备性-性能"权衡与双模式切换。
    bool completeLeakPatch = true;
};

// ---------- 核心求解器 ----------
SolveResult levenbergMarquardt(const NURBSCurve& curve1, const NURBSCurve& curve2,
                               double u0, double v0, const SolverOptions& options = {});

// 统一求解器：返回所有找到的解（可能多个交点），排序按距离升序
std::vector<SolveResult> unifiedSolve(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                      const SolverOptions& options = {});

SolveResult globalOptimization(const NURBSCurve& curve1, const NURBSCurve& curve2,
                               const SolverOptions& options);

// ---------- 辅助函数（对外暴露，对应 JS 导出） ----------
struct ObjectiveResult {
    double f;                // 距离平方
    std::vector<double> gradient; // [df/du, df/dv]（JS 中为 [2g0,2g1]）
    std::vector<std::vector<double>> hessian; // [[h00,h01],[h01,h11]]
    Point2 point1, point2, deriv1, deriv2;
    double distance;
};
ObjectiveResult objectiveFunction(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                  double u, double v);

struct Region { double uStart, uEnd, vStart, vEnd; };
std::vector<Region> findCandidateRegions(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                         int subdivisions = 10);
std::vector<Point2> generateGridSamples(int n1, int n2);

// ---------- 全交点求解 ----------
struct Root { double u, v, distance; Point2 point1, point2; int iterations; };
AABB pieceBox(const NURBSCurve& curve, double a, double b);
std::vector<Root> findAllIntersections(const NURBSCurve& curve1, const NURBSCurve& curve2,
                                       const SolverOptions& options = {});

} // namespace nurbs

#endif // NURBS_UNIFIED_SOLVER_HPP
