// ============================================================================
// nurbs_core.hpp — NURBS 数学核心库（C++ 移植，独立于浏览器/Node）
// ----------------------------------------------------------------------------
// 从 src/nurbs/math_core.js 逐函数移植，语义与数值逐位对齐。
// 算法遵循 The NURBS Book (Piegl & Tiller)：de Boor 求值、解析导数、凸包 AABB。
// 本库无第三方依赖，仅用标准库（vector/cmath/string）。
// ============================================================================
#ifndef NURBS_CORE_HPP
#define NURBS_CORE_HPP

#include <vector>
#include <string>
#include <cmath>

namespace nurbs {

// ---------- 常量（与 JS 一致） ----------
constexpr double kEpsilon = 1e-10;          // 数值计算精度容差
constexpr int    kMaxIterations = 100;      // 默认最大迭代次数
constexpr double kKnotTolerance = 1e-10;    // Knot 向量容差

// ---------- 数据结构 ----------
struct ControlPoint { double x = 0.0; double y = 0.0; double w = 1.0; };
struct Point2       { double x = 0.0; double y = 0.0; };
struct AABB         { double minX = 0.0, minY = 0.0, maxX = 0.0, maxY = 0.0; };

struct NURBSCurve {
    std::vector<ControlPoint> controlPoints;
    std::vector<double> knots;
    int degree = 0;
    double uMin = 0.0, uMax = 1.0;
};

// ---------- 浮点工具 ----------
inline bool floatEqual(double a, double b, double tol = kEpsilon) {
    return std::abs(a - b) < tol;
}

// ---------- 合法性校验（不合法抛 std::runtime_error） ----------
void validateKnotVector(const std::vector<double>& knots, int degree, int nControlPoints);
void validateControlPoints(const std::vector<ControlPoint>& controlPoints);
void validateCurve(const NURBSCurve& curve);

// ---------- Knot 向量生成 ----------
std::vector<double> generateUniformKnotVector(int n, int p);
std::vector<double> generateOpenUniformKnotVector(int n, int p);
std::vector<double> generatePiecewiseBezierKnotVector(int n, int p);
std::vector<double> generateKnotVector(int degree, int numControlPoints,
                                       const std::string& knotType = "clamped");

// ---------- 曲线创建 ----------
NURBSCurve createCurve(const std::vector<Point2>& controlPoints,
                       const std::vector<double>& weights,
                       const std::vector<double>& knots, int degree);

// ---------- 基函数 ----------
double basisFunction(int i, int p, double u, const std::vector<double>& knots);
// de Boor 三角：返回局部非零基函数 N[0..p]，N[i] = N_{span-p+i,p}(u)
std::vector<double> deBoorBasis(int p, double u, const std::vector<double>& knots, int span);

struct BasisWithDerivs {
    std::vector<double> N;      // 基函数（长度 = 控制点数）
    std::vector<double> dN;     // 一阶导数
    std::vector<double> ddN;    // 二阶导数（maxOrder>=2 才有值）
};
BasisWithDerivs basisFunctionsWithDerivatives(int p, double u,
                                              const std::vector<double>& knots,
                                              int maxOrder = 2);
int findSpanIndex(int p, double u, const std::vector<double>& knots);

// ---------- 曲线求值 ----------
Point2 evaluatePoint(const NURBSCurve& curve, double u);
std::vector<Point2> evaluatePoints(const NURBSCurve& curve, const std::vector<double>& uValues);
Point2 evaluateFirstDerivative(const NURBSCurve& curve, double u);
Point2 evaluateSecondDerivative(const NURBSCurve& curve, double u);
Point2 evaluateTangent(const NURBSCurve& curve, double u);
double evaluateCurvature(const NURBSCurve& curve, double u);
// 合并求值：一次扫描同时返回 点(u) 与 一阶导(u)。内部点路径复刻
// evaluateFirstDerivative 的累加（点 = hx/hw、导 = (dhx*hw-hx*dhw)/w2），
// 与分别调用 evaluatePoint + evaluateFirstDerivative 逐位一致；仅省一次
// de Boor 扫描（共享的 hx/hy/hw）。边界(≈uMax)或度0时退化为分别调用保证一致。
std::pair<Point2, Point2> evaluatePointAndDerivative(const NURBSCurve& curve, double u);

// ---------- 几何属性 ----------
AABB computeAABB(const NURBSCurve& curve);
bool aabbIntersect(const AABB& box1, const AABB& box2);

} // namespace nurbs

#endif // NURBS_CORE_HPP
