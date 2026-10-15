// ============================================================================
// nurbs_core.cpp — NURBS 数学核心库实现（C++）
// 逐函数对照 src/nurbs/math_core.js 移植，运算顺序与类型保持一致以保证逐位一致。
// ============================================================================
#include "nurbs/nurbs_core.hpp"
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cctype>

namespace nurbs {

// ============================================================================
// 合法性校验
// ============================================================================
void validateKnotVector(const std::vector<double>& knots, int degree, int nControlPoints) {
    if (knots.empty())
        throw std::invalid_argument("Knot 向量必须是非空数组");
    const auto expectedLength = nControlPoints + degree + 1;
    if (static_cast<int>(knots.size()) != expectedLength)
        throw std::invalid_argument("Knot 向量长度错误：期望 " + std::to_string(expectedLength) +
                                    "，实际 " + std::to_string(knots.size()));
    for (size_t i = 0; i + 1 < knots.size(); ++i) {
        if (knots[i] > knots[i + 1] + kKnotTolerance)
            throw std::invalid_argument("Knot 向量必须非递减");
    }
    // 首尾重节点(clamped)校验由生成阶段完成；此处（逐次求值热路径）不做重复校验。
}

void validateControlPoints(const std::vector<ControlPoint>& controlPoints) {
    if (controlPoints.empty())
        throw std::invalid_argument("控制点数组不能为空");
    for (size_t i = 0; i < controlPoints.size(); ++i) {
        const auto& pt = controlPoints[i];
        if (pt.w <= 0)
            throw std::invalid_argument("控制点 [" + std::to_string(i) + "] 权重必须为正数");
        if (!std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.w))
            throw std::invalid_argument("控制点 [" + std::to_string(i) + "] 包含非有限数值");
    }
}

void validateCurve(const NURBSCurve& curve) {
    validateControlPoints(curve.controlPoints);
    validateKnotVector(curve.knots, curve.degree, static_cast<int>(curve.controlPoints.size()));
    if (!std::isfinite(curve.uMin) || !std::isfinite(curve.uMax))
        throw std::invalid_argument("参数域 uMin/uMax 必须为数值");
    if (curve.uMin >= curve.uMax)
        throw std::invalid_argument("参数域无效：uMin >= uMax");
}

// ============================================================================
// Knot 向量生成
// ============================================================================
std::vector<double> generateUniformKnotVector(int n, int p) {
    const int m = n + p;
    std::vector<double> knots(m + 1);
    for (int i = 0; i <= p; ++i) knots[i] = 0;
    for (int i = m - p; i <= m; ++i) knots[i] = 1;
    const int internalKnotCount = n - p - 1;
    if (internalKnotCount > 0) {
        const double step = 1.0 / (internalKnotCount + 1);
        for (int i = 1; i <= internalKnotCount; ++i)
            knots[p + i] = i * step;
    }
    return knots;
}

std::vector<double> generateOpenUniformKnotVector(int n, int p) {
    return generateUniformKnotVector(n, p);
}

std::vector<double> generatePiecewiseBezierKnotVector(int n, int p) {
    const int segments = n - p;
    if (segments <= 0)
        throw std::runtime_error("控制点数量不足以生成分段贝兹曲线");
    std::vector<double> knots;
    knots.reserve(segments * p + p + 1);
    for (int i = 0; i <= p; ++i) knots.push_back(0);
    for (int i = 1; i < segments; ++i) {
        const double u = static_cast<double>(i) / segments;
        for (int j = 0; j < p; ++j) knots.push_back(u);
    }
    for (int i = 0; i <= p; ++i) knots.push_back(1);
    return knots;
}

std::vector<double> generateKnotVector(int degree, int numControlPoints,
                                       const std::string& knotType) {
    std::string t = knotType;
    std::transform(t.begin(), t.end(), t.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    if (t == "uniform" || t == "clamped")
        return generateUniformKnotVector(numControlPoints, degree);
    if (t == "bezier")
        return generatePiecewiseBezierKnotVector(numControlPoints, degree);
    return generateUniformKnotVector(numControlPoints, degree);
}

// ============================================================================
// 曲线创建
// ============================================================================
NURBSCurve createCurve(const std::vector<Point2>& controlPoints,
                       const std::vector<double>& weights,
                       const std::vector<double>& knots, int degree) {
    if (controlPoints.size() != weights.size())
        throw std::runtime_error("控制点数量必须与权重数量一致");
    NURBSCurve c;
    c.degree = degree;
    c.knots = knots;
    c.uMin = knots[degree];
    c.uMax = knots[knots.size() - 1 - degree];
    c.controlPoints.resize(controlPoints.size());
    for (size_t i = 0; i < controlPoints.size(); ++i) {
        c.controlPoints[i].x = controlPoints[i].x;
        c.controlPoints[i].y = controlPoints[i].y;
        c.controlPoints[i].w = weights[i];
    }
    return c;
}

// ============================================================================
// Cox-de Boor 基函数
// ============================================================================
double basisFunction(int i, int p, double u, const std::vector<double>& knots) {
    if (p == 0) {
        if (u >= knots[i] && u < knots[i + 1]) return 1.0;
        if (i == static_cast<int>(knots.size()) - p - 2 &&
            floatEqual(u, knots[knots.size() - 1], kKnotTolerance)) return 1.0;
        return 0.0;
    }
    double c1 = 0.0, c2 = 0.0;
    const double denom1 = knots[i + p] - knots[i];
    if (!floatEqual(denom1, 0, kKnotTolerance))
        c1 = (u - knots[i]) / denom1 * basisFunction(i, p - 1, u, knots);
    const double denom2 = knots[i + p + 1] - knots[i + 1];
    if (!floatEqual(denom2, 0, kKnotTolerance))
        c2 = (knots[i + p + 1] - u) / denom2 * basisFunction(i + 1, p - 1, u, knots);
    return c1 + c2;
}

std::vector<double> deBoorBasis(int p, double u, const std::vector<double>& knots, int span) {
    std::vector<double> N(p + 1, 0.0);
    std::vector<double> left(p + 1, 0.0), right(p + 1, 0.0);
    N[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double denom = right[r + 1] + left[j - r];
            const double temp = floatEqual(denom, 0, kKnotTolerance) ? 0.0 : N[r] / denom;
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    return N;
}

BasisWithDerivs basisFunctionsWithDerivatives(int p, double u,
                                              const std::vector<double>& knots,
                                              int maxOrder) {
    const int n = static_cast<int>(knots.size()) - p - 2; // 控制点数量 - 1
    BasisWithDerivs out;
    out.N.assign(n + 1, 0.0);
    out.dN.assign(n + 1, 0.0);
    out.ddN.assign(n + 1, 0.0);
    if (p == 0) {
        const int span = findSpanIndex(0, u, knots);
        out.N[span] = 1.0;
        return out;
    }
    const int span = findSpanIndex(p, u, knots);

    // ---- ndu 表 ----
    std::vector<std::vector<double>> ndu(p + 1, std::vector<double>(p + 1, 0.0));
    std::vector<double> left(p + 1, 0.0), right(p + 1, 0.0);
    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = floatEqual(ndu[j][r], 0, kKnotTolerance)
                                    ? 0.0
                                    : ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for (int r = 0; r <= p; ++r)
        out.N[span - p + r] = ndu[r][p];

    // ---- 导数递推（A2.3） ----
    const int nd = std::max(1, std::min(maxOrder, p));
    std::vector<std::vector<double>> ders(nd + 1, std::vector<double>(p + 1, 0.0));
    // a 用 2×(p+1)，s1/s2 交换，与 JS 一致
    std::vector<std::vector<double>> a(2, std::vector<double>(p + 1, 0.0));
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= nd; ++k) {
            double d = 0.0;
            const int rk = r - k;
            const int pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = (rk >= -1) ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            const int t = s1; s1 = s2; s2 = t;
        }
    }
    double factor = p;
    for (int k = 1; k <= nd; ++k) {
        for (int r = 0; r <= p; ++r) ders[k][r] *= factor;
        factor *= (p - k);
    }
    for (int r = 0; r <= p; ++r) {
        out.dN[span - p + r] = ders[1][r];
        if (nd >= 2) out.ddN[span - p + r] = ders[2][r];
    }
    return out;
}

int findSpanIndex(int p, double u, const std::vector<double>& knots) {
    if (floatEqual(u, knots[knots.size() - 1], kKnotTolerance))
        return static_cast<int>(knots.size()) - p - 2;
    int low = p;
    int high = static_cast<int>(knots.size()) - 1 - p;
    while (low < high) {
        const int mid = (low + high + 1) / 2;
        if (u < knots[mid]) high = mid - 1;
        else low = mid;
    }
    return low;
}

// ============================================================================
// 曲线求值
// ============================================================================
Point2 evaluatePoint(const NURBSCurve& curve, double u) {
    validateCurve(curve);
    u = std::max(curve.uMin, std::min(curve.uMax, u));
    const int span = findSpanIndex(curve.degree, u, curve.knots);
    const std::vector<double> N = deBoorBasis(curve.degree, u, curve.knots, span);

    double hx = 0.0, hy = 0.0, hw = 0.0;
    for (int i = 0; i <= curve.degree; ++i) {
        const int idx = span - curve.degree + i;
        if (idx >= 0 && idx < static_cast<int>(curve.controlPoints.size())) {
            const double Nv = N[i];
            const double w = curve.controlPoints[idx].w;
            hx += Nv * w * curve.controlPoints[idx].x;
            hy += Nv * w * curve.controlPoints[idx].y;
            hw += Nv * w;
        }
    }
    if (floatEqual(hw, 0, kEpsilon)) {
        if (floatEqual(u, curve.uMax, kKnotTolerance)) {
            const auto& last = curve.controlPoints[curve.controlPoints.size() - 1];
            return { last.x, last.y };
        }
        throw std::runtime_error("权重和为零，无法投影到笛卡尔坐标");
    }
    const double x = hx / hw;
    const double y = hy / hw;
    if (!std::isfinite(x) || !std::isfinite(y))
        throw std::runtime_error("计算结果包含非有限数值");
    return { x, y };
}

std::vector<Point2> evaluatePoints(const NURBSCurve& curve, const std::vector<double>& uValues) {
    std::vector<Point2> out;
    out.reserve(uValues.size());
    for (double u : uValues) out.push_back(evaluatePoint(curve, u));
    return out;
}

Point2 evaluateFirstDerivative(const NURBSCurve& curve, double u) {
    validateCurve(curve);
    u = std::max(curve.uMin, std::min(curve.uMax, u));
    if (floatEqual(u, curve.uMax, kKnotTolerance))
        u = curve.uMax - 1e-9;
    const int span = findSpanIndex(curve.degree, u, curve.knots);
    const std::vector<double> N = deBoorBasis(curve.degree, u, curve.knots, span);

    if (curve.degree >= 1) {
        const int s2 = findSpanIndex(curve.degree - 1, u, curve.knots);
        const std::vector<double> Nm1 = deBoorBasis(curve.degree - 1, u, curve.knots, s2);
        const int base2 = s2 - (curve.degree - 1);
        std::vector<double> dN(curve.degree + 1, 0.0);
        for (int k = 0; k <= curve.degree; ++k) {
            const int idx = span - curve.degree + k;
            const int j1 = idx - base2;
            const int j2 = idx + 1 - base2;
            const double N1 = (j1 >= 0 && j1 <= curve.degree - 1) ? Nm1[j1] : 0.0;
            const double N2 = (j2 >= 0 && j2 <= curve.degree - 1) ? Nm1[j2] : 0.0;
            const double den1 = curve.knots[idx + curve.degree] - curve.knots[idx];
            const double den2 = curve.knots[idx + curve.degree + 1] - curve.knots[idx + 1];
            dN[k] = curve.degree * ((floatEqual(den1, 0, kKnotTolerance) ? 0.0 : N1 / den1)
                                   - (floatEqual(den2, 0, kKnotTolerance) ? 0.0 : N2 / den2));
        }
        double hx = 0, hy = 0, hw = 0, dhx = 0, dhy = 0, dhw = 0;
        for (int k = 0; k <= curve.degree; ++k) {
            const int idx = span - curve.degree + k;
            if (idx >= 0 && idx < static_cast<int>(curve.controlPoints.size())) {
                const double Nk = N[k], dNk = dN[k];
                const double w = curve.controlPoints[idx].w;
                const double px = curve.controlPoints[idx].x, py = curve.controlPoints[idx].y;
                hx += Nk * w * px; hy += Nk * w * py; hw += Nk * w;
                dhx += dNk * w * px; dhy += dNk * w * py; dhw += dNk * w;
            }
        }
        if (floatEqual(hw, 0, kEpsilon))
            throw std::runtime_error("一阶导数计算：权重和为零");
        const double w2 = hw * hw;
        const double dx = (dhx * hw - hx * dhw) / w2;
        const double dy = (dhy * hw - hy * dhw) / w2;
        if (!std::isfinite(dx) || !std::isfinite(dy))
            throw std::runtime_error("一阶导数计算包含非有限数值");
        return { dx, dy };
    }
    return { 0.0, 0.0 };
}

std::pair<Point2, Point2> evaluatePointAndDerivative(const NURBSCurve& curve, double u) {
    validateCurve(curve);
    const double uClamped = std::max(curve.uMin, std::min(curve.uMax, u));
    const double ud = floatEqual(uClamped, curve.uMax, kKnotTolerance) ? curve.uMax - 1e-9 : uClamped;
    // 边界(≈uMax)或度0：点与导的 u 处理不同，退化为分别调用，保证与
    // evaluatePoint + evaluateFirstDerivative 逐位一致。
    if (uClamped != ud || curve.degree < 1)
        return { evaluatePoint(curve, u), evaluateFirstDerivative(curve, u) };
    // 内部点：一次扫描同时出点+一阶导（公式与分别调用完全相同 → 逐位一致）
    const int span = findSpanIndex(curve.degree, uClamped, curve.knots);
    const std::vector<double> N = deBoorBasis(curve.degree, uClamped, curve.knots, span);
    const int s2 = findSpanIndex(curve.degree - 1, uClamped, curve.knots);
    const std::vector<double> Nm1 = deBoorBasis(curve.degree - 1, uClamped, curve.knots, s2);
    const int base2 = s2 - (curve.degree - 1);
    std::vector<double> dN(curve.degree + 1, 0.0);
    for (int k = 0; k <= curve.degree; ++k) {
        const int idx = span - curve.degree + k;
        const int j1 = idx - base2, j2 = idx + 1 - base2;
        const double N1 = (j1 >= 0 && j1 <= curve.degree - 1) ? Nm1[j1] : 0.0;
        const double N2 = (j2 >= 0 && j2 <= curve.degree - 1) ? Nm1[j2] : 0.0;
        const double den1 = curve.knots[idx + curve.degree] - curve.knots[idx];
        const double den2 = curve.knots[idx + curve.degree + 1] - curve.knots[idx + 1];
        dN[k] = curve.degree * ((floatEqual(den1, 0, kKnotTolerance) ? 0.0 : N1 / den1)
                                - (floatEqual(den2, 0, kKnotTolerance) ? 0.0 : N2 / den2));
    }
    double hx = 0, hy = 0, hw = 0, dhx = 0, dhy = 0, dhw = 0;
    for (int k = 0; k <= curve.degree; ++k) {
        const int idx = span - curve.degree + k;
        if (idx >= 0 && idx < static_cast<int>(curve.controlPoints.size())) {
            const double Nk = N[k], dNk = dN[k];
            const double w = curve.controlPoints[idx].w;
            const double px = curve.controlPoints[idx].x, py = curve.controlPoints[idx].y;
            hx += Nk * w * px;  hy += Nk * w * py;  hw += Nk * w;
            dhx += dNk * w * px; dhy += dNk * w * py; dhw += dNk * w;
        }
    }
    if (floatEqual(hw, 0, kEpsilon))
        throw std::runtime_error("一阶导数计算：权重和为零");
    const double w2 = hw * hw;
    const double dx = (dhx * hw - hx * dhw) / w2;
    const double dy = (dhy * hw - hy * dhw) / w2;
    if (!std::isfinite(dx) || !std::isfinite(dy))
        throw std::runtime_error("一阶导数计算包含非有限数值");
    return { { hx / hw, hy / hw }, { dx, dy } };
}

Point2 evaluateSecondDerivative(const NURBSCurve& curve, double u) {
    validateCurve(curve);
    u = std::max(curve.uMin, std::min(curve.uMax, u));
    const auto bd = basisFunctionsWithDerivatives(curve.degree, u, curve.knots);

    double hx = 0.0, hy = 0.0, hw = 0.0;
    double dhx = 0.0, dhy = 0.0, dhw = 0.0;
    double ddhx = 0.0, ddhy = 0.0, ddhw = 0.0;
    for (size_t i = 0; i < curve.controlPoints.size(); ++i) {
        const double Ni = bd.N[i], dNi = bd.dN[i], ddNi = bd.ddN[i];
        const double w = curve.controlPoints[i].w;
        const double px = curve.controlPoints[i].x, py = curve.controlPoints[i].y;
        hx += Ni * w * px;  hy += Ni * w * py;  hw += Ni * w;
        dhx += dNi * w * px; dhy += dNi * w * py; dhw += dNi * w;
        ddhx += ddNi * w * px; ddhy += ddNi * w * py; ddhw += ddNi * w;
    }
    const double w2 = hw * hw;
    const double w3 = w2 * hw;
    const double x = hx / hw;
    const double y = hy / hw;
    const double dx = (dhx * hw - hx * dhw) / w2;
    const double dy = (dhy * hw - hy * dhw) / w2;
    const double ddx = (ddhx * w2 - 2 * dhx * dhw * hw - hx * ddhw * hw + 2 * hx * dhw * dhw) / w3;
    const double ddy = (ddhy * w2 - 2 * dhy * dhw * hw - hy * ddhw * hw + 2 * hy * dhw * dhw) / w3;
    if (!std::isfinite(ddx) || !std::isfinite(ddy))
        throw std::runtime_error("二阶导数计算包含非有限数值");
    return { ddx, ddy };
}

Point2 evaluateTangent(const NURBSCurve& curve, double u) {
    const Point2 deriv = evaluateFirstDerivative(curve, u);
    const double len = std::sqrt(deriv.x * deriv.x + deriv.y * deriv.y);
    if (len < kEpsilon)
        throw std::runtime_error("切向量为零向量，无法单位化");
    return { deriv.x / len, deriv.y / len };
}

double evaluateCurvature(const NURBSCurve& curve, double u) {
    const Point2 d1 = evaluateFirstDerivative(curve, u);
    const Point2 d2 = evaluateSecondDerivative(curve, u);
    const double d1Len = std::sqrt(d1.x * d1.x + d1.y * d1.y);
    if (d1Len < kEpsilon) return 0.0;
    const double cross = d1.x * d2.y - d1.y * d2.x;
    return std::abs(cross) / std::pow(d1Len, 3);
}

// ============================================================================
// 包围盒
// ============================================================================
AABB computeAABB(const NURBSCurve& curve) {
    validateCurve(curve);
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    for (const auto& pt : curve.controlPoints) {
        minX = std::min(minX, pt.x);
        minY = std::min(minY, pt.y);
        maxX = std::max(maxX, pt.x);
        maxY = std::max(maxY, pt.y);
    }
    return { minX, minY, maxX, maxY };
}

bool aabbIntersect(const AABB& box1, const AABB& box2) {
    return !(box1.maxX < box2.minX || box1.minX > box2.maxX ||
             box1.maxY < box2.minY || box1.minY > box2.maxY);
}

} // namespace nurbs
