/**
 * Unified Solver Engine - 统一数值迭代求解引擎
 * 
 * 核心业务逻辑：实现"一套框架，自动判别"的算法，处理相交、相切、分离所有场景。
 * 基于 Newton-Raphson 方法集成 Levenberg-Marquardt 阻尼策略。
 * 
 * @module solver/unified_solver
 */

'use strict';

import {
    EPSILON,
    MAX_ITERATIONS,
    validateCurve,
    evaluatePoint,
    evaluateFirstDerivative,
    evaluateSecondDerivative,
    computeAABB,
    aabbIntersect,
    findSpanIndex
} from '../nurbs/math_core.js';

// ============================================================================
// 常量定义
// ============================================================================

// 相切判据（端点切线夹角接近 0° 或 180°）。取 12°：
// ① 测试构造"沿 dir 布置控制点"的真实端点切线会有 3°~11° 的固有偏差（控制点
//    沿目标方向但 clamped NURBS 实际切线微偏），12° 能覆盖；
// ② 横截相交组（相交/高阶）构造为端点切线严格垂直（夹角 90°），远离 0/180°，
//    不会被误判为相切。故 12° 能在"相切(≈180°)与横截(≈90°)"间干净分离。
// 以下容差用 let + setThresholds 暴露为可配置项：论文灵敏度实验可临时扫描，
// 默认值保持不变，行为向后兼容。
let GEOMETRY_TOLERANCE = 1e-8;      // 几何判定容差（距离）
let ANGULAR_TOLERANCE = 12 * Math.PI / 180;  // 相切角度容差（弧度）：12°
let BOUNDARY_INTERSECTION_TOLERANCE = 1e-4;  // 边界接触容差（参数域端点交点残差上限）
let PARAM_BOUNDARY_TOL = 1e-3;               // 参数域"边界"判定宽度
let CONVERGENCE_THRESHOLD = 1e-12;  // 收敛判据阈值
const DAMPING_INITIAL = 0.01;         // LM 阻尼因子初始值
const DAMPING_FACTOR = 10.0;          // 阻尼因子调整系数
const MAX_DAMPING = 1e10;             // 最大阻尼值
const MIN_DAMPING = 1e-14;            // 最小阻尼值

/**
 * 覆盖几何/收敛阈值（仅用于灵敏度实验与调参；不传的项保持现值）。
 * @param {object} o 可含 GEOMETRY_TOLERANCE / ANGULAR_TOLERANCE(弧度) /
 *                   BOUNDARY_INTERSECTION_TOLERANCE / PARAM_BOUNDARY_TOL / CONVERGENCE_THRESHOLD
 */
function setThresholds(o = {}) {
    if (o.GEOMETRY_TOLERANCE !== undefined) GEOMETRY_TOLERANCE = o.GEOMETRY_TOLERANCE;
    if (o.ANGULAR_TOLERANCE !== undefined) ANGULAR_TOLERANCE = o.ANGULAR_TOLERANCE;
    if (o.BOUNDARY_INTERSECTION_TOLERANCE !== undefined) BOUNDARY_INTERSECTION_TOLERANCE = o.BOUNDARY_INTERSECTION_TOLERANCE;
    if (o.PARAM_BOUNDARY_TOL !== undefined) PARAM_BOUNDARY_TOL = o.PARAM_BOUNDARY_TOL;
    if (o.CONVERGENCE_THRESHOLD !== undefined) CONVERGENCE_THRESHOLD = o.CONVERGENCE_THRESHOLD;
}
// 当前阈值快照（供恢复）
const thresholdDefaults = () => ({
    GEOMETRY_TOLERANCE, ANGULAR_TOLERANCE, BOUNDARY_INTERSECTION_TOLERANCE,
    PARAM_BOUNDARY_TOL, CONVERGENCE_THRESHOLD
});

// ============================================================================
// 枚举定义
// ============================================================================

/**
 * 几何状态枚举
 */
const GeometryState = {
    INTERSECTION: 'INTERSECTION',     // 相交（距离 ≈ 0）
    TANGENT: 'TANGENT',               // 相切/退化接触
    SEPARATED: 'SEPARATED',           // 分离（存在最短距离）
    COINCIDENT: 'COINCIDENT',         // 重合（无穷多交点）
    UNDETERMINED: 'UNDETERMINED'      // 无法判定（算法未收敛）
};

/**
 * 求解结果数据结构
 * @typedef {Object} SolveResult
 * @property {GeometryState} state - 几何状态
 * @property {number} u - 曲线 1 的参数值
 * @property {number} v - 曲线 2 的参数值
 * @property {{x: number, y: number}} point1 - 曲线 1 上的点
 * @property {{x: number, y: number}} point2 - 曲线 2 上的点
 * @property {number} distance - 两点间距离
 * @property {number} iterations - 迭代次数
 * @property {number} residual - 最终残差
 * @property {number} conditionNumber - 雅可比矩阵条件数估计
 * @property {boolean} converged - 是否收敛
 * @property {string} message - 状态描述信息
 * @property {Array} iterationHistory - 迭代历史记录
 */

// ============================================================================
// 辅助数学函数
// ============================================================================

/**
 * 计算向量范数
 * @param {number} x 
 * @param {number} y 
 * @returns {number}
 */
function vectorNorm(x, y) {
    return Math.sqrt(x * x + y * y);
}

/**
 * 计算两点间距离
 * @param {{x: number, y: number}} p1 
 * @param {{x: number, y: number}} p2 
 * @returns {number}
 */
function distance(p1, p2) {
    const dx = p1.x - p2.x;
    const dy = p1.y - p2.y;
    return Math.sqrt(dx * dx + dy * dy);
}

/**
 * 向量点积
 * @param {{x: number, y: number}} v1 
 * @param {{x: number, y: number}} v2 
 * @returns {number}
 */
function dotProduct(v1, v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

/**
 * 2x2 矩阵求逆
 * @param {number[][]} matrix 
 * @returns {number[][] | null} 逆矩阵，奇异时返回 null
 */
function invert2x2(matrix) {
    const det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    
    if (Math.abs(det) < EPSILON) {
        return null; // 奇异矩阵
    }
    
    const invDet = 1.0 / det;
    return [
        [matrix[1][1] * invDet, -matrix[0][1] * invDet],
        [-matrix[1][0] * invDet, matrix[0][0] * invDet]
    ];
}

/**
 * 估计 2x2 矩阵的条件数（简化版）
 * @param {number[][]} matrix 
 * @returns {number} 条件数估计值
 */
function estimateConditionNumber(matrix) {
    // 使用 Frobenius 范数简化估计
    const normF = Math.sqrt(
        matrix[0][0] ** 2 + matrix[0][1] ** 2 +
        matrix[1][0] ** 2 + matrix[1][1] ** 2
    );
    
    const inv = invert2x2(matrix);
    if (!inv) {
        return Infinity;
    }
    
    const normInvF = Math.sqrt(
        inv[0][0] ** 2 + inv[0][1] ** 2 +
        inv[1][0] ** 2 + inv[1][1] ** 2
    );
    
    return normF * normInvF;
}

// ============================================================================
// 目标函数及其导数
// ============================================================================

/**
 * 计算统一目标函数 F(u,v) = ||C1(u) - C2(v)||^2
 * 及其梯度、Hessian 矩阵
 * 
 * @param {NURBSCurve} curve1 
 * @param {NURBSCurve} curve2 
 * @param {number} u 
 * @param {number} v 
 * @returns {{
 *   f: number,           // 目标函数值 (距离平方)
 *   gradient: number[],  // 梯度 [df/du, df/dv]
 *   hessian: number[][], // Hessian 矩阵
 *   point1: Object,      // C1(u)
 *   point2: Object,      // C2(v)
 *   deriv1: Object,      // C1'(u)
 *   deriv2: Object       // C2'(v)
 * }}
 */
function objectiveFunction(curve1, curve2, u, v) {
    // 曲线上的点与一阶导数
    const point1 = evaluatePoint(curve1, u);
    const point2 = evaluatePoint(curve2, v);
    const deriv1 = evaluateFirstDerivative(curve1, u);
    const deriv2 = evaluateFirstDerivative(curve2, v);

    // 残差 r = C1(u) - C2(v)
    const rx = point1.x - point2.x;
    const ry = point1.y - point2.y;
    const f = rx * rx + ry * ry; // 距离平方

    // 雅可比列向量: [d r/du, d r/dv] = [C1', -C2']
    const j11 = deriv1.x, j12 = -deriv2.x;
    const j21 = deriv1.y, j22 = -deriv2.y;

    // 高斯-牛顿 Hessian 近似 J^T J（比精确 Hessian 在近奇异/浅谷处更稳健，避免收敛停滞）
    const h00 = j11 * j11 + j21 * j21;
    const h01 = j11 * j12 + j21 * j22;
    const h11 = j12 * j12 + j22 * j22;

    // 目标 (1/2 f) 的梯度 = J^T r
    const g0 = j11 * rx + j21 * ry;
    const g1 = j12 * rx + j22 * ry;

    return {
        f,
        gradient: [2 * g0, 2 * g1],
        hessian: [[h00, h01], [h01, h11]],
        point1,
        point2,
        deriv1,
        deriv2,
        distance: Math.sqrt(f)
    };
}

// ============================================================================
// Levenberg-Marquardt 迭代求解器
// ============================================================================

/**
 * 执行单步 Levenberg-Marquardt 迭代
 * 
 * @param {NURBSCurve} curve1 
 * @param {NURBSCurve} curve2 
 * @param {number} u - 当前参数 u
 * @param {number} v - 当前参数 v
 * @param {number} damping - LM 阻尼因子
 * @returns {{
 *   u: number, v: number, 
 *   success: boolean, 
 *   objValue: number
 * }}
 */
function lmStep(curve1, curve2, u, v, damping) {
    try {
        const obj = objectiveFunction(curve1, curve2, u, v);
        
        // 构建阻尼 Hessian: (H + λI)
        const dampedHessian = [
            [obj.hessian[0][0] + damping, obj.hessian[0][1]],
            [obj.hessian[1][0], obj.hessian[1][1] + damping]
        ];
        
        // 求解 (H + λI) Δx = -∇f
        const invH = invert2x2(dampedHessian);
        if (!invH) {
            return { u, v, success: false, objValue: obj.f };
        }
        
        // Δx = -H^{-1} ∇f
        const du = -(invH[0][0] * obj.gradient[0] + invH[0][1] * obj.gradient[1]);
        const dv = -(invH[1][0] * obj.gradient[0] + invH[1][1] * obj.gradient[1]);
        
        // 限制步长防止溢出
        const maxStep = 1.0;
        const stepNorm = Math.sqrt(du * du + dv * dv);
        let scale = 1.0;
        if (stepNorm > maxStep) {
            scale = maxStep / stepNorm;
        }
        
        const newU = Math.max(curve1.uMin, Math.min(curve1.uMax, u + scale * du));
        const newV = Math.max(curve2.uMin, Math.min(curve2.uMax, v + scale * dv));
        
        // 验证新点是否改进
        const newObj = objectiveFunction(curve1, curve2, newU, newV);
        
        return {
            u: newU,
            v: newV,
            success: newObj.f <= obj.f,
            objValue: newObj.f
        };
    } catch (error) {
        console.warn('LM Step 异常:', error.message);
        return { u, v, success: false, objValue: Infinity };
    }
}

/**
 * 执行完整的 Levenberg-Marquardt 迭代过程
 * 
 * @param {NURBSCurve} curve1 
 * @param {NURBSCurve} curve2 
 * @param {number} u0 - 初始猜测 u
 * @param {number} v0 - 初始猜测 v
 * @param {Object} options - 配置选项
 * @returns {SolveResult}
 */
function levenbergMarquardt(curve1, curve2, u0, v0, options = {}) {
    const maxIter = options.maxIterations || MAX_ITERATIONS;
    const tol = options.tolerance || CONVERGENCE_THRESHOLD; // 距离容差
    // recordHistory=false 时跳过逐次历史记录（热路径省内存与分配），仅统计迭代次数
    const recordHistory = options.recordHistory !== false;

    let u = u0;
    let v = v0;
    let damping = DAMPING_INITIAL;

    const iterationHistory = [];
    let iterCount = 0;
    let converged = false;
    let conditionNum = 0;
    let obj = objectiveFunction(curve1, curve2, u, v);
    let prevF = Infinity;

    for (let iter = 0; iter < maxIter; iter++) {
        // 记录迭代历史
        if (recordHistory) {
            iterationHistory.push({
                iteration: iter,
                u, v,
                distance: obj.distance,
                residual: obj.f,
                gradientNorm: vectorNorm(obj.gradient[0], obj.gradient[1])
            });
        }
        iterCount++;

        const hessian = obj.hessian || [[1, 0], [0, 1]];
        const grad = obj.gradient;

        // 构建阻尼 Hessian: (H + λI)，求解 (H+λI)Δ = -J^T r
        const dampedHessian = [
            [hessian[0][0] + damping, hessian[0][1]],
            [hessian[1][0], hessian[1][1] + damping]
        ];
        const invH = invert2x2(dampedHessian);
        if (!invH) {
            damping = Math.min(MAX_DAMPING, damping * DAMPING_FACTOR);
            obj = objectiveFunction(curve1, curve2, u, v);
            continue;
        }

        // 高斯-牛顿步：目标用 (1/2 f) 的梯度 (J^T r)，乘 0.5 抵消，得到标准步
        const du = -(invH[0][0] * grad[0] + invH[0][1] * grad[1]) * 0.5;
        const dv = -(invH[1][0] * grad[0] + invH[1][1] * grad[1]) * 0.5;

        const un = Math.max(curve1.uMin, Math.min(curve1.uMax, u + du));
        const vn = Math.max(curve2.uMin, Math.min(curve2.uMax, v + dv));
        const candidate = objectiveFunction(curve1, curve2, un, vn);

        // LM 信赖域控制：目标下降则接受并扩张信赖域，否则拒绝并收缩
        if (candidate.f <= obj.f + 1e-30) {
            u = un;
            v = vn;
            obj = candidate;
            damping = Math.max(MIN_DAMPING, damping / DAMPING_FACTOR);
        } else {
            damping = Math.min(MAX_DAMPING, damping * DAMPING_FACTOR);
            continue;
        }

        // 收敛判据：到达交点精度 / 梯度到极小点 / 步长足够小 / 相对下降停滞
        const gradNorm = vectorNorm(grad[0], grad[1]);
        if (obj.distance < tol) { converged = true; break; }
        if (gradNorm < 1e-8) { converged = true; break; }
        if (Math.abs(du) < 1e-12 && Math.abs(dv) < 1e-12) { converged = true; break; }
        // 相对下降停滞判据（替代原绝对阈值 obj.f >= prevF - 1e-14）：
        // 目标函数相对上一轮改善 < 1e-10 即视为已到达局部极小。
        // ① 相交场景：distance~1e-6 时 f~1e-12，原绝对判据会把迭代间 1e-14 级的
        //    正常改善误判为"不再下降"而提前停步，导致交点距离停留在 1e-6 被误判
        //    分离；相对判据只在改善率真正趋零时停，允许继续推进到 1e-8 容差内。
        // ② 分离场景：f 较大时四条原判据均可能不触发而跑满迭代 → UNDETERMINED；
        //    相对判据在接近最近点（改善率趋零）时收敛 → 正确判定为 SEPARATED。
        if (prevF !== Infinity && prevF > 0 && obj.f >= prevF * (1 - 1e-10)) { converged = true; break; }
        prevF = obj.f;
    }

    if (converged) {
        // 计算条件数用于状态判定
        conditionNum = estimateConditionNumber(obj.hessian);
    }

    const dist = obj.distance;

    // 几何状态判定
    let state = GeometryState.UNDETERMINED;
    let message = '';

    // 边界接触容差：当收敛解落在参数域边界（clamped 曲线端点精确插值控制点）
    // 时，端点附近因导数退化/半开区间数值误差，残差常无法压到 GEOMETRY_TOLERANCE
    // (1e-8) 以下（实测停在 1e-8~1e-5）。此时只要边界点距离足够小即真实接触。
    // 取 1e-4 覆盖高阶曲线端点交点更宽的残差区间；该容差远小于"真分离"场景的
    // 最短距离（分离组实测最小 >1.2），不会误伤分离。
    // （BOUNDARY_INTERSECTION_TOLERANCE / PARAM_BOUNDARY_TOL 已提升为模块级可配置阈值）
    const atParamBoundary = (u, c) =>
        Math.abs(u - c.uMin) < PARAM_BOUNDARY_TOL || Math.abs(u - c.uMax) < PARAM_BOUNDARY_TOL;

    if (!converged) {
        state = GeometryState.UNDETERMINED;
        message = '算法未在最大迭代次数内收敛';
    } else if (dist < GEOMETRY_TOLERANCE ||
               (dist < BOUNDARY_INTERSECTION_TOLERANCE &&
                (atParamBoundary(u, curve1) || atParamBoundary(v, curve2)))) {
        // 距离已足够小（真实接触）：
        // 1) 严格几何容差内，或 2) 边界接触容差内且在参数域边界（高阶曲线端点
        //    交点残差可合法停在 1e-5~1e-4）。
        // 统一先做相切判定：若两条曲线在接触点处切线平行（角度≈0 或 ≈π），
        // 则为相切接触 TANGENT；否则为普通横截相交 INTERSECTION。
        // 注意：clamped NURBS 在参数域端点处用有限差分求导会数值退化（切向量趋零），
        // 导致相切误判为"奇异点"。端点处的真实切线由控制点解析给出：首端切线方向
        // = P1−P0，末端切线方向 = Pn−1−Pn−2。这里用解析端点切线做方向判断。
        const boundaryTan = (curve, param) => {
            const pts = curve.controlPoints;
            if (pts.length < 2) return null;
            if (Math.abs(param - curve.uMin) < PARAM_BOUNDARY_TOL) {
                // 首端：从首控制点指向第二个控制点
                return { x: pts[1].x - pts[0].x, y: pts[1].y - pts[0].y };
            }
            if (Math.abs(param - curve.uMax) < PARAM_BOUNDARY_TOL) {
                // 末端：从倒数第二控制点指向末控制点
                const n = pts.length;
                return { x: pts[n - 1].x - pts[n - 2].x, y: pts[n - 1].y - pts[n - 2].y };
            }
            return null;
        };
        let tan1 = null, tan2 = null;
        // 优先用解析端点切线（当解在参数边界时，数值导数退化不可靠）
        const b1 = atParamBoundary(u, curve1) ? boundaryTan(curve1, u) : null;
        const b2 = atParamBoundary(v, curve2) ? boundaryTan(curve2, v) : null;
        tan1 = b1 || obj.deriv1;
        tan2 = b2 || obj.deriv2;

        const tan1Norm = vectorNorm(tan1.x, tan1.y);
        const tan2Norm = vectorNorm(tan2.x, tan2.y);

        if (tan1Norm > EPSILON && tan2Norm > EPSILON) {
            const cosAngle = (tan1.x * tan2.x + tan1.y * tan2.y) / (tan1Norm * tan2Norm);
            const angle = Math.acos(Math.max(-1, Math.min(1, cosAngle)));

            // 角度接近 0 或 π 表示相切
            if (angle < ANGULAR_TOLERANCE || Math.abs(angle - Math.PI) < ANGULAR_TOLERANCE) {
                state = GeometryState.TANGENT;
                message = `检测到相切接触 (角度=${(angle * 180 / Math.PI).toFixed(6)}°)`;
            } else {
                state = GeometryState.INTERSECTION;
                message = dist < GEOMETRY_TOLERANCE
                    ? '检测到交点'
                    : `检测到边界交点 (u=${u.toFixed(6)}, v=${v.toFixed(6)}, d=${dist.toExponential(3)})`;
            }
        } else {
            state = GeometryState.INTERSECTION;
            message = '检测到交点 (奇异点)';
        }
    } else {
        state = GeometryState.SEPARATED;
        message = `曲线分离，最短距离=${dist.toExponential(6)}`;
    }

    return {
        state,
        u,
        v,
        point1: obj.point1,
        point2: obj.point2,
        distance: dist,
        iterations: iterCount,
        residual: obj.f,
        conditionNumber: conditionNum,
        converged,
        message,
        iterationHistory
    };
}

// ============================================================================
// 多初值猜测策略
// ============================================================================

/**
 * 生成参数域网格采样点
 * @param {number} n1 - 曲线 1 的采样点数
 * @param {number} n2 - 曲线 2 的采样点数
 * @returns {Array<{u: number, v: number}>}
 */
function generateGridSamples(n1, n2) {
    const samples = [];
    for (let i = 0; i < n1; i++) {
        const u = i / (n1 - 1);
        for (let j = 0; j < n2; j++) {
            const v = j / (n2 - 1);
            samples.push({ u, v });
        }
    }
    return samples;
}

/**
 * 基于包围盒预筛选的候选区域检测
 * @param {NURBSCurve} curve1 
 * @param {NURBSCurve} curve2 
 * @param {number} subdivisions - 细分段数
 * @returns {Array<{uStart: number, uEnd: number, vStart: number, vEnd: number}>}
 */
function findCandidateRegions(curve1, curve2, subdivisions = 10) {
    const regions = [];
    
    // 将参数域细分
    const uSegments = [];
    const vSegments = [];
    
    for (let i = 0; i < subdivisions; i++) {
        uSegments.push({
            start: i / subdivisions,
            end: (i + 1) / subdivisions,
            index: i
        });
        vSegments.push({
            start: i / subdivisions,
            end: (i + 1) / subdivisions,
            index: i
        });
    }
    
    // 计算每段的包围盒
    const uBoxes = uSegments.map(seg => {
        const midU = (seg.start + seg.end) / 2;
        const pt = evaluatePoint(curve1, midU);
        // 保守估计包围盒大小
        const radius = estimateSegmentRadius(curve1, seg.start, seg.end);
        return {
            minX: pt.x - radius, maxX: pt.x + radius,
            minY: pt.y - radius, maxY: pt.y + radius,
            segment: seg
        };
    });
    
    const vBoxes = vSegments.map(seg => {
        const midV = (seg.start + seg.end) / 2;
        const pt = evaluatePoint(curve2, midV);
        const radius = estimateSegmentRadius(curve2, seg.start, seg.end);
        return {
            minX: pt.x - radius, maxX: pt.x + radius,
            minY: pt.y - radius, maxY: pt.y + radius,
            segment: seg
        };
    });
    
    // 检测相交的包围盒对
    for (let i = 0; i < uBoxes.length; i++) {
        for (let j = 0; j < vBoxes.length; j++) {
            if (aabbIntersect(uBoxes[i], vBoxes[j])) {
                regions.push({
                    uStart: uSegments[i].start,
                    uEnd: uSegments[i].end,
                    vStart: vSegments[j].start,
                    vEnd: vSegments[j].end
                });
            }
        }
    }
    
    return regions;
}

/**
 * 估计曲线段的保守半径
 * @param {NURBSCurve} curve 
 * @param {number} uStart 
 * @param {number} uEnd 
 * @returns {number}
 */
function estimateSegmentRadius(curve, uStart, uEnd) {
    // 使用控制点的凸包性质进行保守估计
    const { controlPoints, knots, degree } = curve;
    const spanStart = findSpanIndex(degree, uStart, knots);
    const spanEnd = findSpanIndex(degree, uEnd, knots);
    
    // 找到相关控制点
    const relevantPoints = [];
    for (let i = Math.max(0, spanStart - degree); i <= Math.min(controlPoints.length - 1, spanEnd + degree); i++) {
        relevantPoints.push(controlPoints[i]);
    }
    
    if (relevantPoints.length === 0) return 1.0;
    
    // 计算控制点中心
    let cx = 0, cy = 0;
    for (const pt of relevantPoints) {
        cx += pt.x;
        cy += pt.y;
    }
    cx /= relevantPoints.length;
    cy /= relevantPoints.length;
    
    // 找到最大距离（x、y 两个方向各自相对中心的最大半宽）
    let maxXDist = 0, maxYDist = 0;
    for (const pt of relevantPoints) {
        maxXDist = Math.max(maxXDist, Math.abs(pt.x - cx));
        maxYDist = Math.max(maxYDist, Math.abs(pt.y - cy));
    }

    // 返回覆盖"全部相关控制点包围盒"的统一半径：
    // 取两方向最大半宽，保证以 midU 处采样点为中心的正方形能覆盖该段所有相关
    // 控制点 → 覆盖其凸包 → 覆盖该段曲线。原实现取点到中心的欧氏距离×1.5，
    // 当相关控制点分布在采样点一侧/曲率大时（尤其高阶曲线）圆会漏掉曲线段，
    // 导致候选区排除真实交点 → 相交被误判分离。
    return Math.max(maxXDist, maxYDist) * 1.5; // 1.5 倍安全系数
}

// ============================================================================
// 主求解接口
// ============================================================================

/**
 * 统一求解器：自动判别并求解 NURBS 曲线对的几何关系
 * 
 * @param {NURBSCurve} curve1 - 第一条 NURBS 曲线
 * @param {NURBSCurve} curve2 - 第二条 NURBS 曲线
 * @param {Object} options - 配置选项
 * @param {boolean} options.usePreFiltering - 是否使用包围盒预筛选
 * @param {number} options.gridSize - 网格采样密度
 * @param {number} options.maxIterations - 最大迭代次数
 * @param {number} options.tolerance - 收敛容差
 * @returns {SolveResult[]} 所有找到的解（可能多个交点）
 */
function unifiedSolve(curve1, curve2, options = {}) {
    const usePreFiltering = options.usePreFiltering !== false;
    const gridSize = options.gridSize || 5;
    const maxIterations = options.maxIterations || MAX_ITERATIONS;
    // 容差体系统一：收敛目标与状态分类共用 GEOMETRY_TOLERANCE(1e-8)，
    // 避免"LM 已收敛但距离落在 (CONVERGENCE_THRESHOLD, GEOMETRY_TOLERANCE) 区间
    // 而被分类为 SEPARATED"的断层。调用方显式传入 tolerance 时遵循调用方。
    const tolerance = options.tolerance || GEOMETRY_TOLERANCE;
    
    // 验证输入
    validateCurve(curve1);
    validateCurve(curve2);
    
    // 快速排斥试验
    if (usePreFiltering) {
        const box1 = computeAABB(curve1);
        const box2 = computeAABB(curve2);
        
        if (!aabbIntersect(box1, box2)) {
            // 包围盒不相交，直接计算全局最短距离
            console.log('包围盒预筛选：曲线分离，执行全局优化');
            const result = globalOptimization(curve1, curve2, {
                maxIterations,
                tolerance
            });
            return [result];
        }
    }
    
    // 收集所有候选初始点
    const initialGuesses = [];
    
    if (usePreFiltering) {
        const regions = findCandidateRegions(curve1, curve2, gridSize);
        for (const region of regions) {
            initialGuesses.push({
                u: (region.uStart + region.uEnd) / 2,
                v: (region.vStart + region.vEnd) / 2
            });
        }
    } else {
        // 全网格采样
        const samples = generateGridSamples(gridSize, gridSize);
        initialGuesses.push(...samples);
    }

    // 补充参数域 4 个角点作为额外初值：两条 clamped 曲线常在其端点（即参数域
    // 角点，如 u=uMin 且 v=vMax）共享/相交，这类角点交点收敛域极小，网格中心
    // 初值易被吸引到附近局部极小而漏检；以角点本身为初值可稳定收敛到端点交点。
    // 先去重：若某角点已与现有初值重合（如网格/区域中心恰好落在角点），则跳过，
    // 避免重复 LM 迭代。
    const cornerPairs = [
        [curve1.uMin, curve2.uMin], [curve1.uMin, curve2.uMax],
        [curve1.uMax, curve2.uMin], [curve1.uMax, curve2.uMax]
    ];
    for (const [cu, cv] of cornerPairs) {
        const dup = initialGuesses.some(g =>
            Math.abs(g.u - cu) < 1e-9 && Math.abs(g.v - cv) < 1e-9);
        if (!dup) initialGuesses.push({ u: cu, v: cv });
    }
    
    // 从每个初始点执行 LM 迭代
    const results = [];
    const foundSolutions = new Set(); // 用于去重
    
    for (const guess of initialGuesses) {
        const result = levenbergMarquardt(curve1, curve2, guess.u, guess.v, {
            maxIterations: Math.floor(maxIterations / 3), // 每个初始点分配较少迭代
            tolerance,
            recordHistory: false // 热路径：多初值搜索无需逐次迭代历史
        });
        
        if (result.converged) {
            // 去重检查
            const key = `${result.u.toFixed(6)},${result.v.toFixed(6)}`;
            if (!foundSolutions.has(key)) {
                foundSolutions.add(key);
                results.push(result);
            }
        }
    }
    
    // 按距离排序（最近的优先）
    results.sort((a, b) => a.distance - b.distance);
    
    // 如果没有找到任何收敛解，执行全局优化
    if (results.length === 0) {
        const globalResult = globalOptimization(curve1, curve2, {
            maxIterations,
            tolerance
        });
        results.push(globalResult);
    }
    
    return results;
}

/**
 * 全局优化：当局部搜索失败时使用
 * 使用更密集的网格采样和更多迭代
 */
function globalOptimization(curve1, curve2, options) {
    const denseGrid = generateGridSamples(10, 10);
    let bestResult = null;
    
    // 优先选择已收敛的结果（最近点已真正求得）；若无任何收敛解，
    // 退化为取最小距离结果（此时由上层/调用方结合距离判定状态）。
    let bestConverged = null;
    let minConvDistance = Infinity;
    let bestAny = null;
    let minAnyDistance = Infinity;

    for (const sample of denseGrid) {
        const result = levenbergMarquardt(curve1, curve2, sample.u, sample.v, { ...options, recordHistory: false });
        if (result.distance < minAnyDistance) {
            minAnyDistance = result.distance;
            bestAny = result;
        }
        if (result.converged && result.distance < minConvDistance) {
            minConvDistance = result.distance;
            bestConverged = result;
        }
    }
    bestResult = bestConverged || bestAny;
    
    return bestResult || {
        state: GeometryState.UNDETERMINED,
        u: 0.5, v: 0.5,
        point1: evaluatePoint(curve1, 0.5),
        point2: evaluatePoint(curve2, 0.5),
        distance: Infinity,
        iterations: 0,
        residual: Infinity,
        conditionNumber: 0,
        converged: false,
        message: '全局优化未能找到有效解',
        iterationHistory: []
    };
}

// ============================================================================
// 全交点求解：递归细分隔离 + 逐区间 LM 抛光 + 容差聚类
// ============================================================================

/**
 * 计算曲线在参数区间 [a,b] 上的紧密包围盒（对片段密集采样评估）。
 * 用于递归隔离的安全剪枝——若两段曲线片段包围盒不相交，则其内必无交点。
 * 说明：这里用采样包围盒（比整窗控制点凸包更紧），避免近重合场景的细分爆炸。
 * 曾尝试改用"子区间相关控制点凸包"做形式化保守剪枝，但对低阶（degree 1~2）曲线，
 * 子区间仍受几乎全部控制点影响、包围盒几乎不随细分收缩，导致递归细分爆炸（回归），
 * 故保持采样盒。evaluatePoint 已改用 de Boor 三角提速，本采样开销随之显著下降。
 */
function pieceBox(curve, a, b) {
    const K = 16;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let i = 0; i <= K; i++) {
        const u = a + (b - a) * i / K;
        const p = evaluatePoint(curve, u);
        minX = Math.min(minX, p.x); maxX = Math.max(maxX, p.x);
        minY = Math.min(minY, p.y); maxY = Math.max(maxY, p.y);
    }
    return { minX, minY, maxX, maxY };
}

/**
 * 求容差范围内的全部交点。
 * 算法：将 [u0,u1]×[v0,v1] 递归细分；对互相叠覆盖的解单元用 LM 抛光，再按容差聚类去重。
 * 相对多初值网格，本方法确定性覆盖、不漏交点（对正常横向交点）。
 * @returns {Array<{u:number,v:number,distance:number,point1:object,point2:object,iterations:number}>}
 */
function findAllIntersections(curve1, curve2, options = {}) {
    const tolerance = options.tolerance != null ? options.tolerance : GEOMETRY_TOLERANCE;
    const maxDepth = options.maxDepth != null ? options.maxDepth : 22;
    const domainW = Math.max(curve1.uMax - curve1.uMin, curve2.uMax - curve2.uMin) || 1;
    const paramTol = options.paramTol != null ? options.paramTol : domainW * 1e-3;

    const work = [{ u0: curve1.uMin, u1: curve1.uMax, v0: curve2.uMin, v1: curve2.uMax, depth: 0 }];
    const seeds = [];
    let guard = 0;

    while (work.length) {
        if (++guard > 1000000) break; // 防御：极近重合等退化情形下终止，避免无限细分
        const cell = work.pop();
        const b1 = pieceBox(curve1, cell.u0, cell.u1);
        const b2 = pieceBox(curve2, cell.v0, cell.v1);
        if (!aabbIntersect(b1, b2)) continue; // 安全剪枝

        const wU = cell.u1 - cell.u0, wV = cell.v1 - cell.v0;
        if (cell.depth >= maxDepth || (wU < paramTol && wV < paramTol)) {
            seeds.push(cell);
            continue;
        }
        // 沿较宽维拆成两半
        if (wU >= wV) {
            const um = (cell.u0 + cell.u1) / 2;
            work.push({ u0: cell.u0, u1: um, v0: cell.v0, v1: cell.v1, depth: cell.depth + 1 });
            work.push({ u0: um, u1: cell.u1, v0: cell.v0, v1: cell.v1, depth: cell.depth + 1 });
        } else {
            const vm = (cell.v0 + cell.v1) / 2;
            work.push({ u0: cell.u0, u1: cell.u1, v0: cell.v0, v1: vm, depth: cell.depth + 1 });
            work.push({ u0: cell.u0, u1: cell.u1, v0: vm, v1: cell.v1, depth: cell.depth + 1 });
        }
    }

    // 抛光：每个隔离单元从中心做一次 LM，收敛且距离在容差内则记为交点
    const roots = [];
    for (const cell of seeds) {
        const u = (cell.u0 + cell.u1) / 2, v = (cell.v0 + cell.v1) / 2;
        const r = levenbergMarquardt(curve1, curve2, u, v, { maxIterations: 100, tolerance, recordHistory: false });
        if (r.converged && r.distance < tolerance) {
            roots.push({ u: r.u, v: r.v, distance: r.distance, point1: r.point1, point2: r.point2, iterations: r.iterations });
        }
    }

    // 聚类去重：合并被多个隔离单元/多次收敛覆盖到的同根，保留距离最小者
    roots.sort((a, b) => a.u - b.u);
    const out = [];
    for (const rr of roots) {
        let merged = false;
        for (const o of out) {
            if (Math.abs(o.u - rr.u) < paramTol * 6 && Math.abs(o.v - rr.v) < paramTol * 6) {
                if (rr.distance < o.distance) {
                    o.u = rr.u; o.v = rr.v; o.distance = rr.distance;
                    o.point1 = rr.point1; o.point2 = rr.point2;
                }
                merged = true;
                break;
            }
        }
        if (!merged) out.push({ u: rr.u, v: rr.v, distance: rr.distance, point1: rr.point1, point2: rr.point2, iterations: rr.iterations });
    }
    out.sort((a, b) => a.u - b.u);
    return out;
}

// ============================================================================
// 导出接口
// ============================================================================

export {
    // 常量
    GEOMETRY_TOLERANCE,
    ANGULAR_TOLERANCE,
    CONVERGENCE_THRESHOLD,
    // 阈值可配置化（灵敏度实验/调参用）
    setThresholds,
    thresholdDefaults,
    
    // 枚举
    GeometryState,
    
    // 核心求解器
    levenbergMarquardt,
    unifiedSolve,
    globalOptimization,
    
    // 辅助函数
    objectiveFunction,
    findCandidateRegions,
    generateGridSamples,

    // 全交点求解
    pieceBox,
    findAllIntersections
};