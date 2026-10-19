/**
 * NURBS Math Core - 核心数学基础库
 * 
 * 实现严格的 NURBS 数学定义底层库，不包含业务逻辑，仅提供纯数学计算原语。
 * 遵循《The NURBS Book》(Les Piegl & Wayne Tiller) 标准算法实现。
 * 
 * @module nurbs/math_core
 */

'use strict';

// ============================================================================
// 常量定义
// ============================================================================

const EPSILON = 1e-10;        // 数值计算精度容差
const MAX_ITERATIONS = 100;   // 默认最大迭代次数
const knotTolerance = 1e-10;  // Knot 向量容差

// ============================================================================
// 数据结构定义
// ============================================================================

/**
 * 控制点数据结构
 * @typedef {Object} ControlPoint
 * @property {number} x - X 坐标
 * @property {number} y - Y 坐标
 * @property {number} w - 权重 (weight)，必须 > 0
 */

/**
 * NURBS 曲线数据结构
 * @typedef {Object} NURBSCurve
 * @property {ControlPoint[]} controlPoints - 控制点数组
 * @property {number[]} knots - Knot 向量 (非递减序列)
 * @property {number} degree - 曲线阶数 (degree)
 * @property {number} uMin - 参数域下界
 * @property {number} uMax - 参数域上界
 */

// ============================================================================
// 工具函数
// ============================================================================

/**
 * 浮点数比较（考虑容差）
 * @param {number} a 
 * @param {number} b 
 * @param {number} tol - 容差值
 * @returns {boolean}
 */
function floatEqual(a, b, tol = EPSILON) {
    return Math.abs(a - b) < tol;
}

/**
 * 验证 Knot 向量合法性
 * @param {number[]} knots 
 * @param {number} degree 
 * @param {number} nControlPoints 
 * @throws {Error} 当 Knot 向量不合法时抛出异常
 */
function validateKnotVector(knots, degree, nControlPoints) {
    if (!Array.isArray(knots) || knots.length === 0) {
        throw new Error('Knot 向量必须是非空数组');
    }
    
    const expectedLength = nControlPoints + degree + 1;
    if (knots.length !== expectedLength) {
        throw new Error(`Knot 向量长度错误：期望 ${expectedLength}，实际 ${knots.length}`);
    }
    
    // 检查非递减性
    for (let i = 0; i < knots.length - 1; i++) {
        if (knots[i] > knots[i + 1] + knotTolerance) {
            throw new Error(`Knot 向量必须非递减：knots[${i}]=${knots[i]} > knots[${i+1}]=${knots[i+1]}`);
        }
    }
    
    // 首尾重节点(clamped)校验：位于 generateUniformKnotVector 生成阶段完成，
    // 此处(逐次求值的热路径)不做 console.warn 级别的重复校验，避免刷屏。
    // 合法 clamped 节点向量需 knots[degree]≈knots[0] 且 knots[len-1]≈knots[len-1-degree]。
}

/**
 * 验证控制点与权重合法性
 * @param {ControlPoint[]} controlPoints 
 * @throws {Error} 当控制点或权重不合法时抛出异常
 */
function validateControlPoints(controlPoints) {
    if (!Array.isArray(controlPoints) || controlPoints.length === 0) {
        throw new Error('控制点数组不能为空');
    }
    
    for (let i = 0; i < controlPoints.length; i++) {
        const pt = controlPoints[i];
        if (typeof pt.x !== 'number' || typeof pt.y !== 'number') {
            throw new Error(`控制点 [${i}] 坐标格式错误`);
        }
        // 权重 w 是可选的，如果不存在则默认为 1.0
        if (pt.w !== undefined && (typeof pt.w !== 'number' || pt.w <= 0)) {
            throw new Error(`控制点 [${i}] 权重必须为正数，当前 w=${pt.w}`);
        }
        if (!isFinite(pt.x) || !isFinite(pt.y) || (pt.w !== undefined && !isFinite(pt.w))) {
            throw new Error(`控制点 [${i}] 包含非有限数值 (NaN 或 Infinity)`);
        }
    }
}

/**
 * 验证 NURBS 曲线完整性
 * @param {NURBSCurve} curve 
 * @throws {Error} 当曲线数据不合法时抛出异常
 */
function validateCurve(curve) {
    validateControlPoints(curve.controlPoints);
    validateKnotVector(curve.knots, curve.degree, curve.controlPoints.length);
    
    // 验证参数域
    if (typeof curve.uMin !== 'number' || typeof curve.uMax !== 'number') {
        throw new Error('参数域 uMin/uMax 必须为数值');
    }
    if (curve.uMin >= curve.uMax) {
        throw new Error(`参数域无效：uMin(${curve.uMin}) >= uMax(${curve.uMax})`);
    }
}

// ============================================================================
// Knot 向量生成工具
// ============================================================================

/**
 * 生成均匀 Knot 向量
 * @param {number} n - 控制点数量
 * @param {number} p - 阶数 (degree)
 * @returns {number[]}
 */
function generateUniformKnotVector(n, p) {
    const m = n + p; // Knot 向量索引上限
    const knots = new Array(m + 1);
    
    // 首尾重节点 (clamped)
    for (let i = 0; i <= p; i++) {
        knots[i] = 0;
    }
    for (let i = m - p; i <= m; i++) {
        knots[i] = 1;
    }
    
    // 中间均匀分布（内部简单/单重节点数应为 n-p-1，避免覆盖尾部 clamp 起点 knots[m-p]，
    // 从而让尾部重节点数达到 p+1，曲线末端精确插值到末控制点）
    const internalKnotCount = n - p - 1;
    if (internalKnotCount > 0) {
        const step = 1 / (internalKnotCount + 1);
        for (let i = 1; i <= internalKnotCount; i++) {
            knots[p + i] = i * step;
        }
    }
    
    return knots;
}

/**
 * 生成准均匀 Knot 向量 (Open Uniform)
 * @param {number} n - 控制点数量
 * @param {number} p - 阶数 (degree)
 * @returns {number[]}
 */
function generateOpenUniformKnotVector(n, p) {
    // 准均匀 Knot 向量与均匀 Knot 向量在 clamped 情况下相同
    return generateUniformKnotVector(n, p);
}

/**
 * 生成分段贝塞尔 Knot 向量
 * @param {number} n - 控制点数量
 * @param {number} p - 阶数 (degree)
 * @returns {number[]}
 */
function generatePiecewiseBezierKnotVector(n, p) {
    const segments = n - p; // 贝塞尔段数
    if (segments <= 0) {
        throw new Error('控制点数量不足以生成分段贝兹曲线');
    }
    
    const knots = [];
    
    // 起始重节点
    for (let i = 0; i <= p; i++) {
        knots.push(0);
    }
    
    // 内部节点 (每个内节点重复度为 p)
    for (let i = 1; i < segments; i++) {
        const u = i / segments;
        for (let j = 0; j < p; j++) {
            knots.push(u);
        }
    }
    
    // 结束重节点
    for (let i = 0; i <= p; i++) {
        knots.push(1);
    }
    
    return knots;
}

/**
 * 生成 Knot 向量（统一接口）
 * @param {number} degree - 阶数
 * @param {number} numControlPoints - 控制点数量
 * @param {string} knotType - 类型：'clamped', 'uniform', 'bezier'
 * @returns {number[]}
 */
function generateKnotVector(degree, numControlPoints, knotType = 'clamped') {
    switch (knotType.toLowerCase()) {
        case 'uniform':
        case 'clamped':
            return generateUniformKnotVector(numControlPoints, degree);
        case 'bezier':
            // 方案 A：仅当控制点数 = 阶数+1（合法单段 Bezier）时使用真 Bezier 节点；
            // 多段分片 Bezier 会生成超长节点（不满足 nknots = n+p+1），在 WASM 长度校验
            // 处被静默丢弃 → 曲线沿用旧句柄（Bezier 选择失效）。故多段回退到 clamped 并提示。
            if (numControlPoints === degree + 1)
                return generatePiecewiseBezierKnotVector(numControlPoints, degree);
            console.warn(
                `Knot type "bezier" 需控制点数 = 阶数+1（当前 ${numControlPoints} ≠ ${degree}+1），已回退为 clamped`
            );
            return generateUniformKnotVector(numControlPoints, degree);
        default:
            console.warn(`未知的 Knot 类型 "${knotType}"，使用默认 clamped`);
            return generateUniformKnotVector(numControlPoints, degree);
    }
}

/**
 * 创建 NURBS 曲线对象
 * @param {Array<{x: number, y: number}>} controlPoints - 控制点数组
 * @param {number[]} weights - 权重数组
 * @param {number[]} knots - Knot 向量
 * @param {number} degree - 阶数
 * @returns {NURBSCurve}
 */
function createCurve(controlPoints, weights, knots, degree) {
    if (controlPoints.length !== weights.length) {
        throw new Error('控制点数量必须与权重数量一致');
    }
    
    // 计算参数域
    const uMin = knots[degree];
    const uMax = knots[knots.length - 1 - degree];
    
    // 构建带权重的控制点
    const weightedPoints = controlPoints.map((pt, i) => ({
        x: pt.x,
        y: pt.y,
        w: weights[i]
    }));
    
    return {
        controlPoints: weightedPoints,
        weights,
        knots,
        degree,
        uMin,
        uMax
    };
}

// ============================================================================
// Cox-de Boor 递归公式 - B-Spline 基函数
// ============================================================================

/**
 * 计算 B-Spline 基函数 N_{i,p}(u)
 * 使用 Cox-de Boor 递归公式
 * 
 * @param {number} i - 基函数索引
 * @param {number} p - 阶数 (degree)
 * @param {number} u - 参数值
 * @param {number[]} knots - Knot 向量
 * @returns {number} 基函数值
 */
function basisFunction(i, p, u, knots) {
    // 基础情形：p = 0
    if (p === 0) {
        if (u >= knots[i] && u < knots[i + 1]) {
            return 1.0;
        }
        // 处理最后一个节点的特殊情况
        if (i === knots.length - p - 2 && floatEqual(u, knots[knots.length - 1], knotTolerance)) {
            return 1.0;
        }
        return 0.0;
    }
    
    // 递归情形：p > 0
    let c1 = 0.0, c2 = 0.0;
    
    const denom1 = knots[i + p] - knots[i];
    if (!floatEqual(denom1, 0, knotTolerance)) {
        c1 = (u - knots[i]) / denom1 * basisFunction(i, p - 1, u, knots);
    }
    
    const denom2 = knots[i + p + 1] - knots[i + 1];
    if (!floatEqual(denom2, 0, knotTolerance)) {
        c2 = (knots[i + p + 1] - u) / denom2 * basisFunction(i + 1, p - 1, u, knots);
    }
    
    return c1 + c2;
}

/**
 * de Boor 三角法：一次求出参数 u 处全部 p+1 个非零 B-Spline 基函数。
 * 返回局部数组 N[0..p]，其中 N[i] = N_{span-p+i, p}(u)。
 * 相比递归 Cox-de Boor（无记忆化时最坏 O((p+1)·2^p)），本方法 O(p²) 且数值更稳定。
 * 参考 The NURBS Book Algorithm A2.2/A2.3 的 ndu 表第 p 列。
 *
 * @param {number} p - 阶数
 * @param {number} u - 参数值
 * @param {number[]} knots - Knot 向量
 * @param {number} span - findSpanIndex(p, u, knots) 的跨度
 * @returns {number[]} 局部非零基函数 N[0..p]
 */
function deBoorBasis(p, u, knots, span) {
    const N = new Array(p + 1).fill(0);
    const left = new Array(p + 1).fill(0);
    const right = new Array(p + 1).fill(0);
    N[0] = 1.0;

    for (let j = 1; j <= p; j++) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let saved = 0.0;
        for (let r = 0; r < j; r++) {
            const denom = right[r + 1] + left[j - r];
            // 若 denom≈0（重节点，该局部基函数贡献为 0），temp 记 0 即可
            const temp = floatEqual(denom, 0, knotTolerance) ? 0.0 : N[r] / denom;
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    return N;
}

/**
 * 计算 B-Spline 基函数及其一阶、二阶导数（全部非零项）。
 * 使用 The NURBS Book Algorithm A2.3 的标准递推，一次 O(p²) 同时得到
 * 基函数、一阶导数、二阶导数（取代原先死代码残留 + 公式可疑的实现）。
 * 
 * @param {number} p - 阶数 (degree)
 * @param {number} u - 参数值
 * @param {number[]} knots - Knot 向量
 * @param {number} [maxOrder=2] - 需要的最高导阶（1=只需 N 与一阶导，省去二阶导递推；2=含二阶导）
 * @returns {{N: number[], dN: number[], ddN: number[]}} 
 *          N: 基函数值数组 [N_0, N_1, ..., N_n]
 *          dN: 一阶导数数组（maxOrder>=1）
 *          ddN: 二阶导数数组（maxOrder>=2，否则全 0）
 */
function basisFunctionsWithDerivatives(p, u, knots, maxOrder = 2) {
    const n = knots.length - p - 2; // 控制点数量 - 1
    const N = new Array(n + 1).fill(0);
    const dN = new Array(n + 1).fill(0);
    const ddN = new Array(n + 1).fill(0);
    if (p === 0) {
        const span = findSpanIndex(0, u, knots);
        N[span] = 1.0;
        return { N, dN, ddN };
    }

    const span = findSpanIndex(p, u, knots);

    // ---- ndu 表：ndu[j][r] = N_{span-p+r, j}(u)，第 p 列即所需基函数 ----
    const ndu = [];
    for (let j = 0; j <= p; j++) ndu.push(new Array(p + 1).fill(0));
    const left = new Array(p + 1).fill(0);
    const right = new Array(p + 1).fill(0);
    ndu[0][0] = 1.0;

    for (let j = 1; j <= p; j++) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let saved = 0.0;
        for (let r = 0; r < j; r++) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const temp = floatEqual(ndu[j][r], 0, knotTolerance) ? 0.0 : ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    // 基函数：N_{span-p+r,p} = ndu[r][p]
    for (let r = 0; r <= p; r++) {
        N[span - p + r] = ndu[r][p];
    }

    // ---- 导数递推（A2.3）：ders[k][r] = N^{(k)}_{span-p+r,p}(u)，k=1..nd ----
    const nd = Math.max(1, Math.min(maxOrder, p)); // 需导阶数（至多 maxOrder，且不超 p）
    const ders = [];
    for (let k = 0; k <= nd; k++) ders.push(new Array(p + 1).fill(0));
    const a = [new Array(p + 1).fill(0), new Array(p + 1).fill(0)];

    for (let r = 0; r <= p; r++) {
        let s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (let k = 1; k <= nd; k++) {
            let d = 0.0;
            const rk = r - k;
            const pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            const j1 = (rk >= -1) ? 1 : -rk;
            const j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (let j = j1; j <= j2; j++) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            const t = s1; s1 = s2; s2 = t;
        }
    }

    // 乘以 p!/(p-k)! 因子
    let factor = p;
    for (let k = 1; k <= nd; k++) {
        for (let r = 0; r <= p; r++) ders[k][r] *= factor;
        factor *= (p - k);
    }

    for (let r = 0; r <= p; r++) {
        dN[span - p + r] = ders[1][r];
        if (nd >= 2) ddN[span - p + r] = ders[2][r];
    }

    return { N, dN, ddN };
}

/**
 * 查找参数 u 所在的 Knot 跨度索引
 * 返回 i 使得 knots[i] <= u < knots[i+1]
 * 
 * @param {number} p - 阶数
 * @param {number} u - 参数值
 * @param {number[]} knots - Knot 向量
 * @returns {number} 跨度索引
 */
function findSpanIndex(p, u, knots) {
    // 边界处理
    if (floatEqual(u, knots[knots.length - 1], knotTolerance)) {
        return knots.length - p - 2;
    }
    
    // 二分查找
    let low = p;
    let high = knots.length - 1 - p;
    
    while (low < high) {
        const mid = Math.floor((low + high + 1) / 2);
        if (u < knots[mid]) {
            high = mid - 1;
        } else {
            low = mid;
        }
    }
    
    return low;
}

// ============================================================================
// NURBS 曲线点计算
// ============================================================================

/**
 * 计算 NURBS 曲线上参数 u 处的点 C(u)
 * 使用齐次坐标投影方法
 * 
 * @param {NURBSCurve} curve - NURBS 曲线
 * @param {number} u - 参数值 (必须在 [uMin, uMax] 范围内)
 * @returns {{x: number, y: number}} 笛卡尔坐标点
 * @throws {Error} 当参数超出范围或计算出现 NaN 时抛出异常
 */
function evaluatePoint(curve, u) {
    validateCurve(curve);
    
    // 参数范围限制
    u = Math.max(curve.uMin, Math.min(curve.uMax, u));
    
    const { controlPoints, knots, degree } = curve;
    const span = findSpanIndex(degree, u, knots);
    // de Boor 三角一次求出 p+1 个非零基函数（O(p²)，取代逐控制点递归 O((p+1)·2^p)）
    const N = deBoorBasis(degree, u, knots, span);
    
    // 计算齐次坐标下的点
    let hx = 0.0, hy = 0.0, hw = 0.0;
    
    for (let i = 0; i <= degree; i++) {
        const idx = span - degree + i;
        if (idx >= 0 && idx < controlPoints.length) {
            const Nv = N[i];
            const w = controlPoints[idx].w;
            hx += Nv * w * controlPoints[idx].x;
            hy += Nv * w * controlPoints[idx].y;
            hw += Nv * w;
        }
    }
    
    // 投影到笛卡尔坐标
    if (floatEqual(hw, 0, EPSILON)) {
        // 参数域上界端点：clamped 曲线在 u==uMax 时端点即末控制点。
        // 此处基函数因区间半开判定（u < knots[i+1]）在精确右端点返回 0 而使 hw 归零，
        // 属数值端点情况，直接返回末控制点即可（齐次投影后仍为 P_last）。
        if (floatEqual(u, curve.uMax, knotTolerance)) {
            const last = controlPoints[controlPoints.length - 1];
            return { x: last.x, y: last.y };
        }
        throw new Error(`权重和为零，无法投影到笛卡尔坐标 (u=${u})`);
    }
    
    const x = hx / hw;
    const y = hy / hw;
    
    // 数值稳定性检查
    if (!isFinite(x) || !isFinite(y)) {
        throw new Error(`计算结果包含非有限数值 (u=${u}): x=${x}, y=${y}`);
    }
    
    return { x, y };
}

/**
 * 批量计算 NURBS 曲线上的多个点
 * @param {NURBSCurve} curve 
 * @param {number[]} uValues - 参数值数组
 * @returns {Array<{x: number, y: number}>}
 */
function evaluatePoints(curve, uValues) {
    return uValues.map(u => evaluatePoint(curve, u));
}

// ============================================================================
// NURBS 曲线导数计算
// ============================================================================

/**
 * 计算 NURBS 曲线在参数 u 处的一阶导数 C'(u)
 * 
 * @param {NURBSCurve} curve 
 * @param {number} u 
 * @returns {{x: number, y: number}} 一阶导数向量
 */
function evaluateFirstDerivative(curve, u) {
    validateCurve(curve);
    u = Math.max(curve.uMin, Math.min(curve.uMax, u));

    // 参数域上界端点：基函数在精确右边界(半开区间 [u_i,u_{i+1}))判定为 0 而使求导退化，
    // 轻微内缩到域内一点以正确求得左侧导数。
    if (floatEqual(u, curve.uMax, knotTolerance)) {
        u = curve.uMax - 1e-9;
    }

    const { controlPoints, knots, degree } = curve;
    const span = findSpanIndex(degree, u, knots);
    // 解析一阶导数（取代原先对每项基函数的前向有限差分 h=1e-6：既消除 O(h)≈1e-6
    // 的导数误差，又省去一半的基函数求值开销）。
    // 用 B 样条导数递推：N'_{i,p} = p/(U[i+p]-U[i])·N_{i,p-1} - p/(U[i+p+1]-U[i+1])·N_{i+1,p-1}，
    // 只需两个一维 de Boor 三角形（阶 p 与阶 p-1），避免低阶下 A2.3 二维 ndu 表的分配开销。
    const N = deBoorBasis(degree, u, knots, span);
    let Nm1 = null;
    if (degree >= 1) {
        const s2 = findSpanIndex(degree - 1, u, knots);
        Nm1 = deBoorBasis(degree - 1, u, knots, s2);
        const base2 = s2 - (degree - 1);
        // 一阶导数（局部数组 dN[k] 对应控制点 span-degree+k）
        const dN = new Array(degree + 1).fill(0);
        for (let k = 0; k <= degree; k++) {
            const idx = span - degree + k; // 控制点索引：对应 N_{idx, degree}
            const j1 = idx - base2;        // N_{idx, degree-1} 在 Nm1 中的位置
            const j2 = idx + 1 - base2;    // N_{idx+1, degree-1} 在 Nm1 中的位置
            const N1 = (j1 >= 0 && j1 <= degree - 1) ? Nm1[j1] : 0;
            const N2 = (j2 >= 0 && j2 <= degree - 1) ? Nm1[j2] : 0;
            const den1 = knots[idx + degree] - knots[idx];
            const den2 = knots[idx + degree + 1] - knots[idx + 1];
            dN[k] = degree * ((floatEqual(den1, 0, knotTolerance) ? 0 : N1 / den1)
                            - (floatEqual(den2, 0, knotTolerance) ? 0 : N2 / den2));
        }

        // 齐次坐标下的点与导数
        let hx = 0, hy = 0, hw = 0, dhx = 0, dhy = 0, dhw = 0;
        for (let k = 0; k <= degree; k++) {
            const idx = span - degree + k;
            if (idx >= 0 && idx < controlPoints.length) {
                const Nk = N[k], dNk = dN[k];
                const w = controlPoints[idx].w;
                const px = controlPoints[idx].x, py = controlPoints[idx].y;
                hx += Nk * w * px; hy += Nk * w * py; hw += Nk * w;
                dhx += dNk * w * px; dhy += dNk * w * py; dhw += dNk * w;
            }
        }

        if (floatEqual(hw, 0, EPSILON)) {
            throw new Error(`一阶导数计算：权重和为零 (u=${u})`);
        }
        // 商法则求导：C' = (H'·w - H·w') / w²
        const w2 = hw * hw;
        const dx = (dhx * hw - hx * dhw) / w2;
        const dy = (dhy * hw - hy * dhw) / w2;
        if (!isFinite(dx) || !isFinite(dy)) {
            throw new Error(`一阶导数计算包含非有限数值 (u=${u})`);
        }
        return { x: dx, y: dy };
    }

    // degree === 0：曲线退化为单个控制点，导数为零向量
    return { x: 0, y: 0 };
}

/**
 * 计算 NURBS 曲线在参数 u 处的二阶导数 C''(u)
 * 
 * @param {NURBSCurve} curve 
 * @param {number} u 
 * @returns {{x: number, y: number}} 二阶导数向量
 */
function evaluateSecondDerivative(curve, u) {
    validateCurve(curve);
    u = Math.max(curve.uMin, Math.min(curve.uMax, u));
    
    const { controlPoints, knots, degree } = curve;
    
    // 获取基函数及其导数
    const { N, dN, ddN } = basisFunctionsWithDerivatives(degree, u, knots);
    
    // 计算齐次坐标下的点、一阶导、二阶导
    let hx = 0.0, hy = 0.0, hw = 0.0;
    let dhx = 0.0, dhy = 0.0, dhw = 0.0;
    let ddhx = 0.0, ddhy = 0.0, ddhw = 0.0;
    
    for (let i = 0; i < controlPoints.length; i++) {
        const Ni = N[i];
        const dNi = dN[i];
        const ddNi = ddN[i];
        const w = controlPoints[i].w;
        const px = controlPoints[i].x;
        const py = controlPoints[i].y;
        
        hx += Ni * w * px;
        hy += Ni * w * py;
        hw += Ni * w;
        
        dhx += dNi * w * px;
        dhy += dNi * w * py;
        dhw += dNi * w;
        
        ddhx += ddNi * w * px;
        ddhy += ddNi * w * py;
        ddhw += ddNi * w;
    }
    
    // 二阶商法则求导
    const w2 = hw * hw;
    const w3 = w2 * hw;
    
    const x = hx / hw;
    const y = hy / hw;
    
    const dx = (dhx * hw - hx * dhw) / w2;
    const dy = (dhy * hw - hy * dhw) / w2;
    
    const ddx = (ddhx * w2 - 2 * dhx * dhw * hw - hx * ddhw * hw + 2 * hx * dhw * dhw) / w3;
    const ddy = (ddhy * w2 - 2 * dhy * dhw * hw - hy * ddhw * hw + 2 * hy * dhw * dhw) / w3;
    
    if (!isFinite(ddx) || !isFinite(ddy)) {
        throw new Error(`二阶导数计算包含非有限数值 (u=${u})`);
    }
    
    return { x: ddx, y: ddy };
}

/**
 * 计算 NURBS 曲线在参数 u 处的切向量（单位化）
 * @param {NURBSCurve} curve 
 * @param {number} u 
 * @returns {{x: number, y: number}} 单位切向量
 */
function evaluateTangent(curve, u) {
    const deriv = evaluateFirstDerivative(curve, u);
    const len = Math.sqrt(deriv.x * deriv.x + deriv.y * deriv.y);
    
    if (len < EPSILON) {
        throw new Error(`切向量为零向量，无法单位化 (u=${u})`);
    }
    
    return {
        x: deriv.x / len,
        y: deriv.y / len
    };
}

/**
 * 计算 NURBS 曲线在参数 u 处的曲率
 * @param {NURBSCurve} curve 
 * @param {number} u 
 * @returns {number} 曲率值
 */
function evaluateCurvature(curve, u) {
    const d1 = evaluateFirstDerivative(curve, u);
    const d2 = evaluateSecondDerivative(curve, u);
    
    const d1Len = Math.sqrt(d1.x * d1.x + d1.y * d1.y);
    if (d1Len < EPSILON) {
        return 0; // 奇异点，曲率未定义
    }
    
    // 2D 曲率公式：κ = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
    const cross = d1.x * d2.y - d1.y * d2.x;
    const curvature = Math.abs(cross) / Math.pow(d1Len, 3);
    
    return curvature;
}

// ============================================================================
// 包围盒 (AABB) 计算
// ============================================================================

/**
 * 计算 NURBS 曲线的轴对齐包围盒 (AABB)
 * 基于控制点的凸包性质
 * 
 * @param {NURBSCurve} curve 
 * @returns {{minX: number, minY: number, maxX: number, maxY: number}}
 */
function computeAABB(curve) {
    validateCurve(curve);
    
    const { controlPoints } = curve;
    
    // NURBS 曲线严格包含在其控制点凸包内（所有权重为正时），
    // 故用控制点原始坐标的凸包作为 AABB 是"保守"且数学上正确的：
    // 保证"若 AABB 不相交则曲线必不相交"（用于相交预筛选）。
    // 原实现用 x*w/minW 估计带权范围，权重差异大时会把包围盒错误缩小/偏移，
    // 导致本应相交的曲线被预筛选误判为分离。
    let minX = Infinity, minY = Infinity;
    let maxX = -Infinity, maxY = -Infinity;
    
    for (const pt of controlPoints) {
        minX = Math.min(minX, pt.x);
        minY = Math.min(minY, pt.y);
        maxX = Math.max(maxX, pt.x);
        maxY = Math.max(maxY, pt.y);
    }
    
    return {
        minX,
        minY,
        maxX,
        maxY
    };
}

/**
 * 判断两个 AABB 是否相交
 * @param {Object} box1 
 * @param {Object} box2 
 * @returns {boolean}
 */
function aabbIntersect(box1, box2) {
    return !(box1.maxX < box2.minX || 
             box1.minX > box2.maxX || 
             box1.maxY < box2.minY || 
             box1.minY > box2.maxY);
}

// ============================================================================
// 导出接口
// ============================================================================

export {
    // 常量
    EPSILON,
    MAX_ITERATIONS,
    knotTolerance,
    
    // 数据结构验证
    validateCurve,
    validateControlPoints,
    validateKnotVector,
    
    // Knot 向量生成
    generateKnotVector,
    generateUniformKnotVector,
    generateOpenUniformKnotVector,
    generatePiecewiseBezierKnotVector,
    
    // 曲线创建
    createCurve,
    
    // 基函数计算
    basisFunction,
    basisFunctionsWithDerivatives,
    findSpanIndex,
    
    // 曲线求值
    evaluatePoint,
    evaluatePoints,
    evaluateFirstDerivative,
    evaluateSecondDerivative,
    evaluateTangent,
    evaluateCurvature,
    
    // 几何属性
    computeAABB,
    aabbIntersect
};