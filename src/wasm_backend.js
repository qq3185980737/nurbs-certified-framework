// ============================================================================
// wasm_backend.js — 前端 WASM 数学后端（单源化接入层）
// ----------------------------------------------------------------------------
// 用 C++ 编译的 nurbs.wasm（cpp/wasm/nurbs.js + nurbs.wasm，MODULARIZE 工厂
// createNurbsModule）提供全部数值核，替换旧 JS 版 math_core + unified_solver。
// 暴露与旧数学核完全相同的顶层 API（createCurve / evaluatePoint / evaluatePoints
// / computeAABB / unifiedSolve / findAllIntersections / GeometryState），内部转发
// 到 WASM，从而让浏览器与 C++ 共享同一套算法（单源）。
//
// 环境无关：浏览器经 `<script src="cpp/wasm/nurbs.js">` 暴露全局 createNurbsModule；
// Node 测试可把工厂函数传入 initWasm(factory)。
// 由 build.mjs 内联进 index.html；本文件是唯一事实来源。
// ============================================================================

// ---- UI 工具：Knot 向量生成（纯 JS，非数值核；自 math_core.js 提取保留） ----
export function generateUniformKnotVector(n, p) {
    const m = n + p; // Knot 向量索引上限
    const knots = new Array(m + 1);
    // 首尾重节点 (clamped)
    for (let i = 0; i <= p; i++) knots[i] = 0;
    for (let i = m - p; i <= m; i++) knots[i] = 1;
    // 中间均匀分布（内部节点数 = n-p-1）
    const internalKnotCount = n - p - 1;
    if (internalKnotCount > 0) {
        const step = 1 / (internalKnotCount + 1);
        for (let i = 1; i <= internalKnotCount; i++) knots[p + i] = i * step;
    }
    return knots;
}

export function generatePiecewiseBezierKnotVector(n, p) {
    const segments = n - p; // 贝塞尔段数
    if (segments <= 0) throw new Error('控制点数量不足以生成分段贝兹曲线');
    const knots = [];
    for (let i = 0; i <= p; i++) knots.push(0);
    for (let i = 1; i < segments; i++) {
        const u = i / segments;
        for (let j = 0; j < p; j++) knots.push(u);
    }
    for (let i = 0; i <= p; i++) knots.push(1);
    return knots;
}

export function generateKnotVector(degree, numControlPoints, knotType = 'clamped') {
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

// ---- WASM 模块 / 曲线句柄管理 ----
let mod = null;
let wasmReadyPromise = null;
const curveCache = new WeakMap(); // curve -> {handle, knotsRef, weightsRef, ctrlRef, degree}

export function initWasm(factory) {
    if (wasmReadyPromise) return wasmReadyPromise;
    wasmReadyPromise = (async () => {
        const f = factory || (typeof globalThis !== 'undefined' ? globalThis.createNurbsModule : null);
        if (typeof f !== 'function')
            throw new Error('createNurbsModule 不可用：请先通过 <script src="cpp/wasm/nurbs.js"> 加载 WASM');
        mod = await f();
        if (!mod || typeof mod._nurbs_curve_create !== 'function')
            throw new Error('WASM 模块初始化失败：导出缺失');
        return mod;
    })();
    return wasmReadyPromise;
}

function needMod() {
    if (!mod) throw new Error('WASM 数学核尚未初始化：请在 init() 中先 await initWasm()');
    return mod;
}

// 曲线句柄（懒建；跟随 knots/weights/degree 引用变更自动重建，兼容 UI 原地替换 knots）
function curveHandle(curve) {
    const m = needMod();
    let entry = curveCache.get(curve);
    if (entry && entry.knotsRef === curve.knots && entry.weightsRef === curve.weights
        && entry.ctrlRef === curve.controlPoints && entry.degree === curve.degree) {
        return entry.handle;
    }
    // 数据合法性校验：长度匹配 + 全部控制点/knots 为有限数。
    // 编辑中间态（输入未完成/含 null/NaN）一律沿用旧句柄，不触发 create，
    // 从而杜绝 validateCurve 抛异常 → destroy 悬垂句柄 → OOB 的死循环。
    const n = curve.controlPoints.length;
    const nk = curve.knots.length;
    if (nk !== n + curve.degree + 1) return entry ? entry.handle : null;
    for (let i = 0; i < n; i++) {
        const p = curve.controlPoints[i];
        if (!p || !Number.isFinite(p.x) || !Number.isFinite(p.y) || !Number.isFinite(p.w))
            return entry ? entry.handle : null;
    }
    for (let i = 0; i < nk; i++) {
        if (!Number.isFinite(curve.knots[i])) return entry ? entry.handle : null;
    }
    // 先建新句柄，成功才销毁旧句柄并更新缓存（先建后删，无悬垂窗口）
    const cp = m._malloc(n * 3 * 8), kn = m._malloc(nk * 8);
    const c = cp / 8, k = kn / 8;
    for (let i = 0; i < n; i++) {
        m.HEAPF64[c + i * 3] = curve.controlPoints[i].x;
        m.HEAPF64[c + i * 3 + 1] = curve.controlPoints[i].y;
        m.HEAPF64[c + i * 3 + 2] = curve.controlPoints[i].w;
    }
    for (let i = 0; i < nk; i++) m.HEAPF64[k + i] = curve.knots[i];
    const out = m._malloc(8);
    let rc;
    try {
        rc = m._nurbs_curve_create(cp, n, kn, nk, curve.degree, out);
    } catch (e) {
        m._free(cp); m._free(kn); m._free(out);
        return entry ? entry.handle : null; // create 抛异常：保留旧句柄
    }
    const handle = m.HEAPU32[out >> 2];
    m._free(cp); m._free(kn); m._free(out);
    if (rc !== 0) return entry ? entry.handle : null; // create 失败：保留旧句柄
    if (entry && entry.handle) m._nurbs_curve_destroy(entry.handle);
    curveCache.set(curve, { handle, knotsRef: curve.knots, weightsRef: curve.weights, ctrlRef: curve.controlPoints, degree: curve.degree });
    return handle;
}

// ---- 公共 API（与原 math_core + unified_solver 同名同结构） ----
export function createCurve(controlPoints, weights, knots, degree) {
    if (controlPoints.length !== weights.length)
        throw new Error('控制点数量必须与权重数量一致');
    const uMin = knots[degree];
    const uMax = knots[knots.length - 1 - degree];
    const weightedPoints = controlPoints.map((pt, i) => ({ x: pt.x, y: pt.y, w: weights[i] }));
    return { controlPoints: weightedPoints, weights, knots, degree, uMin, uMax };
}

export function evaluatePoint(curve, u) {
    const m = needMod();
    const h = curveHandle(curve);
    const out = m._malloc(16);
    const rc = m._nurbs_evaluate_point(h, u, out, out + 8);
    if (rc !== 0) { m._free(out); throw new Error('WASM 求值失败 rc=' + rc); }
    const p = { x: m.HEAPF64[out >> 3], y: m.HEAPF64[(out + 8) >> 3] };
    m._free(out);
    return p;
}

export function evaluatePoints(curve, uValues) {
    return uValues.map(u => evaluatePoint(curve, u));
}

export function computeAABB(curve) {
    // NURBS 位于控制点凸包内：用控制点包围盒（含首末点）覆盖整条曲线，供视图适配
    const pts = curve.controlPoints;
    if (pts.length === 0) return { minX: 0, minY: 0, maxX: 0, maxY: 0 };
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const p of pts) {
        if (p.x < minX) minX = p.x; if (p.x > maxX) maxX = p.x;
        if (p.y < minY) minY = p.y; if (p.y > maxY) maxY = p.y;
    }
    return { minX, minY, maxX, maxY };
}

const STATE_STR = ['INTERSECTION', 'TANGENT', 'SEPARATED', 'COINCIDENT', 'UNDETERMINED'];

export const GeometryState = {
    INTERSECTION: 'INTERSECTION',
    TANGENT: 'TANGENT',
    SEPARATED: 'SEPARATED',
    COINCIDENT: 'COINCIDENT',
    UNDETERMINED: 'UNDETERMINED'
};

// NurbsSolveResult 布局（C 结构，88B/条）：
//   state i32@0, pad, u@8, v@16, p1x@24, p1y@32, p2x@40, p2y@48, distance@56, iterations i32@64,
//   converged i32@68, residual@72, condition@80
function readSolveResult(m, base) {
    return {
        state: STATE_STR[m.HEAPU32[base >> 2]] || 'UNDETERMINED',
        u: m.HEAPF64[(base + 8) >> 3],
        v: m.HEAPF64[(base + 16) >> 3],
        point1: { x: m.HEAPF64[(base + 24) >> 3], y: m.HEAPF64[(base + 32) >> 3] },
        point2: { x: m.HEAPF64[(base + 40) >> 3], y: m.HEAPF64[(base + 48) >> 3] },
        distance: m.HEAPF64[(base + 56) >> 3],
        iterations: m.HEAPU32[(base + 64) >> 2],
        converged: m.HEAPU32[(base + 68) >> 2],
        residual: m.HEAPF64[(base + 72) >> 3],
        condition: m.HEAPF64[(base + 80) >> 3]
    };
}

// options: {usePreFiltering, gridSize, maxIterations, tolerance} —— 与 C API
// nurbs_unified_solve_ex 一一对应（equivalence 已证 JS/C++ 同参数逐位一致）
export function unifiedSolve(curve1, curve2, options = {}) {
    const m = needMod();
    const h1 = curveHandle(curve1), h2 = curveHandle(curve2);
    const usePreFiltering = options.usePreFiltering !== undefined ? (options.usePreFiltering ? 1 : 0) : 1;
    const gridSize = options.gridSize !== undefined ? options.gridSize : 6;
    const maxIterations = options.maxIterations !== undefined ? options.maxIterations : 80;
    const tolerance = options.tolerance !== undefined ? options.tolerance : 0;
    const resPtr = m._malloc(8), nPtr = m._malloc(4);
    const rc = m._nurbs_unified_solve_ex(h1, h2, usePreFiltering, gridSize, maxIterations,
        tolerance, resPtr, nPtr);
    const n = m.HEAPU32[nPtr >> 2], res = m.HEAPU32[resPtr >> 2];
    const results = [];
    if (rc === 0 && n >= 1) {
        for (let i = 0; i < n; i++) results.push(readSolveResult(m, res + i * 88));
    } else if (rc !== 0) {
        // rc<0：INVALID(-1)/RUNTIME(-2)。无交点/退化几何属正常无解情形，
        // 按空结果返回，不作为致命错误抛出（避免连续切换场景时误报）。
        m._free(resPtr); m._free(nPtr);
        console.warn('[solve] unifiedSolve rc=', rc, '按无解处理');
        return [];
    }
    if (n >= 1) m._nurbs_solve_results_free(res);
    m._free(resPtr); m._free(nPtr);
    return results;
}

// NurbsRoot 布局（C 结构，64B/条）：
//   u@0, v@8, distance@16, p1x@24, p1y@32, p2x@40, p2y@48, iterations i32@56
// C 侧 find_all 固定容差 1e-8（无 tolerance 入参）；前端 options.tolerance 不生效，以 C++ 为准
export function findAllIntersections(curve1, curve2, options = {}) {
    const m = needMod();
    const h1 = curveHandle(curve1), h2 = curveHandle(curve2);
    const rootsPtr = m._malloc(8), nPtr = m._malloc(4);
    const rc = m._nurbs_find_all_intersections(h1, h2, rootsPtr, nPtr);
    const n = m.HEAPU32[nPtr >> 2], roots = m.HEAPU32[rootsPtr >> 2];
    const out = [];
    if (rc === 0 && n >= 1) {
        for (let i = 0; i < n; i++) {
            const base = roots + i * 64;
            out.push({
                u: m.HEAPF64[base >> 3],
                v: m.HEAPF64[(base + 8) >> 3],
                distance: m.HEAPF64[(base + 16) >> 3],
                point1: { x: m.HEAPF64[(base + 24) >> 3], y: m.HEAPF64[(base + 32) >> 3] },
                point2: { x: m.HEAPF64[(base + 40) >> 3], y: m.HEAPF64[(base + 48) >> 3] },
                iterations: m.HEAPU32[(base + 56) >> 2]
            });
        }
    } else if (rc !== 0) {
        m._free(rootsPtr); m._free(nPtr);
        console.warn('[solve] findAllIntersections rc=', rc, '按无解处理');
        return [];
    }
    if (n >= 1) m._nurbs_roots_free(roots);
    m._free(rootsPtr); m._free(nPtr);
    return out;
}

// NurbsEnumerationResult 布局（C 结构，16B）：proven_count i32@0, tangent_count i32@4,
//   evaluated_boxes i32@8, complete i32@12
// 全局保证性枚举（区间细分 + Krawczyk 证明）：返回证明汇总，供前端显示 certified 状态。
export function enumerateIntersections(curve1, curve2, leafTol = 1e-4) {
    const m = needMod();
    const h1 = curveHandle(curve1), h2 = curveHandle(curve2);
    const out = m._malloc(16);
    const rc = m._nurbs_enumerate_intersections(h1, h2, leafTol, out);
    let res = null;
    if (rc === 0) {
        res = {
            provenCount: m.HEAPU32[out >> 2],
            tangentCount: m.HEAPU32[(out + 4) >> 2],
            evaluatedBoxes: m.HEAPU32[(out + 8) >> 2],
            complete: m.HEAPU32[(out + 12) >> 2] ? true : false
        };
    } else {
        console.warn('[solve] enumerateIntersections rc=', rc);
    }
    m._free(out);
    return res;
}
