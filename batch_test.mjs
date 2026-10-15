// =============================================================================
// NURBS 统一数值框架 —— 批量测试脚本（可复现，扩展覆盖版）
// -----------------------------------------------------------------------------
// 说明：
//   1) 使用固定种子 PRNG（mulberry32），全量输出可精确复现。
//   2) 相交组(INTERSECT)：两条 clamped 曲线共享同一端点控制点，端点精确插值
//      → 必然真实相交。验证 INTERSECTION 识别。
//   3) 相切组(TANGENT)：两条 clamped 曲线共享端点且在端点处切向量平行
//      （共享点 + 相邻控制点沿同一切线方向布置）→ 真实相切退化。
//      验证接触识别（INTERSECTION/TANGENT）及 TANGENT 细分识别率。
//   4) 高阶组(HIGH_ORDER)：更高次(4-5)/更多控制点的共享端点相交曲线对。
//   5) 不相交组(SEPARATED)：第二条曲线包围盒分离平移 → 真实不相交。
//   6) 近接触组(NEAR)：第二条曲线整体微小平移(间距~0.05)，包围盒重叠但曲线
//      分离 → 真实分离且距离很小，验证小间距场景不误判为相交。
//   7) 每个用例调用 unifiedSolve，按返回结果判定几何状态。
//   8) 结果写入 JSON，并在控制台输出统计摘要。
// 运行：node batch_test.mjs
// =============================================================================
import * as mc from './src/nurbs/math_core.js';
import * as sol from './src/solver/unified_solver.js';
import { writeFileSync } from 'fs';
import { fileURLToPath, pathToFileURL } from 'url';
import { dirname, join } from 'path';
const __dirname = dirname(fileURLToPath(import.meta.url));
// 仅当作为主入口直接运行时（node batch_test.mjs）才执行批量测试主流程；
// 被其他脚本 import 时只暴露生成函数/常量，供论文实验脚本复用同一批用例。
const isMain = process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href;

// ---- 可复现 PRNG：mulberry32 ----
function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
        a |= 0; a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
const SEED = 20261012;
const rnd = mulberry32(SEED);
const rand = (lo, hi) => lo + (hi - lo) * rnd();
const randInt = (lo, hi) => Math.floor(rand(lo, hi + 1));

// 各用例组数量
const N_INTERSECT   = 2000;  // 相交
const N_SEPARATED   = 1000;  // 不相交（包围盒分离）
const N_TANGENT     = 400;   // 相切退化
const N_HIGH_ORDER  = 400;   // 高阶曲线相交
const N_NEAR        = 400;   // 近接触分离

// ---- 随机生成一条 NURBS 曲线（默认次 2..3，控制点 4..7）----
function randomCurve(loN = 4, hiN = 7, loDeg = 2, hiDeg = 3) {
    const n = randInt(loN, hiN);
    const degree = randInt(loDeg, Math.min(hiDeg, n - 1));
    const knots = mc.generateUniformKnotVector(n, degree);
    const pts = [];
    const weights = [];
    let x = rand(-1, 1), y = rand(-1, 1);
    for (let i = 0; i < n; i++) {
        x += rand(-0.8, 0.8);
        y += rand(-0.8, 0.8);
        pts.push({ x, y });
        weights.push(rand(0.5, 1.5));
    }
    return mc.createCurve(pts, weights, knots, degree);
}

// ---- 通用：构造两条在 shared 点相接触的曲线 ----
// c1 首控制点 = shared；c2 末控制点 = shared → 两 clamped 曲线端点均精确过 shared 点
function makeSharedPair(opts = {}) {
    const { tangent = false, loN = 4, hiN = 7, loDeg = 2, hiDeg = 3 } = opts;
    const n1 = randInt(loN, hiN), d1 = randInt(loDeg, Math.min(hiDeg, n1 - 1));
    const n2 = randInt(loN, hiN), d2 = randInt(loDeg, Math.min(hiDeg, n2 - 1));
    const shared = { x: rand(-1.5, 1.5), y: rand(-1.5, 1.5) };
    // 相切方向（若需要）：随机单位方向
    const tAngle = rand(0, 2 * Math.PI);
    const dir = { x: Math.cos(tAngle), y: Math.sin(tAngle) };
    // c2 末端切线目标方向：相切组 = dir（与 c1 首端切线平行 → 真·相切）；
    // 横截组(相交/高阶) = dir 逆时针旋转 90°（与 c1 严格垂直 → 干净的横截相交）。
    // 关键：clamped 非有理曲线端点切线方向严格 = P1−P0 / Pn−1−Pn−2，
    // 故只要把对应控制点沿目标方向布置，就能精确控制端点切线，使两类角度分离。
    const c2Dir = tangent
        ? dir
        : { x: -dir.y, y: dir.x };
    const mk = (n, d, isFirstShared) => {
        const knots = mc.generateUniformKnotVector(n, d);
        const pts = [];
        const weights = [];
        let x = shared.x, y = shared.y;
        for (let i = 0; i < n; i++) {
            if (i === 0 && isFirstShared) pts.push({ ...shared });
            else if (i === n - 1 && !isFirstShared) pts.push({ ...shared });
            else if (i === 1 && isFirstShared) {
                // 第二个控制点始终沿 dir 布置 → c1 首端切线方向恒 = dir。
                // 这样共享端点处夹角只由 c2 决定：相切组 c2 沿 dir → 180°(平行)；
                // 横截组 c2 垂直 dir → 90°(严格横截)。两类在角度上干净分离。
                pts.push({ x: shared.x + dir.x * rand(0.3, 1.0), y: shared.y + dir.y * rand(0.3, 1.0) });
                x = pts[pts.length - 1].x; y = pts[pts.length - 1].y;
            } else if (i === n - 2 && !isFirstShared) {
                // 倒数第二控制点沿 c2 目标方向布置 → 端点切线 = c2Dir
                // （相切组 = dir 平行；横截组 = dir 垂直，保证共享端点处横截相交）
                pts.push({ x: shared.x + c2Dir.x * rand(0.3, 1.0), y: shared.y + c2Dir.y * rand(0.3, 1.0) });
                x = pts[pts.length - 1].x; y = pts[pts.length - 1].y;
            } else {
                x += rand(-0.9, 0.9); y += rand(-0.9, 0.9);
                pts.push({ x, y });
            }
            // 权重：非相切组随机（有理曲线，增加构造多样性与鲁棒性考验）；
            // 相切组权重全取 1 —— 只有非有理(权重全等)时端点切线才严格 = P1−P0 方向，
            // 才能构造出"端点切线精确平行于 dir"的**真·相切**。有理曲线权重不等时
            // 端点切线 = (w1·P1 − w0·P0)，即使控制点沿 dir 布置也会偏离 dir 达 ~10°，
            // 那在精确几何下是小角度横截相交而非相切。
            weights.push(tangent ? 1.0 : rand(0.5, 1.5));
        }
        return mc.createCurve(pts, weights, knots, d);
    };
    const c1 = mk(n1, d1, true);
    const c2 = mk(n2, d2, false);
    return { curve1: c1, curve2: c2, sharedPoint: shared, tangent };
}

// ---- 相交组：共享端点 ----
function makeIntersectingPair() {
    const p = makeSharedPair({ tangent: false });
    return { curve1: p.curve1, curve2: p.curve2, groundTruth: 'INTERSECTION', sharedPoint: p.sharedPoint };
}

// ---- 相切组：共享端点 + 端点切线平行 ----
function makeTangentPair() {
    const p = makeSharedPair({ tangent: true });
    return { curve1: p.curve1, curve2: p.curve2, groundTruth: 'TANGENT', sharedPoint: p.sharedPoint };
}

// ---- 高阶组：高次/多控制点相交 ----
function makeHighOrderPair() {
    const p = makeSharedPair({ tangent: false, loN: 7, hiN: 11, loDeg: 3, hiDeg: 5 });
    return { curve1: p.curve1, curve2: p.curve2, groundTruth: 'INTERSECTION', sharedPoint: p.sharedPoint };
}

// ---- 不相交组：包围盒分离平移 ----
function makeSeparatedPair() {
    const c1 = randomCurve();
    const c2 = randomCurve();
    const box1 = mc.computeAABB(c1);
    const box2 = mc.computeAABB(c2);
    const w1 = box1.maxX - box1.minX, h1 = box1.maxY - box1.minY;
    const shiftX = (w1 + (box2.maxX - box2.minX)) / 2 + 3 + rand(0, 2);
    const shiftY = (h1 + (box2.maxY - box2.minY)) / 2 + 3 + rand(0, 2);
    const horizontal = rnd() < 0.5;
    const dx = horizontal ? shiftX : rand(-0.5, 0.5);
    const dy = horizontal ? rand(-0.5, 0.5) : shiftY;
    const shifted = {
        ...c2,
        controlPoints: c2.controlPoints.map(p => ({ x: p.x + dx, y: p.y + dy, w: p.w }))
    };
    return { curve1: c1, curve2: shifted, groundTruth: 'SEPARATED', dx, dy };
}

// ---- 近接触组：间距 ~0.05，包围盒重叠但曲线分离 ----
// 构造要点：控制点 x 严格单调递增（→ B 样条曲线 x 严格单调，曲线自身无回折/自交，
// 平移复制后必不相交），c2 = c1 整体 +y 平移 offset → 真实分离、间距确定且很小。
function makeNearPair() {
    const n = randInt(4, 7);
    const degree = randInt(2, Math.min(3, n - 1));
    const knots = mc.generateUniformKnotVector(n, degree);
    const c1pts = [], c1weights = [];
    let x = rand(-1, 0.2), y = rand(-1, 1);
    for (let i = 0; i < n; i++) {
        x += rand(0.4, 0.9);          // x 严格递增
        y += rand(-0.8, 0.8);
        c1pts.push({ x, y });
        c1weights.push(rand(0.5, 1.5));
    }
    const c1 = mc.createCurve(c1pts, c1weights, knots, degree);
    const offset = rand(0.02, 0.06);  // 微小正间距
    const c2 = {
        ...c1,
        controlPoints: c1.controlPoints.map(p => ({ x: p.x, y: p.y + offset, w: p.w }))
    };
    return { curve1: c1, curve2: c2, groundTruth: 'SEPARATED', offset };
}

// ---- 执行单个求解，返回判定状态 ----
function classify(pair, options = {}) {
    try {
        const results = sol.unifiedSolve(pair.curve1, pair.curve2, options);
        if (!results || results.length === 0) return { state: 'ERROR', detail: 'empty results' };
        const best = results.reduce((a, b) => (a.distance < b.distance ? a : b));
        const r = {
            state: best.state,
            distance: best.distance,
            iterations: best.iterations,
            converged: best.converged,
            nResults: results.length
        };
        return r;
    } catch (e) {
        return { state: 'EXCEPTION', detail: e.message };
    }
}

// 求解器内部有大量 console.log，批量运行时屏蔽
const origLog = console.log.bind(console);
const origError = console.error.bind(console);
let verbose = process.env.VERBOSE === '1';
console.log = function (...a) { if (verbose) origLog(...a); };
console.error = origError;

// ---- 主流程（仅直接运行时执行；被 import 时跳过）----
if (isMain) {
console.log(`[BATCH] seed=${SEED} 相交=${N_INTERSECT} 相切=${N_TANGENT} 高阶=${N_HIGH_ORDER} 分离=${N_SEPARATED} 近接触=${N_NEAR}`);
console.log('[BATCH] 生成用例并执行求解（可能耗时数十秒）...');

const t0 = Date.now();
const intersectCases = [], separatedCases = [], tangentCases = [], highOrderCases = [], nearCases = [];
for (let i = 0; i < N_INTERSECT; i++) intersectCases.push(makeIntersectingPair());
for (let i = 0; i < N_TANGENT; i++) tangentCases.push(makeTangentPair());
for (let i = 0; i < N_HIGH_ORDER; i++) highOrderCases.push(makeHighOrderPair());
for (let i = 0; i < N_SEPARATED; i++) separatedCases.push(makeSeparatedPair());
for (let i = 0; i < N_NEAR; i++) nearCases.push(makeNearPair());

const results = [];

/**
 * @param {string} name
 * @param {Array} cases
 * @param {string[]} expectStates 期望状态集合（任一匹配即"标签正确"）
 * @param {boolean} expectIntersect 几何判据：true 期望 distance 小；false 期望 distance 大
 * @param {string|null} subMetric 附加细分指标名（如相切组报告 TANGENT 识别率）
 */
function runGroup(name, cases, expectStates, expectIntersect, subMetric = null) {
    let correct = 0, wrong = 0, exceptions = 0;
    let geoCorrect = 0, geoWrong = 0;
    let subHit = 0;
    const dists = [];
    const stateCounts = {};
    const geoErrorSamples = [];
    for (const pair of cases) {
        const r = classify(pair, { usePreFiltering: true, gridSize: 6, maxIterations: 80 });
        const state = r.state;
        stateCounts[state] = (stateCounts[state] || 0) + 1;
        if (r.distance !== undefined && isFinite(r.distance)) dists.push(r.distance);
        const isOk = expectStates.includes(state);
        if (r.state === 'EXCEPTION') exceptions++;
        else if (isOk) correct++;
        else wrong++;
        if (subMetric && state === subMetric) subHit++;
        // 几何真值：相交/相切/高阶 → distance 小；分离/近接触 → distance 大。
        // 阈值与求解器边界接触容差 BOUNDARY_INTERSECTION_TOLERANCE=1e-4 对齐：
        // 高阶曲线的端点（边界）交点残差可合法到 ~1e-5，若用更严的 1e-5 会误报。
        const geoOk = expectIntersect
            ? (r.distance !== undefined && isFinite(r.distance) && r.distance < 1e-4)
            : (r.distance !== undefined && isFinite(r.distance) && r.distance > 1e-3);
        if (r.state === 'EXCEPTION') { /* 已计异常 */ }
        else if (geoOk) geoCorrect++;
        else geoWrong++;
        if (!geoOk && r.state !== 'EXCEPTION' && geoErrorSamples.length < 8) {
            geoErrorSamples.push({ group: name, groundTruth: pair.groundTruth, predicted: state, distance: r.distance ?? null });
        }
        results.push({
            group: name, groundTruth: pair.groundTruth,
            predicted: state, distance: r.distance ?? null, iterations: r.iterations ?? null,
            nResults: r.nResults ?? null, error: r.detail ?? null
        });
    }
    dists.sort((a, b) => a - b);
    const nD = dists.length;
    const summary = {
        group: name, cases: cases.length,
        correct, wrong, exceptions,
        accuracy: cases.length ? (correct / cases.length) * 100 : 0,
        geoCorrect, geoWrong,
        geoAccuracy: cases.length ? (geoCorrect / cases.length) * 100 : 0,
        stateCounts,
        distMin: nD ? dists[0] : null,
        distP50: nD ? dists[Math.floor(nD / 2)] : null,
        distMax: nD ? dists[nD - 1] : null,
        avgIterations: cases.length ? (results.filter(r => r.group === name && r.iterations !== null).reduce((s, r) => s + r.iterations, 0) / cases.length) : 0,
        geoErrorSamples
    };
    if (subMetric) {
        summary.subMetric = subMetric;
        summary.subRate = cases.length ? (subHit / cases.length) * 100 : 0;
    }
    return summary;
}

const sInter      = runGroup('INTERSECT',   intersectCases,  ['INTERSECTION', 'TANGENT'], true);
const sTangent    = runGroup('TANGENT',     tangentCases,    ['INTERSECTION', 'TANGENT'], true, 'TANGENT');
const sHigh       = runGroup('HIGH_ORDER',  highOrderCases,  ['INTERSECTION', 'TANGENT'], true);
const sSep        = runGroup('SEPARATED',   separatedCases,  ['SEPARATED'], false);
const sNear       = runGroup('NEAR',        nearCases,       ['SEPARATED'], false);

const totalMs = Date.now() - t0;

// ---- 输出 ----
const LOG = (process.env.VERBOSE === '1') ? console.log : origLog;
LOG('');
LOG('================== 统计摘要 ==================');
origLog(`总耗时: ${(totalMs / 1000).toFixed(2)} s`);
for (const s of [sInter, sTangent, sHigh, sSep, sNear]) {
    LOG('---------------------------------------------');
    LOG(`组: ${s.group}   用例数: ${s.cases}`);
    LOG(`  按状态标签: 正确 ${s.correct}  误判 ${s.wrong}  异常 ${s.exceptions}  (准确率 ${s.accuracy.toFixed(2)}%)`);
    LOG(`  按几何距离: 正确 ${s.geoCorrect}  误判 ${s.geoWrong}  (准确率 ${s.geoAccuracy.toFixed(2)}%)`);
    LOG(`  状态分布: ${JSON.stringify(s.stateCounts)}`);
    if (s.subMetric) LOG(`  细分指标: ${s.subMetric} 识别率 ${s.subRate.toFixed(2)}%`);
    LOG(`  距离 min=${s.distMin === null ? 'N/A' : s.distMin.toExponential(3)}  p50=${s.distP50 === null ? 'N/A' : s.distP50.toExponential(3)}  max=${s.distMax === null ? 'N/A' : s.distMax.toExponential(3)}`);
    LOG(`  平均迭代: ${s.avgIterations.toFixed(1)}`);
}
LOG('=============================================');

// ---- 写入 JSON 结果文件 ----
const out = {
    meta: {
        seed: SEED, generatedAt: new Date().toISOString(),
        runtimeMs: totalMs,
        nIntersect: N_INTERSECT, nTangent: N_TANGENT, nHighOrder: N_HIGH_ORDER,
        nSeparated: N_SEPARATED, nNear: N_NEAR, total: results.length,
        solver: 'unifiedSolve (prefilter+grid6+LM80)',
        note: '相交/相切/高阶经共享端点构造保证真实接触；分离经包围盒分离平移；近接触经微小平移'
    },
    summaries: { intersect: sInter, tangent: sTangent, highOrder: sHigh, separated: sSep, near: sNear },
    sampleResults: results.slice(0, 20),
    allCounts: results.length
};
const DUMP_ALL = process.env.DUMP_ALL === '1';
const payload = DUMP_ALL ? results : results.slice(0, 20);
const jsonPath = join(__dirname, 'batch_results.json');
writeFileSync(jsonPath, JSON.stringify({ ...out, sampleResults: payload, allResults: DUMP_ALL ? results : undefined }, null, 2), 'utf8');
origLog(`[BATCH] 结果已写入 ${jsonPath}${DUMP_ALL ? '（全量）' : '（含全量统计+前20条样例，设置 DUMP_ALL=1 可输出全部）'}`);
} // end isMain

// ---- 供论文实验脚本复用的导出 ----
export {
    SEED, rnd, rand, randInt,
    N_INTERSECT, N_SEPARATED, N_TANGENT, N_HIGH_ORDER, N_NEAR,
    randomCurve, makeSharedPair, makeIntersectingPair, makeTangentPair,
    makeHighOrderPair, makeSeparatedPair, makeNearPair, classify,
    isMain
};
