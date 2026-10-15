// =============================================================================
// NURBS 统一数值框架 —— 对抗 / 退化 / 尺度鲁棒性实验（CAD 方向论文素材）
// -----------------------------------------------------------------------------
// 目的：构造"好解用例"无法区分的对抗场景，验证统一框架(LM+prefilter+多初值)
// 相对纯 Gauss-Newton 基线的鲁棒性优势，并检验已知风险点：
//   1) SHALLOW ：端点浅角横截相交（角度 6°/10°/12°/14°/18°，跨 ANGULAR_TOLERANCE=12°）
//                验证相切/相交分类在阈值附近的稳定性，以及 GN vs LM 的收敛差异
//   2) SCALE_1e3 ：整条曲线坐标放大 1000× → 检验容差尺度不变性（审稿人常问）
//   3) EXTREME_W ：极端权重(1e-2 ~ 1e2)的有理曲线 → 病态权重下鲁棒性
//   4) DEGENERATE：degree=1 线段化退化曲线 → 最简退化
//   5) NEAR_PARALLEL：两曲线内部近平行贴近（不交叉、Hessian 病态）→ GN 最易发散
// 每个族都跑 system(prefilter+LM) 与 gnNoPre 基线，报告族级正确率/收敛率/时间。
// 运行：node experiment_robustness.mjs   [env: CASES=200 缩小规模快速验证]
// =============================================================================
import * as mc from './src/nurbs/math_core.js';
import * as sol from './src/solver/unified_solver.js';
import { writeFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';
import { SEED, rnd, rand, randInt } from './batch_test.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const origLog = (...a) => process.stdout.write(a.join(' ') + '\n');
console.log = function () {};

const SCALE = parseInt(process.env.CASES || '0', 10);
const M = SCALE || 400; // 每族用例数

// ============================================================================
// 复用的求解入口（与 experiment_benchmark 相同）
// ============================================================================
const GEOM_TOL = 1e-8;
function clamp(x, lo, hi) { return Math.max(lo, Math.min(hi, x)); }
function gaussNewton(curve1, curve2, u0, v0, maxIter, tol) {
    let u = u0, v = v0;
    let obj = sol.objectiveFunction(curve1, curve2, u, v);
    let prevF = Infinity, iterations = 0, converged = false;
    for (let iter = 0; iter < maxIter; iter++) {
        iterations++;
        const H = obj.hessian, g = obj.gradient;
        const det = H[0][0] * H[1][1] - H[0][1] * H[1][0];
        if (Math.abs(det) >= 1e-14) {
            const invDet = 1 / det;
            const du = -0.5 * (H[1][1] * g[0] - H[0][1] * g[1]) * invDet;
            const dv = -0.5 * (-H[1][0] * g[0] + H[0][0] * g[1]) * invDet;
            u = clamp(u + du, curve1.uMin, curve1.uMax);
            v = clamp(v + dv, curve2.uMin, curve2.uMax);
        }
        obj = sol.objectiveFunction(curve1, curve2, u, v);
        const gradNorm = Math.hypot(obj.gradient[0], obj.gradient[1]);
        if (obj.distance < tol) { converged = true; break; }
        if (gradNorm < 1e-8) { converged = true; break; }
        if (prevF !== Infinity && prevF > 0 && obj.f >= prevF * (1 - 1e-10)) { converged = true; break; }
        prevF = obj.f;
    }
    return { u, v, distance: obj.distance, iterations, converged };
}
function gnSolve(curve1, curve2, options = {}) {
    const gridSize = options.gridSize || 6, maxIter = options.maxIterations || 50, tol = options.tolerance || GEOM_TOL;
    const guesses = [];
    for (let i = 0; i < gridSize; i++) for (let j = 0; j < gridSize; j++) guesses.push({ u: i / (gridSize - 1), v: j / (gridSize - 1) });
    for (const [cu, cv] of [[curve1.uMin, curve2.uMin], [curve1.uMin, curve2.uMax], [curve1.uMax, curve2.uMin], [curve1.uMax, curve2.uMax]])
        if (!guesses.some(g => Math.abs(g.u - cu) < 1e-9 && Math.abs(g.v - cv) < 1e-9)) guesses.push({ u: cu, v: cv });
    let best = null, bestDist = Infinity, bestConv = null, minConv = Infinity, nConvGuess = 0;
    for (const g of guesses) {
        const r = gaussNewton(curve1, curve2, g.u, g.v, maxIter, tol);
        if (r.converged) nConvGuess++;
        if (r.distance < bestDist) { bestDist = r.distance; best = r; }
        if (r.converged && r.distance < minConv) { minConv = r.distance; bestConv = r; }
    }
    const final = bestConv || best;
    return {
        state: !final.converged ? 'UNDETERMINED' : final.distance < tol ? 'INTERSECTION' : 'SEPARATED',
        distance: final.distance, iterations: final.iterations, converged: final.converged, nConvGuess
    };
}
const system = (pair) => {
    const rs = sol.unifiedSolve(pair.curve1, pair.curve2, { usePreFiltering: true, gridSize: 6, maxIterations: 80 });
    const best = rs.reduce((a, b) => a.distance < b.distance ? a : b);
    return { state: best.state, distance: best.distance, iterations: best.iterations, converged: best.converged };
};

// ============================================================================
// 对抗族生成器
// ============================================================================
// 可控角度浅角横截：两条 clamped 曲线共享端点，端点切线方向夹角 = angleDeg。
// 权重全取 1（非有理）以保证端点切线严格 = P1−P0 方向，从而精确控制夹角。
function makeShallowPair(angleDeg, opts = {}) {
    const { loN = 5, hiN = 8, loDeg = 2, hiDeg = 3 } = opts;
    const n1 = randInt(loN, hiN), d1 = randInt(loDeg, Math.min(hiDeg, n1 - 1));
    const n2 = randInt(loN, hiN), d2 = randInt(loDeg, Math.min(hiDeg, n2 - 1));
    const shared = { x: rand(-1.5, 1.5), y: rand(-1.5, 1.5) };
    const tAngle = rand(0, 2 * Math.PI);
    const dir = { x: Math.cos(tAngle), y: Math.sin(tAngle) };
    const th = angleDeg * Math.PI / 180;
    const c2Dir = { x: dir.x * Math.cos(th) - dir.y * Math.sin(th), y: dir.x * Math.sin(th) + dir.y * Math.cos(th) };
    const mk = (n, d, isFirst) => {
        const knots = mc.generateUniformKnotVector(n, d);
        const pts = [], weights = [];
        let x = shared.x, y = shared.y;
        for (let i = 0; i < n; i++) {
            if (i === 0 && isFirst) pts.push({ ...shared });
            else if (i === n - 1 && !isFirst) pts.push({ ...shared });
            else if (i === 1 && isFirst) pts.push({ x: shared.x + dir.x * rand(0.3, 1.0), y: shared.y + dir.y * rand(0.3, 1.0) });
            else if (i === n - 2 && !isFirst) pts.push({ x: shared.x + c2Dir.x * rand(0.3, 1.0), y: shared.y + c2Dir.y * rand(0.3, 1.0) });
            else { x += rand(-0.9, 0.9); y += rand(-0.9, 0.9); pts.push({ x, y }); }
            weights.push(1.0);
        }
        return mc.createCurve(pts, weights, knots, d);
    };
    return { curve1: mk(n1, d1, true), curve2: mk(n2, d2, false), angleDeg, groundTruth: 'INTERSECTION' };
}

// 随机曲线（可指定权重范围）
function makeCurve(weightLo, weightHi, opts = {}) {
    const { loN = 4, hiN = 7, loDeg = 2, hiDeg = 3 } = opts;
    const n = randInt(loN, hiN), degree = randInt(loDeg, Math.min(hiDeg, n - 1));
    const knots = mc.generateUniformKnotVector(n, degree);
    const pts = [], weights = [];
    let x = rand(-1, 1), y = rand(-1, 1);
    for (let i = 0; i < n; i++) {
        x += rand(-0.8, 0.8); y += rand(-0.8, 0.8);
        pts.push({ x, y }); weights.push(rand(weightLo, weightHi));
    }
    return mc.createCurve(pts, weights, knots, degree);
}
function makeExtremeWeightPair() {
    // 两条权重范围截然不同的有理曲线，制造病态权重；再整体平移让其一端共享 → 真实相交
    const c1 = makeCurve(1e-2, 1e2);
    const c2 = makeCurve(1e-2, 1e2);
    // 把 c2 首控制点移到 c1 首控制点 → 共享端点(权重非1时端点=首控制点*权重归一，仍精确过控制点)
    const p0 = c1.controlPoints[0];
    c2.controlPoints[0] = { ...p0, w: c2.controlPoints[0].w };
    return { curve1: c1, curve2: c2, groundTruth: 'INTERSECTION' };
}
function makeDegenerateLinePair() {
    // degree=1（线性）退化曲线，端点共享 → 线段在端点相交
    const mk = (n) => {
        const knots = mc.generateUniformKnotVector(n, 1);
        const pts = [], weights = [];
        let x = rand(-1, 1), y = rand(-1, 1);
        for (let i = 0; i < n; i++) { x += rand(-1, 1); y += rand(-1, 1); pts.push({ x, y }); weights.push(1.0); }
        return mc.createCurve(pts, weights, knots, 1);
    };
    const c1 = mk(randInt(3, 6));
    const c2 = mk(randInt(3, 6));
    c2.controlPoints[0] = { ...c1.controlPoints[0], w: c2.controlPoints[0].w };
    return { curve1: c1, curve2: c2, groundTruth: 'INTERSECTION' };
}
function makeNearParallelPair() {
    // 两曲线 x 单调递增（不自交/不互交），整体 y 平移 offset(3e-3~1e-2) →
    // 内部近平行贴近、Hessian 病态、真实分离（offset 取 >3e-3 使 >1e-3 判据公平）
    const n = randInt(5, 8), degree = randInt(2, Math.min(3, n - 1));
    const knots = mc.generateUniformKnotVector(n, degree);
    const pts = [], weights = [];
    let x = rand(-1, 0.2), y = rand(-1, 1);
    for (let i = 0; i < n; i++) { x += rand(0.4, 0.9); y += rand(-0.8, 0.8); pts.push({ x, y }); weights.push(1.0); }
    const c1 = mc.createCurve(pts, weights, knots, degree);
    const offset = rand(3e-3, 1e-2);
    const c2 = { ...c1, controlPoints: c1.controlPoints.map(p => ({ x: p.x, y: p.y + offset, w: p.w })) };
    return { curve1: c1, curve2: c2, groundTruth: 'SEPARATED', offset };
}
function makeCoincidentPair() {
    // 近重合曲线：c2 = c1 控制点整体 + 随下标增长的微小法向扰动 →
    // 起始端几乎重合（长段秩亏"山谷"），末端逐渐分离；真实分离、距离很小
    const n = randInt(5, 8), degree = randInt(2, Math.min(3, n - 1));
    const knots = mc.generateUniformKnotVector(n, degree);
    const pts = [], weights = [];
    let x = rand(-1, 0.2), y = rand(-1, 1);
    for (let i = 0; i < n; i++) { x += rand(0.4, 0.9); y += rand(-0.8, 0.8); pts.push({ x, y }); weights.push(1.0); }
    const c1 = mc.createCurve(pts, weights, knots, degree);
    // 法向扰动随 i 从 0 线性增长到 ~0.5 → 起始端近重合
    const grow = rand(0.3, 0.6);
    const c2 = { ...c1, controlPoints: c1.controlPoints.map((p, i) => ({ x: p.x, y: p.y + (i / (c1.controlPoints.length - 1)) * grow, w: p.w })) };
    return { curve1: c1, curve2: c2, groundTruth: 'SEPARATED', grow };
}

// ============================================================================
// 运行与统计
// ============================================================================
function runMethod(cases, expectIntersect, fn) {
    let ok = 0, n = 0, conv = 0, sumIter = 0;
    const times = [], dists = [];
    for (const pair of cases) {
        const t0 = performance.now();
        let r; try { r = fn(pair); } catch (e) { r = { distance: null, converged: false, iterations: null }; }
        times.push(performance.now() - t0); n++;
        if (r.converged) conv++;
        if (r.distance != null && isFinite(r.distance)) dists.push(r.distance);
        const dOk = r.distance != null && isFinite(r.distance) && (expectIntersect ? r.distance < 1e-4 : r.distance > 1e-3);
        if (dOk) ok++;
        if (r.iterations != null) sumIter += r.iterations;
    }
    times.sort((a, b) => a - b); dists.sort((a, b) => a - b);
    const pct = (a, p) => a.length ? a[Math.min(a.length - 1, Math.floor(a.length * p))] : null;
    return {
        n, ok, familyAcc: n ? (ok / n) * 100 : 0, convRate: n ? (conv / n) * 100 : 0,
        avgIter: n ? sumIter / n : 0, p50Time: pct(times, 0.5), p99Time: pct(times, 0.99),
        distP50: pct(dists, 0.5), distMax: pct(dists, 1.0)
    };
}
const fmt = (st) => ({
    family: st.family ?? '', method: st.method ?? '', detail: st.detail ?? '',
    n: st.n ?? '', familyAcc: (st.familyAcc ?? NaN).toFixed(2), convRate: (st.convRate ?? NaN).toFixed(1),
    avgIter: (st.avgIter ?? NaN).toFixed(1),
    p50Time: st.p50Time != null ? st.p50Time.toFixed(2) : 'NA', p99Time: st.p99Time != null ? st.p99Time.toFixed(2) : 'NA',
    distP50: st.distP50 != null ? st.distP50.toExponential(2) : 'NA', distMax: st.distMax != null ? st.distMax.toExponential(2) : 'NA'
});
const H = [
    { label: 'family', key: 'family' }, { label: 'detail', key: 'detail' }, { label: 'method', key: 'method' },
    { label: 'n', key: 'n' }, { label: 'familyAcc%', key: 'familyAcc' }, { label: 'convRate%', key: 'convRate' },
    { label: 'avgIter', key: 'avgIter' }, { label: 'p50Time(ms)', key: 'p50Time' }, { label: 'p99Time(ms)', key: 'p99Time' },
    { label: 'distP50', key: 'distP50' }, { label: 'distMax', key: 'distMax' }
];
const mdTable = (headers, rowOf) =>
    '| ' + headers.map(h => h.label).join(' | ') + ' |\n| ' + headers.map(() => '---').join(' | ') + ' |\n' +
    rowOf.map(r => '| ' + headers.map(h => r[h.key] ?? '').join(' | ') + ' |').join('\n');

// ============================================================================
// 主流程
// ============================================================================
const t0 = Date.now();
const rows = [];
const gen = (name, cases, expectIntersect, detail = '') => {
    for (const [mk, fn] of [['system', system], ['GN-baseline', (pair) => gnSolve(pair.curve1, pair.curve2, { gridSize: 6, maxIterations: 50 }) ]]) {
        const st = runMethod(cases, expectIntersect, fn);
        rows.push({ family: name, detail, method: mk, ...st });
    }
};

// 1) 浅角横截（跨 12° 阈值）—— 关注 system 的 TANGENT/INTERSECTION 状态切换是否稳定 + GN 收敛差异
for (const ang of [6, 10, 12, 14, 18]) {
    const cases = []; for (let i = 0; i < M; i++) cases.push(makeShallowPair(ang));
    gen('SHALLOW', cases, true, `angle=${ang}°`);
}
// 2) 尺度放大 1000×（从浅角 10° 用例放大）
{
    const cases = []; for (let i = 0; i < M; i++) {
        const p = makeShallowPair(10);
        const s = (c) => ({ ...c, controlPoints: c.controlPoints.map(q => ({ x: q.x * 1e3, y: q.y * 1e3, w: q.w })) });
        cases.push({ curve1: s(p.curve1), curve2: s(p.curve2), groundTruth: 'INTERSECTION' });
    }
    gen('SCALE_1e3', cases, true);
}
// 3) 极端权重
{
    const cases = []; for (let i = 0; i < M; i++) cases.push(makeExtremeWeightPair());
    gen('EXTREME_W', cases, true);
}
// 4) degree=1 退化线段
{
    const cases = []; for (let i = 0; i < M; i++) cases.push(makeDegenerateLinePair());
    gen('DEGENERATE', cases, true);
}
// 5) 内部近平行贴近（Hessian 病态，真实分离）
{
    const cases = []; for (let i = 0; i < M; i++) cases.push(makeNearParallelPair());
    gen('NEAR_PARALLEL', cases, false);
}
// 6) 近重合曲线（起始端重合点 → 真实接触；长段秩亏"山谷"考验步控，末端分离）
{
    const cases = []; for (let i = 0; i < M; i++) cases.push(makeCoincidentPair());
    gen('COINCIDENT', cases, true);
}

const totalMs = Date.now() - t0;

origLog('');
origLog('================== 论文实验 · 对抗/退化/尺度鲁棒性 ==================');
origLog(`SEED=${SEED}  M=${M}/族  总耗时=${(totalMs / 1000).toFixed(1)}s`);
origLog('');
origLog(mdTable(H, rows.map(fmt)));
origLog('==============================================================');

const json = { meta: { seed: SEED, generatedAt: new Date().toISOString(), runtimeMs: totalMs }, rows };
writeFileSync(join(__dirname, 'experiment_robustness.json'), JSON.stringify(json, null, 2), 'utf8');
let md = '# NURBS 统一数值框架 —— 对抗/退化/尺度鲁棒性实验结果\n\n';
md += '**seed**: ' + SEED + '  **M/族**: ' + M + '\n\n';
md += mdTable(H, rows.map(fmt)) + '\n\n---\n*由 experiment_robustness.mjs 生成*\n';
writeFileSync(join(__dirname, 'experiment_robustness.md'), md, 'utf8');
origLog('已写入 experiment_robustness.json 与 experiment_robustness.md');
