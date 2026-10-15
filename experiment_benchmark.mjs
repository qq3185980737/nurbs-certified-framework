// =============================================================================
// NURBS 统一数值框架 —— 论文实验脚本（基线对照 + 消融 + 规模扩展）
// -----------------------------------------------------------------------------
// 目的：为期刊（CAD 方向）实验章节产出可复现数据。
//   方法1 system  : unifiedSolve 默认（prefilter + LM 阻尼 + 多初值）—— 论文提出的系统
//   方法2 lm_noPre: 同系统但关闭包围盒预筛选 → 隔离 prefilter 贡献（消融）
//   方法3 gn_noPre: 纯 Gauss-Newton（无 LM 阻尼、无预筛选），相同网格+角点初值 → 隔离 LM 阻尼贡献（基线）
// 所有方法在同一批用例（与 batch_test 相同 SEED=20261012 生成）上运行，保证可比。
// 输出：控制台 Markdown 表 + experiment_results.json + experiment_results.md
// 运行：node experiment_benchmark.mjs  [可设 env: CASES=400 缩小规模快速验证]
// =============================================================================
import * as mc from './src/nurbs/math_core.js';
import * as sol from './src/solver/unified_solver.js';
import { writeFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';
import {
    SEED, N_INTERSECT, N_SEPARATED, N_TANGENT, N_HIGH_ORDER, N_NEAR,
    makeIntersectingPair, makeTangentPair, makeHighOrderPair, makeSeparatedPair, makeNearPair,
    makeSharedPair
} from './batch_test.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));

// 求解器内部有大量 console.log，实验时屏蔽；输出改用 process.stdout.write 直写，
// 避免被 batch_test.mjs import 时的 console.log 覆盖污染（其模块级覆盖先于本模块体执行）。
const origLog = (...a) => process.stdout.write(a.join(' ') + '\n');
console.log = function (...a) {};

const SCALE = parseInt(process.env.CASES || '0', 10);
const nInt = SCALE || N_INTERSECT, nSep = SCALE || N_SEPARATED;
const nTan = SCALE || N_TANGENT, nHigh = SCALE || N_HIGH_ORDER, nNear = SCALE || N_NEAR;

// ============================================================================
// 基线求解器：纯 Gauss-Newton（无 LM 阻尼信赖域）
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

// 纯 GN 多初值搜索（与 unifiedSolve 关闭预筛选时相同的初值集：全网格 + 4 角点）
function gnSolve(curve1, curve2, options = {}) {
    const gridSize = options.gridSize || 6;
    const maxIter = options.maxIterations || 50;
    const tol = options.tolerance || GEOM_TOL;
    const guesses = [];
    for (let i = 0; i < gridSize; i++)
        for (let j = 0; j < gridSize; j++)
            guesses.push({ u: i / (gridSize - 1), v: j / (gridSize - 1) });
    const corners = [[curve1.uMin, curve2.uMin], [curve1.uMin, curve2.uMax],
                     [curve1.uMax, curve2.uMin], [curve1.uMax, curve2.uMax]];
    for (const [cu, cv] of corners)
        if (!guesses.some(g => Math.abs(g.u - cu) < 1e-9 && Math.abs(g.v - cv) < 1e-9))
            guesses.push({ u: cu, v: cv });

    let best = null, bestDist = Infinity, bestConv = null, minConv = Infinity;
    let nConvGuess = 0, totalIter = 0;
    for (const g of guesses) {
        const r = gaussNewton(curve1, curve2, g.u, g.v, maxIter, tol);
        totalIter += r.iterations;
        if (r.converged) nConvGuess++;
        if (r.distance < bestDist) { bestDist = r.distance; best = r; }
        if (r.converged && r.distance < minConv) { minConv = r.distance; bestConv = r; }
    }
    const final = bestConv || best;
    let state;
    if (!final.converged) state = 'UNDETERMINED';
    else if (final.distance < tol) state = 'INTERSECTION';
    else state = 'SEPARATED';
    return {
        state, distance: final.distance,
        iterations: bestConv ? Math.round(totalIter / guesses.length) : final.iterations,
        converged: final.converged, nConvGuess
    };
}

// ============================================================================
// 统一判定：按几何距离给"族级"正确性（与 batch_test 的 geo-correct 对齐）
// ============================================================================
function familyOk(distance, expectIntersect) {
    if (distance === null || distance === undefined || !isFinite(distance)) return false;
    return expectIntersect ? distance < 1e-4 : distance > 1e-3;
}

// ============================================================================
// 运行一个方法在某一族上的统计
// ============================================================================
function runMethod(cases, expectIntersect, fn) {
    let ok = 0, n = 0, conv = 0, sumIter = 0;
    const times = [], dists = [];
    for (const pair of cases) {
        const t0 = performance.now();
        let r;
        try { r = fn(pair); } catch (e) { r = { state: 'EXCEPTION', distance: null, iterations: null, converged: false }; }
        const dt = performance.now() - t0;
        times.push(dt);
        n++;
        if (r.converged) conv++;
        if (r.distance !== null && r.distance !== undefined && isFinite(r.distance)) dists.push(r.distance);
        if (familyOk(r.distance, expectIntersect)) ok++;
        if (r.iterations !== null && r.iterations !== undefined) sumIter += r.iterations;
    }
    times.sort((a, b) => a - b);
    dists.sort((a, b) => a - b);
    const pct = (arr, p) => arr.length ? arr[Math.min(arr.length - 1, Math.floor(arr.length * p))] : null;
    return {
        n, ok, familyAcc: n ? (ok / n) * 100 : 0,
        convRate: n ? (conv / n) * 100 : 0,
        avgIter: n ? sumIter / n : 0,
        p50Time: pct(times, 0.5), p99Time: pct(times, 0.99), maxTime: pct(times, 1.0),
        distP50: pct(dists, 0.5), distMax: pct(dists, 1.0)
    };
}

const system = (pair) => {
    const rs = sol.unifiedSolve(pair.curve1, pair.curve2, { usePreFiltering: true, gridSize: 6, maxIterations: 80 });
    const best = rs.reduce((a, b) => a.distance < b.distance ? a : b);
    return { state: best.state, distance: best.distance, iterations: best.iterations, converged: best.converged };
};
const lmNoPre = (pair) => {
    const rs = sol.unifiedSolve(pair.curve1, pair.curve2, { usePreFiltering: false, gridSize: 6, maxIterations: 80 });
    const best = rs.reduce((a, b) => a.distance < b.distance ? a : b);
    return { state: best.state, distance: best.distance, iterations: best.iterations, converged: best.converged };
};
const gnNoPre = (pair) => gnSolve(pair.curve1, pair.curve2, { gridSize: 6, maxIterations: 50 });

// ============================================================================
// 用例构造（与 batch_test 相同顺序，保证 RNG 状态一致 → 完全相同用例）
// ============================================================================
const t0 = Date.now();
const I = [], T = [], H = [], S = [], N = [];
for (let i = 0; i < nInt; i++) I.push(makeIntersectingPair());
for (let i = 0; i < nTan; i++) T.push(makeTangentPair());
for (let i = 0; i < nHigh; i++) H.push(makeHighOrderPair());
for (let i = 0; i < nSep; i++) S.push(makeSeparatedPair());
for (let i = 0; i < nNear; i++) N.push(makeNearPair());

const families = [
    { name: 'INTERSECT', cases: I, expectIntersect: true },
    { name: 'TANGENT', cases: T, expectIntersect: true },
    { name: 'HIGH_ORDER', cases: H, expectIntersect: true },
    { name: 'SEPARATED', cases: S, expectIntersect: false },
    { name: 'NEAR', cases: N, expectIntersect: false }
];
const methods = { system, lmNoPre, gnNoPre };
const methodLabel = { system: 'system(prefilter+LM)', lmNoPre: 'LM(no prefilter)', gnNoPre: 'GN(no prefilter, baseline)' };

const rows = [];   // 每族×每方法一行
const json = { meta: { seed: SEED, generatedAt: new Date().toISOString(), note: '同 batch_test 用例(SEED=20261012)' }, families: {} };

for (const fam of families) {
    const frow = [];
    for (const [mk, fn] of Object.entries(methods)) {
        const st = runMethod(fam.cases, fam.expectIntersect, fn);
        frow.push({ method: mk, ...st });
        rows.push({ family: fam.name, method: methodLabel[mk], ...st });
    }
    json.families[fam.name] = frow;
}

const totalMs = Date.now() - t0;

// ============================================================================
// 规模扩展（scaling）：固定相交族构造，变化 degree 与控制点数，测系统耗时
// ============================================================================
const scalingRows = [];
{
    const degrees = [2, 3, 5, 7];
    const sizes = [6, 10, 15];
    const M = parseInt(process.env.SCALING_CASES || '40', 10);
    for (const deg of degrees) {
        for (const n of sizes) {
            if (deg > n - 1) continue;
            const cases = [];
            for (let i = 0; i < M; i++)
                cases.push(makeSharedPair({ tangent: false, loN: n, hiN: n, loDeg: deg, hiDeg: deg }));
            const st = runMethod(cases, true, system);
            scalingRows.push({ degree: deg, nControlPoints: n, n: cases.length, avgIter: st.avgIter, p50Time: st.p50Time, maxTime: st.maxTime, familyAcc: st.familyAcc });
        }
    }
}
json.scaling = scalingRows;

// ============================================================================
// 输出：Markdown 表（可直接进论文）+ JSON
// ============================================================================
function fmt(st) {
    return {
        family: st.family ?? '', method: st.method ?? '',
        degree: st.degree ?? '', nControlPoints: st.nControlPoints ?? '',
        n: st.n ?? '', familyAcc: (st.familyAcc ?? NaN).toFixed(2),
        convRate: (st.convRate ?? NaN).toFixed(1), avgIter: (st.avgIter ?? NaN).toFixed(1),
        p50Time: st.p50Time != null ? st.p50Time.toFixed(2) : 'NA',
        p99Time: st.p99Time != null ? st.p99Time.toFixed(2) : 'NA',
        maxTime: st.maxTime != null ? st.maxTime.toFixed(2) : 'NA',
        distP50: st.distP50 != null ? st.distP50.toExponential(2) : 'NA',
        distMax: st.distMax != null ? st.distMax.toExponential(2) : 'NA'
    };
}
function mdTable(headers, rowOf) {
    let s = '| ' + headers.map(h => h.label).join(' | ') + ' |\n| ' + headers.map(() => '---').join(' | ') + ' |\n';
    for (const r of rowOf) s += '| ' + headers.map(h => r[h.key] ?? '').join(' | ') + ' |\n';
    return s;
}

const H1 = [
    { label: 'family', key: 'family' }, { label: 'method', key: 'method' }, { label: 'n', key: 'n' },
    { label: 'familyAcc%', key: 'familyAcc' }, { label: 'convRate%', key: 'convRate' }, { label: 'avgIter', key: 'avgIter' },
    { label: 'p50Time(ms)', key: 'p50Time' }, { label: 'p99Time(ms)', key: 'p99Time' },
    { label: 'distP50', key: 'distP50' }, { label: 'distMax', key: 'distMax' }
];
const H2 = [
    { label: 'degree', key: 'degree' }, { label: 'nControlPoints', key: 'nControlPoints' }, { label: 'n', key: 'n' },
    { label: 'avgIter', key: 'avgIter' }, { label: 'p50Time(ms)', key: 'p50Time' }, { label: 'maxTime(ms)', key: 'maxTime' }
];
const table1 = mdTable(H1, rows.map(fmt));
const table2 = mdTable(H2, scalingRows.map(r => fmt({ ...r, family: '', method: '' })));

origLog('');
origLog('================== 论文实验 · 基线对照与消融 ==================');
origLog(`SEED=${SEED}  总耗时=${(totalMs / 1000).toFixed(1)}s`);
origLog('');
origLog('## 表1：各方法在 5 类用例上的正确率（族级几何判据）');
origLog(table1);
origLog('');
origLog('## 表2：规模扩展（system，相交族）');
origLog(table2);
origLog('==============================================================');

writeFileSync(join(__dirname, 'experiment_results.json'), JSON.stringify({ ...json, runtimeMs: totalMs }, null, 2), 'utf8');
let md = '# NURBS 统一数值框架 —— 论文实验结果（CAD 方向）\n\n';
md += '**seed**: ' + SEED + '   **运行**: ' + new Date().toISOString() + '\n\n';
md += '## 表1：基线对照与消融\n\n' + table1;
md += '\n## 表2：规模扩展（system）\n\n' + table2;
md += '\n---\n*由 experiment_benchmark.mjs 生成，可复现：node experiment_benchmark.mjs*\n';
writeFileSync(join(__dirname, 'experiment_results.md'), md, 'utf8');
origLog('已写入 experiment_results.json 与 experiment_results.md');
