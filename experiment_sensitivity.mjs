// =============================================================================
// NURBS 统一数值框架 —— 阈值灵敏度分析（CAD 论文素材：证明"魔法阈值"不敏感）
// -----------------------------------------------------------------------------
// 目的：回应审稿人"分类阈值(12°/1e-8/1e-4)为何是魔法数"的质疑，量化证明：
//   在合理范围内，族级正确率对阈值选择不敏感（鲁棒），而非刀刃式敏感。
// 关键做法：每个阈值先"一次生成、固定同一批用例"（同一 seed），再对该批用例
//   扫描阈值取值，彻底隔离阈值效应与随机性 → 表中差异纯由阈值引起。
// 扫描对象：
//   1) ANGULAR_TOLERANCE（相切角判据）    ：TANGENT + SHALLOW(10°) 族
//   2) GEOMETRY_TOLERANCE（距离/收敛容差） ：5 个标准族
//   3) BOUNDARY_INTERSECTION_TOLERANCE    ：HIGH_ORDER + TANGENT 族
// 运行：node experiment_sensitivity.mjs   [env: CASES=100 缩小快速验证]
// =============================================================================
import * as sol from './src/solver/unified_solver.js';
import { writeFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';
import { SEED, N_INTERSECT, N_SEPARATED, N_TANGENT, N_HIGH_ORDER, N_NEAR,
    makeIntersectingPair, makeTangentPair, makeHighOrderPair, makeSeparatedPair, makeNearPair, classify } from './batch_test.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const origLog = (...a) => process.stdout.write(a.join(' ') + '\n');
console.log = function () {};

const SCALE = parseInt(process.env.CASES || '0', 10);
const nInt = SCALE || N_INTERSECT, nSep = SCALE || N_SEPARATED;
const nTan = SCALE || N_TANGENT, nHigh = SCALE || N_HIGH_ORDER, nNear = SCALE || N_NEAR;

// ---- 固定用例集（一次生成，供各阈值扫描共用，隔离随机性）----
const I = [], T = [], H = [], S = [], N = [];
for (let i = 0; i < nInt; i++) I.push(makeIntersectingPair());
for (let i = 0; i < nTan; i++) T.push(makeTangentPair());
for (let i = 0; i < nHigh; i++) H.push(makeHighOrderPair());
for (let i = 0; i < nSep; i++) S.push(makeSeparatedPair());
for (let i = 0; i < nNear; i++) N.push(makeNearPair());

// ---- 统计：族级正确率 + 细分指标 ----
// expectIntersect: true → 期望 distance<1e-4；false → 期望 distance>1e-3（与 batch_test geo 判据一致）
function run(cases, expectIntersect, subState = null) {
    let ok = 0, sub = 0;
    for (const pair of cases) {
        const r = classify(pair, { usePreFiltering: true, gridSize: 6, maxIterations: 80 });
        const d = r.distance;
        const dOk = d != null && isFinite(d) && (expectIntersect ? d < 1e-4 : d > 1e-3);
        if (dOk) ok++;
        if (subState && r.state === subState) sub++;
    }
    return {
        n: cases.length,
        familyAcc: cases.length ? (ok / cases.length) * 100 : 0,
        subRate: subState ? (cases.length ? (sub / cases.length) * 100 : 0) : null
    };
}

// ---- 各阈值默认值（restore 用）----
const DEF = sol.thresholdDefaults();

// ---- 本地浅角共享端点对（权重全1，端点切线夹角精确=angleDeg）----
import * as mc from './src/nurbs/math_core.js';
import { rand, randInt } from './batch_test.mjs';
function makeShallow(angleDeg) {
    const n1 = randInt(5, 8), d1 = randInt(2, 3), n2 = randInt(5, 8), d2 = randInt(2, 3);
    const shared = { x: rand(-1.5, 1.5), y: rand(-1.5, 1.5) };
    const ta = rand(0, 2 * Math.PI), dir = { x: Math.cos(ta), y: Math.sin(ta) };
    const th = angleDeg * Math.PI / 180, c2Dir = { x: dir.x * Math.cos(th) - dir.y * Math.sin(th), y: dir.x * Math.sin(th) + dir.y * Math.cos(th) };
    const mk = (n, d, isFirst) => {
        const knots = mc.generateUniformKnotVector(n, d), pts = [], weights = [];
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
    return { curve1: mk(n1, d1, true), curve2: mk(n2, d2, false) };
}
const shallow10 = [];
for (let i = 0; i < nTan; i++) shallow10.push(makeShallow(10));

// ============================================================================
// 1) ANGULAR_TOLERANCE 扫描（弧度→角度）：TANGENT 族 + SHALLOW(10°) 族
// ============================================================================
const angSweep = [6, 8, 10, 12, 14, 16, 20].map(deg => {
    sol.setThresholds({ ANGULAR_TOLERANCE: deg * Math.PI / 180 });
    const tangent = run(T, true, 'TANGENT');
    const shallow = run(shallow10, true);
    return { angleDeg: deg, tangentAcc: tangent.familyAcc, tangentLabelRate: tangent.subRate, shallowAcc: shallow.familyAcc, n: tangent.n };
});
sol.setThresholds(DEF); // 恢复默认

// ============================================================================
// 2) GEOMETRY_TOLERANCE 扫描（距离/收敛容差，指数级）
// ============================================================================
const geoSweep = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10].map(v => {
    sol.setThresholds({ GEOMETRY_TOLERANCE: v, CONVERGENCE_THRESHOLD: v });
    const st = {
        tolerance: v, n: I.length + S.length,
        interAcc: run(I, true).familyAcc, sepAcc: run(S, false).familyAcc,
        nearAcc: run(N, false).familyAcc, highAcc: run(H, true).familyAcc
    };
    sol.setThresholds(DEF);
    return st;
});

// ============================================================================
// 3) BOUNDARY_INTERSECTION_TOLERANCE 扫描
// ============================================================================
const boundSweep = [1e-5, 1e-4, 1e-3].map(v => {
    sol.setThresholds({ BOUNDARY_INTERSECTION_TOLERANCE: v });
    const st = {
        boundaryTol: v, n: H.length,
        highAcc: run(H, true).familyAcc, tangentAcc: run(T, true).familyAcc
    };
    sol.setThresholds(DEF);
    return st;
});

// ============================================================================
// 输出
// ============================================================================
origLog('');
origLog('================== 论文实验 · 阈值灵敏度分析 ==================');
origLog(`SEED=${SEED}  （每阈值扫描共用同一批固定用例，差异纯由阈值引起）`);
origLog('');
origLog('## 表1：ANGULAR_TOLERANCE 扫描（TANGENT 族 n=' + T.length + '，SHALLOW10° n=' + shallow10.length + '）');
origLog('| angle° | TANGENT族Acc% | TANGENT识别率% | SHALLOW10°Acc% |');
origLog('| --- | --- | --- | --- |');
for (const r of angSweep) origLog(`| ${r.angleDeg} | ${r.tangentAcc.toFixed(2)} | ${(r.tangentLabelRate ?? NaN).toFixed(2)} | ${r.shallowAcc.toFixed(2)} |`);
origLog('');
origLog('## 表2：GEOMETRY_TOLERANCE 扫描（各族正确率%）');
origLog('| tolerance | INTERSECT | SEPARATED | NEAR | HIGH_ORDER |');
origLog('| --- | --- | --- | --- | --- |');
for (const r of geoSweep) origLog(`| ${r.tolerance} | ${r.interAcc.toFixed(2)} | ${r.sepAcc.toFixed(2)} | ${r.nearAcc.toFixed(2)} | ${r.highAcc.toFixed(2)} |`);
origLog('');
origLog('## 表3：BOUNDARY_INTERSECTION_TOLERANCE 扫描（正确率%）');
origLog('| boundaryTol | HIGH_ORDER | TANGENT |');
origLog('| --- | --- | --- |');
for (const r of boundSweep) origLog(`| ${r.boundaryTol} | ${r.highAcc.toFixed(2)} | ${r.tangentAcc.toFixed(2)} |`);
origLog('==============================================================');

const json = {
    meta: { seed: SEED, generatedAt: new Date().toISOString(), note: '每阈值扫描共用同一批固定用例，隔离随机性' },
    angular: angSweep, geometry: geoSweep, boundary: boundSweep
};
writeFileSync(join(__dirname, 'experiment_sensitivity.json'), JSON.stringify(json, null, 2), 'utf8');

let md = '# NURBS 统一数值框架 —— 阈值灵敏度分析\n\n';
md += '**seed**: ' + SEED + '  （同一批固定用例，差异纯由阈值引起）\n\n';
md += '## 表1：ANGULAR_TOLERANCE 扫描（TANGENT 族 / SHALLOW10° 族）\n\n| angle° | TANGENT族Acc% | TANGENT识别率% | SHALLOW10°Acc% |\n| --- | --- | --- | --- |\n';
for (const r of angSweep) md += `| ${r.angleDeg} | ${r.tangentAcc.toFixed(2)} | ${(r.tangentLabelRate ?? NaN).toFixed(2)} | ${r.shallowAcc.toFixed(2)} |\n`;
md += '\n## 表2：GEOMETRY_TOLERANCE 扫描\n\n| tolerance | INTERSECT | SEPARATED | NEAR | HIGH_ORDER |\n| --- | --- | --- | --- | --- |\n';
for (const r of geoSweep) md += `| ${r.tolerance} | ${r.interAcc.toFixed(2)} | ${r.sepAcc.toFixed(2)} | ${r.nearAcc.toFixed(2)} | ${r.highAcc.toFixed(2)} |\n`;
md += '\n## 表3：BOUNDARY_INTERSECTION_TOLERANCE 扫描\n\n| boundaryTol | HIGH_ORDER | TANGENT |\n| --- | --- | --- |\n';
for (const r of boundSweep) md += `| ${r.boundaryTol} | ${r.highAcc.toFixed(2)} | ${r.tangentAcc.toFixed(2)} |\n`;
md += '\n---\n*由 experiment_sensitivity.mjs 生成*\n';
writeFileSync(join(__dirname, 'experiment_sensitivity.md'), md, 'utf8');
origLog('已写入 experiment_sensitivity.json / experiment_sensitivity.md');
