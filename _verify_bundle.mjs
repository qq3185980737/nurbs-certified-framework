// _verify_bundle.mjs — 验证 build.mjs 产出的 bundle 与 src 模块完全等价（同输入同输出）
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import * as src from './src/solver/unified_solver.js';
import * as mc from './src/nurbs/math_core.js';
// 复用 batch_test 的种子随机生成器，保证可复现
import { rnd, randomCurve } from './batch_test.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));

// 加载 bundle 到全局（模拟浏览器 window）
globalThis.window = globalThis;
globalThis.document = { createElement: () => ({ getContext: () => null }), body: { appendChild: () => {} } };
const bundleCode = readFileSync(join(__dirname, 'dist', 'nurbs-bundle.js'), 'utf-8');
// 用 vm 在共享全局执行（IIFE 会把导出挂到 globalThis.window 即 globalThis）
import { createRequire } from 'node:module';
const require = createRequire(import.meta.url);
const vm = require('node:vm');
vm.runInThisContext(bundleCode, { filename: 'nurbs-bundle.js' });

const B = globalThis; // bundle 全局
let pass = 0, fail = 0;
const failures = [];

function check(name, cond, extra = '') {
  if (cond) { pass++; }
  else { fail++; failures.push(name + (extra ? ' :: ' + extra : '')); }
}

// 1) 关键全局已暴露
for (const g of ['unifiedSolve','findAllIntersections','evaluatePoint','evaluateFirstDerivative',
                 'createCurve','generateKnotVector','computeAABB','aabbIntersect','GeometryState',
                 'levenbergMarquardt','pieceBox','globalOptimization','setThresholds','thresholdDefaults',
                 'CanvasEngine','Viewport','COLORS']) {
  check(`global.${g}`, typeof B[g] === 'function' || typeof B[g] === 'object');
}
// 渲染器类型冒烟：可实例化（用 stub canvas 构造，仅验证绑定与全局存在）
check('CanvasEngine.prototype.render 存在', typeof B.CanvasEngine?.prototype?.render === 'function');
check('CanvasEngine.prototype.fitCurves 存在', typeof B.CanvasEngine?.prototype?.fitCurves === 'function');
check('COLORS.curve1 存在', typeof B.COLORS?.curve1 === 'string');

// 2) evaluatePoint 逐点等价（多曲线、多 u）
let epDiff = 0, epMax = 0;
for (let t = 0; t < 30; t++) {
  const c = randomCurve();
  for (let i = 0; i <= 8; i++) {
    const u = c.uMin + (c.uMax - c.uMin) * (i / 8);
    const a = mc.evaluatePoint(c, u), b = B.evaluatePoint(c, u);
    const d = Math.hypot(a.x - b.x, a.y - b.y);
    epMax = Math.max(epMax, d);
    if (d > 1e-12) epDiff++;
  }
}
check('evaluatePoint 等价(≤1e-12)', epDiff === 0, `diff=${epDiff}, max=${epMax.toExponential(2)}`);

// 3) unifiedSolve 等价：相同种子生成一批 pair，比对首解 state/u/v/distance
let usDiff = 0, usMax = 0;
let gsCount = 0;
for (let t = 0; t < 80; t++) {
  const c1 = randomCurve(), c2 = randomCurve();
  const o = { gridSize: 5, tolerance: 1e-8, maxIterations: 100, usePreFiltering: (t % 2 === 0) };
  const sa = src.unifiedSolve(c1, c2, o);
  const ba = B.unifiedSolve(c1, c2, o);
  if (!Array.isArray(sa) || !Array.isArray(ba)) { usDiff++; usMax = Infinity; continue; }
  const s0 = sa[0], b0 = ba[0];
  if (!s0 || !b0) { usDiff++; continue; }
  // 比较 状态 + 参数 + 距离
  const stateSame = (s0.state === b0.state);
  const dU = Math.abs(s0.u - b0.u), dV = Math.abs(s0.v - b0.v);
  const dDist = Math.abs(s0.distance - b0.distance);
  usMax = Math.max(usMax, dDist);
  if (!stateSame || dU > 1e-8 || dV > 1e-8 || dDist > 1e-10) {
    usDiff++;
    if (usDiff <= 5) {
      console.log(`  [diff] t=${t} pre=${o.usePreFiltering} src(state=${s0.state},u=${s0.u.toFixed(6)},v=${s0.v.toFixed(6)},d=${s0.distance.toExponential(3)}) vs bdl(state=${b0.state},u=${b0.u?.toFixed(6)},v=${b0.v?.toFixed(6)},d=${b0.distance?.toExponential(3)})`);
    }
  }
}
check('unifiedSolve 等价(状态+u/v+距离)', usDiff === 0, `diff=${usDiff}, maxDistDiff=${usMax.toExponential(2)}`);

// 4) findAllIntersections 等价
let fiDiff = 0, fiMax = 0;
for (let t = 0; t < 40; t++) {
  const c1 = randomCurve(), c2 = randomCurve();
  const o = { tolerance: 1e-8 };
  const sa = src.findAllIntersections(c1, c2, o);
  const ba = B.findAllIntersections(c1, c2, o);
  if (sa.length !== ba.length) { fiDiff++; continue; }
  for (let i = 0; i < sa.length; i++) {
    const d = Math.hypot(sa[i].u - ba[i].u, sa[i].v - ba[i].v);
    fiMax = Math.max(fiMax, d);
    if (d > 1e-8) fiDiff++;
  }
}
check('findAllIntersections 等价(个数+位置)', fiDiff === 0, `diff=${fiDiff}, max=${fiMax.toExponential(2)}`);

const summary =
  `==== 验证结果 ====\npass=${pass}, fail=${fail}\n` +
  `evaluatePoint diff=${epDiff}, max=${epMax.toExponential(2)}\n` +
  `unifiedSolve diff=${usDiff}, maxDistDiff=${usMax.toExponential(2)}\n` +
  `findAllIntersections diff=${fiDiff}, max=${fiMax.toExponential(2)}\n` +
  (failures.length ? '失败项：\n' + failures.map((f) => '  - ' + f).join('\n') : '');
writeFileSync(join(__dirname, '_verify_summary.txt'), summary, 'utf-8');
console.log(summary);
process.exit(fail ? 1 : 0);
