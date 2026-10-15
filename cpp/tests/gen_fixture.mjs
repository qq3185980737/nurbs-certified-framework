// =============================================================================
// gen_fixture.mjs — 生成 JS↔C++ 等价性验证夹具
// -----------------------------------------------------------------------------
// 复用 batch_test.mjs 的同一 mulberry32(seed 20261012) 与用例生成器，
// 产出固定曲线 + JS 求值结果 + unifiedSolve/findAllIntersections 结果，
// 供 cpp/tests/test_equivalence.cpp 重算并逐位对比。
// 运行：node cpp/tests/gen_fixture.mjs
// 输出：cpp/tests/fixture.txt
// =============================================================================
import * as bt from '../../batch_test.mjs';
import * as mc from '../../src/nurbs/math_core.js';
import * as sol from '../../src/solver/unified_solver.js';
import { writeFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const OUT = join(__dirname, 'fixture.txt');

const STATE_INT = { INTERSECTION: 0, TANGENT: 1, SEPARATED: 2, COINCIDENT: 3, UNDETERMINED: 4 };
const N_EVAL = 33;            // 每条曲线 evaluatePoint 采样数（全量规模下控体积）
const SOLVE_OPTS = { usePreFiltering: true, gridSize: 6, maxIterations: 80 };
const FIND_OPTS = { tolerance: 1e-8 };
// 与 batch_test.mjs 主流程完全一致的全量规模
const G = { INTERSECT: 2000, TANGENT: 400, HIGH_ORDER: 400, SEPARATED: 1000, NEAR: 400 };
const N_FIND = 60;            // findAllIntersections 取前 60 对

// ---- 1) 先取前 200 个 rnd 值，校验 mulberry32 双实现一致性 ----
const rndParity = [];
for (let i = 0; i < 200; i++) rndParity.push(bt.rnd());

// ---- 2) 收集曲线 / 求值 / 求解 / 全交点 ----
const curves = [];
const evals = [];
const solves = [];
const finds = [];

function addCurve(c) {
  const id = curves.length;
  curves.push({
    id,
    ncp: c.controlPoints.length,
    deg: c.degree,
    uMin: c.uMin, uMax: c.uMax,
    pts: c.controlPoints.map(p => [p.x, p.y, p.w]),
    knots: c.knots
  });
  return id;
}
function addEvals(id, c) {
  for (let i = 0; i <= N_EVAL; i++) {
    const u = c.uMin + (c.uMax - c.uMin) * i / N_EVAL;
    let p;
    try { p = mc.evaluatePoint(c, u); } catch (e) { continue; }
    evals.push([id, u, p.x, p.y]);
  }
}
function addEvalsFine(id, c) {
  const N = 1000;
  for (let i = 0; i <= N; i++) {
    const u = c.uMin + (c.uMax - c.uMin) * i / N;
    let p;
    try { p = mc.evaluatePoint(c, u); } catch (e) { continue; }
    evals.push([id, u, p.x, p.y]);
  }
}
function addPair(c1, c2, doFind) {
  const id1 = addCurve(c1), id2 = addCurve(c2);
  addEvals(id1, c1);
  addEvals(id2, c2);
  // unifiedSolve
  const res = sol.unifiedSolve(c1, c2, SOLVE_OPTS);
  let best = res[0];
  for (const r of res) if (r.distance < best.distance) best = r;
  solves.push([id1, id2, STATE_INT[best.state], best.u, best.v, best.distance,
               best.iterations, best.converged ? 1 : 0, res.length]);
  // findAllIntersections
  if (doFind) {
    const roots = sol.findAllIntersections(c1, c2, FIND_OPTS);
    finds.push([id1, id2, roots.map(r => [r.u, r.v, r.distance])]);
  }
}

// 全量规模（与 batch_test 主流程一致，RNG 流确定性）
let pairIdx = 0;
const genSet = {
  INTERSECT: () => bt.makeIntersectingPair(),
  TANGENT: () => bt.makeTangentPair(),
  HIGH_ORDER: () => bt.makeHighOrderPair(),
  SEPARATED: () => bt.makeSeparatedPair(),
  NEAR: () => bt.makeNearPair()
};
for (const key of Object.keys(genSet)) {
  for (let i = 0; i < G[key]; i++) {
    const p = genSet[key]();
    addPair(p.curve1, p.curve2, pairIdx++ < N_FIND);
  }
}

// ---- 退化 / 极端几何用例（结构合法、几何退化，覆盖 CAD 鲁棒性）----
function clampedKnots(n, p) {
  const knots = [];
  for (let i = 0; i <= p; i++) knots.push(0);
  const mid = n - p - 1;
  for (let i = 1; i <= mid; i++) knots.push(i / (mid + 1));
  for (let i = 0; i <= p; i++) knots.push(1);
  return knots; // length = n+p+1
}
function mkCurve(pts, weights, degree) {
  return mc.createCurve(pts, weights, clampedKnots(pts.length, degree), degree);
}
const P = (x, y) => ({ x, y });

function addDegeneratePair(c1, c2, label) {
  // 结构合法但可能数值退化：统一用守卫，求解/求值失败则跳过该例
  let id1 = -1, id2 = -1;
  try {
    id1 = addCurve(c1); id2 = addCurve(c2);
    addEvals(id1, c1); addEvals(id2, c2);
    // 细密采样（含多重节点 u=0.5 等敏感点），用于定位 evaluatePoint 逐位差异
    addEvalsFine(id1, c1); addEvalsFine(id2, c2);
  } catch (e) { return; }
  try {
    const res = sol.unifiedSolve(c1, c2, SOLVE_OPTS);
    if (res.length === 0) return;
    let best = res[0];
    for (const r of res) if (r.distance < best.distance) best = r;
    solves.push([id1, id2, STATE_INT[best.state], best.u, best.v, best.distance,
                 best.iterations, best.converged ? 1 : 0, res.length]);
  } catch (e) {}
  try {
    const roots = sol.findAllIntersections(c1, c2, FIND_OPTS);
    finds.push([id1, id2, roots.map(r => [r.u, r.v, r.distance])]);
  } catch (e) {}
  return [id1, id2];
}

{
  // (a) 全重合点曲线（退化为单点）vs 自身 → 重合
  const c = mkCurve([P(2,3), P(2,3), P(2,3), P(2,3)], [1,1,1,1], 3);
  addDegeneratePair(c, c, 'point');
  // (b) 共线重叠段 → 重合/无法判定
  const collA = mkCurve([P(0,0), P(1,0), P(2,0), P(3,0)], [1,1,1,1], 3);
  const collB = mkCurve([P(0,0), P(1,0), P(2,0), P(3,0)], [1,1,1,1], 3);
  addDegeneratePair(collA, collB, 'collinear-overlap');
  // (c) 端点相接（A 终点 = B 起点）→ 边界交点
  const endA = mkCurve([P(0,0), P(2,1), P(4,0)], [1,1,1], 2);
  const endB = mkCurve([P(4,0), P(6,1), P(8,0)], [1,1,1], 2);
  addDegeneratePair(endA, endB, 'endpoint-touch');
  // (d) 折返自重叠曲线 vs 直线
  const fold = mkCurve([P(0,0), P(2,2), P(4,0), P(2,-2), P(0,0)], [1,1,1,1,1], 4);
  const line = mkCurve([P(-1,0), P(5,0)], [1,1], 1);
  addDegeneratePair(fold, line, 'fold-vs-line');
  // (e) 自交"8 字形"曲线 vs 自身偏移
  const fig8 = mkCurve([P(0,0), P(4,4), P(8,0), P(4,-4), P(0,0)], [1,1,1,1,1], 4);
  addDegeneratePair(fig8, fig8, 'self-intersect');
  // (f) 极端尺度：×1e6 与 ×1e-6
  const base = bt.makeIntersectingPair();
  const big1 = mkCurve(base.curve1.controlPoints.map(pt => P(pt.x * 1e6, pt.y * 1e6)), base.curve1.controlPoints.map(pt => pt.w), base.curve1.degree);
  const big2 = mkCurve(base.curve2.controlPoints.map(pt => P(pt.x * 1e6, pt.y * 1e6)), base.curve2.controlPoints.map(pt => pt.w), base.curve2.degree);
  addDegeneratePair(big1, big2, 'scale-1e6');
  const tiny1 = mkCurve(base.curve1.controlPoints.map(pt => P(pt.x * 1e-6, pt.y * 1e-6)), base.curve1.controlPoints.map(pt => pt.w), base.curve1.degree);
  const tiny2 = mkCurve(base.curve2.controlPoints.map(pt => P(pt.x * 1e-6, pt.y * 1e-6)), base.curve2.controlPoints.map(pt => pt.w), base.curve2.degree);
  addDegeneratePair(tiny1, tiny2, 'scale-1e-6');
  // (g) 高重内部节点（knot 0.5 重复 2 次，非均匀但合法）
  const n = 6, p = 3;
  const hk = [0, 0, 0, 0, 0.5, 0.5, 1, 1, 1, 1]; // length = n+p+1 = 10
  const sharpA = mc.createCurve([P(0,0), P(2,4), P(4,-2), P(6,4), P(8,0), P(10,2)], [1,1,1,1,1,1], hk, p);
  const sharpB = mc.createCurve([P(0,2), P(3,-1), P(5,3), P(7,-1), P(9,2), P(11,0)], [1,1,1,1,1,1], hk, p);
  addDegeneratePair(sharpA, sharpB, 'multiplicity-knot');
}

// ---- 3) 写出夹具（空格分隔文本，C++ 侧 std::ifstream 直接读） ----
const lines = [];
const L = (fmt, ...a) => {
  const parts = [];
  if (fmt) parts.push(String(fmt));
  for (const x of a) parts.push(typeof x === 'number' && Number.isFinite(x) ? String(x) : String(x));
  lines.push(parts.join(' '));
};
L('NURBS_FIXTURE');
L('RND', rndParity.length);
for (const v of rndParity) L('', v);
L('NC', curves.length);
for (const c of curves) {
  L('CURVE', c.id, c.ncp, c.deg, c.uMin, c.uMax, c.knots.length);
  for (const p of c.pts) L('', p[0], p[1], p[2]);
  L('K', ...c.knots);
}
L('NE', evals.length);
for (const e of evals) L('EVAL', e[0], e[1], e[2], e[3]);
L('NS', solves.length);
for (const s of solves) L('SOLVE', s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8]);
L('NF', finds.length);
for (const f of finds) {
  L('FIND', f[0], f[1], f[2].length);
  for (const r of f[2]) L('', r[0], r[1], r[2]);
}

const text = lines.join('\n') + '\n';
writeFileSync(OUT, text, 'utf8');
console.log(`[fixture] 写 ${OUT}`);
console.log(`[fixture] 曲线=${curves.length} 求值=${evals.length} 求解=${solves.length} 全交点=${finds.length} (${(text.length/1024).toFixed(1)} KB)`);
