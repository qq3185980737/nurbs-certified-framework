// wasm_front_test.mjs — WASM 前端接入层（wasm_backend.js）端到端验证
// 用 data-URL 以 ESM 加载 wasm_backend（.js 在无 type:module 项目按 CJS 解析，故用 data URL 内联），
// 真调 cpp/wasm/nurbs.wasm，验证：适配器输出结构（evaluatePoint/computeAABB/unifiedSolve/findAll）
// 与数值（对齐 wasm_test 夹具 + 前端默认曲线 + 随机批量自洽）。
import { readFileSync } from 'node:fs';
import createNurbsModule from './cpp/wasm/nurbs.js';

const code = readFileSync(new URL('./src/wasm_backend.js', import.meta.url), 'utf-8');
const wb = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
const { initWasm, createCurve, evaluatePoint, evaluatePoints, computeAABB, unifiedSolve, findAllIntersections, GeometryState } = wb;

await initWasm(createNurbsModule);

const near = (a, b, tol) => Math.abs(a - b) < tol;
const checks = [];
function check(name, ok, detail = '') { checks.push({ name, ok, detail }); }

// ---- 1) 直线夹具：L1 y=x, L2 y=4-x → 交于 (2,2) ----
const makeLine = (x0, y0, x1, y1) => createCurve(
  [{ x: x0, y: y0 }, { x: x1, y: y1 }], [1, 1], [0, 0, 1, 1], 1);
const L1 = makeLine(0, 0, 4, 4), L2 = makeLine(0, 4, 4, 0);

const ev = evaluatePoint(L1, 0.5);
check('evaluatePoint L1@0.5=(2,2)', near(ev.x, 2, 1e-9) && near(ev.y, 2, 1e-9), `got (${ev.x},${ev.y})`);

const evs = evaluatePoints(L1, [0, 0.25, 0.5]);
check('evaluatePoints 返回结构', Array.isArray(evs) && evs.length === 3 && 'x' in evs[0], JSON.stringify(evs[0]));

const box = computeAABB(L1);
check('computeAABB L1 控制点盒', near(box.minX, 0, 1e-12) && near(box.maxX, 4, 1e-12), JSON.stringify(box));

// tolerance 语义 = LM 收敛阈值（非 state 判定；state 用 GEOMETRY_TOLERANCE=1e-8）。
// 与 wasm_test.mjs 对齐：tolerance=0 → 用 C++ 默认收敛阈值(1e-12) → 直线夹具应精确相交。
const sol = unifiedSolve(L1, L2, { usePreFiltering: true, gridSize: 6, maxIterations: 80, tolerance: 0 });
check('unifiedSolve 数组+结构', Array.isArray(sol) && sol.length >= 1, `len=${sol.length}`);
if (sol.length >= 1) {
  const r = sol[0];
  check('unifiedSolve state=INTERSECTION', r.state === GeometryState.INTERSECTION, r.state + ' d=' + r.distance);
  check('unifiedSolve u=v=0.5', near(r.u, 0.5, 1e-4) && near(r.v, 0.5, 1e-4), `u=${r.u} v=${r.v}`);
  check('unifiedSolve point1≈(2,2)', near(r.point1.x, 2, 1e-6) && near(r.point1.y, 2, 1e-6), JSON.stringify(r.point1));
  check('unifiedSolve 字段齐', r.distance != null && r.iterations != null && r.converged != null && 'point2' in r, JSON.stringify({ d: r.distance, it: r.iterations, conv: r.converged }));
}
// 显式 tolerance=1e-6（前端 UI 默认收敛容差量级）→ 距离应落在该容差内（LM 提前收敛，state 可能非相交）
const solUI = unifiedSolve(L1, L2, { usePreFiltering: true, gridSize: 6, maxIterations: 80, tolerance: 1e-6 });
check('unifiedSolve tolerance=1e-6 距离在容差内', solUI.length >= 1 && solUI[0].distance < 1e-6, `d=${solUI[0]?.distance} state=${solUI[0]?.state}`);

const roots = findAllIntersections(L1, L2, { tolerance: 1e-6 });
check('findAll 数组+结构', Array.isArray(roots) && roots.length >= 1, `len=${roots.length}`);
if (roots.length >= 1) {
  const r0 = roots[0];
  check('findAll 交点≈(2,2)', near(r0.point1.x, 2, 1e-6) && near(r0.point1.y, 2, 1e-6), JSON.stringify(r0.point1));
  check('findAll 字段齐', r0.u != null && r0.v != null && r0.distance != null && r0.iterations != null, JSON.stringify({ u: r0.u, v: r0.v, d: r0.distance }));
}

// ---- 2) 前端默认曲线（6 控制点，degree 3，clamped）----
function genKnots(n, p) { const k = []; for (let i = 0; i <= p; i++) k.push(0); for (let i = 1; i <= n - p - 1; i++) k.push(i / (n - p)); for (let i = 0; i <= p; i++) k.push(1); return k; }
const defPts1 = [{x:-3,y:2},{x:-2,y:3},{x:-1,y:1},{x:0,y:-1},{x:1,y:1},{x:2,y:-2}];
const defPts2 = [{x:-2,y:-2},{x:-1,y:1},{x:0,y:-1},{x:1,y:2},{x:2,y:0},{x:3,y:1}];
const C1 = createCurve(defPts1, new Array(6).fill(1), genKnots(6, 3), 3);
const C2 = createCurve(defPts2, new Array(6).fill(1), genKnots(6, 3), 3);
let okDefault = true, detailDefault = '';
try {
  const evC = evaluatePoint(C1, 0.5);
  if (!(evC.x > -1e3 && evC.x < 1e3)) okDefault = false;
  const solC = unifiedSolve(C1, C2, { usePreFiltering: true, gridSize: 6, maxIterations: 80, tolerance: 1e-6 });
  const rootsC = findAllIntersections(C1, C2, { tolerance: 1e-6 });
  detailDefault = `eval@0.5=(${evC.x.toFixed(3)},${evC.y.toFixed(3)}) solve_n=${solC.length} root_n=${rootsC.length} state0=${solC[0]?.state}`;
} catch (e) { okDefault = false; detailDefault = 'ERR ' + e.message; }
check('前端默认曲线（d3×6点）接入', okDefault, detailDefault);

// ---- 3) 随机 50 对自洽（对齐 wasm_test 的 mulberry32 曲线）----
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;var t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
const rnd = mulberry32(20261012);
let pairs = 0, bad = 0;
for (let t = 0; t < 50; t++) {
  const d = 1 + Math.floor(rnd() * 3);
  const n = d + 1 + Math.floor(rnd() * 8);
  const ptsA = [], ptsB = [], kA = [], kB = [];
  const knotOf = (nn, dd) => { const k = []; for (let i = 0; i <= dd; i++) k.push(0); for (let i = 1; i <= nn - dd - 1; i++) k.push(i / (nn - dd)); for (let i = 0; i <= dd; i++) k.push(1); return k; };
  for (let i = 0; i < n; i++) { ptsA.push({ x: rnd() * 10, y: rnd() * 10 }); ptsB.push({ x: rnd() * 10, y: rnd() * 10 }); }
  kA.push(...knotOf(n, d)); kB.push(...knotOf(n, d));
  const A = createCurve(ptsA, new Array(n).fill(1), kA, d);
  const B = createCurve(ptsB, new Array(n).fill(1), kB, d);
  pairs++;
  const solR = unifiedSolve(A, B, { usePreFiltering: true, gridSize: 6, maxIterations: 80, tolerance: 1e-6 });
  let badPair = false;
  for (const r of solR) if (r.state === GeometryState.INTERSECTION && !(r.distance < 1e-6)) badPair = true;
  if (badPair) bad++;
}
check('随机 50 对自洽', bad === 0 && pairs === 50, `${pairs} 对, 失败 ${bad}`);

// ---- 汇总 ----
let fail = 0;
for (const c of checks) { if (!c.ok) fail++; console.log((c.ok ? 'PASS' : 'FAIL') + '  ' + c.name + (c.detail ? '  [' + c.detail + ']' : '')); }
console.log('---');
console.log(`共 ${checks.length} 项，通过 ${checks.length - fail}，失败 ${fail}`);
process.exit(fail === 0 ? 0 : 1);
