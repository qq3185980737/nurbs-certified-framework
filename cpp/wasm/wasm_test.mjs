// wasm_test.mjs — 验证 WASM 版 nurbs 算法库的导出与数值正确性
// 夹具与 cpp/examples/cad_integration.c 一致：
//   L1: (0,0)-(4,4) 直线 y=x；L2: (0,4)-(4,0) 直线 y=4-x → 交于 (2,2)
// 期望：evaluate L1@0.5=(2,2)；unified_solve(L1,L2) → INTERSECTION, u=v=0.5, distance≈0
import createNurbsModule from './nurbs.js';

function makeLine(m, x0, y0, x1, y1) {
  // ctrl: 2 点 (x,y,w=1)，knots:[0,0,1,1]，degree=1，nctrl=2, nknots=4
  const ctrlBytes = 2 * 3 * 8, knotBytes = 4 * 8;
  const cp = m._malloc(ctrlBytes), kn = m._malloc(knotBytes);
  const c = cp / 8, k = kn / 8;
  m.HEAPF64[c] = x0; m.HEAPF64[c+1] = y0; m.HEAPF64[c+2] = 1;
  m.HEAPF64[c+3] = x1; m.HEAPF64[c+4] = y1; m.HEAPF64[c+5] = 1;
  m.HEAPF64[k] = 0; m.HEAPF64[k+1] = 0; m.HEAPF64[k+2] = 1; m.HEAPF64[k+3] = 1;
  const out = m._malloc(8);
  const rc = m._nurbs_curve_create(cp, 2, kn, 4, 1, out);
  const handle = m.HEAPU32[out >> 2];
  m._free(cp); m._free(kn); m._free(out);
  return { handle, rc };
}

function evalPoint(m, curve, u) {
  const out = m._malloc(16);
  const rc = m._nurbs_evaluate_point(curve, u, out, out + 8);
  const x = m.HEAPF64[out >> 3], y = m.HEAPF64[(out + 8) >> 3];
  m._free(out);
  return { rc, x, y };
}

async function main() {
  const m = await createNurbsModule();
  // 检查导出函数存在
  const funcs = ['_nurbs_curve_create','_nurbs_curve_destroy','_nurbs_evaluate_point',
    '_nurbs_evaluate_first_derivative','_nurbs_evaluate_second_derivative',
    '_nurbs_unified_solve','_nurbs_unified_solve_ex','_nurbs_solve_results_free',
    '_nurbs_find_all_intersections','_nurbs_roots_free',
    '_nurbs_bvh_index_create','_nurbs_bvh_index_destroy','_nurbs_bvh_solve'];
  const missing = funcs.filter(f => typeof m[f] !== 'function');
  console.log('导出函数检查: ' + (missing.length ? ('缺失 ' + missing.join(',')) : '全部存在 (' + funcs.length + ')'));
  if (missing.length) { console.log('WASM 导出未生效，需修正 EXPORTED_FUNCTIONS'); return; }

  const L1 = makeLine(m, 0, 0, 4, 4);
  const L2 = makeLine(m, 0, 4, 4, 0);
  console.log('创建曲线 rc: L1=' + L1.rc + ' L2=' + L2.rc + ' (0=OK)');

  // evaluate L1 @ 0.5 -> 期望 (2,2)
  const ev = evalPoint(m, L1.handle, 0.5);
  console.log('evaluate L1@0.5 = (' + ev.x.toFixed(6) + ', ' + ev.y.toFixed(6) + ')  期望 (2,2)  rc=' + ev.rc);
  const evOk = Math.abs(ev.x - 2) < 1e-9 && Math.abs(ev.y - 2) < 1e-9;

  // unified_solve L1 x L2
  const resPtr = m._malloc(8), nPtr = m._malloc(4);
  const rc = m._nurbs_unified_solve(L1.handle, L2.handle, resPtr, nPtr);
  const n = m.HEAPU32[nPtr >> 2];
  const res = m.HEAPU32[resPtr >> 2]; // NurbsSolveResult* 数组
  let state = -99, u = 0, v = 0, dist = -1;
  if (n >= 1) {
    // NurbsSolveResult: int state(off0), pad, double u(8),v(16),p1x(24),p1y(32),p2x(40),p2y(48),distance(56),int iterations(64),converged(68)
    state = m.HEAPU32[res >> 2];
    u = m.HEAPF64[(res + 8) >> 3];
    v = m.HEAPF64[(res + 16) >> 3];
    dist = m.HEAPF64[(res + 56) >> 3];
  }
  console.log('unified_solve rc=' + rc + ' n=' + n + ' state=' + state + ' (0=INTERSECTION)');
  console.log('  u=' + u.toFixed(6) + ' v=' + v.toFixed(6) + ' distance=' + dist.toExponential(3));
  const solveOk = n === 1 && state === 0 && Math.abs(u - 0.5) < 1e-4 && Math.abs(v - 0.5) < 1e-4 && dist < 1e-6;
  if (n >= 1) m._nurbs_solve_results_free(m.HEAPU32[resPtr >> 2]);
  m._free(resPtr); m._free(nPtr);

  // BVH 加速求解对照
  const iaPtr = m._malloc(8), ibPtr = m._malloc(8);
  m._nurbs_bvh_index_create(L1.handle, iaPtr);
  m._nurbs_bvh_index_create(L2.handle, ibPtr);
  const rPtr = m._malloc(8), nPtr2 = m._malloc(4);
  const brc = m._nurbs_bvh_solve(m.HEAPU32[iaPtr >> 2], L1.handle, m.HEAPU32[ibPtr >> 2], L2.handle, rPtr, nPtr2);
  const bn = m.HEAPU32[nPtr2 >> 2];
  const br = m.HEAPU32[rPtr >> 2];
  let bstate = -99, bdist = -1;
  if (bn >= 1) { bstate = m.HEAPU32[br >> 2]; bdist = m.HEAPF64[(br + 56) >> 3]; }
  console.log('BVH solve rc=' + brc + ' n=' + bn + ' state=' + bstate + ' distance=' + bdist.toExponential(3));
  const bvhOk = bn === 1 && bstate === 0 && bdist < 1e-6;
  if (bn >= 1) m._nurbs_solve_results_free(br);
  m._free(rPtr); m._free(nPtr2); m._free(iaPtr); m._free(ibPtr);

  // enumerate_intersections：全局保证性枚举（区间细分 + Krawczyk 证明）汇总
  const ePtr = m._malloc(16);
  const erc = m._nurbs_enumerate_intersections(L1.handle, L2.handle, 1e-4, ePtr);
  let enOk = false;
  if (erc === 0) {
    const pc = m.HEAPU32[ePtr >> 2], tc = m.HEAPU32[(ePtr + 4) >> 2];
    const eb = m.HEAPU32[(ePtr + 8) >> 2], cpl = m.HEAPU32[(ePtr + 12) >> 2];
    enOk = (pc + tc) >= 1 && cpl === 1;  // 检测到解（证明或疑似相切）+ 全域覆盖
    console.log('enumerate_intersections rc=' + erc + ' proven=' + pc + ' tangent=' + tc + ' boxes=' + eb + ' complete=' + cpl);
  } else {
    console.log('enumerate_intersections rc=' + erc + ' (期望 0)');
  }
  m._free(ePtr);

  m._nurbs_curve_destroy(L1.handle);
  m._nurbs_curve_destroy(L2.handle);

  // ---- 随机批量回归：50 对随机曲线，自洽判定 ----
  function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;var t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
  function makeCurve(m, rnd) {
    const d = 1 + Math.floor(rnd() * 3);
    const n = d + 1 + Math.floor(rnd() * 8);
    const ctrlBytes = n * 3 * 8, knotBytes = (n + d + 1) * 8;
    const cp = m._malloc(ctrlBytes), kn = m._malloc(knotBytes);
    const c = cp / 8, k = kn / 8;
    for (let i = 0; i < n; i++) { m.HEAPF64[c + i*3] = rnd() * 10; m.HEAPF64[c + i*3 + 1] = rnd() * 10; m.HEAPF64[c + i*3 + 2] = 1; }
    for (let i = 0; i < n + d + 1; i++) m.HEAPF64[k + i] = i <= d ? 0 : (i > n ? 1 : (i - d) / (n - d));
    const out = m._malloc(8);
    const rc = m._nurbs_curve_create(cp, n, kn, n + d + 1, d, out);
    const handle = m.HEAPU32[out >> 2];
    m._free(cp); m._free(kn); m._free(out);
    return { handle, rc };
  }
  const rnd = mulberry32(20261012);
  let pairs = 0, bad = 0;
  for (let t = 0; t < 50; t++) {
    const A = makeCurve(m, rnd), B = makeCurve(m, rnd);
    if (A.rc !== 0 || B.rc !== 0) { bad++; continue; }
    pairs++;
    const rPtr = m._malloc(8), nPtr = m._malloc(4);
    const rc = m._nurbs_unified_solve(A.handle, B.handle, rPtr, nPtr);
    const n = m.HEAPU32[nPtr >> 2], res = m.HEAPU32[rPtr >> 2];
    let badPair = false;
    for (let i = 0; i < n; i++) {
      const off = res + i * 88;
      const st = m.HEAPU32[off >> 2], dist = m.HEAPF64[(off + 56) >> 3];
      if (st === 0 && !(dist < 1e-6)) badPair = true; // INTERSECTION 必须收敛
    }
    if (badPair) bad++;
    if (n >= 1) m._nurbs_solve_results_free(res);
    m._free(rPtr); m._free(nPtr);
    m._nurbs_curve_destroy(A.handle);
    m._nurbs_curve_destroy(B.handle);
  }
  console.log('随机批量回归: ' + pairs + ' 对, 自洽失败 ' + bad + ' (0=全通过)');
  const randOk = bad === 0 && pairs === 50;

  console.log('---');
  console.log('evaluate: ' + (evOk ? 'PASS' : 'FAIL'));
  console.log('unified_solve: ' + (solveOk ? 'PASS' : 'FAIL'));
  console.log('BVH solve: ' + (bvhOk ? 'PASS' : 'FAIL'));
  console.log('enumerate: ' + (enOk ? 'PASS' : 'FAIL'));
  console.log('random batch: ' + (randOk ? 'PASS' : 'FAIL'));
  const allOk = evOk && solveOk && bvhOk && enOk && randOk;
  console.log('WASM 数值验证: ' + (allOk ? '全部通过' : '存在失败'));
  process.exit(allOk ? 0 : 1);
}

main().catch(e => { console.error('ERR', e); process.exit(1); });
