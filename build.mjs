#!/usr/bin/env node
/**
 * build.mjs — 统一代码库构建脚本（零依赖）
 * ------------------------------------------------------------------
 * 目标：消除 index.html 内联数学/求解器副本与 src/ 的双份漂移。
 * 做法：
 *   1. 读取 src/wasm_backend.js + src/render/canvas_engine.js + src/diagnostics.js（规范单一事实来源）；
 *   2. 去 ESM(import/export)，将模块体合入一个 IIFE，并把全部导出挂到 window 全局
 *      （UI 目前直接以裸全局名调用 evaluatePoint/unifiedSolve/createCurve 等）；
 *      数学核不再内联 JS 版，由 wasm_backend.js 转发到 cpp/wasm/nurbs.wasm（单源 = C++）。
 *   3. 备份当前 index.html，把其内联数学块整体替换为该 bundle（保留 UI/渲染段不动）；
 *   4. 同步产出 dist/nurbs-bundle.js 供 node --check 与参考。
 *
 * 运行：node build.mjs
 * 注意：内联数学块一旦被 bundle 替换后，index.html 的数学部分不得手改；
 *       所有数值/求解逻辑只改 src/，再重新运行本脚本即可同步到浏览器。
 */
import { readFileSync, writeFileSync, copyFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = __dirname;
const SRC = join(ROOT, 'src');
const INDEX = join(ROOT, 'index.html');
const MODULES = [
  join(SRC, 'wasm_backend.js'),
  join(SRC, 'render', 'canvas_engine.js'),
  join(SRC, 'diagnostics.js'),
];

const IDENT = /[A-Za-z_$][\w$]*/g;
const stripComments = (s) => s.replace(/\/\/[^\n]*/g, '').replace(/\/\*[\s\S]*?\*\//g, '');

/** 把单个 ESM 模块体转成"共享作用域内的普通代码"，并收集其导出的名字 */
function transformModule(src) {
  let code = src;
  const exportedNames = [];

  // 1) 去掉 'use strict'（bundle 顶层统一加一次）
  code = code.replace(/['"]use strict['"];?/g, '');

  // 2) 去掉具名 import 块（可能跨行）
  code = code.replace(/import\s*\{[\s\S]*?\}\s*from\s*['"][^'"]+['"]\s*;?/g, '');
  // 3) 去掉默认 import
  code = code.replace(/import\s+[A-Za-z_$][\w$]*\s+from\s*['"][^'"]+['"]\s*;?/g, '');

  // 4) 提取并移除 export { ... }; 块
  code = code.replace(/export\s*\{([\s\S]*?)\};/g, (_m, body) => {
    const clean = stripComments(body);
    for (const m of clean.matchAll(/\b([A-Za-z_$][\w$]*)\b/g)) {
      const name = m[1];
      // 过滤关键字/纯注释噪音
      if (!/^(as|from|const|let|var|function|class)$/.test(name)) exportedNames.push(name);
    }
    return '';
  });

  // 5) export function / export const / export let / export class → 去掉 export 前缀
  code = code.replace(/export\s+(async\s+)?function\s+([A-Za-z_$][\w$]*)/g, (_m, a, n) => {
    exportedNames.push(n);
    return (a ? a + ' ' : '') + 'function ' + n;
  });
  code = code.replace(/export\s+(const|let|var)\s+([A-Za-z_$][\w$]*)/g, (_m, kw, n) => {
    exportedNames.push(n);
    return kw + ' ' + n;
  });
  code = code.replace(/export\s+class\s+([A-Za-z_$][\w$]*)/g, (_m, n) => {
    exportedNames.push(n);
    return 'class ' + n;
  });

  return { code, exportedNames };
}

// 生成 bundle 文本
function buildBundle() {
  const parts = [];
  const allExports = new Set();
  for (const file of MODULES) {
    const raw = readFileSync(file, 'utf-8');
    const { code, exportedNames } = transformModule(raw);
    parts.push(`// ---- module: ${file.replace(/\\/g, '/').split('/').slice(-3).join('/')} ----\n${code}`);
    for (const n of exportedNames) allExports.add(n);
  }

  const body = parts.join('\n\n');
  const regs = [...allExports]
    .map((n) => `        if (typeof ${n} !== 'undefined') __Nurbs.${n} = ${n};`)
    .join('\n');

  return `(function () {
    'use strict';
    const __Nurbs = {};
${body}
    // ---- 导出注册（收集自各模块 export 语句）----
${regs}

    // ---- 暴露到全局：UI 以裸全局名调用这些函数 ----
    for (const k in __Nurbs) {
        if (!(k in window)) window[k] = __Nurbs[k];
    }

    // ---- 诊断（错误覆盖层），调用失败不阻断应用 ----
    try { if (typeof __Nurbs.initDiagnostics === 'function') initDiagnostics(); }
    catch (e) { console.warn('[build] initDiagnostics skipped:', e && e.message); }
})();
`;
}

// 主流程
function main() {
  const bundle = buildBundle();

  // 备份当前 index.html（每次构建前留一份）
  const stamp = new Date().toISOString().slice(0, 10).replace(/-/g, '');
  const backup = join(ROOT, `_备份_index_before_build_${stamp}.html`);
  copyFileSync(INDEX, backup);
  console.log(`[build] backup -> ${backup}`);

  // 注入 index.html：定位内联 <script> 与"应用状态"起点（let curve1），把"数学块 + 内联渲染器"
  // 整段替换为 bundle（bundle 由 src/math_core + src/solver + src/render/canvas_engine + src/diagnostics 生成）。
  let html = readFileSync(INDEX, 'utf-8');
  const scriptTag = '<script>';
  const uiStart = 'let curve1 = null;';
  const sIdx = html.indexOf(scriptTag);
  const uIdx = html.indexOf(uiStart);
  if (sIdx < 0 || uIdx < 0 || uIdx < sIdx) {
    throw new Error('未能定位 index.html 内联 script / 应用状态起点');
  }
  const head = html.slice(0, sIdx + scriptTag.length);
  const tail = html.slice(uIdx);
  const banner =
    '\n' +
    '        // ============================================================================\n' +
    '        // WASM 数学后端 + 渲染引擎 + 诊断 —— 由 build.mjs 从 src/ 自动生成，禁止手改。\n' +
    '        // 数学核单一事实来源：cpp/src/*.cpp（经 cpp/wasm/nurbs.wasm 供浏览器调用）；\n' +
    '        // 接入层：src/wasm_backend.js + src/render/canvas_engine.js + src/diagnostics.js\n' +
    '        // 修改数值/求解逻辑请编辑 cpp/src 后重建 WASM（见 cpp/wasm/README.md）；\n' +
    '        // 修改接入/渲染逻辑请编辑 src/ 后运行 `node build.mjs` 重新生成本段。\n' +
    '        // ============================================================================\n';
  const newHtml = head + banner + bundle + '\n' + tail;
  writeFileSync(INDEX, newHtml, 'utf-8');

  // 同步产出 dist/nurbs-bundle.js 供 node --check / 参考
  const distDir = join(ROOT, 'dist');
  mkdirSync(distDir, { recursive: true });
  writeFileSync(join(distDir, 'nurbs-bundle.js'), bundle, 'utf-8');

  console.log(`[build] index.html updated (${newHtml.length} chars)`);
  console.log(`[build] dist/nurbs-bundle.js written (${bundle.length} chars)`);
  console.log(`[build] exposed globals: ${[...allExportsOf(bundle)].join(', ')}`);
}

function allExportsOf(bundle) {
  const set = new Set();
  const re = /__Nurbs\.([A-Za-z_$][\w$]*)\s*=/g;
  let m;
  while ((m = re.exec(bundle))) set.add(m[1]);
  return set;
}

main();
