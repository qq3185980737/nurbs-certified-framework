# WASM 单源化构建（B2）

把 NURBS 算法库（`cpp/src/*.cpp`）用 Emscripten 编译为 WebAssembly，让浏览器/Node 运行**同一套 C++ 算法**，为消除 JS/C++ 双维护铺路。现有 `index.html` 前端暂不动。

## 产物（`cpp\wasm\`）
- `nurbs.js`    — Emscripten glue，MODULARIZE 工厂函数 `createNurbsModule()`（CommonJS / ESM 均可加载）
- `nurbs.wasm`  — 编译后的 WebAssembly 字节码
- `wasm_test.mjs` — Node 验证脚本（精确夹具 + 随机批量回归）

## 工具链前置
Emscripten SDK 安装于 `D:\emsdk`（本机 SDK 3.1.22，含 LLVM/Node）。克隆源为 gitee openEuler 镜像（github/ghproxy/gitclone 在本机网络均不可达）。

## 导出 API（与 C API 同名，前缀下划线）
`nurbs_curve_create/destroy`、`nurbs_evaluate_point/_first_derivative/_second_derivative`、
`nurbs_unified_solve/_ex`、`nurbs_solve_results_free`、`nurbs_find_all_intersections`、
`nurbs_roots_free`、`nurbs_bvh_index_create/_destroy`、`nurbs_bvh_solve`。

## 构建命令（分两阶段）
> 关键坑：`-s EXPORTED_FUNCTIONS` 的 JS 数组语法经 emcc.bat / PowerShell 传参会被拆词破坏，
> 必须用**响应文件（@.rsp）+ 正斜杠路径**传给 `emcc.py` 才生效（见 `D:\wasmbuild\link.rsp`）。

```powershell
# 1) 编译 4 个源文件为 .o（无 -s 选项，PowerShell 直接传参即可，已含中文路径）
$env:EMSDK="D:\emsdk"; $env:EM_CONFIG="D:\emsdk\.emscripten"
$env:PATH="D:\emsdk\node\14.18.2_64bit\bin;D:\emsdk\upstream\emscripten;D:\emsdk\python\3.9.2-nuget_64bit;"+$env:PATH
$root="D:\qqdownload\NURBS 统一数值框架几何计算系统"; $obj="$root\cpp\wasm\obj"; New-Item -ItemType Directory -Force $obj|Out-Null
emcc -O3 -std=c++17 -I"$root\cpp\include" -I"$root\cpp\include\nurbs" -c "$root\cpp\src\nurbs_core.cpp" -o "$obj\nurbs_core.o"
emcc -O3 -std=c++17 -I"$root\cpp\include" -I"$root\cpp\include\nurbs" -c "$root\cpp\src\unified_solver.cpp" -o "$obj\unified_solver.o"
emcc -O3 -std=c++17 -I"$root\cpp\include" -I"$root\cpp\include\nurbs" -c "$root\cpp\src\accelerated_solver.cpp" -o "$obj\accelerated_solver.o"
emcc -O3 -std=c++17 -I"$root\cpp\include" -I"$root\cpp\include\nurbs" -c "$root\cpp\src\nurbs_c_api.cpp" -o "$obj\nurbs_c_api.o"

# 2) 复制 .o 到 ASCII 路径，用响应文件链接（响应文件必须正斜杠）
Copy-Item "$obj\*.o" "D:\wasmbuild\" -Force
"D:\emsdk\python\3.9.2-nuget_64bit\python.exe" "D:\emsdk\upstream\emscripten\emcc.py" @D:/wasmbuild/link.rsp

# 3) 产物复制回项目 + 验证
Copy-Item "D:\wasmbuild\nurbs.js","D:\wasmbuild\nurbs.wasm" "$root\cpp\wasm\" -Force
node "$root\cpp\wasm\wasm_test.mjs"
```

`D:\wasmbuild\link.rsp` 内容（正斜杠路径 + `-sNAME=` 连写）：
```
-O3
D:/wasmbuild/nurbs_core.o
D:/wasmbuild/unified_solver.o
D:/wasmbuild/accelerated_solver.o
D:/wasmbuild/nurbs_c_api.o
-sWASM=1
-sMODULARIZE=1
-sEXPORT_NAME=createNurbsModule
-sEXPORTED_FUNCTIONS=['_nurbs_curve_create',...,'_malloc','_free']
-sEXPORTED_RUNTIME_METHODS=['HEAPF64','HEAP32','HEAPU8','UTF8ToString','cwrap']
-o
D:/wasmbuild/nurbs.js
```

## 单线程约束
WASM 默认无 pthread，`accelerated_solver.cpp` 的并行块已用 `__EMSCRIPTEN__` 宏编译期为**串行路径**（结果逐位一致，仅损失并行吞吐）；MSVC 构建不受影响。

## 验证结果（2026-10-13）
- 导出函数 13 个全部存在
- evaluate L1@0.5=(2,2) 精确；unified_solve(L1×L2) → INTERSECTION, u=v=0.5, distance≈1e-13
- BVH solve 距离≈1.5e-10
- 随机 50 对曲线自洽回归：0 失败 → **WASM 数值验证全部通过**
- 与 MSVC C++ 语义一致（同为 C++ 编译产物）；与 JS 原实现做容差等价验收，建议在接入前端时按 fixture 对照（见 `cpp\tests\gen_fixture.mjs` + `test_equivalence.cpp`）
