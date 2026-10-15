cmake_minimum_required(VERSION 3.16)
project(nurbs_cpp LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- 独立 NURBS 算法库（静态库，含 C++ 核心 + 纯 C API 包装） ----
add_library(nurbs STATIC
  src/nurbs_core.cpp
  src/unified_solver.cpp
  src/nurbs_c_api.cpp
  src/accelerated_solver.cpp
)
target_include_directories(nurbs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# 数值保真：关闭 fast-math / 允许浮点重组，保证与 JS 逐位一致
# 编码：源文件为 UTF-8，需显式告知 MSVC，避免按 GBK 代码页误读中文注释
if(MSVC)
  target_compile_options(nurbs PRIVATE /fp:precise /utf-8)
else()
  target_compile_options(nurbs PRIVATE -fno-fast-math)
endif()

# ---- 等价性验证程序 ----
enable_testing()
add_executable(test_equivalence tests/test_equivalence.cpp)
target_link_libraries(test_equivalence PRIVATE nurbs)
if(MSVC)
  target_compile_options(test_equivalence PRIVATE /utf-8)
endif()
# 夹具路径作为运行参数传入；默认 tests/fixture.txt
add_test(NAME equivalence COMMAND test_equivalence ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixture.txt)

# ---- 鲁棒性 / 防御性验证程序（不依赖夹具）----
add_executable(test_robust tests/test_robust.cpp)
target_link_libraries(test_robust PRIVATE nurbs)
if(MSVC)
  target_compile_options(test_robust PRIVATE /utf-8)
endif()
add_test(NAME robust COMMAND test_robust)

# ---- 纯 C 接口调用测试（C 编译器编译，链接 C++ 库） ----
add_executable(test_c_api tests/test_c_api.c)
target_link_libraries(test_c_api PRIVATE nurbs)
if(MSVC)
  target_compile_options(test_c_api PRIVATE /utf-8)
endif()
add_test(NAME c_api COMMAND test_c_api)

# ---- 最小 CAD 集成示例（纯 C，经 C API 接入）----
add_executable(cad_integration examples/cad_integration.c)
target_link_libraries(cad_integration PRIVATE nurbs)
if(MSVC)
  target_compile_options(cad_integration PRIVATE /utf-8)
endif()

# ---- 共享库（DLL）：供 Python ctypes / 其它语言 FFI 加载 ----
add_executable(bench_accel tests/bench_accel.cpp)
target_link_libraries(bench_accel PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_accel PRIVATE /fp:precise /utf-8)
endif()
add_executable(bench_batch tests/bench_batch.cpp)
target_link_libraries(bench_batch PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_batch PRIVATE /fp:precise /utf-8)
endif()
# 任务3：与开源库 SISL 对照基准（SISL 库由本机独立 CMake 构建）
add_executable(bench_sisl tests/bench_sisl.cpp)
target_link_libraries(bench_sisl PRIVATE nurbs)
target_include_directories(bench_sisl PRIVATE "D:/sisl_dl/src/SISL-master/include")
target_link_libraries(bench_sisl PRIVATE "D:/sisl_dl/build/sisl.lib")
if(MSVC)
  target_compile_options(bench_sisl PRIVATE /fp:precise /utf-8)
endif()
# A: 保证性求交（区间/Krawczyk 解存在性证明）覆盖率基准
add_executable(bench_prove tests/bench_prove.cpp)
target_link_libraries(bench_prove PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_prove PRIVATE /fp:precise /utf-8)
endif()
# E: 实验矩阵扩展（相切/退化/高次/多重交点/尺度不变性）
add_executable(bench_matrix tests/bench_matrix.cpp)
target_link_libraries(bench_matrix PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_matrix PRIVATE /fp:precise /utf-8)
endif()
# B/C: 解分类 + 条件数/误差界 + 重叠检测
add_executable(bench_classify tests/bench_classify.cpp)
target_link_libraries(bench_classify PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_classify PRIVATE /fp:precise /utf-8)
endif()
# ④: 全局保证性求交（区间细分全局枚举 + Krawczyk 完备证明）
add_executable(bench_enumerate tests/bench_enumerate.cpp)
target_link_libraries(bench_enumerate PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_enumerate PRIVATE /fp:precise /utf-8)
endif()
# ①: GS-Newton vs Krawczyk 区间证明对比
add_executable(bench_prove_gs tests/bench_prove_gs.cpp)
target_link_libraries(bench_prove_gs PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_prove_gs PRIVATE /fp:precise /utf-8)
endif()
# ②: 相切交点精化（阻尼 Gauss-Newton，对标 JCAD 2025 有理插值）
add_executable(bench_tangent tests/bench_tangent.cpp)
target_link_libraries(bench_tangent PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_tangent PRIVATE /fp:precise /utf-8)
endif()
# ③: k-DOP BVH 精剪（14-DOP vs AABB）
add_executable(bench_kdop tests/bench_kdop.cpp)
target_link_libraries(bench_kdop PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_kdop PRIVATE /fp:precise /utf-8)
endif()
# ⑤: 非方系统区间 Krawczyk（三曲线共点 / T 型结点证明）
add_executable(bench_tricurve tests/bench_tricurve.cpp)
target_link_libraries(bench_tricurve PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_tricurve PRIVATE /fp:precise /utf-8)
endif()
# ⑥: 同伦连续追踪 + 参数化区间验证（相切交点保证）
add_executable(bench_homotopy tests/bench_homotopy.cpp)
target_link_libraries(bench_homotopy PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_homotopy PRIVATE /fp:precise /utf-8)
endif()
# ⑦: OpenCL GPU 并行 BVH 叶对预筛（14-DOP 穷举）
add_executable(bench_opencl tests/bench_opencl.cpp)
target_link_libraries(bench_opencl PRIVATE nurbs)
target_include_directories(bench_opencl PRIVATE "D:/Download/opencl/OCL_SDK_Light/include")
target_link_libraries(bench_opencl PRIVATE "D:/Download/opencl/OCL_SDK_Light/lib/x86_64/opencl.lib")
if(MSVC)
  target_compile_options(bench_opencl PRIVATE /fp:precise /utf-8 /EHsc)
endif()
# ⑦ 扩展:更大规模 GPU 吞吐基准（多档 batch 扫描，展示加速比随规模增长）
add_executable(bench_opencl_big tests/bench_opencl_big.cpp)
target_link_libraries(bench_opencl_big PRIVATE nurbs)
target_include_directories(bench_opencl_big PRIVATE "D:/Download/opencl/OCL_SDK_Light/include")
target_link_libraries(bench_opencl_big PRIVATE "D:/Download/opencl/OCL_SDK_Light/lib/x86_64/opencl.lib")
if(MSVC)
  target_compile_options(bench_opencl_big PRIVATE /fp:precise /utf-8 /EHsc)
endif()
# 扩展：3D 空间曲线求交（超定 R²→R³ 广义 Krawczyk，自包含）
add_executable(bench_3d tests/bench_3d.cpp)
if(MSVC)
  target_compile_options(bench_3d PRIVATE /fp:precise /utf-8)
endif()
# 扩展：3D 曲线-曲面求交（方系统 R³→R³ Krawczyk，自包含）
add_executable(bench_cs tests/bench_cs.cpp)
if(MSVC)
  target_compile_options(bench_cs PRIVATE /fp:precise /utf-8)
endif()
# 扩展：3D 曲面-曲面求交（欠定 R⁴→R³ 交线追踪，自包含）
add_executable(bench_ss tests/bench_ss.cpp)
if(MSVC)
  target_compile_options(bench_ss PRIVATE /fp:precise /utf-8)
endif()
# 补充：非均匀 + 有理 NURBS 保证性实测（标准圆弧/椭圆精确表示）
add_executable(bench_nurbs_wide tests/bench_nurbs_wide.cpp)
target_link_libraries(bench_nurbs_wide PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_nurbs_wide PRIVATE /fp:precise /utf-8)
endif()
# 补充：C⁰ 折角（内部重复节点）分段证明
add_executable(bench_kink tests/bench_kink.cpp)
target_link_libraries(bench_kink PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_kink PRIVATE /fp:precise /utf-8)
endif()
# 补充：闭合曲线求交 + wrap-around 同点处理
add_executable(bench_closed tests/bench_closed.cpp)
target_link_libraries(bench_closed PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_closed PRIVATE /fp:precise /utf-8)
endif()
# 综合：NURBS 曲线类型全覆盖 × 跨类型相交（每类型 1000 组）
add_executable(bench_massive tests/bench_massive.cpp)
target_link_libraries(bench_massive PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_massive PRIVATE /fp:precise /utf-8)
endif()
# 综合：曲线×曲面 + 曲面×曲面 跨类型相交（批量证明）
add_executable(bench_massive_surf tests/bench_massive_surf.cpp)
if(MSVC)
  target_compile_options(bench_massive_surf PRIVATE /fp:precise /utf-8)
endif()
# 退化/奇异曲线专项（当前框架行为测试）
add_executable(bench_degenerate tests/bench_degenerate.cpp)
target_link_libraries(bench_degenerate PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_degenerate PRIVATE /fp:precise /utf-8)
endif()
# 退化/奇异曲线专门处理机制（恢复率验证）
add_executable(bench_degenerate_solve tests/bench_degenerate_solve.cpp)
target_link_libraries(bench_degenerate_solve PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_degenerate_solve PRIVATE /fp:precise /utf-8)
endif()
# 退化/奇异曲面专项（检测机制）
add_executable(bench_surf_degenerate tests/bench_surf_degenerate.cpp)
if(MSVC)
  target_compile_options(bench_surf_degenerate PRIVATE /fp:precise /utf-8)
endif()
# 退化曲面 求解/证明恢复（降维 + 点包含）
add_executable(bench_surf_degenerate_solve tests/bench_surf_degenerate_solve.cpp)
if(MSVC)
  target_compile_options(bench_surf_degenerate_solve PRIVATE /fp:precise /utf-8)
endif()
add_executable(bench_unified_degenerate tests/bench_unified_degenerate.cpp)
target_link_libraries(bench_unified_degenerate PRIVATE nurbs)
if(MSVC)
  target_compile_options(bench_unified_degenerate PRIVATE /fp:precise /utf-8)
endif()
add_executable(bench_surf_w_degenerate tests/bench_surf_w_degenerate.cpp)
if(MSVC)
  target_compile_options(bench_surf_w_degenerate PRIVATE /fp:precise /utf-8)
endif()
add_executable(bench_surface_selfintersection tests/bench_surface_selfintersection.cpp)
if(MSVC)
  target_compile_options(bench_surface_selfintersection PRIVATE /fp:precise /utf-8)
endif()
add_executable(bench_surface_overlap tests/bench_surface_overlap.cpp)
if(MSVC)
  target_compile_options(bench_surface_overlap PRIVATE /fp:precise /utf-8)
endif()
add_executable(bench_surface_selfintersection_ext tests/bench_surface_selfintersection_ext.cpp)
if(MSVC)
  target_compile_options(bench_surface_selfintersection_ext PRIVATE /fp:precise /utf-8)
endif()
add_library(nurbs_shared SHARED
  src/nurbs_core.cpp
  src/unified_solver.cpp
  src/nurbs_c_api.cpp
  src/accelerated_solver.cpp
)
set_target_properties(nurbs_shared PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_include_directories(nurbs_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(MSVC)
  target_compile_options(nurbs_shared PRIVATE /fp:precise /utf-8)
endif()
