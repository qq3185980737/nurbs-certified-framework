// ============================================================================
// nurbs_c_api.h — NURBS 算法库的纯 C 接口（extern "C" 包装）
// ----------------------------------------------------------------------------
// 目的：让外部 C / 其它语言（Python/C#/Rust 等经 FFI）或已有 C 内核无需 C++
//       即可链接本库。所有数值均走 C++ 核心，与 JS 逐位一致（见 test_equivalence）。
// 错误约定：每个函数返回 int 错误码，0=成功；-1=参数无效/结构校验失败；
//           -2=求值/数值运行时错误。成功时经输出参数返回结果。
// 所有权：unified_solve / find_all_intersections 返回的数组由 malloc 分配，
//         必须用对应的 _free 函数释放；曲线句柄用 nurbs_curve_destroy 释放。
// ============================================================================
#ifndef NURBS_C_API_H
#define NURBS_C_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- 基础类型 ---- */
typedef struct { double x, y, w; } NurbsControlPoint;   // w 为权重（>0）

/* 不透明曲线句柄 */
typedef struct NurbsCurve NurbsCurve;

/* 不透明 BVH 空间索引句柄（加速求解用，可对多条曲线反复复用） */
typedef struct NurbsBvhIndex NurbsBvhIndex;

/* 几何状态（与 C++/JS GeometryState 一致） */
enum {
    NURBS_STATE_INTERSECTION = 0,
    NURBS_STATE_TANGENT      = 1,
    NURBS_STATE_SEPARATED    = 2,
    NURBS_STATE_COINCIDENT   = 3,
    NURBS_STATE_UNDETERMINED = 4
};

/* unifiedSolve 单解 */
typedef struct {
    int    state;          /* NURBS_STATE_* */
    double u, v;
    double point1_x, point1_y, point2_x, point2_y;
    double distance;
    int    iterations;
    int    converged;      /* 0/1 */
} NurbsSolveResult;

/* findAllIntersections 单根 */
typedef struct {
    double u, v, distance;
    double point1_x, point1_y, point2_x, point2_y;
    int    iterations;
} NurbsRoot;

/* ---- 错误码 ---- */
enum {
    NURBS_OK            = 0,
    NURBS_ERR_INVALID   = -1,   /* 参数无效 / 结构校验失败 */
    NURBS_ERR_RUNTIME   = -2    /* 求值/数值运行时错误 */
};

/* ---- 曲线生命周期 ---- */
/* 控制点数组含权重，knots 长度 = nctrl + degree + 1（clamped 或合法非递减）。
   创建即校验，非法输入返回 NURBS_ERR_INVALID。 */
int nurbs_curve_create(const NurbsControlPoint* ctrl, int nctrl,
                       const double* knots, int nknots, int degree,
                       NurbsCurve** out);
void nurbs_curve_destroy(NurbsCurve* c);

/* ---- 求值 ---- */
int nurbs_evaluate_point(const NurbsCurve* c, double u, double* ox, double* oy);
int nurbs_evaluate_first_derivative(const NurbsCurve* c, double u, double* ox, double* oy);
int nurbs_evaluate_second_derivative(const NurbsCurve* c, double u, double* ox, double* oy);

/* ---- 统一求解 ---- */
/* 默认参数：usePreFiltering=true, gridSize=6, maxIterations=80, tolerance=0(用全局 GEOMETRY_TOLERANCE) */
int nurbs_unified_solve(const NurbsCurve* a, const NurbsCurve* b,
                        NurbsSolveResult** out_results, int* out_n);
/* 显式参数；tolerance<=0 表示用全局默认 */
int nurbs_unified_solve_ex(const NurbsCurve* a, const NurbsCurve* b,
                           int use_prefiltering, int grid_size, int max_iterations, double tolerance,
                           NurbsSolveResult** out_results, int* out_n);
void nurbs_solve_results_free(NurbsSolveResult* r);

/* ---- 全交点 ---- */
/* 默认：tolerance=1e-8, maxDepth=22, paramTol 自动（域宽*1e-3） */
int nurbs_find_all_intersections(const NurbsCurve* a, const NurbsCurve* b,
                                 NurbsRoot** out_roots, int* out_n);
void nurbs_roots_free(NurbsRoot* r);

/* ---- BVH 加速求解（B1：空间索引复用，批量求交吞吐提升） ---- */
/* 曲线索引：建一次、可对任意多次查询复用。索引不持有曲线数据（需与曲线一同存活）。 */
int nurbs_bvh_index_create(const NurbsCurve* curve, NurbsBvhIndex** out);
void nurbs_bvh_index_destroy(NurbsBvhIndex* idx);
/* 用预建索引求解两曲线全部解（结果结构同 unified_solve，取距离最小者为最佳解） */
int nurbs_bvh_solve(const NurbsBvhIndex* ia, const NurbsCurve* ca,
                    const NurbsBvhIndex* ib, const NurbsCurve* cb,
                    NurbsSolveResult** out_results, int* out_n);

/* ---- 阈值配置（全局，影响后续求解；传任何值即覆盖） ---- */
void nurbs_set_geometry_tolerance(double v);        /* 默认 1e-8 */
void nurbs_set_angular_tolerance_deg(double deg);   /* 默认 12° */
void nurbs_set_boundary_tolerance(double v);        /* 默认 1e-4 */
void nurbs_set_param_boundary_tol(double v);        /* 默认 1e-3 */
void nurbs_set_convergence_threshold(double v);     /* 默认 1e-12 */

#ifdef __cplusplus
}
#endif

#endif /* NURBS_C_API_H */
