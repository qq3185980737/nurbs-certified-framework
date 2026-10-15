// ============================================================================
// cad_integration.c — 最小 CAD 集成示例
// ----------------------------------------------------------------------------
// 演示外部 CAD 内核如何仅用纯 C 接口接入 nurbs.lib：
//   1) 构造常见几何曲线（线段 / 二次 / 三次样条）
//   2) 批量求交 + 最近点判定（unifiedSolve）
//   3) 全交点求解（findAllIntersections）
//   4) 错误处理（非法曲线、异常转错误码）
// 编译：cl /TC cad_integration.c nurbs.lib  （或经 CMake example 目标）
// ============================================================================
#include "nurbs/nurbs_c_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const char* state_name(int s) {
    switch (s) {
        case NURBS_STATE_INTERSECTION: return "INTERSECTION 相交";
        case NURBS_STATE_TANGENT:      return "TANGENT 相切";
        case NURBS_STATE_SEPARATED:    return "SEPARATED 分离(最近点)";
        case NURBS_STATE_COINCIDENT:   return "COINCIDENT 重合";
        default:                       return "UNDETERMINED 无法判定";
    }
}

/* 构造 clamped 均匀 knot：长度 = n + p + 1 */
static void clamped_knots(double* k, int n, int p) {
    int i;
    for (i = 0; i <= p; ++i) k[i] = 0.0;
    int mid = n - p - 1;
    for (i = 1; i <= mid; ++i) k[p + i] = i / (double)(mid + 1);
    for (i = 0; i <= p; ++i) k[n + i] = 1.0;
}
static NurbsCurve* make_curve(const NurbsControlPoint* pts, int n, int degree) {
    double* k = (double*)malloc((size_t)(n + degree + 1) * sizeof(double));
    clamped_knots(k, n, degree);
    NurbsCurve* c = NULL;
    int rc = nurbs_curve_create(pts, n, k, n + degree + 1, degree, &c);
    free(k);
    if (rc != NURBS_OK) { printf("[demo] 曲线创建失败，错误码=%d\n", rc); return NULL; }
    return c;
}

int main(void) {
    /* ---- 1) 构造几何 ----
       L1: y=x        (一次)
       L2: y=4-x      (一次)   → 与 L1 交于 (2,2)
       Q : 二次样条   (0,0)(1,2)(2,0)
       C : 三次样条   (0,0)(1,3)(3,-1)(4,1)   */
    const NurbsControlPoint l1[2] = { {0,0,1}, {4,4,1} };
    const NurbsControlPoint l2[2] = { {0,4,1}, {4,0,1} };
    const NurbsControlPoint q3[3] = { {0,0,1}, {1,2,1}, {2,0,1} };
    const NurbsControlPoint c4[4] = { {0,0,1}, {1,3,1}, {3,-1,1}, {4,1,1} };
    NurbsCurve *L1 = make_curve(l1, 2, 1), *L2 = make_curve(l2, 2, 1);
    NurbsCurve *Q = make_curve(q3, 3, 2), *C = make_curve(c4, 4, 3);
    if (!L1 || !L2 || !Q || !C) return 1;

    /* ---- 2) 求交 / 最近点：unifiedSolve 返回所有解，取距离最小 ---- */
    printf("=== unifiedSolve（批量求交 / 最近点）===\n");
    const NurbsCurve* pairs[][2] = { {L1,L2}, {Q,L1}, {C,L2}, {Q,C} };
    const char* name[4] = { "L1×L2", "Q×L1", "C×L2", "Q×C" };
    for (int i = 0; i < 4; ++i) {
        NurbsSolveResult* r = NULL; int n = 0;
        int rc = nurbs_unified_solve(pairs[i][0], pairs[i][1], &r, &n);
        if (rc != NURBS_OK) { printf("  %-6s 求解失败 rc=%d\n", name[i], rc); continue; }
        int b = 0; for (int j = 1; j < n; ++j) if (r[j].distance < r[b].distance) b = j;
        printf("  %-6s 解数=%d | 最佳: %s u=%.6f v=%.6f P1=(%.4f,%.4f) P2=(%.4f,%.4f) d=%.3e iter=%d\n",
               name[i], n, state_name(r[b].state), r[b].u, r[b].v,
               r[b].point1_x, r[b].point1_y, r[b].point2_x, r[b].point2_y,
               r[b].distance, r[b].iterations);
        nurbs_solve_results_free(r);
    }

    /* ---- 3) 全交点：L1×L2 恰一个交点 (2,2) ---- */
    printf("\n=== findAllIntersections（全交点）===\n");
    NurbsRoot* roots = NULL; int nroot = 0;
    int rc = nurbs_find_all_intersections(L1, L2, &roots, &nroot);
    if (rc == NURBS_OK) {
        printf("  L1×L2 根数=%d\n", nroot);
        for (int i = 0; i < nroot; ++i)
            printf("    root#%d u=%.9f v=%.9f d=%.3e P=(%.6f,%.6f)\n",
                   i, roots[i].u, roots[i].v, roots[i].distance, roots[i].point1_x, roots[i].point1_y);
        nurbs_roots_free(roots);
    } else printf("  全交点失败 rc=%d\n", rc);

    /* ---- 4) 求值 + 错误处理 ---- */
    printf("\n=== 求值 & 错误处理 ===\n");
    double x, y;
    if (nurbs_evaluate_point(Q, 0.5, &x, &y) == NURBS_OK) printf("  Q(0.5)=(%.4f,%.4f)\n", x, y);
    const NurbsControlPoint bad[2] = { {0,0,1}, {1,1,-1} };
    const double badk[4] = { 0,0,1,1 };
    NurbsCurve* badC = NULL;
    int badRc = nurbs_curve_create(bad, 2, badk, 4, 1, &badC);
    printf("  非法曲线(负权重)创建 → 错误码=%d %s\n", badRc,
           badRc == NURBS_ERR_INVALID ? "(NURBS_ERR_INVALID，已在入口拦截)" : "(???)");

    nurbs_curve_destroy(L1); nurbs_curve_destroy(L2);
    nurbs_curve_destroy(Q); nurbs_curve_destroy(C);
    printf("\n[demo] 集成示例完成\n");
    return 0;
}
