// ============================================================================
// test_c_api.c — 纯 C 调用测试（C 编译器编译 .c，验证 extern "C" 接口可被
// 外部 C 程序链接调用，且数值与解析解一致）。
// 覆盖：曲线创建/校验、求值/一阶/二阶导、unifiedSolve、findAllIntersections、
//       非法输入错误码、内存释放。
// ============================================================================
#include "nurbs/nurbs_c_api.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static int g_pass = 0, g_fail = 0;
#define CHECK(name, cond) do { if (cond) { g_pass++; } else { g_fail++; printf("[capi] FAIL: %s\n", name); } } while (0)

static int near(double a, double b, double tol) { return fabs(a - b) < tol; }

int main(void) {
    /* ---- 1) 二次 Bezier：P0=(0,0) P1=(1,2) P2=(2,0)，权重全 1，knots clamped [0,0,0,1,1,1] ---- */
    const NurbsControlPoint qp[3] = { {0,0,1}, {1,2,1}, {2,0,1} };
    const double qk[6] = { 0,0,0,1,1,1 };
    NurbsCurve* q = NULL;
    CHECK("create 二次曲线", nurbs_curve_create(qp, 3, qk, 6, 2, &q) == NURBS_OK && q != NULL);

    double x, y;
    CHECK("P(0)", nurbs_evaluate_point(q, 0.0, &x, &y) == NURBS_OK && near(x,0,1e-12) && near(y,0,1e-12));
    CHECK("P(1)", nurbs_evaluate_point(q, 1.0, &x, &y) == NURBS_OK && near(x,2,1e-12) && near(y,0,1e-12));
    CHECK("P(0.5)=(1,1)", nurbs_evaluate_point(q, 0.5, &x, &y) == NURBS_OK && near(x,1,1e-12) && near(y,1,1e-12));
    /* 一阶导 C'(0.5)=(2,0)；二阶导 C''=(0,-8) */
    CHECK("C'(0.5)=(2,0)", nurbs_evaluate_first_derivative(q, 0.5, &x, &y) == NURBS_OK && near(x,2,1e-9) && near(y,0,1e-9));
    CHECK("C''=(0,-8)", nurbs_evaluate_second_derivative(q, 0.5, &x, &y) == NURBS_OK && near(x,0,1e-9) && near(y,-8,1e-9));
    nurbs_curve_destroy(q);

    /* ---- 2) 两条一次直线：y=x 与 y=4-x，交点 (2,2)，u=v=0.5 ---- */
    const NurbsControlPoint la[2] = { {0,0,1}, {4,4,1} };
    const double lak[4] = { 0,0,1,1 };
    const NurbsControlPoint lb[2] = { {0,4,1}, {4,0,1} };
    const double lbk[4] = { 0,0,1,1 };
    NurbsCurve *A = NULL, *B = NULL;
    CHECK("create lineA", nurbs_curve_create(la, 2, lak, 4, 1, &A) == NURBS_OK);
    CHECK("create lineB", nurbs_curve_create(lb, 2, lbk, 4, 1, &B) == NURBS_OK);

    NurbsSolveResult* sol = NULL; int nsol = 0;
    CHECK("unifiedSolve 成功", nurbs_unified_solve(A, B, &sol, &nsol) == NURBS_OK && nsol >= 1);
    if (nsol >= 1) {
        /* 取距离最小的解 */
        int bi = 0; for (int i = 1; i < nsol; ++i) if (sol[i].distance < sol[bi].distance) bi = i;
        CHECK("交点状态=INTERSECTION", sol[bi].state == NURBS_STATE_INTERSECTION);
        CHECK("u≈0.5", near(sol[bi].u, 0.5, 1e-6));
        CHECK("v≈0.5", near(sol[bi].v, 0.5, 1e-6));
        CHECK("distance≈0", sol[bi].distance < 1e-6);
    }
    nurbs_solve_results_free(sol);

    NurbsRoot* roots = NULL; int nroot = 0;
    CHECK("findAllIntersections 成功", nurbs_find_all_intersections(A, B, &roots, &nroot) == NURBS_OK && nroot >= 1);
    if (nroot >= 1) {
        CHECK("根 u≈0.5", near(roots[0].u, 0.5, 1e-6));
        CHECK("根 v≈0.5", near(roots[0].v, 0.5, 1e-6));
        CHECK("根 distance<1e-6", roots[0].distance < 1e-6);
    }
    nurbs_roots_free(roots);
    nurbs_curve_destroy(A); nurbs_curve_destroy(B);

    /* ---- 3) 非法输入：负权重 -> NURBS_ERR_INVALID ---- */
    const NurbsControlPoint bad[2] = { {0,0,1}, {1,1,-1} };
    const double badk[4] = { 0,0,1,1 };
    NurbsCurve* badC = NULL;
    CHECK("负权重报 INVALID", nurbs_curve_create(bad, 2, badk, 4, 1, &badC) == NURBS_ERR_INVALID && badC == NULL);

    /* ---- 4) 阈值配置可调用（不改默认即可返回） ---- */
    nurbs_set_geometry_tolerance(1e-8);
    nurbs_set_angular_tolerance_deg(12.0);
    nurbs_set_convergence_threshold(1e-12);

    /* ---- 5) BVH 加速求解（A1/B1）：索引复用，结果应与 unified_solve 一致 ---- */
    /* 重新建 A、B（上面已 destroy） */
    CHECK("BVH: 重建 lineA", nurbs_curve_create(la, 2, lak, 4, 1, &A) == NURBS_OK);
    CHECK("BVH: 重建 lineB", nurbs_curve_create(lb, 2, lbk, 4, 1, &B) == NURBS_OK);
    NurbsBvhIndex *ia = NULL, *ib = NULL;
    CHECK("BVH: indexA 创建", nurbs_bvh_index_create(A, &ia) == NURBS_OK && ia != NULL);
    CHECK("BVH: indexB 创建", nurbs_bvh_index_create(B, &ib) == NURBS_OK && ib != NULL);
    sol = NULL; nsol = 0;
    CHECK("BVH: solve 成功", nurbs_bvh_solve(ia, A, ib, B, &sol, &nsol) == NURBS_OK && nsol >= 1);
    if (nsol >= 1) {
        int bi = 0; for (int i = 1; i < nsol; ++i) if (sol[i].distance < sol[bi].distance) bi = i;
        CHECK("BVH: 状态=INTERSECTION", sol[bi].state == NURBS_STATE_INTERSECTION);
        CHECK("BVH: u≈0.5", near(sol[bi].u, 0.5, 1e-6));
        CHECK("BVH: v≈0.5", near(sol[bi].v, 0.5, 1e-6));
        CHECK("BVH: distance<1e-6", sol[bi].distance < 1e-6);
    }
    nurbs_solve_results_free(sol);
    nurbs_bvh_index_destroy(ia);
    nurbs_bvh_index_destroy(ib);
    nurbs_curve_destroy(A); nurbs_curve_destroy(B);

    printf("[capi] ==== C API 测试汇总 ==== pass=%d fail=%d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
