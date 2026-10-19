// ============================================================================
// bench_osculate.cpp — 路线 D：相切阶的理论化
//   用"带符号曲率 + 接触阶"判定替代旧的角度启发式（切向量夹角 <12° 即相切），
//   对标 ACIS 分类：
//     INTERSECTION ：切向量不平行（横截）
//     TANGENT      ：一阶接触（切向量平行，但曲率圆不同 → 带符号曲率不相等）
//     OSCULATING   ：二阶接触（切向量平行 且 带符号曲率相等 → 曲率圆贴合/相切）
//   意义：角度启发式只能把"相切"笼统判为一类，无法区分一阶/二阶切触；曲率判据
//   给出接触阶分类，与 CAD 内核（ACIS）的 tangency vs osculation 对齐。
// 自包含：定义多项式 Bezier 曲线 + 解析导数/带符号曲率，不动主库。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

struct P2 { double x, y; };
struct Curve {
    std::vector<P2> P;   // Bezier 控制点（degree = size-1）
    int deg() const { return (int)P.size() - 1; }
    P2 eval(double t) const {
        std::vector<P2> c = P; const int n = (int)c.size();
        for (int r = 1; r < n; ++r) for (int i = 0; i < n - r; ++i) {
            c[i].x = (1 - t) * c[i].x + t * c[i + 1].x;
            c[i].y = (1 - t) * c[i].y + t * c[i + 1].y;
        }
        return c[0];
    }
    P2 deriv(double t, int order) const {   // 解析高阶导数（Bernstein 系数差分）
        std::vector<P2> Q = P; int p = deg();
        for (int o = 0; o < order; ++o) {
            if (p == 0) return { 0, 0 };
            std::vector<P2> R(p);
            for (int i = 0; i < p; ++i) R[i] = { p * (Q[i + 1].x - Q[i].x), p * (Q[i + 1].y - Q[i].y) };
            Q = R; --p;
        }
        const int nn = (int)Q.size(); std::vector<P2> c = Q;
        for (int r = 1; r < nn; ++r) for (int i = 0; i < nn - r; ++i) {
            c[i].x = (1 - t) * c[i].x + t * c[i + 1].x;
            c[i].y = (1 - t) * c[i].y + t * c[i + 1].y;
        }
        return c[0];
    }
    double signedCurvature(double t) const {   // κ = (x'y'' − y'x'')/|r'|³（带符号：弯向）
        const P2 d1 = deriv(t, 1), d2 = deriv(t, 2);
        const double L = std::hypot(d1.x, d1.y);
        if (L < 1e-12) return 0.0;
        return (d1.x * d2.y - d1.y * d2.x) / std::pow(L, 3);
    }
};

// 接触阶分类：0=INTERSECTION, 1=TANGENT(一阶), 2=OSCULATING(二阶)
static int classifyContact(const Curve& c1, double u, const Curve& c2, double v,
                           double angleTolDeg, double curvRelTol) {
    const P2 t1 = c1.deriv(u, 1), t2 = c2.deriv(v, 1);
    const double n1 = std::hypot(t1.x, t1.y), n2 = std::hypot(t2.x, t2.y);
    if (n1 < 1e-12 || n2 < 1e-12) return 1;   // 奇异（顶点）→ 保守判相切
    double ca = (t1.x * t2.x + t1.y * t2.y) / (n1 * n2);
    ca = std::max(-1.0, std::min(1.0, ca));
    const double ang = std::acos(ca) * 180.0 / 3.141592653589793;
    if (ang > angleTolDeg && std::abs(180.0 - ang) > angleTolDeg) return 0;  // 横截
    // 切向量平行 → 相切：用带符号曲率区分接触阶
    const double k1 = c1.signedCurvature(u), k2 = c2.signedCurvature(v);
    const double km = std::max(std::abs(k1), std::abs(k2));
    if (km > 1e-12 && std::abs(k1 - k2) <= curvRelTol * km) return 2;  // 曲率圆相切 → 二阶
    return 1;                                                           // 一阶接触
}

static const char* nm(int r) { return r == 0 ? "INTERSECTION" : r == 1 ? "TANGENT(1st)" : "OSCULATING(2nd)"; }

int main() {
    printf("==== 路线 D：相切阶理论化（带符号曲率 + 接触阶，对标 ACIS）====\n\n");

    // 构造曲线：c(t)=(t, y(t)) 的 Bezier 控制点
    const Curve c_x   = { { {0,0}, {1,0} } };                 // 直线 y=0（x 轴）
    const Curve c_y   = { { {0,0}, {1,1} } };                 // 直线 y=x
    const Curve c_h   = { { {0,0.5}, {1,0.5} } };             // 水平线 y=0.5
    const Curve c_par = { { {0,0}, {0.5,0}, {1,1} } };        // y=t²
    const Curve c_npar= { { {0,0}, {0.5,0}, {1,-1} } };       // y=−t²
    const Curve c_pxp = { { {0,0}, {0.25,0}, {0.5,1.0/6.0}, {0.75,0.5}, {1,2} } }; // y=t²+t⁴

    int pass = 0, fail = 0;

    printf("[1] 横截（角度 45°）y=x 与 y=0.5 → 期望 INTERSECTION\n");
    {
        int r = classifyContact(c_y, 0.5, c_h, 0.5, 12.0, 0.1);
        printf("    分类=%s %s\n", nm(r), r == 0 ? "✓" : "✗");
        r == 0 ? ++pass : ++fail;
    }

    printf("\n[2] 一阶相切（同切向、曲率异号）y=t² 与 y=−t² 在 t=0 → 期望 TANGENT(1st)\n");
    {
        int r = classifyContact(c_par, 0.0, c_npar, 0.0, 12.0, 0.1);
        printf("    带符号曲率 κ1=%+.2f κ2=%+.2f → 分类=%s %s\n",
               c_par.signedCurvature(0.0), c_npar.signedCurvature(0.0), nm(r), r == 1 ? "✓" : "✗");
        r == 1 ? ++pass : ++fail;
    }

    printf("\n[3] 一阶相切（直线曲率 0 vs 抛物线曲率 2）y=0 与 y=t² 在 t=0 → 期望 TANGENT(1st)\n");
    {
        int r = classifyContact(c_x, 0.0, c_par, 0.0, 12.0, 0.1);
        printf("    带符号曲率 κ1=%+.2f κ2=%+.2f → 分类=%s %s\n",
               c_x.signedCurvature(0.0), c_par.signedCurvature(0.0), nm(r), r == 1 ? "✓" : "✗");
        r == 1 ? ++pass : ++fail;
    }

    printf("\n[4] 二阶相切（曲率圆相同）y=t² 与 y=t²+t⁴ 在 t=0 → 期望 OSCULATING(2nd)\n");
    {
        int r = classifyContact(c_par, 0.0, c_pxp, 0.0, 12.0, 0.1);
        printf("    带符号曲率 κ1=%+.2f κ2=%+.2f → 分类=%s %s\n",
               c_par.signedCurvature(0.0), c_pxp.signedCurvature(0.0), nm(r), r == 2 ? "✓" : "✗");
        r == 2 ? ++pass : ++fail;
    }

    printf("\n[5] 对比：旧角度启发式（夹角<12° 即 TANGENT）对 [2][3][4] 全部判 TANGENT，无法区分一阶/二阶\n");
    printf("    曲率判据正确区分：一阶（κ1≠κ2）vs 二阶（κ1=κ2）→ 与 ACIS tangency/osculation 对齐\n");

    printf("\n通过 %d，失败 %d。\n", pass, fail);
    return fail ? 1 : 0;
}
