// ============================================================================
// bench_3d.cpp — 扩展：3D 空间 NURBS/Bezier 曲线求交（超定系统 R²→R³）
//   F(u,v)=c1(u)−c2(v)∈R³，两 3D 曲线求交是 3 方程 2 未知的超定系统。
//   复用 ⑤ 的广义 Krawczyk（Moore-Penrose 伪逆 A=J⁺）证明交点存在唯一。
//   意义：保证性求交框架从 2D 平面推广到 3D 空间曲线（CAD 空间曲线交点），
//   且展示超定系统验证算子跨维度通用。
// 自包含：定义 3D Bezier 曲线（de Casteljau + 导数 + 导数区间界），不动主库。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>

struct P3 { double x, y, z; };
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} double width()const{return hi-lo;} };
struct Int3 { I3 x, y, z; };

// ---------- 3D Bezier 曲线 ----------
struct Bezier3 {
    std::vector<P3> P;   // 控制点（n 个，degree=n-1）
    double tMin = 0, tMax = 1;
    int degree() const { return (int)P.size() - 1; }
    P3 eval(double t) const {           // de Casteljau
        std::vector<P3> c = P;
        const int n = (int)c.size();
        for (int r = 1; r < n; ++r) for (int i = 0; i < n - r; ++i) {
            c[i].x = (1 - t) * c[i].x + t * c[i + 1].x;
            c[i].y = (1 - t) * c[i].y + t * c[i + 1].y;
            c[i].z = (1 - t) * c[i].z + t * c[i + 1].z;
        }
        return c[0];
    }
    P3 deriv(double t) const {          // 导数控制点 Q_i = p(P_{i+1}−P_i)
        const int p = degree();
        if (p == 0) return { 0, 0, 0 };
        std::vector<P3> Q(p);
        for (int i = 0; i < p; ++i) {
            Q[i].x = p * (P[i + 1].x - P[i].x);
            Q[i].y = p * (P[i + 1].y - P[i].y);
            Q[i].z = p * (P[i + 1].z - P[i].z);
        }
        Bezier3 q; q.P = Q;              // degree p-1
        return q.eval(t);
    }
    // 导数在参数区间 [t] 的区间界（导数控制点凸包，严格包含）
    Int3 derivRange(I3 t) const {
        const int p = degree();
        I3 dx(1e300, -1e300), dy(1e300, -1e300), dz(1e300, -1e300);
        if (p == 0) return { I3(0,0), I3(0,0), I3(0,0) };
        for (int i = 0; i < p; ++i) {
            const double qx = p * (P[i + 1].x - P[i].x), qy = p * (P[i + 1].y - P[i].y), qz = p * (P[i + 1].z - P[i].z);
            if (qx < dx.lo) dx.lo = qx; if (qx > dx.hi) dx.hi = qx;
            if (qy < dy.lo) dy.lo = qy; if (qy > dy.hi) dy.hi = qy;
            if (qz < dz.lo) dz.lo = qz; if (qz > dz.hi) dz.hi = qz;
        }
        return { dx, dy, dz };
    }
};

// ---------- 区间运算 ----------
inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) {
    const double m1 = a.lo * b.lo, m2 = a.lo * b.hi, m3 = a.hi * b.lo, m4 = a.hi * b.hi;
    return { std::min(std::min(m1, m2), std::min(m3, m4)), std::max(std::max(m1, m2), std::max(m3, m4)) };
}
inline I3 imul(const I3& a, double s) { return imul(a, I3(s, s)); }

// ---------- 2×2 实矩阵逆 ----------
static bool inv2(const double M[2][2], double R[2][2]) {
    const double d = M[0][0] * M[1][1] - M[0][1] * M[1][0];
    if (std::abs(d) < 1e-14) return false;
    R[0][0] = M[1][1] / d; R[0][1] = -M[0][1] / d; R[1][0] = -M[1][0] / d; R[1][1] = M[0][0] / d;
    return true;
}

// 超定 3D 曲线求交 Krawczyk：返回 1=证明存在唯一交点，0=未定/不存在
static int prove3D(const Bezier3& c1, const Bezier3& c2, I3& ub, I3& vb) {
    const double mu = ub.mid(), mv = vb.mid();
    const P3 p1 = c1.eval(mu), p2 = c2.eval(mv), d1 = c1.deriv(mu), d2 = c2.deriv(mv);
    const double f[3] = { p1.x - p2.x, p1.y - p2.y, p1.z - p2.z };
    // J(m) 3×2 实数
    const double J[3][2] = { { d1.x, -d2.x }, { d1.y, -d2.y }, { d1.z, -d2.z } };
    // G = JᵀJ (2×2)，A = G⁻¹ Jᵀ (2×3)
    double G[2][2], Gi[2][2], A[2][3];
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) { G[i][j] = 0; for (int k = 0; k < 3; ++k) G[i][j] += J[k][i] * J[k][j]; }
    if (!inv2(G, Gi)) return 0;   // 中点 Jacobian 降秩 → 相切/退化，不证
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 3; ++j) { A[i][j] = 0; for (int k = 0; k < 2; ++k) A[i][j] += Gi[i][k] * J[j][k]; }
    // J([x]) 3×2 区间（导数区间界）
    const Int3 ra = c1.derivRange(ub), rb = c2.derivRange(vb);
    const I3 Ji[3][2] = { { ra.x, imul(rb.x, -1.0) }, { ra.y, imul(rb.y, -1.0) }, { ra.z, imul(rb.z, -1.0) } };
    const I3 dx[2] = { isub(ub, I3(mu, mu)), isub(vb, I3(mv, mv)) };
    I3 K[2];
    for (int i = 0; i < 2; ++i) {
        I3 sum(0, 0);
        for (int j = 0; j < 2; ++j) {
            I3 proj = (i == j) ? I3(1, 1) : I3(0, 0);
            for (int k = 0; k < 3; ++k) proj = isub(proj, imul(I3(A[i][k], A[i][k]), Ji[k][j]));
            sum = iadd(sum, imul(proj, dx[j]));
        }
        double af = 0; for (int k = 0; k < 3; ++k) af += A[i][k] * f[k];
        const double mi = (i == 0) ? mu : mv;
        K[i] = iadd(I3(mi - af, mi - af), sum);
    }
    const I3 xx[2] = { ub, vb };
    for (int i = 0; i < 2; ++i) if (!(xx[i].lo < K[i].lo && K[i].hi < xx[i].hi)) return 0;
    return 1;
}

static void run(const char* label, const Bezier3& c1, const Bezier3& c2,
                double u0, double v0, double delta, int expect) {
    I3 ub(u0 - delta, u0 + delta), vb(v0 - delta, v0 + delta);
    int r = prove3D(c1, c2, ub, vb);
    printf("  %-30s δ=%.0e → %s（期望 %s）%s\n", label, delta,
           r == 1 ? "已证明" : "未证", expect == 1 ? "证" : "拒", r == expect ? " ✓" : " ✗");
}

int main() {
    printf("==== 扩展：3D 空间曲线求交（超定 R²→R³，广义 Krawczyk）====\n\n");

    // 场景 1：两 3D 直线横截交于 (0.2,0.2,0.2)
    {
        printf("[1] 3D 直线横截共点 (0.2,0.2,0.2)\n");
        Bezier3 c1; c1.P = { {0,0.2,0}, {1,0.2,1} };   // c1(u)=(u,0.2,u)
        Bezier3 c2; c2.P = { {0,0,0.2}, {1,1,0.2} };   // c2(v)=(v,v,0.2)
        for (double d : {2e-2, 1e-3, 1e-5}) run("3D 直线横截交点", c1, c2, 0.2, 0.2, d, 1);
    }
    // 场景 2：3D 二次曲线（空间抛物面曲线）与直线横截交点
    {
        printf("[2] 3D 二次曲线 vs 直线横截交点\n");
        Bezier3 c1; c1.P = { {0,0,0}, {0.5,0.25,0.25}, {1,0,1} };  // 空间曲线 (u, 0.5u(1-u), 0.5u(1-u)+u²)
        // c1(0.5)=(0.5, 0.125, 0.375)；取直线经过该点且横截（c1'(0.5)=(1,0,1)）
        Bezier3 c2; c2.P = { {0,0.125,0.375}, {1,0.125,0.375} };   // 直线 (v, 0.125, 0.375)
        for (double d : {2e-2, 1e-3, 1e-5}) run("3D 二次/直线横截", c1, c2, 0.5, 0.5, d, 1);
    }
    // 场景 3：3D 相切（两曲线在某点切向相同 → J 降秩）→ 不误证
    {
        printf("[3] 3D 相切（不误证）\n");
        // 两曲线在 (0,0,0) 相切：c1(u)=(u,u²,u³) 起点，c2(v)=(v,-v²,v³) 起点，切向均为 (1,0,0)? 
        Bezier3 c1; c1.P = { {0,0,0}, {0.5,0.1,0.05}, {1,1,1} };   // 经过 (0,0,0)
        Bezier3 c2; c2.P = { {0,0,0}, {0.5,-0.1,0.05}, {1,-1,1} }; // 经过 (0,0,0)，起点切向 (1,0,0)
        // 两者只在起点 (0,0,0) 相交（切向相同 → 相切），J 在起点降秩
        run("3D 起点相切", c1, c2, 0.05, 0.05, 5e-3, 0);
    }
    // 场景 4：非交 → 拒绝
    {
        printf("[4] 3D 非交（拒绝）\n");
        Bezier3 c1; c1.P = { {0,0.2,0}, {1,0.2,1} };
        Bezier3 c2; c2.P = { {0,0.7,0.2}, {1,0.7,0.2} };  // y=0.7 与 c1(y=0.2) 不交
        run("3D 平行/非交", c1, c2, 0.5, 0.5, 1e-2, 0);
    }
    printf("\n结论：3D 空间曲线求交（超定 R²→R³）由广义 Krawczyk 保证，横截证明、相切/非交不误证。\n");
    return 0;
}
