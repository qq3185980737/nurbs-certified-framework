// ============================================================================
// bench_cs.cpp — 扩展：3D 曲线-曲面求交（方系统 R³→R³，Krawczyk 证明，迁移到统一 NurbsSurface）
//   曲线 C(u)⊂R³ 与双线性 Bezier 曲面 S(s,t)⊂R³ 求交：F(u,s,t)=C(u)−S(s,t)=0，
//   3 方程 3 未知 → 方系统，用 3×3 区间 Krawczyk 证明交点存在唯一。
//   这是 CAD 最经典原语（曲线与曲面交点），把保证性框架从曲线-曲线推到曲线-曲面。
//   曲面统一存储结构：nurbs/NurbsSurface（双线性 makeBilinear）。
// 自包含（不动主库）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <chrono>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
using SurfacePoint = nurbs::SurfacePoint;
static int g_pass = 0, g_fail = 0;
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} double width()const{return hi-lo;} };
struct Int3 { I3 x, y, z; };

inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) {
    const double m1 = a.lo * b.lo, m2 = a.lo * b.hi, m3 = a.hi * b.lo, m4 = a.hi * b.hi;
    return { std::min(std::min(m1, m2), std::min(m3, m4)), std::max(std::max(m1, m2), std::max(m3, m4)) };
}
inline I3 idiv(const I3& a, const I3& b) {
    if (b.lo <= 0.0 && b.hi >= 0.0) return { 1e300, -1e300 };   // 含 0 → 空区间
    double mn = 1e300, mx = -1e300;
    const double c[4] = { a.lo / b.lo, a.lo / b.hi, a.hi / b.lo, a.hi / b.hi };
    for (int i = 0; i < 4; ++i) { if (c[i] < mn) mn = c[i]; if (c[i] > mx) mx = c[i]; }
    return { mn, mx };
}

// ---------- 3D Bezier 曲线 ----------
struct Bezier3 {
    std::vector<P3> P;
    P3 eval(double t) const {
        std::vector<P3> c = P; const int n = (int)c.size();
        for (int r = 1; r < n; ++r) for (int i = 0; i < n - r; ++i) {
            c[i].x = (1 - t) * c[i].x + t * c[i + 1].x;
            c[i].y = (1 - t) * c[i].y + t * c[i + 1].y;
            c[i].z = (1 - t) * c[i].z + t * c[i + 1].z;
        }
        return c[0];
    }
    P3 deriv(double t) const {
        const int p = (int)P.size() - 1; if (p == 0) return { 0, 0, 0 };
        std::vector<P3> Q(p);
        for (int i = 0; i < p; ++i) Q[i] = { p * (P[i + 1].x - P[i].x), p * (P[i + 1].y - P[i].y), p * (P[i + 1].z - P[i].z) };
        Bezier3 q; q.P = Q; return q.eval(t);
    }
    Int3 derivRange() const {   // 全定义域导数界（导数控制点凸包）
        const int p = (int)P.size() - 1;
        I3 dx(1e300, -1e300), dy(1e300, -1e300), dz(1e300, -1e300);
        for (int i = 0; i < p; ++i) {
            const double qx = p * (P[i + 1].x - P[i].x), qy = p * (P[i + 1].y - P[i].y), qz = p * (P[i + 1].z - P[i].z);
            if (qx < dx.lo) dx.lo = qx; if (qx > dx.hi) dx.hi = qx;
            if (qy < dy.lo) dy.lo = qy; if (qy > dy.hi) dy.hi = qy;
            if (qz < dz.lo) dz.lo = qz; if (qz > dz.hi) dz.hi = qz;
        }
        return { dx, dy, dz };
    }
};

// ---------- 双线性 NurbsSurface 的偏导全定义域界（由控制点凸包给出） ----------
static Int3 derivSRange(const nurbs::NurbsSurface& S) { // ∂S/∂s 界
    // 双线性 ∂S/∂s=(1-t)(P10−P00)+t(P11−P01)；界 = (P10−P00),(P11−P01) 凸包
    P3 A{ S.C[1].x-S.C[0].x, S.C[1].y-S.C[0].y, S.C[1].z-S.C[0].z },
       B{ S.C[3].x-S.C[2].x, S.C[3].y-S.C[2].y, S.C[3].z-S.C[2].z };
    return { { std::min(A.x,B.x), std::max(A.x,B.x) }, { std::min(A.y,B.y), std::max(A.y,B.y) }, { std::min(A.z,B.z), std::max(A.z,B.z) } };
}
static Int3 derivTRange(const nurbs::NurbsSurface& S) { // ∂S/∂t 界
    P3 C{ S.C[2].x-S.C[0].x, S.C[2].y-S.C[0].y, S.C[2].z-S.C[0].z },
       D{ S.C[3].x-S.C[1].x, S.C[3].y-S.C[1].y, S.C[3].z-S.C[1].z };
    return { { std::min(C.x,D.x), std::max(C.x,D.x) }, { std::min(C.y,D.y), std::max(C.y,D.y) }, { std::min(C.z,D.z), std::max(C.z,D.z) } };
}

// 平面双线性 patch 构造（统一 NurbsSurface，P00→C[0], P10→C[1], P01→C[2], P11→C[3]）
static nurbs::NurbsSurface makePatch(P3 p00, P3 p10, P3 p01, P3 p11) {
    SurfacePoint cp4[4] = { p00, p10, p01, p11 };
    return nurbs::NurbsSurface::makeBilinear(cp4);
}

// ---------- 3×3 实矩阵逆 ----------
static bool inv3(const double M[3][3], double R[3][3]) {
    const double d0 = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                    - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                    + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (std::abs(d0) < 1e-14) return false;
    const double inv = 1.0 / d0;
    R[0][0] =  (M[1][1]*M[2][2]-M[1][2]*M[2][1])*inv;
    R[0][1] = -(M[0][1]*M[2][2]-M[0][2]*M[2][1])*inv;
    R[0][2] =  (M[0][1]*M[1][2]-M[0][2]*M[1][1])*inv;
    R[1][0] = -(M[1][0]*M[2][2]-M[1][2]*M[2][0])*inv;
    R[1][1] =  (M[0][0]*M[2][2]-M[0][2]*M[2][0])*inv;
    R[1][2] = -(M[0][0]*M[1][2]-M[0][2]*M[1][0])*inv;
    R[2][0] =  (M[1][0]*M[2][1]-M[1][1]*M[2][0])*inv;
    R[2][1] = -(M[0][0]*M[2][1]-M[0][1]*M[2][0])*inv;
    R[2][2] =  (M[0][0]*M[1][1]-M[0][1]*M[1][0])*inv;
    return true;
}

// 方系统 3D 曲线-曲面求交 Krawczyk：返回 1=证明存在唯一交点，0=未定/不存在
static int proveCS(const Bezier3& C, const nurbs::NurbsSurface& S, I3& ub, I3& sb, I3& tb) {
    const double mu = ub.mid(), ms = sb.mid(), mt = tb.mid();
    const P3 cp = C.eval(mu), sp = S.eval(ms, mt);
    const double f[3] = { cp.x - sp.x, cp.y - sp.y, cp.z - sp.z };
    const P3 cd = C.deriv(mu), ss = S.du(ms, mt), st = S.dv(ms, mt);
    const double Jm[3][3] = { { cd.x, -ss.x, -st.x }, { cd.y, -ss.y, -st.y }, { cd.z, -ss.z, -st.z } };
    double Ji[3][3];
    if (!inv3(Jm, Ji)) return 0;   // 中点 Jacobian 奇异 → 相切/退化，不证
    // J([x]) 3×3 区间
    const Int3 ca = C.derivRange(), sa = derivSRange(S), ta = derivTRange(S);
    const I3 Jint[3][3] = { { ca.x, imul(sa.x, I3(-1,-1)), imul(ta.x, I3(-1,-1)) },
                            { ca.y, imul(sa.y, I3(-1,-1)), imul(ta.y, I3(-1,-1)) },
                            { ca.z, imul(sa.z, I3(-1,-1)), imul(ta.z, I3(-1,-1)) } };
    const I3 dx[3] = { isub(ub, I3(mu, mu)), isub(sb, I3(ms, ms)), isub(tb, I3(mt, mt)) };
    const double m[3] = { mu, ms, mt };
    I3 K[3];
    for (int i = 0; i < 3; ++i) {
        I3 sum(0, 0);
        for (int j = 0; j < 3; ++j) {
            I3 proj = (i == j) ? I3(1, 1) : I3(0, 0);
            for (int k = 0; k < 3; ++k) proj = isub(proj, imul(I3(Ji[i][k], Ji[i][k]), Jint[k][j]));
            sum = iadd(sum, imul(proj, dx[j]));
        }
        double jf = 0; for (int k = 0; k < 3; ++k) jf += Ji[i][k] * f[k];
        K[i] = iadd(I3(m[i] - jf, m[i] - jf), sum);
    }
    const I3 xx[3] = { ub, sb, tb };
    for (int i = 0; i < 3; ++i) if (!(xx[i].lo < K[i].lo && K[i].hi < xx[i].hi)) return 0;
    return 1;
}

// 方系统 3D 曲线-曲面求交 GS-Newton（Hansen-Sengupta 区间 Gauss-Seidel 消元，3×3）：
// 线性化 A(x−m)=−Yf(m)（A=Y·J([x])，Y=Jm⁻¹），逐变量区间消元，利用最新区间，比
// 一步 Krawczyk 更紧。返回 1=证明存在唯一，0=未定，-1=区间变空（无解）。
static int proveCSGS(const Bezier3& C, const nurbs::NurbsSurface& S, I3& ub, I3& sb, I3& tb) {
    const double mu = ub.mid(), ms = sb.mid(), mt = tb.mid();
    const P3 cp = C.eval(mu), sp = S.eval(ms, mt);
    const double f[3] = { cp.x - sp.x, cp.y - sp.y, cp.z - sp.z };
    const P3 cd = C.deriv(mu), ss = S.du(ms, mt), st = S.dv(ms, mt);
    const double Jm[3][3] = { { cd.x, -ss.x, -st.x }, { cd.y, -ss.y, -st.y }, { cd.z, -ss.z, -st.z } };
    double Ji[3][3];
    if (!inv3(Jm, Ji)) return 0;   // 中点 Jacobian 奇异 → 相切/退化，不证
    const Int3 ca = C.derivRange(), sa = derivSRange(S), ta = derivTRange(S);
    const I3 Jint[3][3] = { { ca.x, imul(sa.x, I3(-1,-1)), imul(ta.x, I3(-1,-1)) },
                            { ca.y, imul(sa.y, I3(-1,-1)), imul(ta.y, I3(-1,-1)) },
                            { ca.z, imul(sa.z, I3(-1,-1)), imul(ta.z, I3(-1,-1)) } };
    // A = Y·J([x])（3×3 区间），b = −Y f(m)
    I3 A[3][3]; double b[3];
    for (int i = 0; i < 3; ++i) {
        b[i] = 0;
        for (int k = 0; k < 3; ++k) b[i] -= Ji[i][k] * f[k];
        for (int j = 0; j < 3; ++j) {
            I3 s(0, 0);
            for (int k = 0; k < 3; ++k) s = iadd(s, imul(I3(Ji[i][k], Ji[i][k]), Jint[k][j]));
            A[i][j] = s;
        }
    }
    const double m[3] = { mu, ms, mt };
    I3 x[3] = { ub, sb, tb };       // 当前盒
    I3 y[3];
    for (int i = 0; i < 3; ++i) y[i] = isub(x[i], I3(m[i], m[i]));
    // Gauss-Seidel：逐变量用最新区间
    for (int i = 0; i < 3; ++i) {
        if (A[i][i].lo <= 0.0 && A[i][i].hi >= 0.0) return 0;   // 主元含 0，无法区间消元
        I3 rhs = I3(b[i], b[i]);
        for (int j = 0; j < 3; ++j) if (j != i) rhs = isub(rhs, imul(A[i][j], y[j]));
        I3 xn = iadd(I3(m[i], m[i]), idiv(rhs, A[i][i]));
        xn.lo = std::max(x[i].lo, xn.lo); xn.hi = std::min(x[i].hi, xn.hi);
        if (xn.lo > xn.hi) return -1;   // 无解
        x[i] = xn; y[i] = isub(x[i], I3(m[i], m[i]));
    }
    // 证明：新盒严格内缩到旧盒内部（存在唯一）
    const bool inU = x[0].lo > ub.lo && x[0].hi < ub.hi;
    const bool inS = x[1].lo > sb.lo && x[1].hi < sb.hi;
    const bool inT = x[2].lo > tb.lo && x[2].hi < tb.hi;
    return (inU && inS && inT) ? 1 : 0;
}

// 修正牛顿 + 误差界（interval Newton-corrector，3×3）：
// 在 Krawczyk/GS 证明存在唯一后，对解盒 x 反复施加区间 Newton 步
//   x ← x ∩ { m − Yf(m) + (I − Y J([x]))(x−m) }
// 直到盒不再收缩，所得窄盒 x 即解的严格误差界（比"证明存在"更强）。
static bool correctNewt3(const Bezier3& C, const nurbs::NurbsSurface& S,
                         I3& ub, I3& sb, I3& tb) {
    I3 x[3] = { ub, sb, tb };
    for (int iter = 0; iter < 40; ++iter) {
        const double m[3] = { x[0].mid(), x[1].mid(), x[2].mid() };
        const P3 cp = C.eval(m[0]), sp = S.eval(m[1], m[2]);
        const double f[3] = { cp.x - sp.x, cp.y - sp.y, cp.z - sp.z };
        const P3 cd = C.deriv(m[0]), ss = S.du(m[1], m[2]), st = S.dv(m[1], m[2]);
        const double Jm[3][3] = { { cd.x, -ss.x, -st.x }, { cd.y, -ss.y, -st.y }, { cd.z, -ss.z, -st.z } };
        double Ji[3][3];
        if (!inv3(Jm, Ji)) return false;
        const Int3 ca = C.derivRange(), sa = derivSRange(S), ta = derivTRange(S);
        const I3 Jint[3][3] = { { ca.x, imul(sa.x, I3(-1,-1)), imul(ta.x, I3(-1,-1)) },
                                { ca.y, imul(sa.y, I3(-1,-1)), imul(ta.y, I3(-1,-1)) },
                                { ca.z, imul(sa.z, I3(-1,-1)), imul(ta.z, I3(-1,-1)) } };
        I3 xn[3];
        for (int i = 0; i < 3; ++i) {
            double af = 0; for (int k = 0; k < 3; ++k) af += Ji[i][k] * f[k];
            I3 sum(0, 0);
            for (int j = 0; j < 3; ++j) {
                I3 proj = (i == j) ? I3(1, 1) : I3(0, 0);
                for (int k = 0; k < 3; ++k) proj = isub(proj, imul(I3(Ji[i][k], Ji[i][k]), Jint[k][j]));
                sum = iadd(sum, imul(proj, isub(x[j], I3(m[j], m[j]))));
            }
            xn[i] = iadd(I3(m[i] - af, m[i] - af), sum);
            xn[i].lo = std::max(x[i].lo, xn[i].lo); xn[i].hi = std::min(x[i].hi, xn[i].hi);
            if (xn[i].lo > xn[i].hi) return false;   // 盒变空 → 校正失败
        }
        double shrink = 0;
        for (int i = 0; i < 3; ++i) shrink = std::max(shrink, (x[i].width() - xn[i].width()) / x[i].width());
        x[0] = xn[0]; x[1] = xn[1]; x[2] = xn[2];
        if (shrink < 1e-12) break;
    }
    ub = x[0]; sb = x[1]; tb = x[2];
    return true;
}

static void run(const char* label, const Bezier3& C, const nurbs::NurbsSurface& S,
                double u0, double s0, double t0, double delta, int expect) {
    I3 ub(u0 - delta, u0 + delta), sb(s0 - delta, s0 + delta), tb(t0 - delta, t0 + delta);
    int r = proveCS(C, S, ub, sb, tb);
    printf("  %-30s δ=%.0e → %s（期望 %s）%s\n", label, delta,
           r == 1 ? "已证明" : "未证", expect == 1 ? "证" : "拒", r == expect ? " ✓" : " ✗");
    if (r == expect) ++g_pass; else ++g_fail;
}

int main() {
    printf("==== 扩展：3D 曲线-曲面求交（方系统 R³→R³，Krawczyk 证明，统一 NurbsSurface）====\n\n");

    // 曲面：xz 平面 y=0 的双线性 patch（s∈x∈[-0.5,0.5]，t∈z∈[0,0.5]）
    nurbs::NurbsSurface S = makePatch({ -0.5, 0, 0 }, { 0.5, 0, 0 }, { -0.5, 0, 0.5 }, { 0.5, 0, 0.5 });

    // 场景 1：曲线穿过曲面（横截交点），C(u)=(0, u-0.5, 0.25)，交点 (u,s,t)=(0.5,0.5,0.5)
    {
        printf("[1] 曲线穿过曲面（横截）\n");
        Bezier3 C; C.P = { {0,-0.5,0.25}, {0,0.5,0.25} };
        for (double d : {2e-2, 1e-3, 1e-5}) run("曲线穿面横截", C, S, 0.5, 0.5, 0.5, d, 1);
    }
    // 场景 2：曲线顶点与曲面相切（仅接触一点，C'(0.5)=0 → 奇异）→ 不误证
    {
        printf("[2] 曲线与曲面相切（不误证）\n");
        Bezier3 C; C.P = { {0,0.25,0.25}, {0,-0.25,0.25}, {0,0.25,0.25} };
        run("曲线顶点切面", C, S, 0.5, 0.5, 0.5, 1e-2, 0);
        run("曲线顶点切面(细盒)", C, S, 0.5, 0.5, 0.5, 1e-4, 0);
    }
    // 场景 3：曲线不与曲面相交（y≥0.1，不触 y=0 平面）→ 拒绝
    {
        printf("[3] 曲线不与曲面相交（拒绝）\n");
        Bezier3 C; C.P = { {0,0.1,0.25}, {0,1.0,0.25} };
        run("曲线在面上方", C, S, 0.5, 0.5, 0.5, 1e-2, 0);
    }
    // 场景 4：曲线在曲面内部一段（退化重叠：曲线落在平面 y=0 上）→ 非孤立横截，Krawczyk 不应误证唯一
    {
        printf("[4] 曲线落在曲面内（退化重叠）\n");
        Bezier3 C; C.P = { {-0.4,0,0.25}, {0.4,0,0.25} };   // C=(u',0,0.25) ⊂ 平面 y=0
        run("曲线在面内", C, S, 0.5, 0.5, 0.5, 1e-2, 0);   // 重叠非孤立解 → 不证
    }
    // 场景 5：批量 93 组竖直曲线 C(u)=(0,u-0.5,z0) 穿过 xz 平面（z0∈[0.01,0.49]，横截）
    {
        printf("\n[5] 批量 93 组竖直曲线横截穿面（z0∈[0.01,0.49]）\n");
        const int N = 93; int okN = 0;
        for (int kk = 0; kk < N; ++kk) {
            const double z0 = 0.01 + 0.48 * kk / (N - 1);
            const double t0 = 2 * z0;                 // patch z=0.5t → 曲线 z=z0 交点 t=2z0
            Bezier3 C; C.P = { {0,-0.5,z0}, {0,0.5,z0} };
            I3 ub(0.5-1e-3,0.5+1e-3), sb(0.5-1e-3,0.5+1e-3), tb(t0-1e-3,t0+1e-3);
            if (proveCS(C, S, ub, sb, tb)) ++okN;
        }
        printf("      批量 93/93 曲线穿面横截证明: %d\n", okN);
        g_pass += okN; g_fail += (N - okN);
    }
    // ③ GS-Newton(3×3) vs Krawczyk(3×3) 区间证明对比（δ 扫描）
    {
        printf("\n[6] 3×3 GS-Newton vs Krawczyk 区间证明对比（δ 扫描）\n");
        struct CS { Bezier3 C; nurbs::NurbsSurface S; double u0, s0, t0; };
        std::vector<CS> cases;
        { Bezier3 C; C.P = { {0,-0.5,0.25}, {0,0.5,0.25} }; cases.push_back({ C, S, 0.5, 0.5, 0.5 }); }
        for (int kk = 0; kk < 93; ++kk) {
            const double z0 = 0.01 + 0.48 * kk / (93 - 1);
            const double t0 = 2 * z0;
            Bezier3 C; C.P = { {0,-0.5,z0}, {0,0.5,z0} };
            cases.push_back({ C, S, 0.5, 0.5, t0 });
        }
        printf("横截解总数 %zu\n", cases.size());
        printf("δ      | Krawczyk 覆盖率 | GS 覆盖率 | Krawczyk 耗时(µs) | GS 耗时(µs)\n");
        printf("-------+----------------+-----------+-------------------+-----------\n");
        using clk = std::chrono::steady_clock;
        for (double d : { 2e-2, 1e-3, 1e-4, 1e-5 }) {
            int okK = 0, okG = 0; long long tK = 0, tG = 0;
            for (auto& cs : cases) {
                I3 ub(cs.u0-d, cs.u0+d), sb(cs.s0-d, cs.s0+d), tb(cs.t0-d, cs.t0+d);
                auto a = clk::now(); int pk = proveCS(cs.C, cs.S, ub, sb, tb); tK += (clk::now() - a).count();
                I3 ub2(cs.u0-d, cs.u0+d), sb2(cs.s0-d, cs.s0+d), tb2(cs.t0-d, cs.t0+d);
                auto b2 = clk::now(); int pg = proveCSGS(cs.C, cs.S, ub2, sb2, tb2); tG += (clk::now() - b2).count();
                if (pk == 1) ++okK; if (pg == 1) ++okG;
            }
            printf("%-6.0e | %14.1f%% | %9.1f%% | %17.2f | %10.2f\n",
                   d, 100.0 * okK / cases.size(), 100.0 * okG / cases.size(),
                   (double)tK / cases.size(), (double)tG / cases.size());
        }
    }
    // ④ 修正牛顿 + 误差界（interval Newton-corrector）
    {
        printf("\n[7] 修正牛顿 + 误差界（interval Newton-corrector）：证明存在唯一后给出解的严格窄误差界\n");
        for (int kk = 0; kk < 93; kk += 20) {
            const double z0 = 0.01 + 0.48 * kk / (93 - 1);
            const double t0 = 2 * z0;
            Bezier3 C; C.P = { {0,-0.5,z0}, {0,0.5,z0} };
            I3 ub(0.5-1e-3,0.5+1e-3), sb(0.5-1e-3,0.5+1e-3), tb(t0-1e-3,t0+1e-3);
            if (!proveCS(C, S, ub, sb, tb)) { printf("  z0=%.4f 证明失败，跳过\n", z0); continue; }
            const double w0 = std::max(ub.width(), std::max(sb.width(), tb.width()));
            I3 ub2(0.5-1e-3,0.5+1e-3), sb2(0.5-1e-3,0.5+1e-3), tb2(t0-1e-3,t0+1e-3);
            const bool ok = correctNewt3(C, S, ub2, sb2, tb2);
            if (!ok) { printf("  z0=%.4f 校正失败\n", z0); ++g_fail; continue; }
            const double w1 = std::max(ub2.width(), std::max(sb2.width(), tb2.width()));
            printf("  z0=%.4f  初始盒宽=%.1e → 误差界=%.2e  (收缩 %.1e 倍)\n", z0, w0, w1, w0 / w1);
            printf("         u∈[%.6f,%.6f] s∈[%.6f,%.6f] t∈[%.8f,%.8f]（真值 u=0.5 s=0.5 t=%.4f）\n",
                   ub2.lo, ub2.hi, sb2.lo, sb2.hi, tb2.lo, tb2.hi, t0);
            ++g_pass;
        }
    }
    printf("\n结论：曲线-曲面求交（方系统 R³→R³）由 3×3 Krawczyk 保证；横截证明、相切/退化重叠/非交不误证。\n");
    printf("      曲面表示统一为 nurbs::NurbsSurface（双线性 makeBilinear）。\n");
    printf("      通过 %d，失败 %d。\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
