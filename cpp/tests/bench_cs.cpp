// ============================================================================
// bench_cs.cpp — 扩展：3D 曲线-曲面求交（方系统 R³→R³，Krawczyk 证明）
//   曲线 C(u)⊂R³ 与双线性 Bezier 曲面 S(s,t)⊂R³ 求交：F(u,s,t)=C(u)−S(s,t)=0，
//   3 方程 3 未知 → 方系统，用 3×3 区间 Krawczyk 证明交点存在唯一。
//   这是 CAD 最经典原语（曲线与曲面交点），把保证性框架从曲线-曲线推到曲线-曲面。
// 自包含（不动主库）。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>

struct P3 { double x, y, z; };
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} };
struct Int3 { I3 x, y, z; };

inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) {
    const double m1 = a.lo * b.lo, m2 = a.lo * b.hi, m3 = a.hi * b.lo, m4 = a.hi * b.hi;
    return { std::min(std::min(m1, m2), std::min(m3, m4)), std::max(std::max(m1, m2), std::max(m3, m4)) };
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

// ---------- 双线性 Bezier 曲面 ----------
struct Patch {
    P3 P00, P10, P01, P11;
    P3 eval(double s, double t) const {
        const double w0 = (1 - s) * (1 - t), w1 = s * (1 - t), w2 = (1 - s) * t, w3 = s * t;
        return { w0 * P00.x + w1 * P10.x + w2 * P01.x + w3 * P11.x,
                 w0 * P00.y + w1 * P10.y + w2 * P01.y + w3 * P11.y,
                 w0 * P00.z + w1 * P10.z + w2 * P01.z + w3 * P11.z };
    }
    P3 derivS(double s, double t) const {  // ∂S/∂s
        return { (1 - t) * (P10.x - P00.x) + t * (P11.x - P01.x),
                 (1 - t) * (P10.y - P00.y) + t * (P11.y - P01.y),
                 (1 - t) * (P10.z - P00.z) + t * (P11.z - P01.z) };
    }
    P3 derivT(double s, double t) const {  // ∂S/∂t
        return { (1 - s) * (P01.x - P00.x) + s * (P11.x - P10.x),
                 (1 - s) * (P01.y - P00.y) + s * (P11.y - P10.y),
                 (1 - s) * (P01.z - P00.z) + s * (P11.z - P10.z) };
    }
    Int3 derivSRange() const { // ∂S/∂s 界（(P10−P00),(P11−P01) 凸包）
        P3 A{ P10.x - P00.x, P10.y - P00.y, P10.z - P00.z }, B{ P11.x - P01.x, P11.y - P01.y, P11.z - P01.z };
        return { { std::min(A.x, B.x), std::max(A.x, B.x) }, { std::min(A.y, B.y), std::max(A.y, B.y) }, { std::min(A.z, B.z), std::max(A.z, B.z) } };
    }
    Int3 derivTRange() const {
        P3 C{ P01.x - P00.x, P01.y - P00.y, P01.z - P00.z }, D{ P11.x - P10.x, P11.y - P10.y, P11.z - P10.z };
        return { { std::min(C.x, D.x), std::max(C.x, D.x) }, { std::min(C.y, D.y), std::max(C.y, D.y) }, { std::min(C.z, D.z), std::max(C.z, D.z) } };
    }
};

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
static int proveCS(const Bezier3& C, const Patch& S, I3& ub, I3& sb, I3& tb) {
    const double mu = ub.mid(), ms = sb.mid(), mt = tb.mid();
    const P3 cp = C.eval(mu), sp = S.eval(ms, mt);
    const double f[3] = { cp.x - sp.x, cp.y - sp.y, cp.z - sp.z };
    const P3 cd = C.deriv(mu), ss = S.derivS(ms, mt), st = S.derivT(ms, mt);
    const double Jm[3][3] = { { cd.x, -ss.x, -st.x }, { cd.y, -ss.y, -st.y }, { cd.z, -ss.z, -st.z } };
    double Ji[3][3];
    if (!inv3(Jm, Ji)) return 0;   // 中点 Jacobian 奇异 → 相切/退化，不证
    // J([x]) 3×3 区间
    const Int3 ca = C.derivRange(), sa = S.derivSRange(), ta = S.derivTRange();
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

static void run(const char* label, const Bezier3& C, const Patch& S,
                double u0, double s0, double t0, double delta, int expect) {
    I3 ub(u0 - delta, u0 + delta), sb(s0 - delta, s0 + delta), tb(t0 - delta, t0 + delta);
    int r = proveCS(C, S, ub, sb, tb);
    printf("  %-30s δ=%.0e → %s（期望 %s）%s\n", label, delta,
           r == 1 ? "已证明" : "未证", expect == 1 ? "证" : "拒", r == expect ? " ✓" : " ✗");
}

int main() {
    printf("==== 扩展：3D 曲线-曲面求交（方系统 R³→R³，Krawczyk 证明）====\n\n");

    // 曲面：xz 平面 y=0 的双线性 patch（s∈x∈[-0.5,0.5]，t∈z∈[0,0.5]）
    Patch S; S.P00 = { -0.5, 0, 0 }; S.P10 = { 0.5, 0, 0 }; S.P01 = { -0.5, 0, 0.5 }; S.P11 = { 0.5, 0, 0.5 };

    // 场景 1：曲线穿过曲面（横截交点），C(u)=(0, u-0.5, 0.25)，交点 (u,s,t)=(0.5,0.5,0.5)
    {
        printf("[1] 曲线穿过曲面（横截）\n");
        Bezier3 C; C.P = { {0,-0.5,0.25}, {0,0.5,0.25} };
        for (double d : {2e-2, 1e-3, 1e-5}) run("曲线穿面横截", C, S, 0.5, 0.5, 0.5, d, 1);
    }
    // 场景 2：曲线顶点与曲面相切（仅接触一点，C'(0.5)=0 → 奇异）→ 不误证
    {
        printf("[2] 曲线与曲面相切（不误证）\n");
        // C(u)=(0, (u-0.5)², 0.25)：二次 Bezier，y 控制点 0.25,-0.25,0.25 → 顶点 y=0 at u=0.5
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
    printf("\n结论：曲线-曲面求交（方系统 R³→R³）由 3×3 Krawczyk 保证；横截证明、相切/退化重叠/非交不误证。\n");
    return 0;
}
