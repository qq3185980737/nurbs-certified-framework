// ============================================================================
// bench_massive_surf.cpp — 综合：曲线×曲面 + 曲面×曲面 跨类型相交（批量证明，迁移到统一 NurbsSurface）
//   主库仅支持 NURBSCurve；曲面用统一 NurbsSurface 证明算子（proveCS/proveSS）。
//   [A] 曲线×曲面：3 类 3D 曲线（直线/二次/三次）竖直穿过随机平面 patch，
//       交点参数解析已知 → proveCS（方系统）批量证明，各 200 组。
//   [B] 曲面×曲面：两平面 patch 交线（欠定），随机倾斜 → proveSS 逐点证明。
//   展示不同曲线/曲面类型之间相交的保证性覆盖。
//   曲面统一存储结构：nurbs/NurbsSurface（双线性 makeBilinear）。自包含。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
using SurfacePoint = nurbs::SurfacePoint;
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} };
struct Int3 { I3 x, y, z; };
inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) {
    const double m1 = a.lo * b.lo, m2 = a.lo * b.hi, m3 = a.hi * b.lo, m4 = a.hi * b.hi;
    return { std::min(std::min(m1, m2), std::min(m3, m4)), std::max(std::max(m1, m2), std::max(m3, m4)) };
}

struct Bezier3 {
    std::vector<P3> P;
    P3 eval(double t) const { std::vector<P3> c = P; const int n = (int)c.size();
        for (int r = 1; r < n; ++r) for (int i = 0; i < n - r; ++i) { c[i].x = (1-t)*c[i].x + t*c[i+1].x; c[i].y = (1-t)*c[i].y + t*c[i+1].y; c[i].z = (1-t)*c[i].z + t*c[i+1].z; } return c[0]; }
    P3 deriv(double t) const { const int p = (int)P.size() - 1; if (p == 0) return {0,0,0};
        std::vector<P3> Q(p); for (int i = 0; i < p; ++i) Q[i] = { p*(P[i+1].x-P[i].x), p*(P[i+1].y-P[i].y), p*(P[i+1].z-P[i].z) }; Bezier3 q; q.P = Q; return q.eval(t); }
    Int3 derivRange() const { const int p = (int)P.size() - 1; I3 dx(1e300,-1e300), dy(1e300,-1e300), dz(1e300,-1e300);
        for (int i = 0; i < p; ++i) { const double qx = p*(P[i+1].x-P[i].x), qy = p*(P[i+1].y-P[i].y), qz = p*(P[i+1].z-P[i].z);
            if (qx < dx.lo) dx.lo = qx; if (qx > dx.hi) dx.hi = qx; if (qy < dy.lo) dy.lo = qy; if (qy > dy.hi) dy.hi = qy; if (qz < dz.lo) dz.lo = qz; if (qz > dz.hi) dz.hi = qz; } return { dx, dy, dz }; }
};

// ---------- 双线性 NurbsSurface 的偏导全定义域界（由控制点凸包给出） ----------
static Int3 derivSRange(const nurbs::NurbsSurface& S) {
    P3 A{ S.C[1].x-S.C[0].x, S.C[1].y-S.C[0].y, S.C[1].z-S.C[0].z },
       B{ S.C[3].x-S.C[2].x, S.C[3].y-S.C[2].y, S.C[3].z-S.C[2].z };
    return { { std::min(A.x,B.x), std::max(A.x,B.x) }, { std::min(A.y,B.y), std::max(A.y,B.y) }, { std::min(A.z,B.z), std::max(A.z,B.z) } };
}
static Int3 derivTRange(const nurbs::NurbsSurface& S) {
    P3 C{ S.C[2].x-S.C[0].x, S.C[2].y-S.C[0].y, S.C[2].z-S.C[0].z },
       D{ S.C[3].x-S.C[1].x, S.C[3].y-S.C[1].y, S.C[3].z-S.C[1].z };
    return { { std::min(C.x,D.x), std::max(C.x,D.x) }, { std::min(C.y,D.y), std::max(C.y,D.y) }, { std::min(C.z,D.z), std::max(C.z,D.z) } };
}
static nurbs::NurbsSurface makePatch(P3 p00, P3 p10, P3 p01, P3 p11) {
    SurfacePoint cp4[4] = { p00, p10, p01, p11 };
    return nurbs::NurbsSurface::makeBilinear(cp4);
}

static bool inv3(const double M[3][3], double R[3][3]) {
    const double d0 = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (std::abs(d0) < 1e-14) return false; const double inv = 1.0 / d0;
    R[0][0]=(M[1][1]*M[2][2]-M[1][2]*M[2][1])*inv; R[0][1]=-(M[0][1]*M[2][2]-M[0][2]*M[2][1])*inv; R[0][2]=(M[0][1]*M[1][2]-M[0][2]*M[1][1])*inv;
    R[1][0]=-(M[1][0]*M[2][2]-M[1][2]*M[2][0])*inv; R[1][1]=(M[0][0]*M[2][2]-M[0][2]*M[2][0])*inv; R[1][2]=-(M[0][0]*M[1][2]-M[0][2]*M[1][0])*inv;
    R[2][0]=(M[1][0]*M[2][1]-M[1][1]*M[2][0])*inv; R[2][1]=-(M[0][0]*M[2][1]-M[0][1]*M[2][0])*inv; R[2][2]=(M[0][0]*M[1][1]-M[0][1]*M[1][0])*inv;
    return true;
}

// 曲线×曲面 方系统 R³→R³ Krawczyk
static int proveCS(const Bezier3& C, const nurbs::NurbsSurface& S, I3& ub, I3& sb, I3& tb) {
    const double mu=ub.mid(), ms=sb.mid(), mt=tb.mid();
    const P3 cp=C.eval(mu), sp=S.eval(ms,mt);
    const double f[3] = { cp.x-sp.x, cp.y-sp.y, cp.z-sp.z };
    const P3 cd=C.deriv(mu), ss=S.du(ms,mt), st=S.dv(ms,mt);
    const double Jm[3][3] = { {cd.x,-ss.x,-st.x}, {cd.y,-ss.y,-st.y}, {cd.z,-ss.z,-st.z} };
    double Ji[3][3]; if (!inv3(Jm, Ji)) return 0;
    const Int3 ca=C.derivRange(), sa=derivSRange(S), ta=derivTRange(S);
    const I3 Jint[3][3] = { {ca.x, imul(sa.x,I3(-1,-1)), imul(ta.x,I3(-1,-1))}, {ca.y, imul(sa.y,I3(-1,-1)), imul(ta.y,I3(-1,-1))}, {ca.z, imul(sa.z,I3(-1,-1)), imul(ta.z,I3(-1,-1))} };
    const I3 dx[3] = { isub(ub,I3(mu,mu)), isub(sb,I3(ms,ms)), isub(tb,I3(mt,mt)) };
    const double m[3] = { mu, ms, mt };
    I3 K[3];
    for (int i = 0; i < 3; ++i) { I3 sum(0,0);
        for (int j = 0; j < 3; ++j) { I3 proj = (i==j)?I3(1,1):I3(0,0);
            for (int k = 0; k < 3; ++k) proj = isub(proj, imul(I3(Ji[i][k],Ji[i][k]), Jint[k][j]));
            sum = iadd(sum, imul(proj, dx[j])); }
        double jf = 0; for (int k = 0; k < 3; ++k) jf += Ji[i][k]*f[k];
        K[i] = iadd(I3(m[i]-jf, m[i]-jf), sum); }
    const I3 xx[3] = { ub, sb, tb };
    for (int i = 0; i < 3; ++i) if (!(xx[i].lo < K[i].lo && K[i].hi < xx[i].hi)) return 0;
    return 1;
}
// 曲面×曲面 交线逐点（固定 t1，方系统 R³→R³）
static int proveSS(const nurbs::NurbsSurface& S1, const nurbs::NurbsSurface& S2, double t0, I3& s1b, I3& s2b, I3& t2b) {
    const double m1=s1b.mid(), m2=s2b.mid(), mt=t2b.mid();
    const P3 sp=S1.eval(m1,t0), qp=S2.eval(m2,mt);
    const double f[3] = { sp.x-qp.x, sp.y-qp.y, sp.z-qp.z };
    const P3 s1s=S1.du(m1,t0), s2s=S2.du(m2,mt), s2t=S2.dv(m2,mt);
    const double Jm[3][3] = { {s1s.x,-s2s.x,-s2t.x}, {s1s.y,-s2s.y,-s2t.y}, {s1s.z,-s2s.z,-s2t.z} };
    double Ji[3][3]; if (!inv3(Jm, Ji)) return 0;
    const P3 a1 = S1.du(0.0, t0);
    const Int3 s2a=derivSRange(S2), t2a=derivTRange(S2);
    const I3 Jint[3][3] = { {I3(a1.x,a1.x), imul(s2a.x,I3(-1,-1)), imul(t2a.x,I3(-1,-1))}, {I3(a1.y,a1.y), imul(s2a.y,I3(-1,-1)), imul(t2a.y,I3(-1,-1))}, {I3(a1.z,a1.z), imul(s2a.z,I3(-1,-1)), imul(t2a.z,I3(-1,-1))} };
    const I3 dx[3] = { isub(s1b,I3(m1,m1)), isub(s2b,I3(m2,m2)), isub(t2b,I3(mt,mt)) };
    const double m[3] = { m1, m2, mt };
    I3 K[3];
    for (int i = 0; i < 3; ++i) { I3 sum(0,0);
        for (int j = 0; j < 3; ++j) { I3 proj = (i==j)?I3(1,1):I3(0,0);
            for (int k = 0; k < 3; ++k) proj = isub(proj, imul(I3(Ji[i][k],Ji[i][k]), Jint[k][j]));
            sum = iadd(sum, imul(proj, dx[j])); }
        double jf = 0; for (int k = 0; k < 3; ++k) jf += Ji[i][k]*f[k];
        K[i] = iadd(I3(m[i]-jf, m[i]-jf), sum); }
    const I3 xx[3] = { s1b, s2b, t2b };
    for (int i = 0; i < 3; ++i) if (!(xx[i].lo < K[i].lo && K[i].hi < xx[i].hi)) return 0;
    return 1;
}

uint32_t seed = 20261015u;
double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

int main() {
    printf("==== 综合：曲线×曲面 + 曲面×曲面 跨类型相交（批量证明，统一 NurbsSurface）====\n\n");

    // [A] 曲线×曲面：3 类 3D 曲线竖直穿过平面 patch（y 控制点等差 → 交点 u=0.5）
    {
        printf("[A] 曲线 × 平面 patch（方系统）各类型 200 组\n");
        const char* cname[] = { "3D直线", "3D二次", "3D三次" };
        for (int ct = 0; ct < 3; ++ct) {
            int n = ct + 2;                 // 控制点数 2/3/4
            long long sols = 0, proven = 0;
            for (int k = 0; k < 200; ++k) {
                nurbs::NurbsSurface S = makePatch({-1,0,0},{1,0,0},{-1,0,1},{1,0,1});
                double x0 = -0.7 + 1.4 * rnd(), z0 = 0.1 + 0.8 * rnd();
                double A = 0.3 + 1.2 * rnd(), dx = 0.2 + 0.6 * rnd(), dz = 0.1 + 0.4 * rnd();
                Bezier3 C; C.P.resize(n);
                for (int i = 0; i < n; ++i) {
                    double u1 = (double)i / (n - 1);
                    double yc = A * (u1 - 0.5);
                    double xc = x0 + dx * (2.0 * u1 - 1.0);
                    double zc = z0 + dz * (2.0 * u1 - 1.0);
                    C.P[i] = { xc, yc, zc };
                }
                double s = (x0 + 1.0) / 2.0, t = z0;
                I3 ub(0.5 - 1e-3, 0.5 + 1e-3), sb(s - 1e-3, s + 1e-3), tb(t - 1e-3, t + 1e-3);
                if (proveCS(C, S, ub, sb, tb)) proven++;
                sols++;
            }
            printf("  %-8s : 200 组，证明 %3lld → 覆盖率 %.1f%%\n", cname[ct], proven, 100.0 * proven / sols);
        }
    }

    // [B] 曲面×曲面：两平面 patch 交线（欠定），随机倾斜 S2，逐 t1 证明
    {
        printf("\n[B] 曲面 × 曲面（交线逐点证明，欠定系统）随机 20 组\n");
        long long pts = 0, proven = 0;
        for (int k = 0; k < 20; ++k) {
            nurbs::NurbsSurface S1 = makePatch({-0.5,0,0},{0.5,0,0},{-0.5,0,0.5},{0.5,0,0.5});
            double m = 0.3 + 1.0 * rnd(), b;
            double s2star = 0.2 + 0.6 * rnd(); b = -m * s2star;
            double xb = -0.5, xt = 0.5;
            nurbs::NurbsSurface S2 = makePatch({ xb, m*0 + b, 0 }, { xt, m*1 + b, 0 }, { xb, m*0 + b, 0.5 }, { xt, m*1 + b, 0.5 });
            for (double t1 : {0.2, 0.4}) {
                I3 s1b(s2star-1e-3, s2star+1e-3), s2b(s2star-1e-3, s2star+1e-3), t2b(t1-1e-3, t1+1e-3);
                if (proveSS(S1, S2, t1, s1b, s2b, t2b)) proven++;
                pts++;
            }
        }
        printf("  交线采样点证明 %lld/%lld → 覆盖率 %.1f%%\n", proven, pts, 100.0 * proven / (pts ? pts : 1));
    }

    printf("\n结论：不同曲线/曲面类型之间相交（曲线×曲面、曲面×曲面）由统一 NurbsSurface 证明算子批量保证，\n");
    printf("      覆盖 3D 直线/二次/三次 × 平面 patch 及两平面交线。\n");
    return 0;
}
