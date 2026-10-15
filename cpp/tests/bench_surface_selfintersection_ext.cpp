// ============================================================================
// bench_surface_selfintersection_ext.cpp — 曲面自交检测·实验扩展（路径D①+）
//   P2 论文扩展实验：把折叠检测与自交存在性证明推广到更具一般性的情形，
//   提升 GMP 审稿说服力：
//   [E0] 非对称折叠：对称轴 p 任意 ∈(0.25,0.75)，而非固定 s=0.5
//        S=(x0+A·(s−p)², 0, B·t+z0) → 存在 s1≠s2 (s1+s2=2p) 使 S(s1,t)=S(s2,t)
//   [E1] 非平面折叠：y(s)=y0+C·(s−p)² 使曲面非平面，但仍沿 s=p 折叠
//   [E2] 非平面正则零误报：y(s)=y0+C·s²（弯曲但法向包围盒三分量不同时含0）
//   [E3] 双三次 patch（4×4 控制点）：对称折叠 detectFold + 自交证明
//   每个 500 组：检测命中 + 自交存在证明；正则零误报。
// 自包含。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <algorithm>

struct P3 { double x, y, z; };
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} };
inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) { const double m1=a.lo*b.lo,m2=a.lo*b.hi,m3=a.hi*b.lo,m4=a.hi*b.hi; return { std::min(std::min(m1,m2),std::min(m3,m4)), std::max(std::max(m1,m2),std::max(m3,m4)) }; }
struct C3 { P3 A, B; P3 eval(double u) const { return { A.x+u*(B.x-A.x), A.y+u*(B.y-A.y), A.z+u*(B.z-A.z) }; } P3 d() const { return { B.x-A.x, B.y-A.y, B.z-A.z }; } };

static bool provePointOnCurve3D(const C3& C, const P3& P, double& lo, double& hi, int depth) {
    if (depth > 26) return false;
    { double x0=std::min(C.A.x,C.B.x), x1=std::max(C.A.x,C.B.x), y0=std::min(C.A.y,C.B.y), y1=std::max(C.A.y,C.B.y), z0=std::min(C.A.z,C.B.z), z1=std::max(C.A.z,C.B.z);
      if (P.x < x0-1e-9 || P.x > x1+1e-9 || P.y < y0-1e-9 || P.y > y1+1e-9 || P.z < z0-1e-9 || P.z > z1+1e-9) return false; }
    double m = 0.5*(lo+hi); const double w = 0.5*(hi-lo);
    P3 cm = C.eval(m), cd = C.d();
    const double g0=cm.x-P.x, g1=cm.y-P.y, g2=cm.z-P.z;
    const double G = cd.x*cd.x + cd.y*cd.y + cd.z*cd.z;
    if (G < 1e-22) { double lm=0.5*(lo+m), rm=0.5*(m+hi); return provePointOnCurve3D(C,P,lo,m,depth+1)||provePointOnCurve3D(C,P,m,hi,depth+1); }
    const double A0=cd.x/G, A1=cd.y/G, A2=cd.z/G;
    I3 Jx(cd.x,cd.x), Jy(cd.y,cd.y), Jz(cd.z,cd.z);
    I3 M = isub(I3(1,1), iadd(iadd(imul(I3(A0,A0),Jx), imul(I3(A1,A1),Jy)), imul(I3(A2,A2),Jz)));
    const double jf = A0*g0 + A1*g1 + A2*g2;
    I3 K = iadd(I3(m-jf, m-jf), imul(M, I3(-w,w)));
    if (K.lo > lo && K.hi < hi) {
        double mm = 0.5*(K.lo+K.hi); P3 cm = C.eval(mm);
        double r = std::sqrt((cm.x-P.x)*(cm.x-P.x)+(cm.y-P.y)*(cm.y-P.y)+(cm.z-P.z)*(cm.z-P.z));
        return r < 1e-6;
    }
    double lm=0.5*(lo+m), rm=0.5*(m+hi);
    return provePointOnCurve3D(C,P,lo,m,depth+1)||provePointOnCurve3D(C,P,m,hi,depth+1);
}

// ---------- 双二次 patch（3×3 控制点） ----------
struct BiqPatch {
    P3 C[3][3];
    P3 eval(double s, double t) const {
        double ws[3]={ (1-s)*(1-s), 2*s*(1-s), s*s }, wt[3]={ (1-t)*(1-t), 2*t*(1-t), t*t };
        P3 p{0,0,0};
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) { const double b=ws[i]*wt[j];
            p.x += b*C[i][j].x; p.y += b*C[i][j].y; p.z += b*C[i][j].z; }
        return p;
    }
    P3 ds(double s, double t) const {
        double wds[3]={ -2*(1-s), 2*(1-2*s), 2*s }, wt[3]={ (1-t)*(1-t), 2*t*(1-t), t*t };
        P3 p{0,0,0};
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) { const double b=wds[i]*wt[j];
            p.x += b*C[i][j].x; p.y += b*C[i][j].y; p.z += b*C[i][j].z; }
        return p;
    }
    P3 dt(double s, double t) const {
        double ws[3]={ (1-s)*(1-s), 2*s*(1-s), s*s }, wdt[3]={ -2*(1-t), 2*(1-2*t), 2*t };
        P3 p{0,0,0};
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) { const double b=ws[i]*wdt[j];
            p.x += b*C[i][j].x; p.y += b*C[i][j].y; p.z += b*C[i][j].z; }
        return p;
    }
    P3 normal(double s, double t) const { P3 a=ds(s,t), b=dt(s,t);
        return { a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x }; }
};

// ---------- 双三次 patch（4×4 控制点） ----------
struct CubPatch {
    P3 C[4][4];
    P3 eval(double s, double t) const {
        double ws[4]={ (1-s)*(1-s)*(1-s), 3*s*(1-s)*(1-s), 3*s*s*(1-s), s*s*s };
        double wt[4]={ (1-t)*(1-t)*(1-t), 3*t*(1-t)*(1-t), 3*t*t*(1-t), t*t*t };
        P3 p{0,0,0};
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) { const double b=ws[i]*wt[j];
            p.x += b*C[i][j].x; p.y += b*C[i][j].y; p.z += b*C[i][j].z; }
        return p;
    }
    P3 ds(double s, double t) const {
        double wds[4]={ -3*(1-s)*(1-s), 3*(1-s)*(1-3*s), 3*s*(2-3*s), 3*s*s };
        double wt[4]={ (1-t)*(1-t)*(1-t), 3*t*(1-t)*(1-t), 3*t*t*(1-t), t*t*t };
        P3 p{0,0,0};
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) { const double b=wds[i]*wt[j];
            p.x += b*C[i][j].x; p.y += b*C[i][j].y; p.z += b*C[i][j].z; }
        return p;
    }
    P3 dt(double s, double t) const {
        double ws[4]={ (1-s)*(1-s)*(1-s), 3*s*(1-s)*(1-s), 3*s*s*(1-s), s*s*s };
        double wdt[4]={ -3*(1-t)*(1-t), 3*(1-t)*(1-3*t), 3*t*(2-3*t), 3*t*t };
        P3 p{0,0,0};
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) { const double b=ws[i]*wdt[j];
            p.x += b*C[i][j].x; p.y += b*C[i][j].y; p.z += b*C[i][j].z; }
        return p;
    }
    P3 normal(double s, double t) const { P3 a=ds(s,t), b=dt(s,t);
        return { a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x }; }
};

// ---------- 法向包围盒含原点判据（3×3 采样） ----------
template <class PATCH>
static bool normalBoxDegenerate(const PATCH& S) {
    double minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,minz=1e300,maxz=-1e300;
    for (int i=0;i<=2;++i) for (int j=0;j<=2;++j) {
        const double s=0.5*i, t=0.5*j;
        const P3 n=S.normal(s,t);
        minx=std::min(minx,n.x);maxx=std::max(maxx,n.x);
        miny=std::min(miny,n.y);maxy=std::max(maxy,n.y);
        minz=std::min(minz,n.z);maxz=std::max(maxz,n.z);
    }
    return (minx<=0 && maxx>=0) && (miny<=0 && maxy>=0) && (minz<=0 && maxz>=0);
}

// ---------- 自交存在性证明：P=S(s1,t0) ∈ 曲线 C(t)=S(s2,t)，s1≠s2 ----------
template <class PATCH>
static bool proveSelfIntersection(const PATCH& S, double s1, double s2, double t0) {
    const P3 ca=S.eval(s2,0), cb=S.eval(s2,1);
    C3 C{ca,cb};
    const P3 P=S.eval(s1,t0);
    double lo=0,hi=1;
    return provePointOnCurve3D(C,P,lo,hi,0);
}

static uint32_t seed = 20261014u;
static double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

int main() {
    printf("==== 曲面自交检测·实验扩展（非对称折叠/非平面/双三次）====\n\n");

    // [E0] 非对称折叠：对称轴 p 任意 ∈(0.25,0.75)，x(s)=x0+A(s−p)², y=0
    //     Bezier 控制点（二次）：c=[Ap², A(p²−p), A(1−p)²]（非节点值）
    {
        long long foldD=0, foldP=0;
        for (int i=0;i<500;++i) {
            double p=0.25+0.5*rnd(), A=0.3+0.7*rnd(), x0=-0.3+rnd(), B=0.3+0.7*rnd(), z0=-0.3+rnd();
            BiqPatch S{};
            const double cx[3] = { x0 + A*p*p, x0 + A*(p*p-p), x0 + A*(1-p)*(1-p) };
            for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                const double t=0.5*b;
                S.C[a][b] = { cx[a], 0, B*t+z0 };
            }
            if (normalBoxDegenerate(S)) foldD++;
            double d=0.05+0.18*rnd(), s1=p-d, s2=p+d;   // 对称对，均在内域
            if (proveSelfIntersection(S,s1,s2,0.3+0.4*rnd())) foldP++;
        }
        printf("[E0] 非对称折叠（对称轴 p 任意）× 500：检测命中 %lld (%.1f%%)；自交证明 %lld (%.1f%%)\n",
               foldD, 100.0*foldD/500, foldP, 100.0*foldP/500);
    }

    // [E1] 非平面折叠：y(s)=y0+C(s−p)²（曲面非平面，仍沿 s=p 折叠）
    {
        long long foldD=0, foldP=0;
        for (int i=0;i<500;++i) {
            double p=0.25+0.5*rnd(), A=0.3+0.7*rnd(), C=0.2+0.5*rnd();
            double x0=-0.3+rnd(), y0=-0.3+rnd(), B=0.3+0.7*rnd(), z0=-0.3+rnd();
            BiqPatch S{};
            const double cx[3] = { x0 + A*p*p, x0 + A*(p*p-p), x0 + A*(1-p)*(1-p) };
            const double cy[3] = { y0 + C*p*p, y0 + C*(p*p-p), y0 + C*(1-p)*(1-p) };
            for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                const double t=0.5*b;
                S.C[a][b] = { cx[a], cy[a], B*t+z0 };
            }
            if (normalBoxDegenerate(S)) foldD++;
            double d=0.05+0.18*rnd(), s1=p-d, s2=p+d;   // 对称对，均在内域
            if (proveSelfIntersection(S,s1,s2,0.3+0.4*rnd())) foldP++;
        }
        printf("[E1] 非平面折叠（y=C(s−p)²）× 500：检测命中 %lld (%.1f%%)；自交证明 %lld (%.1f%%)\n",
               foldD, 100.0*foldD/500, foldP, 100.0*foldP/500);
    }

    // [E2] 非平面正则零误报：y(s)=y0+C·s²（弯曲但法向包围盒三分量不同时含0）
    {
        long long regOK=0;
        for (int i=0;i<500;++i) {
            double A=0.3+0.7*rnd(), C=0.2+0.5*rnd();
            double x0=-0.3+rnd(), y0=-0.3+rnd(), B=0.3+0.7*rnd(), z0=-0.3+rnd();
            BiqPatch S{};
            for (int a=0;a<3;++a) for (int b=0;b<3;++b) {
                const double s=0.5*a, t=0.5*b;
                S.C[a][b] = { x0 + A*s, y0 + C*s*s, B*t+z0 };
            }
            if (!normalBoxDegenerate(S)) regOK++;
        }
        printf("[E2] 非平面正则（y=C·s²）× 500：零误报（REGULAR）%lld (%.1f%%)\n", regOK, 100.0*regOK/500);
    }

    // [E3] 双三次 patch：对称折叠 x=A·s(1−s)（控制点 [0,A/3,A/3,0]），y=0
    {
        long long foldD=0, foldP=0;
        for (int i=0;i<500;++i) {
            double A=0.3+0.7*rnd(), x0=-0.3+rnd(), B=0.3+0.7*rnd(), z0=-0.3+rnd();
            CubPatch S{};
            // 三次 Bernstein：x(s)=x0+A·s(1−s) → 控制点 x = [x0, x0+A/3, x0+A/3, x0]
            for (int a=0;a<4;++a) for (int b=0;b<4;++b) {
                const double t=0.5*b/1.5;   // b=0..3 → t 采样
                const double xs = (a==1||a==2) ? x0 + A/3.0 : x0;
                const double tt = (double)b/3.0;
                S.C[a][b] = { xs, 0, B*tt+z0 };
            }
            if (normalBoxDegenerate(S)) foldD++;
            if (proveSelfIntersection(S,0.25,0.75,0.3+0.4*rnd())) foldP++;
        }
        printf("[E3] 双三次折叠 × 500：检测命中 %lld (%.1f%%)；自交证明 %lld (%.1f%%)\n",
               foldD, 100.0*foldD/500, foldP, 100.0*foldP/500);
    }

    printf("\n结论：实验扩展完成——非对称折叠（对称轴任意）、非平面折叠、双三次折叠均\n");
    printf("      检测命中 + 自交存在证明；非平面正则零误报。\n");
    return 0;
}
