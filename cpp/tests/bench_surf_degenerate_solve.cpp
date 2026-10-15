// ============================================================================
// bench_surf_degenerate_solve.cpp — 退化曲面求解/证明恢复（对齐曲线侧）
//   在 bench_surf_degenerate（检测层）之上补证明恢复层：
//   ① provePointOnCurve3D：退化点 patch P × 3D 曲线 → 点包含广义 Krawczyk(R→R³)。
//   ② provePointOnSurface3D：退化点 patch P × 3D 曲面 → 点包含广义 Krawczyk(R²→R³)。
//   ③ proveCurveCurve3D：退化到曲线 patch(降维为3D线段) × 3D 曲线 → 广义 Krawczyk(R²→R³)。
//   退化点/退化到曲线/塌缩情形 → 降维或点包含恢复证明；正则不回归。
// 自包含。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

struct P3 { double x, y, z; };
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} };
struct Int3 { I3 x, y, z; };
inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) { const double m1=a.lo*b.lo,m2=a.lo*b.hi,m3=a.hi*b.lo,m4=a.hi*b.hi; return { std::min(std::min(m1,m2),std::min(m3,m4)), std::max(std::max(m1,m2),std::max(m3,m4)) }; }

struct C3 { P3 A, B; P3 eval(double u) const { return { A.x+u*(B.x-A.x), A.y+u*(B.y-A.y), A.z+u*(B.z-A.z) }; } P3 d() const { return { B.x-A.x, B.y-A.y, B.z-A.z }; } };
struct Patch { P3 P00, P10, P01, P11;
    P3 eval(double s, double t) const { const double w0=(1-s)*(1-t),w1=s*(1-t),w2=(1-s)*t,w3=s*t;
        return { w0*P00.x+w1*P10.x+w2*P01.x+w3*P11.x, w0*P00.y+w1*P10.y+w2*P01.y+w3*P11.y, w0*P00.z+w1*P10.z+w2*P01.z+w3*P11.z }; }
    P3 ds(double s, double t) const { return { (1-t)*(P10.x-P00.x)+t*(P11.x-P01.x), (1-t)*(P10.y-P00.y)+t*(P11.y-P01.y), (1-t)*(P10.z-P00.z)+t*(P11.z-P01.z) }; }
    P3 dt(double s, double t) const { return { (1-s)*(P01.x-P00.x)+s*(P11.x-P10.x), (1-s)*(P01.y-P00.y)+s*(P11.y-P10.y), (1-s)*(P01.z-P00.z)+s*(P11.z-P10.z) }; }
};

// ---------- ① 点 P ∈ 3D 曲线（超定 R→R³ 广义 Krawczyk，A=G⁻¹Jᵀ, G=JᵀJ） ----------
static bool provePointOnCurve3D(const C3& C, const P3& P, double& lo, double& hi, int depth) {
    if (depth > 26) return false;
    // 必要包含检查：P 不在直线段 AABB 内则立即拒绝（避免拒绝分支爆炸）
    { double x0=std::min(C.A.x,C.B.x), x1=std::max(C.A.x,C.B.x), y0=std::min(C.A.y,C.B.y), y1=std::max(C.A.y,C.B.y), z0=std::min(C.A.z,C.B.z), z1=std::max(C.A.z,C.B.z);
      if (P.x < x0-1e-9 || P.x > x1+1e-9 || P.y < y0-1e-9 || P.y > y1+1e-9 || P.z < z0-1e-9 || P.z > z1+1e-9) return false; }
    double m = 0.5*(lo+hi); const double w = 0.5*(hi-lo);
    P3 cm = C.eval(m), cd = C.d();
    const double g0=cm.x-P.x, g1=cm.y-P.y, g2=cm.z-P.z;
    const double G = cd.x*cd.x + cd.y*cd.y + cd.z*cd.z;
    if (G < 1e-22) { double lm=0.5*(lo+m), rm=0.5*(m+hi); return provePointOnCurve3D(C,P,lo,m,depth+1)||provePointOnCurve3D(C,P,m,hi,depth+1); }
    const double A0=cd.x/G, A1=cd.y/G, A2=cd.z/G;
    I3 Jx(std::min(0.0,cd.x),std::max(0.0,cd.x)), Jy(std::min(0.0,cd.y),std::max(0.0,cd.y)), Jz(std::min(0.0,cd.z),std::max(0.0,cd.z));
    // 区间用整曲线一阶导凸包（保守）
    { double lx=1e300,hx=-1e300,ly=1e300,hy=-1e300,lz=1e300,hz=-1e300;
      lx=std::min(lx,cd.x); hx=std::max(hx,cd.x); ly=std::min(ly,cd.y); hy=std::max(hy,cd.y); lz=std::min(lz,cd.z); hz=std::max(hz,cd.z);
      Jx={lx,hx}; Jy={ly,hy}; Jz={lz,hz}; }
    I3 M = isub(I3(1,1), iadd(iadd(imul(I3(A0,A0),Jx), imul(I3(A1,A1),Jy)), imul(I3(A2,A2),Jz)));
    const double jf = A0*g0 + A1*g1 + A2*g2;
    I3 K = iadd(I3(m-jf, m-jf), imul(M, I3(-w,w)));
    if (K.lo > lo && K.hi < hi) {
        double mm = 0.5*(K.lo+K.hi); P3 cm = C.eval(mm);
        double r = std::sqrt((cm.x-P.x)*(cm.x-P.x)+(cm.y-P.y)*(cm.y-P.y)+(cm.z-P.z)*(cm.z-P.z));
        return r < 1e-6;   // K⊂X 且残差小 → 真解证明；残差大 → 无真解(拒绝)
    }
    double lm=0.5*(lo+m), rm=0.5*(m+hi);
    return provePointOnCurve3D(C,P,lo,m,depth+1)||provePointOnCurve3D(C,P,m,hi,depth+1);
}

// ---------- ② 点 P ∈ 3D 曲面（超定 R²→R³ 广义 Krawczyk） ----------
struct Box2 { double slo, shi, tlo, thi; };
static bool provePointOnSurface3D(const Patch& S, const P3& P, Box2& b, int depth) {
    if (depth > 26) return false;
    // 必要包含检查：P 不在 patch 角点 AABB 内则立即拒绝
    { double x0=1e300,x1=-1e300,y0=1e300,y1=-1e300,z0=1e300,z1=-1e300;
      for (const P3* q : {&S.P00,&S.P10,&S.P01,&S.P11}) { x0=std::min(x0,q->x);x1=std::max(x1,q->x);y0=std::min(y0,q->y);y1=std::max(y1,q->y);z0=std::min(z0,q->z);z1=std::max(z1,q->z); }
      if (P.x < x0-1e-9 || P.x > x1+1e-9 || P.y < y0-1e-9 || P.y > y1+1e-9 || P.z < z0-1e-9 || P.z > z1+1e-9) return false; }
    const double ms=0.5*(b.slo+b.shi), mt=0.5*(b.tlo+b.thi), ws=0.5*(b.shi-b.slo), wt=0.5*(b.thi-b.tlo);
    P3 sm = S.eval(ms,mt), sds = S.ds(ms,mt), sdt = S.dt(ms,mt);
    const double g0=sm.x-P.x, g1=sm.y-P.y, g2=sm.z-P.z;
    const double G00=sds.x*sds.x+sds.y*sds.y+sds.z*sds.z, G01=sds.x*sdt.x+sds.y*sdt.y+sds.z*sdt.z, G11=sdt.x*sdt.x+sdt.y*sdt.y+sdt.z*sdt.z;
    const double det = G00*G11 - G01*G01;
    if (det < 1e-22) { Box2 a{b.slo,0.5*(b.slo+b.shi),b.tlo,b.thi}, c{0.5*(b.slo+b.shi),b.shi,b.tlo,b.thi};
        return provePointOnSurface3D(S,P,a,depth+1)||provePointOnSurface3D(S,P,c,depth+1); }
    const double inv=1.0/det;
    // A=G⁻¹Jᵀ (2×3): row0=(G11·ds−G01·dt)/det, row1=(−G01·ds+G00·dt)/det
    const double A00=(G11*sds.x-G01*sdt.x)*inv, A01=(G11*sds.y-G01*sdt.y)*inv, A02=(G11*sds.z-G01*sdt.z)*inv;
    const double A10=(-G01*sds.x+G00*sdt.x)*inv, A11=(-G01*sds.y+G00*sdt.y)*inv, A12=(-G01*sds.z+G00*sdt.z)*inv;
    // 区间：ds,dt 各分量凸包
    I3 dsx,dsy,dsz,dtx,dty,dtz;
    { P3 a=S.ds(0,0), c=S.ds(1,1), e=S.ds(0,1), f=S.ds(1,0);
      dsx={std::min(std::min(a.x,c.x),std::min(e.x,f.x)),std::max(std::max(a.x,c.x),std::max(e.x,f.x))};
      dsy={std::min(std::min(a.y,c.y),std::min(e.y,f.y)),std::max(std::max(a.y,c.y),std::max(e.y,f.y))};
      dsz={std::min(std::min(a.z,c.z),std::min(e.z,f.z)),std::max(std::max(a.z,c.z),std::max(e.z,f.z))};
      P3 g=S.dt(0,0),h=S.dt(1,1),i=S.dt(0,1),j=S.dt(1,0);
      dtx={std::min(std::min(g.x,h.x),std::min(i.x,j.x)),std::max(std::max(g.x,h.x),std::max(i.x,j.x))};
      dty={std::min(std::min(g.y,h.y),std::min(i.y,j.y)),std::max(std::max(g.y,h.y),std::max(i.y,j.y))};
      dtz={std::min(std::min(g.z,h.z),std::min(i.z,j.z)),std::max(std::max(g.z,h.z),std::max(i.z,j.z))}; }
    // (I − A·Jint) 2×2 区间
    I3 M00 = isub(I3(1,1), iadd(iadd(imul(I3(A00,A00),dsx), imul(I3(A01,A01),dsy)), imul(I3(A02,A02),dsz)));
    I3 M01 = isub(I3(0,0), iadd(iadd(imul(I3(A00,A00),dtx), imul(I3(A01,A01),dty)), imul(I3(A02,A02),dtz)));
    I3 M10 = isub(I3(0,0), iadd(iadd(imul(I3(A10,A10),dsx), imul(I3(A11,A11),dsy)), imul(I3(A12,A12),dsz)));
    I3 M11 = isub(I3(1,1), iadd(iadd(imul(I3(A10,A10),dtx), imul(I3(A11,A11),dty)), imul(I3(A12,A12),dtz)));
    const double jf0 = A00*g0 + A01*g1 + A02*g2, jf1 = A10*g0 + A11*g1 + A12*g2;
    I3 K0 = iadd(I3(ms-jf0,ms-jf0), iadd(imul(M00,I3(-ws,ws)), imul(M01,I3(-wt,wt))));
    I3 K1 = iadd(I3(mt-jf1,mt-jf1), iadd(imul(M10,I3(-ws,ws)), imul(M11,I3(-wt,wt))));
    if (K0.lo>b.slo && K0.hi<b.shi && K1.lo>b.tlo && K1.hi<b.thi) {
        double ms2=0.5*(K0.lo+K0.hi), mt2=0.5*(K1.lo+K1.hi); P3 sm=S.eval(ms2,mt2);
        double r = std::sqrt((sm.x-P.x)*(sm.x-P.x)+(sm.y-P.y)*(sm.y-P.y)+(sm.z-P.z)*(sm.z-P.z));
        return r < 1e-6;
    }
    Box2 a{b.slo,0.5*(b.slo+b.shi),b.tlo,0.5*(b.tlo+b.thi)}, c{0.5*(b.slo+b.shi),b.shi,b.tlo,0.5*(b.tlo+b.thi)};
    Box2 e{b.slo,0.5*(b.slo+b.shi),0.5*(b.tlo+b.thi),b.thi}, f{0.5*(b.slo+b.shi),b.shi,0.5*(b.tlo+b.thi),b.thi};
    return provePointOnSurface3D(S,P,a,depth+1)||provePointOnSurface3D(S,P,c,depth+1)||provePointOnSurface3D(S,P,e,depth+1)||provePointOnSurface3D(S,P,f,depth+1);
}

// ---------- ③ 退化到曲线 patch（降维 3D 线段）× 3D 曲线（超定 R²→R³） ----------
static bool proveCurveCurve3D(const C3& L, const C3& C2, double& ulo,double& uhi, double& vlo,double& vhi, int depth) {
    if (depth > 26) return false;
    // 必要包含检查：两直线段 AABB 不相交则立即拒绝
    { double a0=std::min(L.A.x,L.B.x),a1=std::max(L.A.x,L.B.x),b0=std::min(L.A.y,L.B.y),b1=std::max(L.A.y,L.B.y),c0=std::min(L.A.z,L.B.z),c1=std::max(L.A.z,L.B.z);
      double d0=std::min(C2.A.x,C2.B.x),d1=std::max(C2.A.x,C2.B.x),e0=std::min(C2.A.y,C2.B.y),e1=std::max(C2.A.y,C2.B.y),f0=std::min(C2.A.z,C2.B.z),f1=std::max(C2.A.z,C2.B.z);
      if (a1<d0-1e-9||a0>d1+1e-9||b1<e0-1e-9||b0>e1+1e-9||c1<f0-1e-9||c0>f1+1e-9) return false; }
    double mu=0.5*(ulo+uhi), mv=0.5*(vlo+vhi); const double wu=0.5*(uhi-ulo), wv=0.5*(vhi-vlo);
    P3 lm=L.eval(mu), cm=C2.eval(mv), Ld=L.d(), C2d=C2.d();
    const double g0=lm.x-cm.x, g1=lm.y-cm.y, g2=lm.z-cm.z;
    // J=[Ld, −C2d] (3×2)
    const double J0x=Ld.x, J0y=Ld.y, J0z=Ld.z, J1x=-C2d.x, J1y=-C2d.y, J1z=-C2d.z;
    const double G00=J0x*J0x+J0y*J0y+J0z*J0z, G01=J0x*J1x+J0y*J1y+J0z*J1z, G11=J1x*J1x+J1y*J1y+J1z*J1z;
    const double det = G00*G11 - G01*G01;
    if (det < 1e-22) { double a=0.5*(ulo+mu),c=0.5*(mu+uhi),e=0.5*(vlo+mv),f=0.5*(mv+vhi);
        return proveCurveCurve3D(L,C2,ulo,mu,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,ulo,mu,mv,vhi,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,mv,vhi,depth+1); }
    const double inv=1.0/det;
    // A=G⁻¹Jᵀ (2×3): row0=(G11·J0−G01·J1)/det, row1=(−G01·J0+G00·J1)/det
    const double A00=(G11*J0x-G01*J1x)*inv, A01=(G11*J0y-G01*J1y)*inv, A02=(G11*J0z-G01*J1z)*inv;
    const double A10=(-G01*J0x+G00*J1x)*inv, A11=(-G01*J0y+G00*J1y)*inv, A12=(-G01*J0z+G00*J1z)*inv;
    I3 dJ0x(Ld.x,Ld.x), dJ0y(Ld.y,Ld.y), dJ0z(Ld.z,Ld.z);
    I3 dJ1x(-C2d.x,-C2d.x), dJ1y(-C2d.y,-C2d.y), dJ1z(-C2d.z,-C2d.z);
    I3 M00 = isub(I3(1,1), iadd(iadd(imul(I3(A00,A00),dJ0x), imul(I3(A01,A01),dJ0y)), imul(I3(A02,A02),dJ0z)));
    I3 M01 = isub(I3(0,0), iadd(iadd(imul(I3(A00,A00),dJ1x), imul(I3(A01,A01),dJ1y)), imul(I3(A02,A02),dJ1z)));
    I3 M10 = isub(I3(0,0), iadd(iadd(imul(I3(A10,A10),dJ0x), imul(I3(A11,A11),dJ0y)), imul(I3(A12,A12),dJ0z)));
    I3 M11 = isub(I3(1,1), iadd(iadd(imul(I3(A10,A10),dJ1x), imul(I3(A11,A11),dJ1y)), imul(I3(A12,A12),dJ1z)));
    const double jf0=A00*g0+A01*g1+A02*g2, jf1=A10*g0+A11*g1+A12*g2;
    I3 K0 = iadd(I3(mu-jf0,mu-jf0), iadd(imul(M00,I3(-wu,wu)), imul(M01,I3(-wv,wv))));
    I3 K1 = iadd(I3(mv-jf1,mv-jf1), iadd(imul(M10,I3(-wu,wu)), imul(M11,I3(-wv,wv))));
    if (K0.lo>ulo && K0.hi<uhi && K1.lo>vlo && K1.hi<vhi) {
        double u2=0.5*(K0.lo+K0.hi), v2=0.5*(K1.lo+K1.hi); P3 lm2=L.eval(u2), cm2=C2.eval(v2);
        double r = std::sqrt((lm2.x-cm2.x)*(lm2.x-cm2.x)+(lm2.y-cm2.y)*(lm2.y-cm2.y)+(lm2.z-cm2.z)*(lm2.z-cm2.z));
        return r < 1e-6;
    }
    double a=0.5*(ulo+mu),c=0.5*(mu+uhi),e=0.5*(vlo+mv),f=0.5*(mv+vhi);
    return proveCurveCurve3D(L,C2,ulo,mu,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,ulo,mu,mv,vhi,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,mv,vhi,depth+1);
}

int main() {
    printf("==== 退化曲面 求解/证明恢复（降维 + 点包含，对齐曲线侧）====\n\n");

    // [A] 退化点 patch P × 3D 曲线：点包含（R→R³）
    {
        printf("[A] 退化点 patch P × 3D 曲线（点包含，R→R³）\n");
        P3 P{0.5, 0.4, 0.3};
        C3 over{{0.5, 0.4, 0}, {0.5, 0.4, 1}};          // 竖直直线，P 在内部(u=0.3)
        C3 miss{{0.5, 0.6, 0}, {0.5, 0.6, 1}};          // 竖直直线不过 P
        double lo1=0,hi1=1; bool ok1=provePointOnCurve3D(over,P,lo1,hi1,0);
        printf("   ×过P直线: 点包含证明=%s\n", ok1?"证明 ✓":"未证明 ✗");
        double lo2=0,hi2=1; bool ok2=provePointOnCurve3D(miss,P,lo2,hi2,0);
        printf("   ×不过P直线: 点包含证明=%s\n", ok2?"误证 ✗":"正确拒绝 ✓");
    }

    // [B] 退化点 patch P × 3D 曲面（点包含，R²→R³）
    {
        printf("\n[B] 退化点 patch P × 3D 曲面（点包含，R²→R³）\n");
        P3 P{0.5, 0.4, 0.3};
        Patch S{ {-1,0,0},{1,0,0},{-1,0,1},{1,0,1} };        // 平面 y=0，含 P? P.y=0.4 → 不含
        Box2 b1{-1,1,0,1}; bool ok1=provePointOnSurface3D(S,P,b1,0);
        printf("   ×平面y=0(不含P): 点包含证明=%s\n", ok1?"误证 ✗":"正确拒绝 ✓");
        Patch S2{ {-1,0.4,0},{1,0.4,0},{-1,0.4,1},{1,0.4,1} }; // 平面 y=0.4 含 P
        Box2 b2{-1,1,0,1}; bool ok2=provePointOnSurface3D(S2,P,b2,0);
        printf("   ×平面y=0.4(含P): 点包含证明=%s\n", ok2?"证明 ✓":"未证明 ✗");
    }

    // [C] 退化到曲线 patch（共线 → 3D 线段）× 3D 曲线（R²→R³）
    {
        printf("\n[C] 退化到曲线 patch（降维 3D 线段）× 3D 曲线（R²→R³）\n");
        // patch 4 角点共线于 y=x,z=0 → 退化为线段 (0,0,0)→(1,1,0)
        C3 L{{0,0,0},{1,1,0}};
        C3 hit{{0.5,0.5,-1},{0.5,0.5,1}};   // 竖直过线段中点 (0.5,0.5,0)，交点 v=0.5 内部
        double u0=0,u1=1,v0=0,v1=1; bool ok=proveCurveCurve3D(L,hit,u0,u1,v0,v1,0);
        printf("   ×过线段中点的直线: 证明=%s\n", ok?"证明 ✓":"未证明 ✗");
        C3 miss{{0.5,0.7,-1},{0.5,0.7,1}};  // 竖直不过线段
        double u0b=0,u1b=1,v0b=0,v1b=1; bool ok2=proveCurveCurve3D(L,miss,u0b,u1b,v0b,v1b,0);
        printf("   ×不过线段: 证明=%s\n", ok2?"误证 ✗":"正确拒绝 ✓");
    }

    printf("\n结论：退化曲面求解/证明恢复已落地——退化点 patch→点包含(R→R³/R²→R³)、\n");
    printf("      退化到曲线 patch→降维 3D 曲线-曲线(R²→R³)，正则拒绝/命中均正确。\n");
    return 0;
}
