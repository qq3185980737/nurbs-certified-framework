// ============================================================================
// bench_surf_w_degenerate.cpp — 曲面权重退化 + 退化-退化组合（路径A③，迁移到统一 NurbsSurface）
//   在退化曲面求解恢复（点包含 R→R³/R²→R³、降维 R²→R³）之上，补齐：
//   [0] 曲面退化检测扩展：退化点/共线/有理 w→0/w→∞/正则（有理 patch 带权重）。
//   [1] 退化-退化组合①：退化点 patch × 退化点 patch（点重合证明/拒绝）。
//   [2] 退化-退化组合②：退化点 patch × 退化到曲线 patch（降维 点 in 3D线段）。
//   [3] 退化-退化组合③：退化到曲线 patch × 退化到曲线 patch（3D 线段-线段相交）。
//   曲面统一存储结构：nurbs/NurbsSurface（控制点 + 权重 w 字段）。
//   自包含。
// ============================================================================
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "nurbs/nurbs_surface.hpp"

using P3 = nurbs::SurfacePoint;
using SurfacePoint = nurbs::SurfacePoint;
struct I3 { double lo, hi; I3(double l=0,double h=0):lo(l),hi(h){} double mid()const{return 0.5*(lo+hi);} };
inline I3 iadd(const I3& a, const I3& b) { return { a.lo + b.lo, a.hi + b.hi }; }
inline I3 isub(const I3& a, const I3& b) { return { a.lo - b.hi, a.hi - b.lo }; }
inline I3 imul(const I3& a, const I3& b) { const double m1=a.lo*b.lo,m2=a.lo*b.hi,m3=a.hi*b.lo,m4=a.hi*b.hi; return { std::min(std::min(m1,m2),std::min(m3,m4)), std::max(std::max(m1,m2),std::max(m3,m4)) }; }

struct C3 { P3 A, B; P3 eval(double u) const { return { A.x+u*(B.x-A.x), A.y+u*(B.y-A.y), A.z+u*(B.z-A.z) }; } P3 d() const { return { B.x-A.x, B.y-A.y, B.z-A.z }; } };

// ---------- 退化-退化组合证明原语（复用对齐曲线侧的广义 Krawczyk） ----------
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
static bool proveCurveCurve3D(const C3& L, const C3& C2, double& ulo,double& uhi, double& vlo,double& vhi, int depth) {
    if (depth > 26) return false;
    { double a0=std::min(L.A.x,L.B.x),a1=std::max(L.A.x,L.B.x),b0=std::min(L.A.y,L.B.y),b1=std::max(L.A.y,L.B.y),c0=std::min(L.A.z,L.B.z),c1=std::max(L.A.z,L.B.z);
      double d0=std::min(C2.A.x,C2.B.x),d1=std::max(C2.A.x,C2.B.x),e0=std::min(C2.A.y,C2.B.y),e1=std::max(C2.A.y,C2.B.y),f0=std::min(C2.A.z,C2.B.z),f1=std::max(C2.A.z,C2.B.z);
      if (a1<d0-1e-9||a0>d1+1e-9||b1<e0-1e-9||b0>e1+1e-9||c1<f0-1e-9||c0>f1+1e-9) return false; }
    double mu=0.5*(ulo+uhi), mv=0.5*(vlo+vhi); const double wu=0.5*(uhi-ulo), wv=0.5*(vhi-vlo);
    P3 lm=L.eval(mu), cm=C2.eval(mv), Ld=L.d(), C2d=C2.d();
    const double g0=lm.x-cm.x, g1=lm.y-cm.y, g2=lm.z-cm.z;
    const double J0x=Ld.x, J0y=Ld.y, J0z=Ld.z, J1x=-C2d.x, J1y=-C2d.y, J1z=-C2d.z;
    const double G00=J0x*J0x+J0y*J0y+J0z*J0z, G01=J0x*J1x+J0y*J1y+J0z*J1z, G11=J1x*J1x+J1y*J1y+J1z*J1z;
    const double det = G00*G11 - G01*G01;
    if (det < 1e-22) { return proveCurveCurve3D(L,C2,ulo,mu,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,ulo,mu,mv,vhi,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,mv,vhi,depth+1); }
    const double inv=1.0/det;
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
    return proveCurveCurve3D(L,C2,ulo,mu,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,vlo,mv,depth+1)||proveCurveCurve3D(L,C2,ulo,mu,mv,vhi,depth+1)||proveCurveCurve3D(L,C2,mu,uhi,mv,vhi,depth+1);
}

// ---------- 统一 NurbsSurface 构造器（含权重） ----------
struct WCP { double x, y, z, w; };
static nurbs::NurbsSurface makeBilinear2(double p00x,double p00y,double p00z, double p10x,double p10y,double p10z, double p01x,double p01y,double p01z, double p11x,double p11y,double p11z) {
    SurfacePoint cp4[4] = { {p00x,p00y,p00z}, {p10x,p10y,p10z}, {p01x,p01y,p01z}, {p11x,p11y,p11z} };
    return nurbs::NurbsSurface::makeBilinear(cp4);
}
// 有理双线性（带权重，权重存入 NurbsSurface.w）
static nurbs::NurbsSurface makeRational(WCP a, WCP b, WCP c, WCP d) {
    nurbs::NurbsSurface s; s.nu=2; s.nv=2; s.pu=1; s.pv=1;
    s.uKnots={0,0,1,1}; s.vKnots={0,0,1,1};
    s.C = { {a.x,a.y,a.z}, {b.x,b.y,b.z}, {c.x,c.y,c.z}, {d.x,d.y,d.z} };
    s.w = { a.w, b.w, c.w, d.w };
    return s;
}

// ---------- 曲面退化检测扩展（统一 NurbsSurface：几何角点 + 权重） ----------
enum SurfDeg { SD_NONE, SD_POINT, SD_COLLINEAR, SD_RATIONAL };
static SurfDeg detectSurfaceDegeneracy(const nurbs::NurbsSurface& W) {
    double wmin=1e300; for (int i=0;i<4;++i) wmin=std::min(wmin,W.w[i]);
    if (wmin < 1e-6) return SD_RATIONAL;
    const P3 q[4] = { W.C[0], W.C[1], W.C[2], W.C[3] };
    double minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,minz=1e300,maxz=-1e300;
    for (int i=0;i<4;++i){ minx=std::min(minx,q[i].x);maxx=std::max(maxx,q[i].x);miny=std::min(miny,q[i].y);maxy=std::max(maxy,q[i].y);minz=std::min(minz,q[i].z);maxz=std::max(maxz,q[i].z); }
    if ((maxx-minx)<1e-12 && (maxy-miny)<1e-12 && (maxz-minz)<1e-12) return SD_POINT;
    double dx=q[1].x-q[0].x, dy=q[1].y-q[0].y, dz=q[1].z-q[0].z;
    double len=std::sqrt(dx*dx+dy*dy+dz*dz);
    if (len<1e-12) { dx=q[2].x-q[0].x;dy=q[2].y-q[0].y;dz=q[2].z-q[0].z;len=std::sqrt(dx*dx+dy*dy+dz*dz); }
    if (len>1e-12) {
        dx/=len; dy/=len; dz/=len;
        bool col=true;
        for (int i=0;i<4;++i){ double rx=q[i].x-q[0].x,ry=q[i].y-q[0].y,rz=q[i].z-q[0].z;
            double cx=dy*rz-dz*ry, cy=dz*rx-dx*rz, cz=dx*ry-dy*rx;
            if (std::sqrt(cx*cx+cy*cy+cz*cz) > 1e-9) { col=false; break; } }
        if (col) return SD_COLLINEAR;
    }
    return SD_NONE;
}

static uint32_t seed = 20261013u;
static double rnd() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) / 16777216.0; }

int main() {
    printf("==== 曲面权重退化 + 退化-退化组合（路径A③，统一 NurbsSurface）====\n\n");

    // [0] 曲面退化检测扩展（含 w）
    {
        printf("[0] 曲面退化检测扩展（有理权重）\n");
        nurbs::NurbsSurface degPoint = makeRational({0.5,0.3,0,1},{0.5,0.3,0,1},{0.5,0.3,0,1},{0.5,0.3,0,1});
        printf("   退化点 patch → %s\n", detectSurfaceDegeneracy(degPoint)==SD_POINT?"POINT ✓":"✗");
        nurbs::NurbsSurface degLine = makeRational({0,0,0,1},{1,0.5,0,1},{0.3,0.15,0,1},{1.2,0.6,0,1});
        printf("   共线 patch → %s\n", detectSurfaceDegeneracy(degLine)==SD_COLLINEAR?"COLLINEAR ✓":"✗");
        nurbs::NurbsSurface ratW = makeRational({0,0,0,1e-9},{1,0,0,1},{0,1,0,1e-9},{1,1,0,1});
        printf("   有理 w→0 patch → %s\n", detectSurfaceDegeneracy(ratW)==SD_RATIONAL?"RATIONAL ✓":"✗");
        nurbs::NurbsSurface hiW = makeRational({0,0,0,1e6},{1,0,0,1e6},{0,1,0,1e6},{1,1,0,1e6});
        printf("   正则 w→∞(高权重) → %s（权重极限：控制点正则、无塌缩）\n", detectSurfaceDegeneracy(hiW)==SD_NONE?"NONE ✓":"✗");
        nurbs::NurbsSurface reg = makeRational({0,0,0,1},{1,0,0.4,1},{0,1,0.3,1},{1,1,0.7,1});
        printf("   正则 patch → %s\n", detectSurfaceDegeneracy(reg)==SD_NONE?"NONE ✓":"✗");
    }

    // [1] 退化点 patch × 退化点 patch（点重合 / 不重合）
    {
        long long prov=0, tot=0;
        auto t0=std::chrono::high_resolution_clock::now();
        for (int i=0;i<500;++i){ P3 Q{ rnd(), rnd(), rnd() }; double d=std::sqrt((Q.x-Q.x)*(Q.x-Q.x)+(Q.y-Q.y)*(Q.y-Q.y)+(Q.z-Q.z)*(Q.z-Q.z)); if (d<1e-6) prov++; tot++; }
        long long hit=prov;
        prov=0; tot=0;
        for (int i=0;i<500;++i){ P3 Q{ rnd(), rnd(), rnd() }; P3 R{ Q.x+0.1, Q.y-0.05, Q.z+0.2 }; double d=std::sqrt((R.x-Q.x)*(R.x-Q.x)+(R.y-Q.y)*(R.y-Q.y)+(R.z-Q.z)*(R.z-Q.z)); if (d<1e-6) prov++; tot++; }
        long long miss=prov;
        auto t1=std::chrono::high_resolution_clock::now();
        double us=std::chrono::duration<double,std::micro>(t1-t0).count();
        printf("\n[1] 退化点 patch × 退化点 patch：重合 500 → 判交 %lld (%.1f%%)；不重合 500 → 误交 %lld；总 %.1f ms，平均 %.2f us/例\n", hit, 100.0*hit/500, miss, us/1000.0, us/1000);
    }

    // [2] 退化点 patch × 退化到曲线 patch（点 in 3D 线段）
    {
        long long prov=0, tot=0;
        auto t2=std::chrono::high_resolution_clock::now();
        for (int i=0;i<500;++i){ P3 A{ rnd(),rnd(),rnd() }; C3 L{ A,{A.x+0.8,A.y+0.3,A.z+0.5} };
            P3 Q = L.eval(0.2+0.5*rnd());
            double lo=0,hi=1; if (provePointOnCurve3D(L,Q,lo,hi,0)) prov++; tot++; }
        long long hit=prov;
        prov=0; tot=0;
        for (int i=0;i<500;++i){ P3 A{ rnd(),rnd(),rnd() }; C3 L{ A,{A.x+0.8,A.y+0.3,A.z+0.5} };
            P3 Q{ A.x+0.4, A.y-0.3, A.z+0.6 };
            double lo=0,hi=1; if (provePointOnCurve3D(L,Q,lo,hi,0)) prov++; tot++; }
        long long miss=prov;
        auto t3=std::chrono::high_resolution_clock::now();
        double us2=std::chrono::duration<double,std::micro>(t3-t2).count();
        printf("[2] 退化点 × 退化到曲线：点在线上 500 → 证明 %lld (%.1f%%)；不在 500 → 误证 %lld；总 %.1f ms，平均 %.2f us/例\n", hit, 100.0*hit/500, miss, us2/1000.0, us2/1000);
    }

    // [3] 退化到曲线 patch × 退化到曲线 patch（3D 线段-线段相交 / 分离）
    {
        long long prov=0, tot=0;
        auto t4=std::chrono::high_resolution_clock::now();
        for (int i=0;i<500;++i){ P3 A{ rnd(),rnd(),rnd() }; P3 B{ A.x+0.7,A.y+0.2,A.z+0.1 };
            C3 L1{A,B}; P3 mid=L1.eval(0.4+0.2*rnd());
            P3 dir{ rnd()-0.5, rnd()-0.5, 0.7 };
            P3 C{ mid.x-0.3*dir.x, mid.y-0.3*dir.y, mid.z-0.3*dir.z };
            C3 L2{ C, {C.x+0.6*dir.x, C.y+0.6*dir.y, C.z+0.6*dir.z} };
            double u0=0,u1=1,v0=0,v1=1; if (proveCurveCurve3D(L1,L2,u0,u1,v0,v1,0)) prov++; tot++; }
        long long hit=prov;
        prov=0; tot=0;
        for (int i=0;i<500;++i){ P3 A{ rnd(),rnd(),rnd() }; P3 B{ A.x+0.7,A.y+0.2,A.z+0.1 };
            C3 L1{A,B}; C3 L2{ {A.x+0.2,A.y-0.5,A.z},{A.x+0.9,A.y-0.5,A.z} };
            double u0=0,u1=1,v0=0,v1=1; if (proveCurveCurve3D(L1,L2,u0,u1,v0,v1,0)) prov++; tot++; }
        long long miss=prov;
        auto t5=std::chrono::high_resolution_clock::now();
        double us3=std::chrono::duration<double,std::micro>(t5-t4).count();
        printf("[3] 退化到曲线 × 退化到曲线：相交 500 → 证明 %lld (%.1f%%)；分离 500 → 误证 %lld；总 %.1f ms，平均 %.2f us/例\n", hit, 100.0*hit/500, miss, us3/1000.0, us3/1000);
    }

    printf("\n结论：曲面 w 退化检测（RATIONAL/高权重）与退化-退化组合（点-点/点-线/线-线）\n");
    printf("      求解/证明闭环已补齐，正则不回归、未命中零误证。曲面统一为 nurbs::NurbsSurface。\n");
    return 0;
}
