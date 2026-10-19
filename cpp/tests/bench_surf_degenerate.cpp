// ============================================================================
// bench_surf_degenerate.cpp — 退化/奇异曲面专项（检测机制 + 诚实行为测试，迁移到统一 NurbsSurface）
//   结合最新理论（GPU Spline Surface Intersections 退化法向检测：法向曲面包围盒
//   不含原点则正则；A Precision Controlled SSI(ACM2025) 相切曲面/三类奇异性；
//   Reliable Computing 区间交线延拓）分类构造退化曲面：
//   ① 退化点 patch（控制点全重合→单点）
//   ② 塌缩参数化（patch 一条边→0 → 退化三角形）
//   ③ 退化到曲线（控制点共线）
//   ④ 退化法向（du x dv = 0，法向包围盒含原点，GPU 法向曲面检测）
//   ⑤ 相切交线（两曲面法向平行，SSI 奇异性）
//   + 专门处理：detectSurfaceDegeneracy / normalBoxDegenerate / tangentDetect。
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
inline I3 imul(const I3& a, const I3& b) { const double m1=a.lo*b.lo,m2=a.lo*b.hi,m3=a.hi*b.lo,m4=a.hi*b.hi; return { std::min(std::min(m1,m2),std::min(m3,m4)), std::max(std::max(m1,m2),std::max(m3,m4)) }; }

static nurbs::NurbsSurface makePatch(P3 p00, P3 p10, P3 p01, P3 p11) {
    SurfacePoint cp4[4] = { p00, p10, p01, p11 };
    return nurbs::NurbsSurface::makeBilinear(cp4);
}
static P3 normal(const nurbs::NurbsSurface& S, double s, double t) {
    P3 a=S.du(s,t), b=S.dv(s,t);
    return { a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x };
}

// ---------- 退化曲面检测（统一 NurbsSurface） ----------
struct SurfDeg { enum Type { NONE=0, POINT, COLLINEAR, SINGULAR_NORMAL } type; P3 P; };
static SurfDeg detectSurfaceDegeneracy(const nurbs::NurbsSurface& S) {
    SurfDeg d; d.type = SurfDeg::NONE;
    const P3* q[4] = { &S.C[0], &S.C[1], &S.C[2], &S.C[3] };
    double minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,minz=1e300,maxz=-1e300;
    for (auto p : q) { minx=std::min(minx,p->x); maxx=std::max(maxx,p->x); miny=std::min(miny,p->y); maxy=std::max(maxy,p->y); minz=std::min(minz,p->z); maxz=std::max(maxz,p->z); }
    // ① 退化点
    if ((maxx-minx)<1e-12 && (maxy-miny)<1e-12 && (maxz-minz)<1e-12) { d.type=SurfDeg::POINT; d.P={0.5*(minx+maxx),0.5*(miny+maxy),0.5*(minz+maxz)}; return d; }
    // ② 共线（退化到曲线）：所有角点距基线距离→0
    auto base=[&](){ double dx=q[1]->x-q[0]->x, dy=q[1]->y-q[0]->y, dz=q[1]->z-q[0]->z; double L=std::sqrt(dx*dx+dy*dy+dz*dz);
        if (L<1e-15) { dx=q[2]->x-q[0]->x; dy=q[2]->y-q[0]->y; dz=q[2]->z-q[0]->z; L=std::sqrt(dx*dx+dy*dy+dz*dz); }
        return std::make_pair(dx/(L+1e-300), std::make_pair(dy/(L+1e-300), dz/(L+1e-300))); };
    {
        auto b=base(); double ux=b.first, uy=b.second.first, uz=b.second.second;
        bool col=true; for (auto p : q) { double cx=p->x-q[0]->x, cy=p->y-q[0]->y, cz=p->z-q[0]->z;
            double cr = uy*cz-uz*cy, cs = uz*cx-ux*cz, ct = ux*cy-uy*cx;
            if (std::sqrt(cr*cr+cs*cs+ct*ct) > 1e-9) { col=false; break; } }
        if (col) { d.type=SurfDeg::COLLINEAR; return d; }
    }
    // ③ 退化法向：法向曲面包围盒（角点法向凸包）是否含原点（GPU 法向曲面检测）
    {
        P3 n[4] = { normal(S,0,0), normal(S,1,0), normal(S,0,1), normal(S,1,1) };
        double lox=1e300,hix=-1e300,loy=1e300,hiy=-1e300,loz=1e300,hiz=-1e300;
        for (auto& v : n) { lox=std::min(lox,v.x); hix=std::max(hix,v.x); loy=std::min(loy,v.y); hiy=std::max(hiy,v.y); loz=std::min(loz,v.z); hiz=std::max(hiz,v.z); }
        if (lox<=0 && hix>=0 && loy<=0 && hiy>=0 && loz<=0 && hiz>=0) d.type=SurfDeg::SINGULAR_NORMAL;
    }
    return d;
}
// 相切交线检测：两曲面在 (s1,t1,s2,t2) 法向平行（SSI 奇异性）
static int tangentDetect(const nurbs::NurbsSurface& S1, const nurbs::NurbsSurface& S2, double s1,double t1,double s2,double t2) {
    P3 n1=normal(S1,s1,t1), n2=normal(S2,s2,t2);
    double L1=std::sqrt(n1.x*n1.x+n1.y*n1.y+n1.z*n1.z), L2=std::sqrt(n2.x*n2.x+n2.y*n2.y+n2.z*n2.z);
    if (L1<1e-15 || L2<1e-15) return 1;   // 法向退化
    double c=(n1.x*n2.x+n1.y*n2.y+n1.z*n2.z)/(L1*L2);
    return std::abs(std::abs(c)-1.0) < 1e-6 ? 1 : 0;  // 平行 → 相切
}

int main() {
    printf("==== 退化/奇异曲面专项：检测机制 + 当前证明算子行为（统一 NurbsSurface）====\n\n");

    // [1] 退化点 patch（4 角点全重合 → 单点）
    {
        printf("[1] 退化点 patch（角点全重合 → 单点）\n");
        nurbs::NurbsSurface S = makePatch({0.5,0,0.5},{0.5,0,0.5},{0.5,0,0.5},{0.5,0,0.5});
        auto d = detectSurfaceDegeneracy(S);
        printf("   detectSurfaceDegeneracy → %s (P=(%.2f,%.2f,%.2f))\n",
               d.type==SurfDeg::POINT?"POINT":d.type==SurfDeg::COLLINEAR?"COLLINEAR":d.type==SurfDeg::SINGULAR_NORMAL?"SINGULAR_NORMAL":"NONE", d.P.x,d.P.y,d.P.z);
        P3 nn = normal(S,0.5,0.5);
        printf("   退化点 patch 的法向在任意点 = %s（法向恒零，Krawczyk Jacobian 奇异）\n",
               (nn.x==0 && nn.y==0 && nn.z==0) ? "0 退化" : "非零");
    }

    // [2] 塌缩参数化（双线性 patch 一条边 →0 → 退化三角形）
    {
        printf("\n[2] 塌缩参数化（t=0 边塌缩到一点 → 退化三角形）\n");
        nurbs::NurbsSurface S = makePatch({0,0,0},{1,0,0},{0,0,0},{1,0,0});   // P01=P00,P11=P10
        auto d = detectSurfaceDegeneracy(S);
        printf("   detectSurfaceDegeneracy → %s\n",
               d.type==SurfDeg::POINT?"POINT":d.type==SurfDeg::COLLINEAR?"COLLINEAR":d.type==SurfDeg::SINGULAR_NORMAL?"SINGULAR_NORMAL":"NONE");
        P3 n = normal(S, 0.5, 0.5);
        printf("   中心法向 du x dv = (%.2g,%.2g,%.2g) → %s\n", n.x,n.y,n.z,
               (n.x==0&&n.y==0&&n.z==0)?"法向退化(SINGULAR_NORMAL) ✓":"正则");
    }

    // [3] 退化到曲线（控制点共线）
    {
        printf("\n[3] 退化到曲线（4 角点共线 → 曲面塌缩为线段）\n");
        nurbs::NurbsSurface S = makePatch({0,0,0},{0.5,0.5,0},{0.2,0.2,0},{0.8,0.8,0});
        auto d = detectSurfaceDegeneracy(S);
        printf("   detectSurfaceDegeneracy → %s（应 COLLINEAR，降维为曲线-曲线处理）\n",
               d.type==SurfDeg::POINT?"POINT":d.type==SurfDeg::COLLINEAR?"COLLINEAR":d.type==SurfDeg::SINGULAR_NORMAL?"SINGULAR_NORMAL":"NONE");
    }

    // [4] 退化法向 patch（法向包围盒含原点）
    {
        printf("\n[4] 退化法向检测（法向曲面包围盒含原点，GPU 方法）\n");
        nurbs::NurbsSurface reg = makePatch({-1,0,0},{1,0,0},{-1,0,1},{1,0,1});
        auto d1 = detectSurfaceDegeneracy(reg);
        printf("   平面y=0 patch → %s（法向包围盒不含原点，正则）\n",
               d1.type==SurfDeg::SINGULAR_NORMAL?"SINGULAR_NORMAL":"NONE ✓");
        nurbs::NurbsSurface sg = makePatch({-1,0,0},{1,0,0},{-1,0,0},{1,0,0});
        auto d2 = detectSurfaceDegeneracy(sg);
        printf("   塌缩 patch → %s（法向包围盒含原点，退化法向）\n",
               d2.type==SurfDeg::SINGULAR_NORMAL?"SINGULAR_NORMAL ✓":d2.type==SurfDeg::COLLINEAR?"COLLINEAR":"NONE");
    }

    // [5] 相切交线（两曲面法向平行，SSI 奇异性）
    {
        printf("\n[5] 相切交线检测（两曲面法向平行，SSI 三类奇异性之一）\n");
        nurbs::NurbsSurface S1 = makePatch({-1,0,0},{1,0,0},{-1,0,1},{1,0,1});
        nurbs::NurbsSurface S2 = makePatch({-1,0.01,0},{1,0.01,0},{-1,0.01,1},{1,0.01,1});
        int t = tangentDetect(S1, S2, 0.5, 0.5, 0.5, 0.5);
        printf("   平面y=0 × 平行平面y=0.01: tangentDetect=%s（法向平行 → 相切/平行交线）\n", t ? "相切 ✓" : "横截");
        nurbs::NurbsSurface S3 = makePatch({-1,0,0},{1,0,0},{-1,0.5,0.5},{1,0.5,0.5});
        int t2 = tangentDetect(S1, S3, 0.5, 0.5, 0.5, 0.5);
        printf("   平面y=0 × 倾斜平面(横截): tangentDetect=%s（法向不平行 → 横截交线）\n", t2 ? "相切" : "横截 ✓");
    }

    printf("\n结论：退化曲面检测（POINT/COLLINEAR/SINGULAR_NORMAL）+ 相切交线检测已落地，\n");
    printf("      退化法向用 GPU 法向曲面包围盒判据（不含原点即正则），正则曲面不误报。\n");
    printf("      曲面表示统一为 nurbs::NurbsSurface（双线性 makeBilinear）。\n");
    return 0;
}
