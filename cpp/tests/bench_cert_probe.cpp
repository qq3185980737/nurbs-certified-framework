// probe: test kernel slice Krawczyk on the plane-plane demo pair
#include <cstdio>
#include <array>
#include <vector>
#include "nurbs/cert_ssi.hpp"

using namespace nurbs::cert;
using nurbs::Interval;

static RBezierSurface plane_z0() {
    std::vector<std::array<double,4>> P(4);
    P[0]={-1,-1,0,1}; P[1]={-1,1,0,1}; P[2]={1,-1,0,1}; P[3]={1,1,0,1};
    return make_surface(2,2,P);
}
static RBezierSurface plane_x(double x0) {
    std::vector<std::array<double,4>> P(4);
    P[0]={x0,-1,-1,1}; P[1]={x0,-1,1,1}; P[2]={x0,1,-1,1}; P[3]={x0,1,1,1};
    return make_surface(2,2,P);
}
int main() {
    auto A = plane_z0();
    auto B = plane_x(0.5);
    // 交线：u=0.75, t=0.5, v=s
    for (double hv : {0.01, 0.05, 0.1, 0.2}) {
        int code = probe_box(A, B, 0.75, 0.4, 0.4, 0.5, hv);
        std::printf("probe hv=%g at (0.75,0.4,0.4,0.5): code=%d\n", hv, code);
    }
    // 一个覆盖盒：u[0.65,0.85] v[0.2,0.4] s[0.2,0.4] t[0.4,0.6]
    Box4 box;
    box.u = Interval(0.65,0.85); box.v = Interval(0.2,0.4);
    box.s = Interval(0.2,0.4); box.t = Interval(0.4,0.6);
    std::printf("krawczyk_slice on cover box: %d\n", (int)krawczyk_slice(A,B,box));
    Box4 lb;
    lb.u = Interval(0.625, 0.75); lb.v = Interval(0.0, 0.125);
    lb.s = Interval(0.0, 0.125); lb.t = Interval(0.375, 0.5);
    std::printf("leaf box krawczyk_slice: %d\n", (int)krawczyk_slice(A,B,lb));
    for (int k=0;k<4;++k) std::printf("  leaf force_k=%d code=%d\n", k, krawczyk_slice_code(A,B,lb,k));
    for (int k=0;k<4;++k) std::printf("  force_k=%d code=%d\n", k, krawczyk_slice_code(A,B,box,k));
    return 0;
}
