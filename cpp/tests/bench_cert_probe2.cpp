// diagnostic: replicate kernel slice Krawczyk for one leaf box
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
    Box4 X;
    X.u = Interval(0.625, 0.75); X.v = Interval(0.0, 0.125);
    X.s = Interval(0.0, 0.125); X.t = Interval(0.375, 0.5);
    Interval S1[3], S1u[3], S1v[3], S2[3], S2u[3], S2v[3];
    int ok1 = eval_rational_duv(A, X.u, X.v, S1, S1u, S1v);
    int ok2 = eval_rational_duv(B, X.s, X.t, S2, S2u, S2v);
    std::printf("eval A ok=%d S1=(%.17g,%.17g,%.17g) S1u=(%.17g..%.17g,%.17g..%.17g,%.17g..%.17g)\n",
        ok1, S1[0].mid(), S1[1].mid(), S1[2].mid(),
        S1u[0].lo, S1u[0].hi, S1u[1].lo, S1u[1].hi, S1u[2].lo, S1u[2].hi);
    std::printf("eval B ok=%d S2=(%.17g,%.17g,%.17g) S2u=(%.17g..%.17g,%.17g..%.17g,%.17g..%.17g)\n",
        ok2, S2[0].mid(), S2[1].mid(), S2[2].mid(),
        S2u[0].lo, S2u[0].hi, S2u[1].lo, S2u[1].hi, S2u[2].lo, S2u[2].hi);
    // J3 for kfix=v (idx u,s,t)
    double J3m[3][3] = {
        {S1u[0].mid(), -S2u[0].mid(), -S2v[0].mid()},
        {S1u[1].mid(), -S2u[1].mid(), -S2v[1].mid()},
        {S1u[2].mid(), -S2u[2].mid(), -S2v[2].mid()}};
    double inv[3][3]; bool invok = invert3(J3m, inv);
    std::printf("invert3 ok=%d inv diag=(%.6g,%.6g,%.6g)\n", (int)invok, inv[0][0], inv[1][1], inv[2][2]);
    // K
    double m3[3] = {X.u.mid(), X.s.mid(), X.t.mid()};
    double Fm[3] = {S1[0].mid()-S2[0].mid(), S1[1].mid()-S2[1].mid(), S1[2].mid()-S2[2].mid()};
    for (int r=0;r<3;++r) {
        double k = m3[r];
        for (int c=0;c<3;++c) k -= inv[r][c]*Fm[c];
        std::printf("K[%d]=%.17g  X3=[%.17g,%.17g] inside=%d\n", r, k,
            (r==0?X.u.lo:(r==1?X.s.lo:X.t.lo)), (r==0?X.u.hi:(r==1?X.s.hi:X.t.hi)),
            (int)(k >= (r==0?X.u.lo:(r==1?X.s.lo:X.t.lo)) && k <= (r==0?X.u.hi:(r==1?X.s.hi:X.t.hi))));
    }
    return 0;
}
