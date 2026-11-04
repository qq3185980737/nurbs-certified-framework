// debug_pair.cpp — 复现 probe_pairs 中某一个 (fidA,piA,fidB,piB) 的 bbox 与 classify 结果
#include <cstdio>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include "nurbs/cert_ssi.hpp"
using nurbs::cert::RBezierSurface;
using nurbs::Interval;
struct AABB { double xlo=1e300, ylo=1e300, zlo=1e300, xhi=-1e300, yhi=-1e300, zhi=-1e300; };
static AABB patch_bb(const RBezierSurface& s) {
    AABB b;
    double wmin = s.wmin > 0 ? s.wmin : 1e-9;
    double wmax = s.wmax > 0 ? s.wmax : 1.0;
    for (const auto& p : s.P) {
        double xl = (p[0].lo >= 0) ? p[0].lo / wmax : p[0].lo / wmin;
        double xh = (p[0].hi <= 0) ? p[0].hi / wmax : p[0].hi / wmin;
        double yl = (p[1].lo >= 0) ? p[1].lo / wmax : p[1].lo / wmin;
        double yh = (p[1].hi <= 0) ? p[1].hi / wmax : p[1].hi / wmin;
        double zl = (p[2].lo >= 0) ? p[2].lo / wmax : p[2].lo / wmin;
        double zh = (p[2].hi <= 0) ? p[2].hi / wmax : p[2].hi / wmin;
        b.xlo = std::min(b.xlo, xl); b.xhi = std::max(b.xhi, xh);
        b.ylo = std::min(b.ylo, yl); b.yhi = std::max(b.yhi, yh);
        b.zlo = std::min(b.zlo, zl); b.zhi = std::max(b.zhi, zh);
    }
    return b;
}
static RBezierSurface make_patch(int nu, int nv, const std::vector<std::array<double,4>>& cw) {
    RBezierSurface s; s.nu = nu; s.nv = nv; s.P.resize((std::size_t)nu*nv);
    for (std::size_t k = 0; k < s.P.size(); ++k) {
        double w = cw[k][3];
        s.P[k][0] = Interval::point(cw[k][0]*w); s.P[k][1] = Interval::point(cw[k][1]*w);
        s.P[k][2] = Interval::point(cw[k][2]*w); s.P[k][3] = Interval::point(w);
    }
    s.update_weight_range(); return s;
}
int main(int argc, char** argv) {
    // usage: data.txt fidA piA fidB piB
    std::ifstream f(argv[1]);
    std::string tok; std::vector<std::vector<RBezierSurface>> faces;
    while (f >> tok) {
        if (tok == "NSURF") { int n; f >> n; }
        else if (tok == "SURF") {
            int fid, pi, du, dv, nu, nv; f >> fid >> pi >> du >> dv >> nu >> nv;
            std::vector<std::array<double,4>> pts((std::size_t)nu*nv);
            for (int i=0;i<nu;++i) for(int j=0;j<nv;++j){ double x,y,z,w; f>>x>>y>>z>>w; pts[(std::size_t)i*nv+j]={x,y,z,w}; }
            std::string e; f >> e;
            if (fid >= (int)faces.size()) faces.resize(fid+1);
            faces[fid].push_back(make_patch(nu,nv,pts));
        } else { std::string l; std::getline(f,l); }
    }
    int a=atoi(argv[2]), ia=atoi(argv[3]), b=atoi(argv[4]), ib=atoi(argv[5]);
    AABB ba=patch_bb(faces[a][ia]), bb=patch_bb(faces[b][ib]);
    auto ov=[&](const AABB& x,const AABB& y){ return x.xhi>=y.xlo&&y.xhi>=x.xlo&&x.yhi>=y.ylo&&y.yhi>=x.ylo&&x.zhi>=y.zlo&&y.zhi>=x.zlo; };
    std::printf("A[%d][%d] w=[%.9g,%.9g] bb=[%.6f,%.6f]x[%.6f,%.6f]x[%.6f,%.6f]\n", a,ia,faces[a][ia].wmin,faces[a][ia].wmax,ba.xlo,ba.xhi,ba.ylo,ba.yhi,ba.zlo,ba.zhi);
    std::printf("B[%d][%d] w=[%.9g,%.9g] bb=[%.6f,%.6f]x[%.6f,%.6f]x[%.6f,%.6f]\n", b,ib,faces[b][ib].wmin,faces[b][ib].wmax,bb.xlo,bb.xhi,bb.ylo,bb.yhi,bb.zlo,bb.zhi);
    std::printf("aabb_overlap=%d\n", (int)ov(ba,bb));
    auto pc = nurbs::cert::classify_pair(faces[a][ia], faces[b][ib], 0.02, 4000000L);
    std::printf("classify type=%d cover=%ld certified=%ld undecided=%ld comps=%ld cop=%d\n", pc.type, pc.cover, pc.certified, pc.undecided, pc.components, pc.coplanar_overlap);
    return 0;
}
