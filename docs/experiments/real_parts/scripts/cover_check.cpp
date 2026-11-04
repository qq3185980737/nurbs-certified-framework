// ============================================================================
// cover_check.cpp — 覆盖检验（需求 4）：OCCT 交线采样点是否落在认证 cover 盒的
//   三维像包围盒内。
//   原理：cover(A,B,tau) 返回参数域 4D 盒 (u,v,s,t) 的集合，覆盖全部真交点参数。
//   若真交点 p = A(u,v) = B(s,t)，则 p 必属于某个盒在 A（及 B）上的区间求值像。
//   故"p 不在任何 A-盒像内" ⇒ cover 漏掉该交点（不健全）。本工具用区间 de Casteljau
//   （nurbs::cert::eval_net）对每个盒求保守三维像 bbox，检查 OCCT 采样点的包含性。
//   注意：这是**必要条件**检验（保守），通过则支持覆盖正确性，失败则说明漏覆盖。
// 用法: cover_check <data.txt> <pairs.csv> <points.csv> <out.csv> [tau] [node_cap]
// ============================================================================
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::Interval;

struct AABB3 { double xlo=1e300, ylo=1e300, zlo=1e300, xhi=-1e300, yhi=-1e300, zhi=-1e300; };

static RBezierSurface make_patch(int nu, int nv, const std::vector<std::array<double,4>>& cart_w) {
    RBezierSurface s; s.nu = nu; s.nv = nv;
    s.P.resize((std::size_t)nu * nv);
    for (std::size_t k = 0; k < s.P.size(); ++k) {
        double w = cart_w[k][3];
        s.P[k][0] = Interval::point(cart_w[k][0] * w);
        s.P[k][1] = Interval::point(cart_w[k][1] * w);
        s.P[k][2] = Interval::point(cart_w[k][2] * w);
        s.P[k][3] = Interval::point(w);
    }
    s.update_weight_range();
    return s;
}

// 用区间 de Casteljau 对给定参数盒求三维保守像 bbox
static AABB3 image_bb(const RBezierSurface& S, const nurbs::cert::Box4& X, bool useA) {
    nurbs::cert::HP out;
    nurbs::cert::eval_net(S.P, S.nu, S.nv, useA ? X.u : X.s, useA ? X.v : X.t, out);
    AABB3 b;
    double wlo = out[3].lo, whi = out[3].hi;
    if (whi <= 0) return b;              // 退化：忽略
    if (wlo < 1e-12) wlo = 1e-12;
    // 保守像区间：x = wx/w，w∈[wlo,whi]>0。
    //   wx>=0 → [wx/whi, wx/wlo]；wx<0 → [wx/wlo, wx/whi]
    // （此前直接写 lo=wx.lo/whi、hi=wx.hi/wlo 对负坐标得到反向区间，导致覆盖检验假 miss）
    double lo[3], hi[3];
    for (int c = 0; c < 3; ++c) {
        lo[c] = (out[c].lo >= 0) ? out[c].lo / whi : out[c].lo / wlo;
        hi[c] = (out[c].hi <= 0) ? out[c].hi / whi : out[c].hi / wlo;
    }
    b.xlo = lo[0]; b.xhi = hi[0];
    b.ylo = lo[1]; b.yhi = hi[1];
    b.zlo = lo[2]; b.zhi = hi[2];
    return b;
}

static bool inside(const AABB3& b, double x, double y, double z, double tol) {
    return x >= b.xlo - tol && x <= b.xhi + tol &&
           y >= b.ylo - tol && y <= b.yhi + tol &&
           z >= b.zlo - tol && z <= b.zhi + tol;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::printf("usage: %s <data.txt> <pairs.csv> <points.csv> <out.csv> [tau] [node_cap]\n", argv[0]);
        return 1;
    }
    const char* datap = argv[1];
    const char* pairsp = argv[2];
    const char* pointsp = argv[3];
    const char* outp = argv[4];
    double tau = (argc >= 6) ? std::atof(argv[5]) : 0.02;
    long node_cap = (argc >= 7) ? std::atol(argv[6]) : 4000000L;

    std::ifstream f(datap);
    if (!f) { std::printf("cannot open %s\n", datap); return 1; }
    std::map<std::pair<int,int>, RBezierSurface> patches;
    std::string tok;
    while (f >> tok) {
        if (tok == "NSURF") { int n; f >> n; }
        else if (tok == "SURF") {
            int fid, pi, du, dv, nu, nv;
            f >> fid >> pi >> du >> dv >> nu >> nv;
            std::vector<std::array<double,4>> pts((std::size_t)nu * nv);
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j) {
                    double x,y,z,w; f >> x >> y >> z >> w;
                    pts[(std::size_t)i*nv+j] = {x,y,z,w};
                }
            std::string e; f >> e;
            patches[std::make_pair(fid, pi)] = make_patch(nu, nv, pts);
        } else { std::string line; std::getline(f, line); }
    }

    // 读 pairs.csv 的表头，定位列
    std::ifstream pf(pairsp);
    if (!pf) { std::printf("cannot open %s\n", pairsp); return 1; }
    std::string line;
    std::getline(pf, line);
    std::vector<std::array<int,4>> pairList;
    while (std::getline(pf, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell; std::vector<std::string> cells;
        while (std::getline(ss, cell, ',')) cells.push_back(cell);
        if (cells.size() < 5) continue;
        pairList.push_back({std::atoi(cells[1].c_str()), std::atoi(cells[2].c_str()),
                            std::atoi(cells[3].c_str()), std::atoi(cells[4].c_str())});
    }
    pf.close();

    // 读 points.csv：part,fidA,piA,fidB,piB,line,k,x,y,z
    std::ifstream qf(pointsp);
    if (!qf) { std::printf("cannot open %s\n", pointsp); return 1; }
    std::getline(qf, line);
    std::map<std::array<int,4>, std::vector<std::array<double,3>>> pts;
    while (std::getline(qf, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell; std::vector<std::string> cells;
        while (std::getline(ss, cell, ',')) cells.push_back(cell);
        if (cells.size() < 10) continue;
        std::array<int,4> key = {std::atoi(cells[1].c_str()), std::atoi(cells[2].c_str()),
                                 std::atoi(cells[3].c_str()), std::atoi(cells[4].c_str())};
        pts[key].push_back({std::atof(cells[7].c_str()), std::atof(cells[8].c_str()),
                            std::atof(cells[9].c_str())});
    }
    qf.close();

    FILE* out = std::fopen(outp, "w");
    std::fprintf(out, "fidA,piA,fidB,piB,ncover,npoints,inside_A,inside_B,inside_both,min_miss_Linf\n");
    long nPairsDone = 0, nPtsTotal = 0, nInsideBoth = 0, nPairsAllInside = 0, nPairsMiss = 0;
    for (auto& key : pairList) {
        auto ia = patches.find({key[0], key[1]});
        auto ib = patches.find({key[2], key[3]});
        if (ia == patches.end() || ib == patches.end()) continue;
        auto pit = pts.find(key);
        if (pit == pts.end() || pit->second.empty()) continue;
        auto r = nurbs::cert::cover(ia->second, ib->second, tau, 60, node_cap);
        std::vector<AABB3> imA, imB;
        imA.reserve(r.boxes.size()); imB.reserve(r.boxes.size());
        for (const auto& X : r.boxes) { imA.push_back(image_bb(ia->second, X, true)); imB.push_back(image_bb(ib->second, X, false)); }
        long inA = 0, inB = 0, inBoth = 0;
        double worst = 0.0;
        for (const auto& p : pit->second) {
            bool fa = false, fb = false;
            double bestMiss = 1e300;
            for (std::size_t k = 0; k < r.boxes.size(); ++k) {
                if (inside(imA[k], p[0], p[1], p[2], 1e-7)) fa = true;
                if (inside(imB[k], p[0], p[1], p[2], 1e-7)) fb = true;
                // 到 A 盒像的 L∞ 距离（仅用于诊断）
                double dx = std::max({imA[k].xlo - p[0], 0.0, p[0] - imA[k].xhi});
                double dy = std::max({imA[k].ylo - p[1], 0.0, p[1] - imA[k].yhi});
                double dz = std::max({imA[k].zlo - p[2], 0.0, p[2] - imA[k].zhi});
                bestMiss = std::min(bestMiss, std::max({dx, dy, dz}));
            }
            if (fa) ++inA;
            if (fb) ++inB;
            if (fa && fb) ++inBoth;
            worst = std::max(worst, bestMiss);
            ++nPtsTotal;
        }
        nInsideBoth += inBoth;
        if (inBoth == (long)pit->second.size()) ++nPairsAllInside;
        else ++nPairsMiss;
        ++nPairsDone;
        std::fprintf(out, "%d,%d,%d,%d,%ld,%zu,%ld,%ld,%ld,%.6g\n",
                     key[0], key[1], key[2], key[3], (long)r.boxes.size(), pit->second.size(),
                     inA, inB, inBoth, worst);
        std::fflush(out);
    }
    std::fclose(out);
    std::printf("COVERCHECK pairs=%ld points=%ld inside_both=%ld pairs_all_inside=%ld pairs_with_miss=%ld "
                "tau=%.4f cap=%ld\n", nPairsDone, nPtsTotal, nInsideBoth, nPairsAllInside, nPairsMiss, tau, node_cap);
    return 0;
}
