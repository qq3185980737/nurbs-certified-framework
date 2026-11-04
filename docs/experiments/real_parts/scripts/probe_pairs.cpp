// ============================================================================
// probe_pairs.cpp — 真实零件 patch-pair 抽样探针
//   复用 bench_real_batch.cpp 的解析 / bbox / classify_pair 逻辑，但：
//     * 统计精确的 face-pair / patch-pair / bbox 预筛数量（count 模式，总是执行）
//     * 从"face bbox 相交且至少一对 patch bbox 相交"的 face-pair 中抽样 npos 对，
//       对每个 face-pair 取第一对 bbox 相交 patch-pair 调 classify_pair，逐对计时
//     * 另抽 nneg 对 bbox **不相交** 的 patch-pair（同 face-pair 内）作阴性对照
//   输出 CSV（每行一对），供 occt_check.py 用 GeomAPI_IntSS 做交叉校验。
// 用法: probe_pairs <data.txt> <out.csv> [tau] [node_cap] [npos] [nneg] [seed]
// ============================================================================
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <random>
#include <algorithm>
#include <cstdlib>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::Interval;

struct AABB { double xlo=1e300, ylo=1e300, zlo=1e300, xhi=-1e300, yhi=-1e300, zhi=-1e300; };
static AABB patch_bb(const RBezierSurface& s) {
    AABB b;
    double wmin = s.wmin > 0 ? s.wmin : 1e-9;
    double wmax = s.wmax > 0 ? s.wmax : 1.0;
    for (const auto& p : s.P) {
        // 正确保守 AABB：w>0，x=wx/w。
        //   wx>=0: x∈[wx/wmax, wx/wmin]；wx<0: x∈[wx/wmin, wx/wmax]
        // （原 bench 版直接写 xl=lo/wmax, xh=hi/wmin，负坐标会得到反向/空区间 —— 见 summary.md）
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
static bool aabb_overlap(const AABB& a, const AABB& b) {
    return a.xhi >= b.xlo && b.xhi >= a.xlo && a.yhi >= b.ylo && b.yhi >= a.ylo && a.zhi >= b.zlo && b.zhi >= a.zlo;
}
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

struct Cand { int a, b, pa, pb, xa, xb; };  // x*= 第一对 bbox 相交 / 不相交 patch-pair

int main(int argc, char** argv) {
    if (argc < 3) { std::printf("usage: %s <data.txt> <out.csv> [tau] [node_cap] [npos] [nneg] [seed]\n", argv[0]); return 1; }
    const char* path = argv[1];
    const char* outcsv = argv[2];
    double tau = (argc >= 4) ? std::atof(argv[3]) : 0.02;
    long node_cap = (argc >= 5) ? std::atol(argv[4]) : 4000000L;
    int npos = (argc >= 6) ? std::atoi(argv[5]) : 50;
    int nneg = (argc >= 7) ? std::atoi(argv[6]) : 50;
    unsigned long long seed = (argc >= 8) ? std::strtoull(argv[7], nullptr, 10) : 12345ULL;

    std::ifstream f(path);
    if (!f) { std::printf("cannot open %s\n", path); return 1; }
    std::vector<std::vector<RBezierSurface>> faces;
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
            if (fid >= (int)faces.size()) faces.resize(fid + 1);
            faces[fid].push_back(make_patch(nu, nv, pts));
        } else { std::string line; std::getline(f, line); }
    }
    int nf = (int)faces.size();
    // patch 索引 = 数据文件中的 (fid, pi) 顺序（bench 解析顺序）
    std::vector<std::vector<AABB>> pbb(nf);
    long totalPatches = 0;
    std::vector<AABB> fbb(nf);
    for (int a = 0; a < nf; ++a) {
        pbb[a].resize(faces[a].size());
        AABB u;
        for (std::size_t k = 0; k < faces[a].size(); ++k) {
            pbb[a][k] = patch_bb(faces[a][k]);
            const AABB& ba = pbb[a][k];
            u.xlo = std::min(u.xlo, ba.xlo); u.xhi = std::max(u.xhi, ba.xhi);
            u.ylo = std::min(u.ylo, ba.ylo); u.yhi = std::max(u.yhi, ba.yhi);
            u.zlo = std::min(u.zlo, ba.zlo); u.zhi = std::max(u.zhi, ba.zhi);
        }
        fbb[a] = u;
        totalPatches += (long)faces[a].size();
    }

    long facePairs = 0, overlapFacePairs = 0, patchPairs = 0, overlapPatchPairs = 0;
    std::vector<Cand> cands;
    std::mt19937_64 rng(seed);
    for (int a = 0; a < nf; ++a) {
        for (int b = a + 1; b < nf; ++b) {
            ++facePairs;
            if (!aabb_overlap(fbb[a], fbb[b])) continue;
            ++overlapFacePairs;
            Cand c{-1,-1,-1,-1,-1,-1};
            c.a = a; c.b = b;
            for (std::size_t ia = 0; ia < faces[a].size(); ++ia) {
                for (std::size_t ib = 0; ib < faces[b].size(); ++ib) {
                    ++patchPairs;
                    bool ov = aabb_overlap(pbb[a][ia], pbb[b][ib]);
                    if (ov) { ++overlapPatchPairs; if (c.xa < 0) { c.xa = (int)ia; c.xb = (int)ib; } }
                    else      { if (c.pa < 0) { c.pa = (int)ia; c.pb = (int)ib; } }
                }
            }
            if (c.xa >= 0 || c.pa >= 0) cands.push_back(c);
        }
    }

    // 抽样：正样本只从"存在 bbox 相交 patch-pair"的 face-pair 中取；阴性从存在 bbox 不相交 patch-pair 中取
    std::vector<int> posIdx, negIdx;
    for (std::size_t i = 0; i < cands.size(); ++i) {
        if (cands[i].xa >= 0) posIdx.push_back((int)i);
        if (cands[i].pa >= 0) negIdx.push_back((int)i);
    }
    std::shuffle(posIdx.begin(), posIdx.end(), rng);
    std::shuffle(negIdx.begin(), negIdx.end(), rng);
    if ((int)posIdx.size() > npos) posIdx.resize(npos);
    if ((int)negIdx.size() > nneg) negIdx.resize(nneg);

    FILE* out = std::fopen(outcsv, "w");
    if (!out) { std::printf("cannot write %s\n", outcsv); return 1; }
    std::fprintf(out, "part,fidA,piA,fidB,piB,src,bbox_overlap,type,cop_overlap,coincidence,cover,certified,undecided,components,ms\n");
    std::fflush(out);

    auto emit = [&](int a, int ia, int b, int ib, const char* src, int bbox_ov) {
        const RBezierSurface& A = faces[a][ia];
        const RBezierSurface& B = faces[b][ib];
        auto t0 = std::chrono::high_resolution_clock::now();
        auto pc = nurbs::cert::classify_pair(A, B, tau, node_cap);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::fprintf(out, "%s,%d,%d,%d,%d,%s,%d,%d,%d,%.6g,%ld,%ld,%ld,%ld,%.3f\n",
                     path, a, ia, b, ib, src, bbox_ov, pc.type, pc.coplanar_overlap,
                     pc.coincidence, pc.cover, pc.certified, pc.undecided, pc.components, ms);
        std::fflush(out);
    };
    for (int idx : posIdx) emit(cands[idx].a, cands[idx].xa, cands[idx].b, cands[idx].xb, "pos", 1);
    for (int idx : negIdx) emit(cands[idx].a, cands[idx].pa, cands[idx].b, cands[idx].pb, "neg", 0);
    std::fclose(out);

    std::printf("PROBE data=%s faces=%d patches=%ld facepairs=%ld overlap_facepairs=%ld "
                "patchpairs=%ld overlap_patchpairs=%ld cand=%zu pos=%d neg=%d tau=%.4f cap=%ld\n",
                path, nf, totalPatches, facePairs, overlapFacePairs, patchPairs,
                overlapPatchPairs, cands.size(), (int)posIdx.size(), (int)negIdx.size(), tau, node_cap);
    return 0;
}
