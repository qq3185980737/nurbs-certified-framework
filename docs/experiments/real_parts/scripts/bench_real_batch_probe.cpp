// ============================================================================
// bench_real_batch_probe.cpp — bench_real_batch.cpp 的**插桩副本**（放在
// docs/experiments/real_parts/scripts/ 下，不改动 cpp/ 源文件）。
//   与 cpp/tests/bench_real_batch.cpp 输出完全相同的 parsed / RESULT / certified_boxes 行，
//   另外输出一行 PAIRSTAT：
//     - patch-pair 单次 classify_pair 耗时的 median / p95 / p99 / max（微秒）
//     - 触顶次数 cap_hit（cover >= node_cap/2）与 cover_gt2500（cover > 2500）
//     - patch-pair 级 type 分布与各 type 的总耗时
//   用途：给"真实零件基准"报告每件的中位/最大单对耗时与预算触顶次数。
// 用法: bench_real_batch_probe <data.txt> [tau] [node_cap] [detail.out]
// ============================================================================
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <algorithm>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::Interval;

struct AABB { double xlo=1e300, ylo=1e300, zlo=1e300, xhi=-1e300, yhi=-1e300, zhi=-1e300; };
static AABB patch_bb(const RBezierSurface& s) {
    AABB b;
    double wmin = s.wmin > 0 ? s.wmin : 1e-9;
    double wmax = s.wmax > 0 ? s.wmax : 1.0;
    for (const auto& p : s.P) {
        // 正确保守 AABB（w>0）：负坐标时 wx/wmin 才是下界、wx/wmax 才是上界。
        // 原 bench 版 xl=lo/wmax、xh=hi/wmin 对负坐标产生反向区间 → 空 AABB → 误跳过真实相交对。
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

static double pct(std::vector<double>& v, double q) {
    if (v.empty()) return 0.0;
    std::size_t k = (std::size_t)((v.size() - 1) * q);
    return v[k];
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: %s <data.txt> [tau] [node_cap] [detail.out]\n", argv[0]); return 1; }
    setvbuf(stdout, nullptr, _IOLBF, 0);   // 行缓冲：被 timeout 杀掉时已输出行不丢
    const char* path = argv[1];
    double tau = (argc >= 3) ? std::atof(argv[2]) : 0.02;
    long node_cap = (argc >= 4) ? std::atol(argv[3]) : 4000000L;
    const char* detail = (argc >= 5) ? argv[4] : nullptr;
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
    long totalPatches = 0; for (auto& v : faces) totalPatches += (int)v.size();
    std::printf("parsed %d faces (%ld patches), tau=%.4f\n", nf, totalPatches, tau);

    FILE* det = detail ? std::fopen(detail, "w") : nullptr;

    std::vector<AABB> faces_bb(nf);
    for (int a = 0; a < nf; ++a) {
        AABB u;
        for (const auto& A : faces[a]) {
            AABB ba = patch_bb(A);
            u.xlo = std::min(u.xlo, ba.xlo); u.xhi = std::max(u.xhi, ba.xhi);
            u.ylo = std::min(u.ylo, ba.ylo); u.yhi = std::max(u.yhi, ba.yhi);
            u.zlo = std::min(u.zlo, ba.zlo); u.zhi = std::max(u.zhi, ba.zhi);
        }
        faces_bb[a] = u;
    }

    long facePairs = 0, disjoint = 0, t1d = 0, t2d = 0, tcop = 0, tund = 0;
    long patchPairs = 0, bboxSkip = 0, certBoxes = 0, undecBoxes = 0, branches = 0;
    long pairType[5] = {0,0,0,0,0};
    double pairMsType[5] = {0,0,0,0,0};
    long capHit = 0, coverGt2500 = 0;
    std::vector<double> pairMs;
    double msTotal = 0.0;

    // 定期输出累计 PAIRSTAT（超时被杀时保留部分统计）
    auto emit_stats = [&](bool finalFlag) {
        std::vector<double> v = pairMs;
        std::sort(v.begin(), v.end());
        long np = (long)v.size();
        double med = np ? (np % 2 ? v[np/2] : 0.5*(v[np/2-1]+v[np/2])) : 0.0;
        std::printf("PAIRSTAT tau=%.4f cap=%ld n=%ld med_ms=%.3f p95_ms=%.3f p99_ms=%.3f max_ms=%.3f "
                    "cap_hit=%ld cover_gt2500=%ld pair_t0=%ld pair_t1=%ld pair_t2=%ld pair_t3=%ld pair_t4=%ld "
                    "und_pair_ratio=%.4f t0_ms=%.1f t1_ms=%.1f t2_ms=%.1f t3_ms=%.1f t4_ms=%.1f%s\n",
                    tau, node_cap, np, med, pct(v, 0.95), pct(v, 0.99), np ? v[np-1] : 0.0,
                    capHit, coverGt2500, pairType[0], pairType[1], pairType[2], pairType[3], pairType[4],
                    np ? (double)pairType[3] / np : 0.0,
                    pairMsType[0], pairMsType[1], pairMsType[2], pairMsType[3], pairMsType[4],
                    finalFlag ? "" : " partial=1");
    };
    auto tStart = std::chrono::high_resolution_clock::now();
    double nextEmit = 15.0;   // 秒
    for (int a = 0; a < nf; ++a) {
        for (int b = a + 1; b < nf; ++b) {
            ++facePairs;
            if (!aabb_overlap(faces_bb[a], faces_bb[b])) { ++disjoint; continue; }
            int fp_type = -1; long fp_cert = 0; int fp_comps = 0; int hit = 0;
            for (const auto& A : faces[a]) {
                AABB ba = patch_bb(A);
                for (const auto& B : faces[b]) {
                    ++patchPairs;
                    AABB bb = patch_bb(B);
                    if (!aabb_overlap(ba, bb)) { ++bboxSkip; continue; }
                    auto t0 = std::chrono::high_resolution_clock::now();
                    auto pc = nurbs::cert::classify_pair(A, B, tau, node_cap);
                    auto t1 = std::chrono::high_resolution_clock::now();
                    double ms = std::chrono::duration<double,std::milli>(t1-t0).count();
                    msTotal += ms;
                    pairMs.push_back(ms);
                    ++hit;
                    if (pc.cover >= node_cap / 2) ++capHit;
                    if (pc.cover > 2500) ++coverGt2500;
                    if (pc.type >= 0 && pc.type <= 4) { ++pairType[pc.type]; pairMsType[pc.type] += ms; }
                    if (pc.type != 0) {
                        fp_type = std::max(fp_type, pc.type);
                        fp_cert += pc.certified; fp_comps += (int)pc.components;
                    }
                    certBoxes += pc.certified; undecBoxes += pc.undecided;
                }
            }
            if (fp_type == 1) { ++t1d; branches += fp_comps; }
            else if (fp_type == 2) ++t2d;
            else if (fp_type == 4) ++tcop;
            else if (fp_type == 3) ++tund;
            else ++disjoint;
            if (det) { std::fprintf(det, "face %4d-%4d hit=%d type=%d cert=%ld comps=%d\n", a, b, hit, fp_type, fp_cert, fp_comps); std::fflush(det); }
            double el = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tStart).count();
            if (el >= nextEmit) { emit_stats(false); nextEmit = el + 15.0; }
        }
    }
    if (det) std::fclose(det);
    long pairTry = patchPairs - bboxSkip;
    std::printf("RESULT faces=%d facepairs=%ld | patchpairs=%ld bboxskip=%ld try=%ld | disjoint=%ld | 1D=%ld (comps=%ld) | 2D=%ld | coplanar=%ld | undecided=%ld\n",
                nf, facePairs, patchPairs, bboxSkip, pairTry, disjoint, t1d, branches, t2d, tcop, tund);
    std::printf("       certified_boxes=%ld undecided_boxes=%ld time=%.1f ms (%.3f ms/facepair, %.3f ms/try)\n",
                certBoxes, undecBoxes, msTotal, msTotal / std::max(1L, facePairs), msTotal / std::max(1L, pairTry));

    emit_stats(true);
    return 0;
}
