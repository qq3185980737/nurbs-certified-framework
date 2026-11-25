// ============================================================================
// bench_cert_export_nurbs.cpp 鈥?璺嚎 C锛歅7锛圢URBS鈫払茅zier锛夎瘉涔﹂摼婕旂ず
//   杈撳叆锛氫袱涓甫鍐呴儴鑺傜偣鐨?NURBS锛堝弻绾挎€э紝U=V=[0,0,0.5,1,1]锛?//     A: 骞抽潰 z=0.5
//     B: 骞抽潰 z = x + 0.375
//   鍐呮牳鐢?P7 绮剧‘杞崲锛圔oehm 鑺傜偣鎻掑叆锛夊緱鍒?4脳4 B茅zier 鍧楋紱
//   鍙栧潡 (su=0, sv=0) 涓€瀵瑰鍑鸿瘉涔︼細nurbs 杈撳叆 + 鍧楁帶鍒剁綉 + 瑕嗙洊鏍?+ 鍒囩墖璁よ瘉鐩掋€?//   楠岃瘉鍣ㄧ嫭绔嬬敤绮剧‘鏈夌悊鏁伴噸鍋?Boehm 杞崲骞舵瘮瀵瑰潡鎺у埗缃戯紙P7 璇佷功锛夛紝
//   鍐嶉獙璇佽鐩栨爲涓?Krawczyk 璁よ瘉銆?// ============================================================================
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <array>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::NurbsPatch;
using nurbs::cert::Box4;
using nurbs::Interval;

// 3x3 鎺у埗鐐圭綉鏍硷細x_i = 0,0.5,1; y_j = 0,0.5,1
static NurbsPatch make_nurbs_plane(double z0, double slope_x) {
    NurbsPatch S;
    S.du = S.dv = 1;
    S.nu = S.nv = 3;
    S.U = {0.0, 0.0, 0.5, 1.0, 1.0};
    S.V = {0.0, 0.0, 0.5, 1.0, 1.0};
    S.P.resize(9);
    double xs[3] = {0.0, 0.5, 1.0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double x = xs[i], y = xs[j], z = z0 + slope_x * x;
            S.P[(std::size_t)i * 3 + j] = { Interval::point(x), Interval::point(y),
                           Interval::point(z), Interval::point(1.0) };
        }
    return S;
}

static void put_hex(std::string& s, double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a", v);
    s += "\"";
    s += buf;
    s += "\"";
}

static void put_box(std::string& s, const Box4& b) {
    const Interval* iv[4] = { &b.u, &b.v, &b.s, &b.t };
    s += "[";
    for (int k = 0; k < 4; ++k) {
        put_hex(s, iv[k]->lo); s += ","; put_hex(s, iv[k]->hi);
        if (k < 3) s += ",";
    }
    s += "]";
}

static void emit_tree(std::string& s, const RBezierSurface& A, const RBezierSurface& B,
                      const Box4& box, int depth, double tau, int max_depth) {
    double wmax = std::max(std::max(box.u.width(), box.v.width()),
                           std::max(box.s.width(), box.t.width()));
    bool leaf = (wmax <= tau) || (depth >= max_depth);
    bool excl = !leaf && nurbs::cert::disjoint(nurbs::cert::hull_aabb(A), nurbs::cert::hull_aabb(B));
    s += "{\"box\":";
    put_box(s, box);
    if (leaf) { s += ",\"kind\":\"leaf\",\"status\":\"output\"}"; return; }
    if (excl) { s += ",\"kind\":\"leaf\",\"status\":\"excluded\"}"; return; }
    int axis = 0;
    double w = box.u.width();
    if (box.v.width() > w) { w = box.v.width(); axis = 1; }
    if (box.s.width() > w) { w = box.s.width(); axis = 2; }
    if (box.t.width() > w) { axis = 3; }
    RBezierSurface L, R;
    Box4 b1 = box, b2 = box;
    if (axis == 0) {
        double mid = 0.5 * (box.u.lo + box.u.hi);
        nurbs::cert::split_u(A, L, R);
        b1.u = Interval(box.u.lo, mid); b2.u = Interval(mid, box.u.hi);
        s += ",\"kind\":\"node\",\"axis\":0,\"children\":[";
        emit_tree(s, L, B, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, R, B, b2, depth + 1, tau, max_depth); s += "]}";
    } else if (axis == 1) {
        double mid = 0.5 * (box.v.lo + box.v.hi);
        nurbs::cert::split_v(A, L, R);
        b1.v = Interval(box.v.lo, mid); b2.v = Interval(mid, box.v.hi);
        s += ",\"kind\":\"node\",\"axis\":1,\"children\":[";
        emit_tree(s, L, B, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, R, B, b2, depth + 1, tau, max_depth); s += "]}";
    } else if (axis == 2) {
        double mid = 0.5 * (box.s.lo + box.s.hi);
        nurbs::cert::split_u(B, L, R);
        b1.s = Interval(box.s.lo, mid); b2.s = Interval(mid, box.s.hi);
        s += ",\"kind\":\"node\",\"axis\":2,\"children\":[";
        emit_tree(s, A, L, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, A, R, b2, depth + 1, tau, max_depth); s += "]}";
    } else {
        double mid = 0.5 * (box.t.lo + box.t.hi);
        nurbs::cert::split_v(B, L, R);
        b1.t = Interval(box.t.lo, mid); b2.t = Interval(mid, box.t.hi);
        s += ",\"kind\":\"node\",\"axis\":3,\"children\":[";
        emit_tree(s, A, L, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, A, R, b2, depth + 1, tau, max_depth); s += "]}";
    }
}

static void emit_patch(std::string& s, const RBezierSurface& P) {
    char nb[32];
    s += "{\"nu\":";
    std::snprintf(nb, sizeof nb, "%d", P.nu); s += nb;
    s += ",\"nv\":";
    std::snprintf(nb, sizeof nb, "%d", P.nv); s += nb;
    s += ",\"poles\":[";
    for (std::size_t k = 0; k < P.P.size(); ++k) {
        s += "[";
        for (int c = 0; c < 4; ++c) {
            put_hex(s, P.P[k][c].mid());
            if (c < 3) s += ",";
        }
        s += "]";
        if (k + 1 < P.P.size()) s += ",";
    }
    s += "]}";
}

static void emit_nurbs(std::string& s, const NurbsPatch& S) {
    char nb[32];
    s += "{\"du\":";
    std::snprintf(nb, sizeof nb, "%d", S.du); s += nb;
    s += ",\"dv\":";
    std::snprintf(nb, sizeof nb, "%d", S.dv); s += nb;
    s += ",\"nu\":";
    std::snprintf(nb, sizeof nb, "%d", S.nu); s += nb;
    s += ",\"nv\":";
    std::snprintf(nb, sizeof nb, "%d", S.nv); s += nb;
    s += ",\"poles\":[";
    for (std::size_t k = 0; k < S.P.size(); ++k) {
        s += "[";
        for (int c = 0; c < 4; ++c) {
            put_hex(s, S.P[k][c].mid());
            if (c < 3) s += ",";
        }
        s += "]";
        if (k + 1 < S.P.size()) s += ",";
    }
    s += "],\"U\":[";
    for (std::size_t i = 0; i < S.U.size(); ++i) { put_hex(s, S.U[i]); if (i + 1 < S.U.size()) s += ","; }
    s += "],\"V\":[";
    for (std::size_t i = 0; i < S.V.size(); ++i) { put_hex(s, S.V[i]); if (i + 1 < S.V.size()) s += ","; }
    s += "]}";
}

int main(int argc, char** argv) {
    const char* out = (argc >= 2) ? argv[1] : "cert_nurbs.json";
    double tau = (argc >= 3) ? std::atof(argv[2]) : 0.2;
    auto NA = make_nurbs_plane(0.5, 0.0);      // z = 0.5
    auto NB = make_nurbs_plane(0.375, 1.0);    // z = x + 0.375
    auto blocksA = nurbs::cert::to_bezier_patches(NA);
    auto blocksB = nurbs::cert::to_bezier_patches(NB);
    // 鍙?(su=0,sv=0) 鍧楋細block index 0
    const RBezierSurface& A = blocksA[0];
    const RBezierSurface& B = blocksB[0];
    std::string s = "{\"tau\":";
    put_hex(s, tau);
    s += ",\"nurbs\":[";
    emit_nurbs(s, NA); s += ","; emit_nurbs(s, NB);
    s += "],\"patches\":[";
    emit_patch(s, A); s += ","; emit_patch(s, B);
    s += "],\"tree\":";
    Box4 root{ Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0) };
    emit_tree(s, A, B, root, 0, tau, 40);
    s += ",\"certified\":[";
    // 浜ょ嚎锛堝潡鍧愭爣锛夛細u=s=0.25, v=t
    double t0s[4] = {0.125, 0.375, 0.625, 0.875};
    for (int k = 0; k < 4; ++k) {
        double t0 = t0s[k];
        Box4 cb;
        cb.u = Interval(0.25, 0.25);
        cb.v = Interval(t0 - 0.1, t0 + 0.1);
        cb.s = Interval(0.25, 0.25);
        cb.t = Interval(t0 - 0.1, t0 + 0.1);
        if (k) s += ",";
        s += "{\"box\":";
        put_box(s, cb);
        s += ",\"kind\":\"slice_existence\"}";
    }
    s += "]}\n";
    FILE* f = std::fopen(out, "wb");
    if (!f) { std::printf("cannot write %s\n", out); return 1; }
    std::fwrite(s.data(), 1, s.size(), f);
    std::fclose(f);
    std::printf("nurbs certificate written: %s (%zu bytes, tau=%a, blocksA=%zu blocksB=%zu)\n",
                out, s.size(), tau, blocksA.size(), blocksB.size());
    return 0;
}

