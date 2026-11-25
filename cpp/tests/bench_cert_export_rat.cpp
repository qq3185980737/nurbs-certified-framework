// ============================================================================
// bench_cert_export_rat.cpp 鈥?璺嚎 C 鏈夌悊鏉冮噸婕旂ず璇佷功锛?//   鏈夌悊鍦嗘煴闈㈢墖锛堝洓鍒嗕箣涓€鍦?脳 绾挎€э級 脳 骞抽潰 z=0.5
//   浜ょ嚎锛歷=0.5, (s,t)=(x(u),y(u))锛寀鈭圼0,1]
//   瀵煎嚭锛氬嚑浣曪紙榻愭鎺у埗缃戯紝hex float锛? 4D 瑕嗙洊鏍?+ 鍒囩墖璁よ瘉鐩?// ============================================================================
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <array>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;

// 鍥涘垎涔嬩竴鍦嗭紙鏈夌悊浜屾锛夛細鎺у埗鐐?(1,0),(1,1),(0,1)锛屾潈閲?1,鈭?/2,1锛泎 绾挎€?0..1
static RBezierSurface cylinder_quarter() {
    const double w1 = std::sqrt(2.0) / 2.0;
    std::vector<std::array<double,4>> P(6);
    double xs[3] = {1.0, 1.0, 0.0};
    double ys[3] = {0.0, 1.0, 1.0};
    double ws[3] = {1.0, w1, 1.0};
    for (int i = 0; i < 3; ++i) {
        P[i * 2 + 0] = { xs[i] * ws[i], ys[i] * ws[i], 0.0, ws[i] };   // v=0 鈫?z=0
        P[i * 2 + 1] = { xs[i] * ws[i], ys[i] * ws[i], ws[i], ws[i] }; // v=1 鈫?z=1
    }
    return nurbs::cert::make_surface(3, 2, P);
}

static RBezierSurface plane_z05() {
    std::vector<std::array<double,4>> P(4);
    P[0] = { 0.0, 0.0, 0.5, 1.0 };
    P[1] = { 0.0, 1.0, 0.5, 1.0 };
    P[2] = { 1.0, 0.0, 0.5, 1.0 };
    P[3] = { 1.0, 1.0, 0.5, 1.0 };
    return nurbs::cert::make_surface(2, 2, P);
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
    if (leaf) {
        s += ",\"kind\":\"leaf\",\"status\":\"output\"}";
        return;
    }
    if (excl) {
        s += ",\"kind\":\"leaf\",\"status\":\"excluded\"}";
        return;
    }
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

int main(int argc, char** argv) {
    const char* out = (argc >= 2) ? argv[1] : "cert_rat.json";
    double tau = (argc >= 3) ? std::atof(argv[2]) : 0.2;
    auto A = cylinder_quarter();
    auto B = plane_z05();
    std::string s = "{\"tau\":";
    put_hex(s, tau);
    s += ",\"patches\":[";
    emit_patch(s, A); s += ","; emit_patch(s, B);
    s += "],\"tree\":";
    Box4 root{ Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0) };
    emit_tree(s, A, B, root, 0, tau, 40);
    s += ",\"certified\":[";
    {
        double u0s[4] = {0.125, 0.375, 0.625, 0.875};
        for (int k = 0; k < 4; ++k) {
            double u0 = u0s[k];
            Interval U(u0, u0), V(0.5, 0.5);
            Interval S1[3], S1u[3], S1v[3];
            nurbs::cert::eval_rational_duv(A, U, V, S1, S1u, S1v);
            double x0 = S1[0].mid(), y0 = S1[1].mid();
            Box4 cb;
            cb.u = Interval(u0, u0);
            cb.v = Interval(0.4, 0.6);
            cb.s = Interval(x0 - 0.1, x0 + 0.1);
            cb.t = Interval(y0 - 0.1, y0 + 0.1);
            auto clamp = [](Interval& iv) {
                if (iv.lo < 0.0) iv = Interval(0.0, iv.hi);
                if (iv.hi > 1.0) iv = Interval(iv.lo, 1.0);
            };
            clamp(cb.u); clamp(cb.v); clamp(cb.s); clamp(cb.t);
            if (k) s += ",";
            s += "{\"box\":";
            put_box(s, cb);
            s += ",\"kind\":\"slice_existence\"}";
        }
    }
    s += "]}\n";
    FILE* f = std::fopen(out, "wb");
    if (!f) { std::printf("cannot write %s\n", out); return 1; }
    std::fwrite(s.data(), 1, s.size(), f);
    std::fclose(f);
    std::printf("rational certificate written: %s (%zu bytes, tau=%a)\n", out, s.size(), tau);
    return 0;
}

