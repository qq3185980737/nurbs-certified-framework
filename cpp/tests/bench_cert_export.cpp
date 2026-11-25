// ============================================================================
// bench_cert_export.cpp — 路线 C 原型：导出 SSI 证书（几何 + 覆盖树）为 JSON
//   证书内容：
//     patches: 两个有理 Bezier 面片的齐次控制网（十六进制浮点，精确可逆）
//     tree:    4D 参数盒细分树；叶子 = excluded（附 AABB 不相交证书）或 output
//   供独立验证器（Python + Fraction 精确有理数）重放与检查。
//   演示用解析案例：平面 z=0 × 平面 x=0.5（交线：x=0.5, z=0 的直线）
// ============================================================================
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

static RBezierSurface plane_z0() {
    std::vector<std::array<double, 4>> P(4);
    // 齐次 (wx, wy, wz, w)：z=0 平面，u,v ∈ [0,1] → (2u-1, 2v-1, 0)
    P[0] = { -1.0, -1.0, 0.0, 1.0 };
    P[1] = { -1.0,  1.0, 0.0, 1.0 };
    P[2] = {  1.0, -1.0, 0.0, 1.0 };
    P[3] = {  1.0,  1.0, 0.0, 1.0 };
    return nurbs::cert::make_surface(2, 2, P);
}

static RBezierSurface plane_x(double x0) {
    std::vector<std::array<double, 4>> P(4);
    // x=x0 平面，s,t ∈ [0,1] → (x0, 2s-1, 2t-1)
    P[0] = { x0, -1.0, -1.0, 1.0 };
    P[1] = { x0, -1.0,  1.0, 1.0 };
    P[2] = { x0,  1.0, -1.0, 1.0 };
    P[3] = { x0,  1.0,  1.0, 1.0 };
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
                      const RBezierSurface& rootA, const RBezierSurface& rootB,
                      const Box4& box, int depth, double tau, int max_depth) {
    double wmax = std::max(std::max(box.u.width(), box.v.width()),
                           std::max(box.s.width(), box.t.width()));
    bool leaf = (wmax <= tau) || (depth >= max_depth);
    bool excl = !leaf && nurbs::cert::disjoint(nurbs::cert::hull_aabb(A), nurbs::cert::hull_aabb(B));
    s += "{\"box\":";
    put_box(s, box);
    if (leaf) {
        s += ",\"kind\":\"leaf\",\"status\":\"output\"";
        bool k4 = (nurbs::cert::krawczyk_4d_code(rootA, rootB, box, false, nullptr, nullptr) == 0);
        s += k4 ? ",\"claimed_4d\":true}" : ",\"claimed_4d\":false}";
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
        emit_tree(s, L, B, rootA, rootB, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, R, B, rootA, rootB, b2, depth + 1, tau, max_depth); s += "]}";
    } else if (axis == 1) {
        double mid = 0.5 * (box.v.lo + box.v.hi);
        nurbs::cert::split_v(A, L, R);
        b1.v = Interval(box.v.lo, mid); b2.v = Interval(mid, box.v.hi);
        s += ",\"kind\":\"node\",\"axis\":1,\"children\":[";
        emit_tree(s, L, B, rootA, rootB, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, R, B, rootA, rootB, b2, depth + 1, tau, max_depth); s += "]}";
    } else if (axis == 2) {
        double mid = 0.5 * (box.s.lo + box.s.hi);
        nurbs::cert::split_u(B, L, R);
        b1.s = Interval(box.s.lo, mid); b2.s = Interval(mid, box.s.hi);
        s += ",\"kind\":\"node\",\"axis\":2,\"children\":[";
        emit_tree(s, A, L, rootA, rootB, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, A, R, rootA, rootB, b2, depth + 1, tau, max_depth); s += "]}";
    } else {
        double mid = 0.5 * (box.t.lo + box.t.hi);
        nurbs::cert::split_v(B, L, R);
        b1.t = Interval(box.t.lo, mid); b2.t = Interval(mid, box.t.hi);
        s += ",\"kind\":\"node\",\"axis\":3,\"children\":[";
        emit_tree(s, A, L, rootA, rootB, b1, depth + 1, tau, max_depth); s += ",";
        emit_tree(s, A, R, rootA, rootB, b2, depth + 1, tau, max_depth); s += "]}";
    }
}

int main(int argc, char** argv) {
    const char* out = (argc >= 2) ? argv[1] : "cert.json";
    double tau = (argc >= 3) ? std::atof(argv[2]) : 0.25;
    auto A = plane_z0();
    auto B = plane_x(0.5);
    std::string s = "{\"tau\":";
    put_hex(s, tau);
    s += ",\"patches\":[";
    {
        auto put_patch = [&](const RBezierSurface& P) {
            s += "{\"nu\":";
            char nb[32]; std::snprintf(nb, sizeof nb, "%d", P.nu); s += nb;
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
        };
        put_patch(A); s += ","; put_patch(B);
    }
    s += "],\"tree\":";
    Box4 root{ Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0) };
    emit_tree(s, A, B, A, B, root, 0, tau, 40);
    s += ",\"certified\":[";
    {
        // 演示：4 个切片解（u=0.75, t=0.5, v=s=s0），hv=0.1 的盒
        const double s0s[4] = {0.125, 0.375, 0.625, 0.875};
        for (int k = 0; k < 4; ++k) {
            double s0 = s0s[k];
            Box4 cb;
            cb.u = Interval(0.75 - 0.1, 0.75 + 0.1);
            cb.v = Interval(s0, s0);
            cb.s = Interval(s0 - 0.1, s0 + 0.1);
            cb.t = Interval(0.5 - 0.1, 0.5 + 0.1);
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
    std::printf("certificate written: %s (%zu bytes, tau=%a)\n", out, s.size(), tau);
    return 0;
}
