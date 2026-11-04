// ============================================================================
// bench_cert_batch.cpp — 路线 C 批量证书生成（基准实验）
//   mode=acis   <data.txt> <outdir> <tau> : 解析 SURF/NURB 块，逐 face-pair 生成证书
//   mode=random <N> <outdir> <tau>        : 生成 N 个随机相交对证书
//     - 60% 随机双线性平面片 × 随机双线性平面片（AABB 相交的保留）
//     - 40% 有理四分之一圆柱面片 × 平面 z=c
//   证书: {"tau","patches":[A,B],"tree":...,"certified":[]}
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <random>
#include <chrono>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::NurbsPatch;
using nurbs::cert::Box4;
using nurbs::Interval;

static RBezierSurface make_patch(int nu, int nv, const std::vector<std::array<double, 4>>& cart_w) {
    RBezierSurface s;
    s.nu = nu; s.nv = nv;
    s.P.resize((std::size_t)nu * nv);
    for (std::size_t k = 0; k < s.P.size(); ++k) {
        double w = cart_w[k][3];
        s.P[k][0] = Interval::point(cart_w[k][0] * w);
        s.P[k][1] = Interval::point(cart_w[k][1] * w);
        s.P[k][2] = Interval::point(cart_w[k][2] * w);
        s.P[k][3] = Interval::point(w);
    }
    return s;
}

static void put_hex(std::string& s, double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a", v);
    s += "\""; s += buf; s += "\"";
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
        bool k4 = (nurbs::cert::krawczyk_4d_code(rootA, rootB, box, false, nullptr, nullptr) == 0);
        bool ks = (nurbs::cert::krawczyk_slice(rootA, rootB, box) == true);
        s += ",\"kind\":\"leaf\",\"status\":\"output\"";
        s += k4 ? ",\"claimed_4d\":true" : ",\"claimed_4d\":false";
        s += ks ? ",\"claimed_slice\":true" : ",\"claimed_slice\":false";
        s += "}";
        return;
    }
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

static bool write_cert(const RBezierSurface& A, const RBezierSurface& B, double tau,
                       int max_depth, const std::string& path) {
    for (const auto& P : { A.P, B.P })
        for (const auto& p : P)
            if (p[3].mid() <= 0.0) return false;   // 原型只支持正权
    std::string s = "{\"tau\":";
    put_hex(s, tau);
    s += ",\"patches\":[";
    emit_patch(s, A); s += ","; emit_patch(s, B);
    s += "],\"tree\":";
    Box4 root{ Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(0.0, 1.0) };
    emit_tree(s, A, B, A, B, root, 0, tau, max_depth);
    s += ",\"certified\":[]}";
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fwrite(s.data(), 1, s.size(), f);
    std::fclose(f);
    return true;
}

static RBezierSurface plane_patch(double cx, double cy, double cz,
                                  const double U[3], const double V[3]) {
    // A(i,j) = c + (i-0.5)U + (j-0.5)V
    std::vector<std::array<double, 4>> P(4);
    int k = 0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            double a = (i - 0.5), b = (j - 0.5);
            P[k++] = { cx + a * U[0] + b * V[0], cy + a * U[1] + b * V[1], cz + a * U[2] + b * V[2], 1.0 };
        }
    return make_patch(2, 2, P);
}

static RBezierSurface quarter_cylinder() {
    const double w1 = std::sqrt(2.0) / 2.0;
    double xs[3] = {1.0, 1.0, 0.0}, ys[3] = {0.0, 1.0, 1.0}, ws[3] = {1.0, w1, 1.0};
    std::vector<std::array<double, 4>> P(6);
    for (int i = 0; i < 3; ++i) {
        P[i * 2 + 0] = { xs[i] * ws[i], ys[i] * ws[i], 0.0, ws[i] };
        P[i * 2 + 1] = { xs[i] * ws[i], ys[i] * ws[i], ws[i], ws[i] };
    }
    return make_patch(3, 2, P);
}

static RBezierSurface plane_z(double z) {
    std::vector<std::array<double, 4>> P(4);
    P[0] = { 0.0, 0.0, z, 1.0 };
    P[1] = { 0.0, 1.0, z, 1.0 };
    P[2] = { 1.0, 0.0, z, 1.0 };
    P[3] = { 1.0, 1.0, z, 1.0 };
    return make_patch(2, 2, P);
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::printf("usage: %s <acis|random> <data|N> <outdir> <tau> [max_depth]\n", argv[0]);
        return 1;
    }
    std::string mode = argv[1];
    std::string outdir = argv[3];
    double tau = std::atof(argv[4]);
    int max_depth = (argc >= 6) ? std::atoi(argv[5]) : 40;
    int made = 0;
    if (mode == "acis") {
        std::ifstream f(argv[2]);
        if (!f) { std::printf("cannot open %s\n", argv[2]); return 1; }
        std::vector<std::vector<RBezierSurface>> surfs;
        std::string tok;
        while (f >> tok) {
            if (tok == "NSURF") { int n; f >> n; surfs.resize(n); }
            else if (tok == "NURB") {
                int fid, conv, du, dv, nu, nv;
                f >> fid >> conv >> du >> dv >> nu >> nv;
                (void)conv;
                if (fid >= (int)surfs.size()) surfs.resize(fid + 1);
                NurbsPatch N; N.du = du; N.dv = dv; N.nu = nu; N.nv = nv;
                N.U.resize(nu + du + 1); for (auto& x : N.U) f >> x;
                N.V.resize(nv + dv + 1); for (auto& x : N.V) f >> x;
                N.P.resize((std::size_t)nu * nv);
                for (int i = 0; i < nu; ++i)
                    for (int j = 0; j < nv; ++j) {
                        double x, y, z, w; f >> x >> y >> z >> w;
                        nurbs::cert::HP hp;
                        hp[0] = Interval::point(x * w);
                        hp[1] = Interval::point(y * w);
                        hp[2] = Interval::point(z * w);
                        hp[3] = Interval::point(w);
                        N.P[(std::size_t)i * nv + j] = hp;
                    }
                std::string e; f >> e;
                for (auto& p : nurbs::cert::to_bezier_patches(N)) surfs[fid].push_back(p);
            } else if (tok == "SURF") {
                int fid, pi, du, dv, nu, nv;
                f >> fid >> pi >> du >> dv >> nu >> nv;
                std::vector<std::array<double, 4>> pts((std::size_t)nu * nv);
                for (int i = 0; i < nu; ++i)
                    for (int j = 0; j < nv; ++j) {
                        double x, y, z, w; f >> x >> y >> z >> w;
                        pts[(std::size_t)i * nv + j] = { x, y, z, w };
                    }
                std::string e; f >> e;
                surfs[fid].push_back(make_patch(nu, nv, pts));
            } else { std::string line; std::getline(f, line); }
        }
        int ns = (int)surfs.size();
        for (int a = 0; a < ns; ++a)
            for (int b = a + 1; b < ns; ++b)
                for (const auto& A : surfs[a])
                    for (const auto& B : surfs[b]) {
                        char name[512];
                        std::snprintf(name, sizeof name, "%s/acis_%02d_%02d_%d.json", outdir.c_str(), a, b, made);
                        if (write_cert(A, B, tau, max_depth, name)) ++made;
                    }
    } else {
        int N = std::atoi(argv[2]);
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        std::uniform_real_distribution<double> c01(0.05, 0.95);
        int tries = 0;
        while (made < N && tries < N * 20) {
            ++tries;
            RBezierSurface A, B;
            if (made % 5 == 0) {           // 20%: 有理圆柱 × 平面
                A = quarter_cylinder();
                B = plane_z(c01(rng));
            } else {                        // 80%: 随机平面 × 随机平面
                auto rnd_dir = [&]() {
                    double a = std::acos(uni(rng)), t = 3.141592653589793 * uni(rng);
                    double s = std::sin(a);
                    return std::array<double, 3>{ s * std::cos(t), s * std::sin(t), std::cos(a) };
                };
                double U1[3] = { 0, 0, 0 }, V1[3] = { 0, 0, 0 }, U2[3] = { 0, 0, 0 }, V2[3] = { 0, 0, 0 };
                auto d1 = rnd_dir(), d2 = rnd_dir();
                // 正交化
                auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
                    return std::array<double, 3>{ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
                };
                auto e1 = d1, e2 = cross(d1, d2);
                double n = std::sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
                if (n < 1e-9) continue;
                for (int k = 0; k < 3; ++k) { e2[k] /= n; U1[k] = e1[k]; V1[k] = e2[k]; }
                double e3[3];
                e3[0] = U1[1] * V1[2] - U1[2] * V1[1];
                e3[1] = U1[2] * V1[0] - U1[0] * V1[2];
                e3[2] = U1[0] * V1[1] - U1[1] * V1[0];
                auto d3 = rnd_dir(), d4 = cross(d3, d1);
                n = std::sqrt(d4[0] * d4[0] + d4[1] * d4[1] + d4[2] * d4[2]);
                if (n < 1e-9) continue;
                for (int k = 0; k < 3; ++k) { U2[k] = d3[k]; V2[k] = d4[k] / n; }
                A = plane_patch(uni(rng), uni(rng), uni(rng), U1, V1);
                B = plane_patch(uni(rng), uni(rng), uni(rng), U2, V2);
                (void)e3;
            }
            if (!nurbs::cert::disjoint(nurbs::cert::hull_aabb(A), nurbs::cert::hull_aabb(B))) {
                char name[512];
                std::snprintf(name, sizeof name, "%s/rand_%05d.json", outdir.c_str(), made);
                if (write_cert(A, B, tau, max_depth, name)) ++made;
            }
        }
    }
    std::printf("wrote %d certificates to %s (mode=%s tau=%.3f)\n", made, outdir.c_str(), mode.c_str(), tau);
    return 0;
}



