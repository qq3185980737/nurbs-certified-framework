// ============================================================================
// bench_cert_export_types.cpp — 路线 C：重合 / 相切证书类型
//   mode = coinc  : 两片共面双线性参数块，重合区域用 4 个 4D 盒 + 仿射参数对应
//                   A(u,v)=(u,v,0); B(s,t)=(0.25+0.5s, 0.25+0.5t, 0)
//                   重合对应: s = 2u-0.5, t = 2v-0.5 (alpha=2,beta=0,gamma=-0.5,
//                             delta=0,epsilon=2,zeta=-0.5)
//   mode = tangent : A 平面 z=0；B 双二次抛物面 B(s,t)=(s,t,(s-0.5)^2+(t-0.5)^2)
//                   相切点 (u,v,s,t)=(0.5,0.5,0.5,0.5)
//   证书均含 patches + tree(root leaf output) + 额外字段
// ============================================================================
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <cstring>

#include "nurbs/cert_ssi.hpp"

using nurbs::cert::RBezierSurface;
using nurbs::cert::Box4;
using nurbs::Interval;

static void put_hex(std::string& s, double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a", v);
    s += "\""; s += buf; s += "\"";
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

static void put_box(std::string& s, const Box4& b) {
    const Interval* iv[4] = { &b.u, &b.v, &b.s, &b.t };
    s += "[";
    for (int k = 0; k < 4; ++k) {
        put_hex(s, iv[k]->lo); s += ","; put_hex(s, iv[k]->hi);
        if (k < 3) s += ",";
    }
    s += "]";
}

static RBezierSurface plane_uv() {
    std::vector<std::array<double,4>> P(4);
    P[0] = { 0.0, 0.0, 0.0, 1.0 };
    P[1] = { 0.0, 1.0, 0.0, 1.0 };
    P[2] = { 1.0, 0.0, 0.0, 1.0 };
    P[3] = { 1.0, 1.0, 0.0, 1.0 };
    return nurbs::cert::make_surface(2, 2, P);
}

static RBezierSurface plane_uv_25() {
    // B(s,t) = (0.25+0.5s, 0.25+0.5t, 0)
    std::vector<std::array<double,4>> P(4);
    P[0] = { 0.25, 0.25, 0.0, 1.0 };
    P[1] = { 0.25, 0.75, 0.0, 1.0 };
    P[2] = { 0.75, 0.25, 0.0, 1.0 };
    P[3] = { 0.75, 0.75, 0.0, 1.0 };
    return nurbs::cert::make_surface(2, 2, P);
}

static RBezierSurface paraboloid() {
    // B(s,t) = (s, t, (s-0.5)^2 + (t-0.5)^2)，双二次
    // x_ij = i/2, y_ij = j/2, z_ij = cf[i] + cg[j], cf=cg=[0.25,-0.25,0.25]
    std::vector<std::array<double,4>> P(9);
    double cf[3] = {0.25, -0.25, 0.25};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double x = 0.5 * i, y = 0.5 * j, z = cf[i] + cf[j];
            P[i * 3 + j] = { x, y, z, 1.0 };
        }
    return nurbs::cert::make_surface(3, 3, P);
}



static void emit_root_leaf(std::string& s) {
    s += "{\"box\":[";
    for (int k = 0; k < 4; ++k) {
        put_hex(s, 0.0); s += ","; put_hex(s, 1.0);
        if (k < 3) s += ",";
    }
    s += "],\"kind\":\"leaf\",\"status\":\"output\"}";
}

int main(int argc, char** argv) {
    const char* out = (argc >= 2) ? argv[1] : "cert_types.json";
    const char* mode = (argc >= 3) ? argv[2] : "coinc";
    std::string s;
    if (std::strcmp(mode, "coinc") == 0) {
        auto A = plane_uv();
        auto B = plane_uv_25();
        s = "{\"kind\":\"coincidence\",\"patches\":[";
        emit_patch(s, A); s += ","; emit_patch(s, B);
        s += "],\"tree\":";
        emit_root_leaf(s);
        s += ",\"coincidence\":{\"alpha\":";
        put_hex(s, 2.0); s += ",\"beta\":"; put_hex(s, 0.0);
        s += ",\"gamma\":"; put_hex(s, -0.5);
        s += ",\"delta\":"; put_hex(s, 0.0);
        s += ",\"epsilon\":"; put_hex(s, 2.0);
        s += ",\"zeta\":"; put_hex(s, -0.5);
        s += ",\"boxes\":[";
        double grid[2] = {0.0, 0.5};
        bool first = true;
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) {
                double s0 = grid[a], s1 = s0 + 0.5;
                double t0 = grid[b], t1 = t0 + 0.5;
                Box4 bx;
                bx.u = Interval(0.25 + 0.5 * s0, 0.25 + 0.5 * s1);
                bx.v = Interval(0.25 + 0.5 * t0, 0.25 + 0.5 * t1);
                bx.s = Interval(s0, s1);
                bx.t = Interval(t0, t1);
                if (!first) s += ",";
                first = false;
                s += "{\"box\":";
                put_box(s, bx);
                s += "}";
            }
        s += "]}}";
    } else {
        auto A = plane_uv();
        auto B = paraboloid();
        s = "{\"kind\":\"tangency\",\"patches\":[";
        emit_patch(s, A); s += ","; emit_patch(s, B);
        s += "],\"tree\":";
        emit_root_leaf(s);
        s += ",\"tangent_points\":[{\"point\":[";
        put_hex(s, 0.5); s += ","; put_hex(s, 0.5); s += ","; put_hex(s, 0.5); s += ","; put_hex(s, 0.5);
        s += "],\"box\":";
        Box4 bx;
        bx.u = Interval(0.4, 0.6); bx.v = Interval(0.4, 0.6);
        bx.s = Interval(0.4, 0.6); bx.t = Interval(0.4, 0.6);
        put_box(s, bx);
        s += "}]}";
    }
    FILE* f = std::fopen(out, "wb");
    if (!f) { std::printf("cannot write %s\n", out); return 1; }
    std::fwrite(s.data(), 1, s.size(), f);
    std::fclose(f);
    std::printf("typed certificate written: %s (%zu bytes, mode=%s)\n", out, s.size(), mode);
    return 0;
}
