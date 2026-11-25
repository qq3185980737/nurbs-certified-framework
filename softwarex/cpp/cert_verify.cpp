// ============================================================================
// cert_verify.cpp -- Route C: independent C++ exact-rational certificate verifier
//   Build (WSL): g++ -O3 -march=native -std=c++17 cert_verify.cpp -o cert_verify -lgmpxx -lgmp
//   Usage: cert_verify <cert.json> [--chunked] [--jobs N]
//          cert_verify --batch <dir>
// ============================================================================
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <gmpxx.h>

using Q = mpq_class;

// ---------------------------------------------------------------- JSON
struct J {
    enum Type { NUL, BOOL, NUM, STR, ARR, OBJ } t = NUL;
    bool b = false; double num = 0.0; std::string s;
    std::vector<J> a;
    std::vector<std::pair<std::string, J>> o;
    const J* get(const std::string& k) const {
        for (auto& kv : o) if (kv.first == k) return &kv.second;
        return nullptr;
    }
    bool has(const std::string& k) const { return get(k) != nullptr; }
};

struct Parser {
    const std::string& s; size_t i = 0;
    explicit Parser(const std::string& src) : s(src) {}
    void ws() { while (i < s.size() && std::isspace((unsigned char)s[i])) ++i; }
    char peek() { ws(); return i < s.size() ? s[i] : '\0'; }
    void expect(char c) {
        ws();
        if (i >= s.size() || s[i] != c) throw std::runtime_error("json parse error");
        ++i;
    }
    std::string str() {
        expect('"');
        std::string out;
        while (i < s.size() && s[i] != '"') {
            char c = s[i++];
            if (c != '\\') { out.push_back(c); continue; }
            char e = s[i++];
            switch (e) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'n': out.push_back('\n'); break;
                case 't': out.push_back('\t'); break;
                case 'r': out.push_back('\r'); break;
                default: out.push_back('?'); break;
            }
        }
        ++i; return out;
    }
    J parse() {
        char c = peek(); J j;
        if (c == '{') {
            j.t = J::OBJ; ++i;
            if (peek() == '}') { ++i; return j; }
            while (true) {
                std::string k = str(); expect(':');
                j.o.emplace_back(k, parse());
                if (peek() == ',') { ++i; continue; }
                expect('}'); break;
            }
        } else if (c == '[') {
            j.t = J::ARR; ++i;
            if (peek() == ']') { ++i; return j; }
            while (true) {
                j.a.push_back(parse());
                if (peek() == ',') { ++i; continue; }
                expect(']'); break;
            }
        } else if (c == '"') { j.t = J::STR; j.s = str(); }
        else if (c == 't' || c == 'f') {
            j.t = J::BOOL;
            if (s.compare(i, 4, "true") == 0) { j.b = true; i += 4; } else { j.b = false; i += 5; }
        } else if (c == 'n') { j.t = J::NUL; i += 4; }
        else {
            j.t = J::NUM; size_t st = i;
            if (s[i] == '-') ++i;
            while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i] == '.' || s[i] == 'e' ||
                                    s[i] == 'E' || s[i] == '+' || s[i] == '-')) ++i;
            j.num = std::stod(s.substr(st, i - st));
        }
        return j;
    }
};

static J parse_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    std::string src = ss.str(); Parser p(src); return p.parse();
}

// ------------------------------------------------------------- rationals
static Q q_from_double(double d) {
    if (d == 0.0) return Q(0);
    uint64_t bits; std::memcpy(&bits, &d, 8);
    bool neg = (bits >> 63) != 0;
    int exp = (int)((bits >> 52) & 0x7FF);
    uint64_t frac = bits & ((1ULL << 52) - 1);
    mpz_class mant, den = 1;
    if (exp == 0) { mant = mpz_class((unsigned long)frac); exp = -1074; }
    else { mant = mpz_class((unsigned long)(frac | (1ULL << 52))); exp = exp - 1075; }
    if (exp >= 0) mant <<= exp; else den <<= (-exp);
    Q q(mant, den); q.canonicalize();
    return neg ? -q : q;
}
static Q parse_rat(const std::string& s) {
    if (s.find('/') != std::string::npos) {
        Q q; mpq_set_str(q.get_mpq_t(), s.c_str(), 10); q.canonicalize(); return q;
    }
    return q_from_double(std::stod(s));
}
static double qd(const Q& q) { return q.get_d(); }

// -------------------------------------------------------------- intervals
struct I {
    Q lo, hi;
    I() : lo(0), hi(0) {}
    I(const Q& a, const Q& b) : lo(a), hi(b) {}
    static I pt(const Q& a) { return I(a, a); }
    Q width() const { return hi - lo; }
};
static I iadd(const I& a, const I& b) { return I(a.lo + b.lo, a.hi + b.hi); }
static I isub(const I& a, const I& b) { return I(a.lo - b.hi, a.hi - b.lo); }
static I imul(const I& a, const I& b) {
    Q p[4] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    Q lo = p[0], hi = p[0];
    for (int k = 1; k < 4; ++k) { lo = std::min(lo, p[k]); hi = std::max(hi, p[k]); }
    return I(lo, hi);
}
static I idiv(const I& a, const I& b) {
    if (b.lo <= 0 && 0 <= b.hi) throw std::runtime_error("interval division by zero");
    return imul(a, I(1 / b.hi, 1 / b.lo));
}

// ---------------------------------------------------------------- patch
using HP = std::array<Q, 4>;
struct Patch {
    int nu = 0, nv = 0;
    std::vector<HP> P;
    Q wmin, wmax;
    const HP& at(int i, int j) const { return P[(size_t)i * nv + j]; }
};
static Patch load_patch(const J& obj) {
    Patch p;
    p.nu = (int)obj.get("nu")->num;
    p.nv = (int)obj.get("nv")->num;
    p.P.resize((size_t)p.nu * p.nv);
    size_t idx = 0;
    for (auto& row : obj.get("poles")->a) {
        HP hp;
        for (int c = 0; c < 4; ++c) hp[c] = parse_rat(row.a[c].s);
        p.P[idx++] = hp;
    }
    p.wmin = p.wmax = p.P[0][3];
    for (auto& hp : p.P) { p.wmin = std::min(p.wmin, hp[3]); p.wmax = std::max(p.wmax, hp[3]); }
    return p;
}
static void split1d(const std::vector<HP>& pts, std::vector<HP>& L, std::vector<HP>& R) {
    int n = (int)pts.size();
    std::vector<std::vector<HP>> tri; tri.push_back(pts);
    Q h(1, 2);
    for (int r = 1; r < n; ++r) {
        std::vector<HP> row(n - r);
        for (int i = 0; i < n - r; ++i)
            for (int c = 0; c < 4; ++c) row[i][c] = tri[r - 1][i][c] * h + tri[r - 1][i + 1][c] * h;
        tri.push_back(std::move(row));
    }
    L.resize(n); R.resize(n);
    for (int r = 0; r < n; ++r) { L[r] = tri[r][0]; R[r] = tri[n - 1 - r][r]; }
}
static Patch split_u(const Patch& s, bool right) {
    Patch r; r.nu = s.nu; r.nv = s.nv; r.P.resize(s.P.size()); r.wmin = s.wmin; r.wmax = s.wmax;
    for (int j = 0; j < s.nv; ++j) {
        std::vector<HP> col(s.nu), l, rr;
        for (int i = 0; i < s.nu; ++i) col[i] = s.at(i, j);
        split1d(col, l, rr);
        const std::vector<HP>& src = right ? rr : l;
        for (int i = 0; i < s.nu; ++i) r.P[(size_t)i * s.nv + j] = src[i];
    }
    return r;
}
static Patch split_v(const Patch& s, bool right) {
    Patch r; r.nu = s.nu; r.nv = s.nv; r.P.resize(s.P.size()); r.wmin = s.wmin; r.wmax = s.wmax;
    for (int i = 0; i < s.nu; ++i) {
        std::vector<HP> row(s.nv), l, rr;
        for (int j = 0; j < s.nv; ++j) row[j] = s.at(i, j);
        split1d(row, l, rr);
        const std::vector<HP>& src = right ? rr : l;
        for (int j = 0; j < s.nv; ++j) r.P[(size_t)i * s.nv + j] = src[j];
    }
    return r;
}
static std::array<Q, 4> point_h(const Patch& s, const Q& u, const Q& v) {
    std::vector<std::vector<HP>> net(s.nu, std::vector<HP>(s.nv));
    for (int i = 0; i < s.nu; ++i)
        for (int j = 0; j < s.nv; ++j) net[i][j] = s.at(i, j);
    for (int r = 1; r < s.nu; ++r)
        for (int i = 0; i + r < s.nu; ++i)
            for (int j = 0; j < s.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    net[i][j][c] = net[i][j][c] * (1 - u) + net[i + 1][j][c] * u;
    for (int r = 1; r < s.nv; ++r)
        for (int i = 0; i < s.nu; ++i)
            for (int j = 0; j + r < s.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    net[i][j][c] = net[i][j][c] * (1 - v) + net[i][j + 1][c] * v;
    return net[0][0];
}
static std::array<Q, 3> point_eval(const Patch& s, const Q& u, const Q& v) {
    auto h = point_h(s, u, v);
    return {h[0] / h[3], h[1] / h[3], h[2] / h[3]};
}
static std::array<I, 4> eval_net(const Patch& s, const I& U, const I& V) {
    std::vector<std::vector<std::array<I, 4>>> net(s.nu, std::vector<std::array<I, 4>>(s.nv));
    for (int i = 0; i < s.nu; ++i)
        for (int j = 0; j < s.nv; ++j)
            for (int c = 0; c < 4; ++c) net[i][j][c] = I::pt(s.at(i, j)[c]);
    for (int r = 1; r < s.nu; ++r)
        for (int i = 0; i + r < s.nu; ++i)
            for (int j = 0; j < s.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    net[i][j][c] = iadd(net[i][j][c], imul(isub(net[i + 1][j][c], net[i][j][c]), U));
    for (int r = 1; r < s.nv; ++r)
        for (int i = 0; i < s.nu; ++i)
            for (int j = 0; j + r < s.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    net[i][j][c] = iadd(net[i][j][c], imul(isub(net[i][j + 1][c], net[i][j][c]), V));
    std::array<I, 4> out;
    for (int c = 0; c < 4; ++c) out[c] = net[0][0][c];
    return out;
}
static I clampw(const Patch& s, const I& w) { return I(std::max(w.lo, s.wmin), std::min(w.hi, s.wmax)); }
static std::array<I, 3> cart_iv(const Patch& s, const I& U, const I& V) {
    auto h = eval_net(s, U, V); I w = clampw(s, h[3]);
    std::array<I, 3> o;
    for (int c = 0; c < 3; ++c) o[c] = idiv(h[c], w);
    return o;
}
static Patch deriv(const Patch& s, int axis) {
    Patch d;
    if (axis == 0) {
        if (s.nu < 2) return d;
        d.nu = s.nu - 1; d.nv = s.nv; d.P.resize((size_t)d.nu * d.nv); d.wmin = s.wmin; d.wmax = s.wmax;
        for (int i = 0; i < d.nu; ++i)
            for (int j = 0; j < d.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    d.P[(size_t)i * d.nv + j][c] = (s.at(i + 1, j)[c] - s.at(i, j)[c]) * (s.nu - 1);
    } else {
        if (s.nv < 2) return d;
        d.nu = s.nu; d.nv = s.nv - 1; d.P.resize((size_t)d.nu * d.nv); d.wmin = s.wmin; d.wmax = s.wmax;
        for (int i = 0; i < d.nu; ++i)
            for (int j = 0; j < d.nv; ++j)
                for (int c = 0; c < 4; ++c)
                    d.P[(size_t)i * d.nv + j][c] = (s.at(i, j + 1)[c] - s.at(i, j)[c]) * (s.nv - 1);
    }
    return d;
}
static std::array<I, 3> du_iv(const Patch& s, const Patch& d, const I& U, const I& V) {
    if (d.nu == 0 || d.nv == 0) return {I::pt(0), I::pt(0), I::pt(0)};
    auto h = eval_net(s, U, V); auto hu = eval_net(d, U, V); I w = clampw(s, h[3]);
    std::array<I, 3> o;
    for (int c = 0; c < 3; ++c) o[c] = idiv(isub(imul(hu[c], w), imul(h[c], hu[3])), imul(w, w));
    return o;
}
static std::array<Q, 3> deriv_point(const Patch& s, const Patch& d, const Q& u, const Q& v, bool u_axis) {
    if (d.nu == 0 || d.nv == 0) return {Q(0), Q(0), Q(0)};
    auto h = point_h(s, u, v); auto hd = point_h(d, u, v);
    std::array<Q, 3> o;
    for (int c = 0; c < 3; ++c) o[c] = (hd[c] * h[3] - h[c] * hd[3]) / (h[3] * h[3]);
    (void)u_axis;
    return o;
}

// 3x3 inverse over Q
static bool inv3(const Q m[3][3], Q out[3][3]) {
    Q det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if (det == 0) return false;
    Q cof[3][3] = {
        {m[1][1] * m[2][2] - m[1][2] * m[2][1], -(m[1][0] * m[2][2] - m[1][2] * m[2][0]), m[1][0] * m[2][1] - m[1][1] * m[2][0]},
        {-(m[0][1] * m[2][2] - m[0][2] * m[2][1]), m[0][0] * m[2][2] - m[0][2] * m[2][0], -(m[0][0] * m[2][1] - m[0][1] * m[2][0])},
        {m[0][1] * m[1][2] - m[0][2] * m[1][1], -(m[0][0] * m[1][2] - m[0][2] * m[1][0]), m[0][0] * m[1][1] - m[0][1] * m[1][0]}};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) out[c][r] = cof[r][c] / det;
    return true;
}

// ------------------------------------------------------------ AABB / dist
static bool disjoint(const std::array<I, 3>& a, const std::array<I, 3>& b) {
    for (int c = 0; c < 3; ++c) if (a[c].hi < b[c].lo || b[c].hi < a[c].lo) return true;
    return false;
}
static std::array<I, 3> hull_aabb(const Patch& s) {
    std::array<I, 3> a;
    bool first = true;
    for (auto& hp : s.P) {
        std::array<I, 3> c;
        for (int k = 0; k < 3; ++k) c[k] = idiv(I::pt(hp[k]), I::pt(hp[3]));
        if (first) { a = c; first = false; }
        else for (int k = 0; k < 3; ++k) {
            a[k].lo = std::min(a[k].lo, c[k].lo);
            a[k].hi = std::max(a[k].hi, c[k].hi);
        }
    }
    return a;
}
static Q seg_dist2(const Q p[3], const Q a[3], const Q b[3]) {
    Q dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
    Q L2 = dx * dx + dy * dy + dz * dz;
    if (L2 == 0) return (p[0] - a[0]) * (p[0] - a[0]) + (p[1] - a[1]) * (p[1] - a[1]) + (p[2] - a[2]) * (p[2] - a[2]);
    Q dot = (p[0] - a[0]) * dx + (p[1] - a[1]) * dy + (p[2] - a[2]) * dz;
    if (dot <= 0) return (p[0] - a[0]) * (p[0] - a[0]) + (p[1] - a[1]) * (p[1] - a[1]) + (p[2] - a[2]) * (p[2] - a[2]);
    if (dot >= L2) return (p[0] - b[0]) * (p[0] - b[0]) + (p[1] - b[1]) * (p[1] - b[1]) + (p[2] - b[2]) * (p[2] - b[2]);
    Q base = (p[0] - a[0]) * (p[0] - a[0]) + (p[1] - a[1]) * (p[1] - a[1]) + (p[2] - a[2]) * (p[2] - a[2]);
    return base - dot * dot / L2;
}

// ------------------------------------------------------------- Krawczyk
static bool krawczyk_slice(const Patch& A, const Patch& B, const Patch& dAu, const Patch& dAv,
                           const Patch& dBu, const Patch& dBv, const I& U, const I& V, const I& S,
                           const I& T) {
    I X[4] = {U, V, S, T};
    for (int kfix = 0; kfix < 4; ++kfix) {
        int idx[3], t = 0;
        for (int c = 0; c < 4; ++c) if (c != kfix) idx[t++] = c;
        try {
            auto Au = du_iv(A, dAu, U, V), Av = du_iv(A, dAv, U, V);
            auto Bu = du_iv(B, dBu, S, T), Bv = du_iv(B, dBv, S, T);
            I J[3][4];
            for (int c = 0; c < 3; ++c) {
                J[c][0] = Au[c]; J[c][1] = Av[c];
                J[c][2] = I(-Bu[c].hi, -Bu[c].lo);
                J[c][3] = I(-Bv[c].hi, -Bv[c].lo);
            }
            Q m[4] = {U.lo + U.width() / 2, V.lo + V.width() / 2, S.lo + S.width() / 2, T.lo + T.width() / 2};
            Q Fm[3];
            {
                auto pa = point_eval(A, m[0], m[1]);
                auto pb = point_eval(B, m[2], m[3]);
                for (int c = 0; c < 3; ++c) Fm[c] = pa[c] - pb[c];
            }
            Q J3m[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) J3m[r][c] = (J[r][idx[c]].lo + J[r][idx[c]].hi) / 2;
            Q Y[3][3];
            if (!inv3(J3m, Y)) continue;
            I X3[3];
            for (int i = 0; i < 3; ++i) X3[i] = X[idx[i]];
            I W[3][3];
            for (int r = 0; r < 3; ++r)
                for (int d = 0; d < 3; ++d) {
                    I acc = (r == d) ? I::pt(1) : I::pt(0);
                    for (int c = 0; c < 3; ++c) acc = isub(acc, imul(I::pt(Y[r][c]), J[c][idx[d]]));
                    W[r][d] = acc;
                }
            I K[3];
            for (int r = 0; r < 3; ++r) {
                I acc = I::pt(m[idx[r]]);
                for (int c = 0; c < 3; ++c) acc = isub(acc, imul(I::pt(Y[r][c]), I::pt(Fm[c])));
                for (int d = 0; d < 3; ++d)
                    acc = iadd(acc, imul(W[r][d], I(X3[d].lo - m[idx[d]], X3[d].hi - m[idx[d]])));
                K[r] = acc;
            }
            bool inside = true;
            for (int r = 0; r < 3; ++r) {
                if (!(X3[r].lo <= K[r].lo && K[r].hi <= X3[r].hi)) { inside = false; break; }
                Q w = K[r].width(), xw = X3[r].width();
                if (!(w <= Q(9, 10) * xw || xw < Q(1, 1000000000000))) { inside = false; break; }
            }
            if (inside) return true;
        } catch (...) {
            continue;
        }
    }
    return false;
}

static bool krawczyk4d(const Patch& A, const Patch& B, const Patch& dAu, const Patch& dAv,
                       const Patch& dBu, const Patch& dBv, const I& U, const I& V, const I& S,
                       const I& T) {
    try {
        auto Au = du_iv(A, dAu, U, V), Av = du_iv(A, dAv, U, V);
        auto Bu = du_iv(B, dBu, S, T), Bv = du_iv(B, dBv, S, T);
        I J[3][4];
        for (int c = 0; c < 3; ++c) {
            J[c][0] = Au[c]; J[c][1] = Av[c];
            J[c][2] = I(-Bu[c].hi, -Bu[c].lo);
            J[c][3] = I(-Bv[c].hi, -Bv[c].lo);
        }
        I X[4] = {U, V, S, T};
        Q m[4] = {U.lo + U.width() / 2, V.lo + V.width() / 2, S.lo + S.width() / 2, T.lo + T.width() / 2};
        auto pa = point_eval(A, m[0], m[1]);
        auto pb = point_eval(B, m[2], m[3]);
        Q Fm[3];
        for (int c = 0; c < 3; ++c) Fm[c] = pa[c] - pb[c];
        Q Jm[3][4];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c) Jm[r][c] = (J[r][c].lo + J[r][c].hi) / 2;
        Q M[3][3];
        for (int r = 0; r < 3; ++r)
            for (int s2 = 0; s2 < 3; ++s2) {
                Q acc = 0;
                for (int d = 0; d < 4; ++d) acc += Jm[r][d] * Jm[s2][d];
                M[r][s2] = acc;
            }
        Q Minv[3][3];
        if (!inv3(M, Minv)) return false;
        Q Y[4][3];
        for (int a = 0; a < 4; ++a)
            for (int c = 0; c < 3; ++c) {
                Q acc = 0;
                for (int d = 0; d < 3; ++d) acc += Jm[d][a] * Minv[d][c];
                Y[a][c] = acc;
            }
        I W[4][4];
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                I acc = (a == b) ? I::pt(1) : I::pt(0);
                for (int c = 0; c < 3; ++c) acc = isub(acc, imul(I::pt(Y[a][c]), J[c][b]));
                W[a][b] = acc;
            }
        I K[4];
        for (int a = 0; a < 4; ++a) {
            I acc = I::pt(m[a]);
            for (int c = 0; c < 3; ++c) acc = isub(acc, imul(I::pt(Y[a][c]), I::pt(Fm[c])));
            for (int b = 0; b < 4; ++b) acc = iadd(acc, imul(W[a][b], I(X[b].lo - m[b], X[b].hi - m[b])));
            K[a] = acc;
        }
        for (int a = 0; a < 4; ++a)
            if (!(X[a].lo <= K[a].lo && K[a].hi <= X[a].hi)) return false;
        return true;
    } catch (...) {
        return false;
    }
}

// Newton polish + 小盒精确 slice Krawczyk
static bool slice_claim_ok(const Patch& A, const Patch& B, const Patch& dAu, const Patch& dAv,
                           const Patch& dBu, const Patch& dBv,
                           const std::array<std::array<Q, 2>, 4>& box) {
    auto Fd = [&](const Q& u, const Q& v, const Q& s, const Q& t, int c) {
        auto a = point_eval(A, u, v);
        auto b = point_eval(B, s, t);
        return qd(a[c] - b[c]);
    };
    for (int kfix = 0; kfix < 4; ++kfix) {
        int idx[3], ti = 0;
        for (int c = 0; c < 4; ++c) if (c != kfix) idx[ti++] = c;
        double p[4];
        for (int c = 0; c < 4; ++c) p[c] = qd((box[c][0] + box[c][1]) / 2);
        bool ok = true;
        for (int it = 0; it < 40; ++it) {
            Q u = q_from_double(p[0]), v = q_from_double(p[1]), s = q_from_double(p[2]), t = q_from_double(p[3]);
            double Fv[3];
            for (int c = 0; c < 3; ++c) Fv[c] = qd(point_eval(A, u, v)[c] - point_eval(B, s, t)[c]);
            if (std::max({std::fabs(Fv[0]), std::fabs(Fv[1]), std::fabs(Fv[2])}) < 1e-14) break;
            auto Au = deriv_point(A, dAu, u, v, true);
            auto Av = deriv_point(A, dAv, u, v, false);
            auto Bu = deriv_point(B, dBu, s, t, true);
            auto Bv = deriv_point(B, dBv, s, t, false);
            double J[3][4];
            for (int c = 0; c < 3; ++c) {
                J[c][0] = qd(Au[c]); J[c][1] = qd(Av[c]);
                J[c][2] = -qd(Bu[c]); J[c][3] = -qd(Bv[c]);
            }
            double J3[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) J3[r][c] = J[r][idx[c]];
            // 3x3 double inverse
            double det = J3[0][0] * (J3[1][1] * J3[2][2] - J3[1][2] * J3[2][1]) -
                         J3[0][1] * (J3[1][0] * J3[2][2] - J3[1][2] * J3[2][0]) +
                         J3[0][2] * (J3[1][0] * J3[2][1] - J3[1][1] * J3[2][0]);
            if (std::fabs(det) < 1e-300) { ok = false; break; }
            double inv[3][3] = {
                {(J3[1][1] * J3[2][2] - J3[1][2] * J3[2][1]) / det,
                 -(J3[1][0] * J3[2][2] - J3[1][2] * J3[2][0]) / det,
                 (J3[1][0] * J3[2][1] - J3[1][1] * J3[2][0]) / det},
                {-(J3[0][1] * J3[2][2] - J3[0][2] * J3[2][1]) / det,
                 (J3[0][0] * J3[2][2] - J3[0][2] * J3[2][0]) / det,
                 -(J3[0][0] * J3[2][1] - J3[0][1] * J3[2][0]) / det},
                {(J3[0][1] * J3[1][2] - J3[0][2] * J3[1][1]) / det,
                 -(J3[0][0] * J3[1][2] - J3[0][2] * J3[1][0]) / det,
                 (J3[0][0] * J3[1][1] - J3[0][1] * J3[1][0]) / det}};
            // 上面按 cofactor 行填写 => 需转置得到 C^T/det（真正的逆）
            for (int r = 0; r < 3; ++r)
                for (int c = r + 1; c < 3; ++c) std::swap(inv[r][c], inv[c][r]);
            double delta[3];
            for (int r = 0; r < 3; ++r) {
                double acc = 0;
                for (int c = 0; c < 3; ++c) acc += inv[r][c] * Fv[c];
                delta[r] = -acc;
            }
            if (std::getenv("CERT_DEBUG") && it == 0 && kfix == 0) {
                std::printf("[newton] p0=(%.17g,%.17g,%.17g,%.17g) F=(%.3g,%.3g,%.3g) delta=(%.3g,%.3g,%.3g)\\n",
                            p[0], p[1], p[2], p[3], Fv[0], Fv[1], Fv[2], delta[0], delta[1], delta[2]);
                std::printf("[newton] J3=[[%.6g,%.6g,%.6g],[%.6g,%.6g,%.6g],[%.6g,%.6g,%.6g]]\\n",
                            J3[0][0], J3[0][1], J3[0][2], J3[1][0], J3[1][1], J3[1][2], J3[2][0], J3[2][1], J3[2][2]);
            }
            for (int r = 0; r < 3; ++r) p[idx[r]] += delta[r];
            for (int c = 0; c < 4; ++c)
                if (p[c] < 0 || p[c] > 1) { if (std::getenv("CERT_DEBUG")) std::printf("[slice] kfix=%d out of domain at it=%d coord=%d val=%.17g\\n", kfix, it, c, p[c]); ok = false; break; }
            if (!ok) break;
        }
        if (!ok) { if (std::getenv("CERT_DEBUG")) std::printf("[slice] kfix=%d newton failed\\n", kfix); continue; }
        Q h = q_from_double(1e-9 * std::max(1.0, std::max({std::fabs(p[0]), std::fabs(p[1]), std::fabs(p[2]), std::fabs(p[3])})));
        I sb[4];
        bool good = true;
        for (int c = 0; c < 4; ++c) {
            Q lo = std::max(Q(q_from_double(p[c]) - h), Q(box[c][0]));
            Q hi = std::min(Q(q_from_double(p[c]) + h), Q(box[c][1]));
            if (lo > hi) { good = false; break; }
            sb[c] = I(lo, hi);
        }
        if (!good) { if (std::getenv("CERT_DEBUG")) std::printf("[slice] kfix=%d small box invalid\\n", kfix); continue; }
        bool kk = krawczyk_slice(A, B, dAu, dAv, dBu, dBv, sb[0], sb[1], sb[2], sb[3]);
        if (std::getenv("CERT_DEBUG"))
            std::printf("[slice] kfix=%d p=(%.17g,%.17g,%.17g,%.17g) h=%.3g kraw=%d\\n", kfix, p[0], p[1], p[2], p[3], qd(h), (int)kk);
        if (kk) return true;
    }
    return false;
}

// ------------------------------------------------------------ SHA-256
static const uint32_t K256[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
static void sha256(const std::string& msg, uint8_t out[32]) {
    std::vector<uint8_t> data(msg.begin(), msg.end());
    uint64_t bitlen = (uint64_t)data.size() * 8;
    data.push_back(0x80);
    while (data.size() % 64 != 56) data.push_back(0);
    for (int i = 7; i >= 0; --i) data.push_back((uint8_t)(bitlen >> (i * 8)));
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    for (size_t off = 0; off < data.size(); off += 64) {
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (uint32_t)data[off + 4 * i] << 24 | (uint32_t)data[off + 4 * i + 1] << 16 |
                   (uint32_t)data[off + 4 * i + 2] << 8 | data[off + 4 * i + 3];
        for (int i = 16; i < 64; ++i) {
            uint32_t s0 = (w[i - 15] >> 7 | w[i - 15] << 25) ^ (w[i - 15] >> 18 | w[i - 15] << 14) ^ (w[i - 15] >> 3);
            uint32_t s1 = (w[i - 2] >> 17 | w[i - 2] << 15) ^ (w[i - 2] >> 19 | w[i - 2] << 13) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            uint32_t S1 = (e >> 6 | e << 26) ^ (e >> 11 | e << 21) ^ (e >> 25 | e << 7);
            uint32_t ch = (e & f) ^ (~e & g);
            uint32_t t1 = hh + S1 + ch + K256[i] + w[i];
            uint32_t S0 = (a >> 2 | a << 30) ^ (a >> 13 | a << 19) ^ (a >> 22 | a << 10);
            uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t t2 = S0 + maj;
            hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }
    for (int i = 0; i < 8; ++i) {
        out[4 * i] = (uint8_t)(h[i] >> 24); out[4 * i + 1] = (uint8_t)(h[i] >> 16);
        out[4 * i + 2] = (uint8_t)(h[i] >> 8); out[4 * i + 3] = (uint8_t)h[i];
    }
}
static std::string sha256_hex(const std::string& s) {
    uint8_t d[32]; sha256(s, d);
    char buf[65];
    for (int i = 0; i < 32; ++i) std::snprintf(buf + 2 * i, 3, "%02x", d[i]);
    return std::string(buf, 64);
}

// canonical JSON serialization (Python json.dumps(sort_keys=True,separators=(",",":")))
static std::string esc(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '"') o += "\\\"";
        else if (c == '\\') o += "\\\\";
        else o += c;
    }
    o += "\"";
    return o;
}
static std::string canon(const J& j) {
    if (j.t == J::OBJ) {
        std::vector<std::pair<std::string, const J*>> kv;
        for (auto& p : j.o) kv.emplace_back(p.first, &p.second);
        std::sort(kv.begin(), kv.end());
        std::string o = "{";
        for (size_t i = 0; i < kv.size(); ++i) {
            if (i) o += ",";
            o += esc(kv[i].first) + ":" + canon(*kv[i].second);
        }
        return o + "}";
    }
    if (j.t == J::ARR) {
        std::string o = "[";
        for (size_t i = 0; i < j.a.size(); ++i) { if (i) o += ","; o += canon(j.a[i]); }
        return o + "]";
    }
    if (j.t == J::STR) return esc(j.s);
    if (j.t == J::BOOL) return j.b ? "true" : "false";
    if (j.t == J::NUL) return "null";
    double d = j.num;
    if (d == (long long)d) { char b[32]; std::snprintf(b, sizeof b, "%lld", (long long)d); return b; }
    char b[64]; std::snprintf(b, sizeof b, "%.17g", d); return b;
}

// ---------------------------------------------------------------- P7 Boehm
static std::vector<HP> insert_u(const std::vector<HP>& P, int nu, int nv, const std::vector<Q>& U,
                                int du, const Q& t, std::vector<Q>& Uout) {
    int k = du;
    for (int i = du; i < nu; ++i) if (U[i] <= t && t < U[i + 1]) { k = i; break; }
    std::vector<HP> Pout((size_t)(nu + 1) * nv);
    for (int j = 0; j < nv; ++j) {
        for (int i = 0; i <= k - du; ++i) Pout[(size_t)i * nv + j] = P[(size_t)i * nv + j];
        for (int i = k - du + 1; i <= k; ++i) {
            Q den = U[i + du] - U[i];
            Q alpha = (den != 0) ? (t - U[i]) / den : Q(0);
            for (int c = 0; c < 4; ++c)
                Pout[(size_t)i * nv + j][c] = alpha * P[(size_t)i * nv + j][c] +
                                              (1 - alpha) * P[(size_t)(i - 1) * nv + j][c];
        }
        for (int i = k + 1; i <= nu; ++i) Pout[(size_t)i * nv + j] = P[(size_t)(i - 1) * nv + j];
    }
    Uout = U; Uout.insert(Uout.begin() + k + 1, t);
    return Pout;
}
static std::vector<HP> insert_v(const std::vector<HP>& P, int nu, int nv, const std::vector<Q>& V,
                                int dv, const Q& t, std::vector<Q>& Vout) {
    int k = dv;
    for (int i = dv; i < nv; ++i) if (V[i] <= t && t < V[i + 1]) { k = i; break; }
    std::vector<HP> Pout((size_t)nu * (nv + 1));
    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j <= k - dv; ++j) Pout[(size_t)i * (nv + 1) + j] = P[(size_t)i * nv + j];
        for (int j = k - dv + 1; j <= k; ++j) {
            Q den = V[j + dv] - V[j];
            Q alpha = (den != 0) ? (t - V[j]) / den : Q(0);
            for (int c = 0; c < 4; ++c)
                Pout[(size_t)i * (nv + 1) + j][c] = alpha * P[(size_t)i * nv + j][c] +
                                                    (1 - alpha) * P[(size_t)i * nv + j - 1][c];
        }
        for (int j = k + 1; j <= nv; ++j) Pout[(size_t)i * (nv + 1) + j] = P[(size_t)i * nv + j - 1];
    }
    Vout = V; Vout.insert(Vout.begin() + k + 1, t);
    return Pout;
}
struct Blocks { int du, dv, nu, nv; std::vector<std::vector<HP>> blocks; };
static Blocks p7_convert(const J& nur) {
    int du = (int)nur.get("du")->num, dv = (int)nur.get("dv")->num;
    int nu = (int)nur.get("nu")->num, nv = (int)nur.get("nv")->num;
    std::vector<HP> P((size_t)nu * nv);
    size_t idx = 0;
    for (auto& row : nur.get("poles")->a) {
        HP hp;
        for (int c = 0; c < 4; ++c) hp[c] = parse_rat(row.a[c].s);
        P[idx++] = hp;
    }
    std::vector<Q> U, V;
    for (auto& x : nur.get("U")->a) U.push_back(parse_rat(x.s));
    for (auto& x : nur.get("V")->a) V.push_back(parse_rat(x.s));
    std::vector<Q> distinctU;
    for (std::size_t i = 1; i + 1 < U.size(); ++i)
        if (U[i] > U[0] && U[i] < U.back()) distinctU.push_back(U[i]);
    std::sort(distinctU.begin(), distinctU.end());
    distinctU.erase(std::unique(distinctU.begin(), distinctU.end()), distinctU.end());
    for (auto& t : distinctU) {
        int mult = 0;
        for (auto& x : U) if (x == t) ++mult;
        while (mult < du) {
            std::vector<Q> Un;
            P = insert_u(P, nu, nv, U, du, t, Un);
            U = Un; ++nu; ++mult;
        }
    }
    std::vector<Q> distinctV;
    for (std::size_t i = 1; i + 1 < V.size(); ++i)
        if (V[i] > V[0] && V[i] < V.back()) distinctV.push_back(V[i]);
    std::sort(distinctV.begin(), distinctV.end());
    distinctV.erase(std::unique(distinctV.begin(), distinctV.end()), distinctV.end());
    for (auto& t : distinctV) {
        int mult = 0;
        for (auto& x : V) if (x == t) ++mult;
        while (mult < dv) {
            std::vector<Q> Vn;
            P = insert_v(P, nu, nv, V, dv, t, Vn);
            V = Vn; ++nv; ++mult;
        }
    }
    std::vector<Q> uk, vk;
    for (std::size_t i = 0; i + 1 < U.size(); ++i) if (U[i] < U[i + 1]) uk.push_back(U[i]);
    uk.push_back(U.back());
    for (std::size_t i = 0; i + 1 < V.size(); ++i) if (V[i] < V[i + 1]) vk.push_back(V[i]);
    vk.push_back(V.back());
    Blocks out; out.du = du; out.dv = dv; out.nu = nu; out.nv = nv;
    for (std::size_t su = 0; su + 1 < uk.size(); ++su)
        for (std::size_t sv = 0; sv + 1 < vk.size(); ++sv) {
            std::vector<HP> blk((size_t)(du + 1) * (dv + 1));
            for (int i = 0; i <= du; ++i)
                for (int j = 0; j <= dv; ++j)
                    blk[(size_t)i * (dv + 1) + j] =
                        P[(size_t)(su * du + i) * nv + (sv * dv + j)];
            out.blocks.push_back(blk);
        }
    return out;
}
static bool patch_equals(const Patch& p, const std::vector<HP>& blk, int du, int dv) {
    if (p.nu != du + 1 || p.nv != dv + 1) return false;
    for (int i = 0; i <= du; ++i)
        for (int j = 0; j <= dv; ++j)
            for (int c = 0; c < 4; ++c)
                if (p.at(i, j)[c] != blk[(size_t)i * (dv + 1) + j][c]) return false;
    return true;
}

// ------------------------------------------------------------------ stats
struct Stats {
    long nodes = 0, excluded = 0, output = 0, cert_ok = 0, cert_fail = 0;
    long p7_ok = 0, p7_fail = 0, k4_ok = 0, k4_fail = 0;
    long coinc_ok = 0, tan_ok = 0, miss_ok = 0, cover_ok = 0;
    long branch_ok = 0, branch_fail = 0;
};

static std::array<std::array<Q, 2>, 4> parse_box(const J& b) {
    std::array<std::array<Q, 2>, 4> box;
    for (int k = 0; k < 4; ++k) {
        box[k][0] = parse_rat(b.a[2 * k].s);
        box[k][1] = parse_rat(b.a[2 * k + 1].s);
    }
    return box;
}

// --------------------------------------------------------- tree verify
struct Ctx {
    Patch A, B, dAu, dAv, dBu, dBv;
    Stats* st;
    std::vector<std::string>* errs;
};
static bool unit_box(const std::array<std::array<Q, 2>, 4>& b) {
    for (int k = 0; k < 4; ++k)
        if (b[k][0] != 0 || b[k][1] != 1) return false;
    return true;
}
static bool positive_weights(const Patch& p) {
    for (auto& hp : p.P) if (hp[3] <= 0) return false;
    return true;
}
static void verify_tree(Ctx& c, const J& node, const Patch& A, const Patch& B) {
    c.st->nodes++;
    if (node.get("kind") && node.get("kind")->s == "leaf") {
        const J* stj = node.get("status");
        if (!stj || (stj->s != "excluded" && stj->s != "output")) {
            c.errs->push_back("unknown leaf status");
            return;
        }
        if (stj->s == "excluded") {
            c.st->excluded++;
            if (!disjoint(hull_aabb(A), hull_aabb(B)))
                c.errs->push_back("excluded leaf is NOT disjoint");
            return;
        }
        c.st->output++;
        auto box = parse_box(*node.get("box"));
        I U = I(box[0][0], box[0][1]), V = I(box[1][0], box[1][1]);
        I S = I(box[2][0], box[2][1]), T = I(box[3][0], box[3][1]);
        auto cl = node.get("claimed_slice");
        auto c4 = node.get("claimed_4d");
        if (cl && cl->b) {
            if (slice_claim_ok(c.A, c.B, c.dAu, c.dAv, c.dBu, c.dBv, box))
                c.st->cert_ok++;
            else { c.st->cert_fail++; c.errs->push_back("claimed_slice fails"); }
        }
        if (c4 && c4->b) {
            if (krawczyk4d(c.A, c.B, c.dAu, c.dAv, c.dBu, c.dBv, U, V, S, T))
                c.st->k4_ok++;
            else { c.st->k4_fail++; c.errs->push_back("claimed_4d fails"); }
        }
        return;
    }
    int axis = (int)node.get("axis")->num;
    const J& ch = *node.get("children");
    auto pb = parse_box(*node.get("box"));
    for (int side = 0; side < 2; ++side) {
        auto cb = parse_box(*ch.a[side].get("box"));
        for (int k = 0; k < 4; ++k) {
            Q mid = (pb[k][0] + pb[k][1]) / 2;
            Q lo = pb[k][0], hi = pb[k][1];
            if (k == axis) {
                lo = (side == 0) ? pb[k][0] : mid;
                hi = (side == 0) ? mid : pb[k][1];
            }
            if (cb[k][0] != lo || cb[k][1] != hi) {
                c.errs->push_back("children do not partition parent");
                break;
            }
        }
    }
    if (axis == 0) {
        verify_tree(c, ch.a[0], split_u(A, false), B);
        verify_tree(c, ch.a[1], split_u(A, true), B);
    } else if (axis == 1) {
        verify_tree(c, ch.a[0], split_v(A, false), B);
        verify_tree(c, ch.a[1], split_v(A, true), B);
    } else if (axis == 2) {
        verify_tree(c, ch.a[0], A, split_u(B, false));
        verify_tree(c, ch.a[1], A, split_u(B, true));
    } else {
        verify_tree(c, ch.a[0], A, split_v(B, false));
        verify_tree(c, ch.a[1], A, split_v(B, true));
    }
}

// ---------------------------------------------------- coincidence/tangency
static bool is_bilinear(const Patch& p) { return p.nu == 2 && p.nv == 2; }
static void check_coincidence(const J& co, const Patch& A, const Patch& B, Stats& st,
                              std::vector<std::string>& errs) {
    if (!is_bilinear(A) || !is_bilinear(B)) { errs.push_back("coincidence needs bilinear patches"); return; }
    Q al = parse_rat(co.get("alpha")->s), be = parse_rat(co.get("beta")->s),
      ga = parse_rat(co.get("gamma")->s), de = parse_rat(co.get("delta")->s),
      ep = parse_rat(co.get("epsilon")->s), ze = parse_rat(co.get("zeta")->s);
    auto weighted = [](const HP& p) {
        HP h;
        const Q& w = p[3];
        h[0] = p[0] * w; h[1] = p[1] * w; h[2] = p[2] * w; h[3] = w;
        return h;
    };
    HP A0 = weighted(A.at(0, 0)), A10 = weighted(A.at(1, 0));
    HP A01 = weighted(A.at(0, 1)), A11 = weighted(A.at(1, 1));
    HP B0 = weighted(B.at(0, 0)), B10 = weighted(B.at(1, 0));
    HP B01 = weighted(B.at(0, 1)), B11 = weighted(B.at(1, 1));
    HP Au, Av, Auv, Bs, Bt, Bst;
    for (int c = 0; c < 4; ++c) {
        Au[c] = A10[c] - A0[c];
        Av[c] = A01[c] - A0[c];
        Auv[c] = A11[c] - A10[c] - A01[c] + A0[c];
        Bs[c] = B10[c] - B0[c];
        Bt[c] = B01[c] - B0[c];
        Bst[c] = B11[c] - B10[c] - B01[c] + B0[c];
    }
    HP c00, cu, cv, cuv;
    for (int c = 0; c < 4; ++c) {
        c00[c] = B0[c] + ga * Bs[c] + ze * Bt[c] + ga * ze * Bst[c];
        cu[c] = al * Bs[c] + de * Bt[c] + (al * ze + ga * de) * Bst[c];
        cv[c] = be * Bs[c] + ep * Bt[c] + (be * ze + ga * ep) * Bst[c];
        cuv[c] = (al * ep + be * de) * Bst[c];
    }
    bool ok = A0[3] != 0 && c00[3] != 0;
    if (ok) {
        Q lam = c00[3] / A0[3];
        ok = (lam != 0);
        for (int c = 0; ok && c < 4; ++c) {
            if (c00[c] != lam * A0[c]) ok = false;
            if (cu[c] != lam * Au[c]) ok = false;
            if (cv[c] != lam * Av[c]) ok = false;
            if (cuv[c] != lam * Auv[c]) ok = false;
        }
    }
    if (!ok) { errs.push_back("affine coincidence identity fails"); return; }
    long cnt = 0;
    for (auto& item : co.get("boxes")->a) {
        auto b = parse_box(*item.get("box"));
        bool good = true;
        for (int k = 0; k < 4; ++k) if (b[k][0] < 0 || b[k][1] > 1 || b[k][0] > b[k][1]) good = false;
        Q smin = std::min({al * b[0][0] + be * b[1][0] + ga, al * b[0][0] + be * b[1][1] + ga,
                           al * b[0][1] + be * b[1][0] + ga, al * b[0][1] + be * b[1][1] + ga});
        Q smax = std::max({al * b[0][0] + be * b[1][0] + ga, al * b[0][0] + be * b[1][1] + ga,
                           al * b[0][1] + be * b[1][0] + ga, al * b[0][1] + be * b[1][1] + ga});
        Q tmin = std::min({de * b[0][0] + ep * b[1][0] + ze, de * b[0][0] + ep * b[1][1] + ze,
                           de * b[0][1] + ep * b[1][0] + ze, de * b[0][1] + ep * b[1][1] + ze});
        Q tmax = std::max({de * b[0][0] + ep * b[1][0] + ze, de * b[0][0] + ep * b[1][1] + ze,
                           de * b[0][1] + ep * b[1][0] + ze, de * b[0][1] + ep * b[1][1] + ze});
        if (smin < b[2][0] || smax > b[2][1] || tmin < b[3][0] || tmax > b[3][1]) good = false;
        if (good) ++cnt;
        else errs.push_back("coincidence box mapping fails");
    }
    st.coinc_ok += cnt;
}
static bool cross_zero(const Q a[3], const Q b[3]) {
    Q c[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    return c[0] == 0 && c[1] == 0 && c[2] == 0;
}
static bool nonzero(const Q a[3]) { return a[0] != 0 || a[1] != 0 || a[2] != 0; }
static void check_tangency(const J& tp, const Patch& A, const Patch& B, const Patch& dAu,
                           const Patch& dAv, const Patch& dBu, const Patch& dBv, Stats& st,
                           std::vector<std::string>& errs) {
    long cnt = 0;
    for (auto& item : tp.a) {
        const J& p = *item.get("point");
        Q u = parse_rat(p.a[0].s), v = parse_rat(p.a[1].s), s = parse_rat(p.a[2].s), t = parse_rat(p.a[3].s);
        bool ok = u >= 0 && u <= 1 && v >= 0 && v <= 1 && s >= 0 && s <= 1 && t >= 0 && t <= 1;
        auto a = point_eval(A, u, v);
        auto b = point_eval(B, s, t);
        for (int c = 0; c < 3; ++c) if (a[c] != b[c]) ok = false;
        Q nA[3], nB[3];
        {
            auto au = deriv_point(A, dAu, u, v, true), av = deriv_point(A, dAv, u, v, false);
            auto bu = deriv_point(B, dBu, s, t, true), bv = deriv_point(B, dBv, s, t, false);
            nA[0] = au[1] * av[2] - au[2] * av[1];
            nA[1] = au[2] * av[0] - au[0] * av[2];
            nA[2] = au[0] * av[1] - au[1] * av[0];
            nB[0] = bu[1] * bv[2] - bu[2] * bv[1];
            nB[1] = bu[2] * bv[0] - bu[0] * bv[2];
            nB[2] = bu[0] * bv[1] - bu[1] * bv[0];
        }
        if (!nonzero(nA) || !nonzero(nB) || !cross_zero(nA, nB)) ok = false;
        if (ok) ++cnt;
        else errs.push_back("tangent point certificate fails");
    }
    st.tan_ok += cnt;
}
// --------------------------------------------------- polyline / distances
struct Poly { std::vector<std::array<Q, 3>> pts; };
static Poly parse_poly(const J& arr) {
    Poly p;
    for (auto& q : arr.a) {
        std::array<Q, 3> pt;
        for (int c = 0; c < 3; ++c) pt[c] = parse_rat(q.a[c].s);
        p.pts.push_back(pt);
    }
    return p;
}
static Q poly_dist2(const std::array<Q, 3>& p, const std::vector<Poly>& polys) {
    Q best = -1;
    for (auto& pl : polys)
        for (size_t i = 0; i + 1 < pl.pts.size(); ++i) {
            Q d2 = seg_dist2(p.data(), pl.pts[i].data(), pl.pts[i + 1].data());
            if (best < 0 || d2 < best) best = d2;
        }
    return best;
}
static void check_miss(const J& cert, const Patch& A, const Patch& B, Stats& st,
                       std::vector<std::string>& errs) {
    Q delta2 = parse_rat(cert.get("delta2")->s);
    std::vector<Poly> polys;
    for (auto& pl : cert.get("occt_polyline")->a) polys.push_back(parse_poly(pl));
    long cnt = 0;
    for (auto& item : cert.get("solution_points")->a) {
        auto& pt = *item.get("point");
        Q u = parse_rat(pt.a[0].s), v = parse_rat(pt.a[1].s), s = parse_rat(pt.a[2].s), t = parse_rat(pt.a[3].s);
        bool ok = (u >= 0 && u <= 1 && v >= 0 && v <= 1 && s >= 0 && s <= 1 && t >= 0 && t <= 1);
        auto a = point_eval(A, u, v);
        auto b = point_eval(B, s, t);
        for (int c = 0; c < 3; ++c) if (a[c] != b[c]) ok = false;
        Q d2 = poly_dist2(a, polys);
        if (!(d2 > delta2)) ok = false;
        if (ok) ++cnt;
        else errs.push_back("miss certificate fails");
    }
    st.miss_ok += cnt;
}
static void check_cover(const J& cert, const Patch& A, const Patch& B, Stats& st,
                        std::vector<std::string>& errs) {
    Q eps2 = parse_rat(cert.get("eps2")->s);
    std::vector<Poly> polys;
    for (auto& pl : cert.get("occt_polyline")->a) polys.push_back(parse_poly(pl));
    long cnt = 0;
    for (auto& item : cert.get("solution_points")->a) {
        auto& pt = *item.get("point");
        Q u = parse_rat(pt.a[0].s), v = parse_rat(pt.a[1].s), s = parse_rat(pt.a[2].s), t = parse_rat(pt.a[3].s);
        bool ok = (u >= 0 && u <= 1 && v >= 0 && v <= 1 && s >= 0 && s <= 1 && t >= 0 && t <= 1);
        auto a = point_eval(A, u, v);
        auto b = point_eval(B, s, t);
        for (int c = 0; c < 3; ++c) if (a[c] != b[c]) ok = false;
        Q d2 = poly_dist2(a, polys);
        if (d2 > eps2) ok = false;
        if (ok) ++cnt;
        else errs.push_back("coverage certificate fails");
    }
    st.cover_ok += cnt;
}

// ---------------------------------------------------- branch coverage
static double seg_dist2_d(const Q p[3], const Q a[3], const Q b[3]) {
    double ax = qd(a[0]), ay = qd(a[1]), az = qd(a[2]);
    double dx = qd(b[0]) - ax, dy = qd(b[1]) - ay, dz = qd(b[2]) - az;
    double px = qd(p[0]) - ax, py = qd(p[1]) - ay, pz = qd(p[2]) - az;
    double L2 = dx * dx + dy * dy + dz * dz;
    if (L2 == 0) return px * px + py * py + pz * pz;
    double dot = px * dx + py * dy + pz * dz;
    if (dot <= 0) return px * px + py * py + pz * pz;
    if (dot >= L2) {
        double qx = qd(p[0]) - qd(b[0]), qy = qd(p[1]) - qd(b[1]), qz = qd(p[2]) - qd(b[2]);
        return qx * qx + qy * qy + qz * qz;
    }
    return px * px + py * py + pz * pz - dot * dot / L2;
}
static Q emax_corner(const std::array<I, 3>& e, const std::array<Q, 3>& a, const std::array<Q, 3>& b) {
    Q best = -1;
    for (int sx = 0; sx < 2; ++sx)
        for (int sy = 0; sy < 2; ++sy)
            for (int sz = 0; sz < 2; ++sz) {
                Q p[3] = {sx ? e[0].hi : e[0].lo, sy ? e[1].hi : e[1].lo, sz ? e[2].hi : e[2].lo};
                Q d2 = seg_dist2(p, a.data(), b.data());
                if (best < 0 || d2 > best) best = d2;
            }
    return best;
}
static void check_branch(const J& cert, const Patch& A, Stats& st, std::vector<std::string>& errs) {
    Q eps2 = parse_rat(cert.get("eps2")->s);
    std::vector<Poly> polys;
    for (auto& pl : cert.get("occt_polyline")->a) polys.push_back(parse_poly(pl));
    std::vector<std::pair<std::array<Q, 3>, std::array<Q, 3>>> segs;
    for (auto& pl : polys)
        for (size_t i = 0; i + 1 < pl.pts.size(); ++i) segs.emplace_back(pl.pts[i], pl.pts[i + 1]);
    std::function<void(const J&)> walk = [&](const J& n) {
        if (n.get("kind") && n.get("kind")->s == "leaf") {
            if (n.get("status")->s != "output") return;
            auto box = parse_box(*n.get("box"));
            I U(box[0][0], box[0][1]), V(box[1][0], box[1][1]);
            std::array<I, 3> enc;
            try { enc = cart_iv(A, U, V); } catch (...) { st.branch_fail++; errs.push_back("branch enclosure undefined"); return; }
            double bestv = 1e300; int besti = -1;
            for (size_t i = 0; i < segs.size(); ++i) {
                Q pp[3] = {enc[0].lo, enc[1].lo, enc[2].lo};
                double worst = 0;
                for (int sx = 0; sx < 2; ++sx)
                    for (int sy = 0; sy < 2; ++sy)
                        for (int sz = 0; sz < 2; ++sz) {
                            Q q[3] = {sx ? enc[0].hi : enc[0].lo, sy ? enc[1].hi : enc[1].lo,
                                      sz ? enc[2].hi : enc[2].lo};
                            worst = std::max(worst, seg_dist2_d(q, segs[i].first.data(), segs[i].second.data()));
                        }
                if (worst < bestv) { bestv = worst; besti = (int)i; }
            }
            if (besti < 0) { st.branch_fail++; errs.push_back("branch: empty polyline"); return; }
            Q d2 = emax_corner(enc, segs[besti].first, segs[besti].second);
            if (d2 <= eps2) st.branch_ok++;
            else { st.branch_fail++; errs.push_back("branch coverage fails in a leaf"); }
            return;
        }
        for (auto& ch : n.get("children")->a) walk(ch);
    };
    walk(*cert.get("tree"));
}

// ------------------------------------------------------------- verify one
static void check_claims(const J& cert, std::vector<std::string>& errs) {
    if (!cert.has("claims")) return;
    const J& cl = *cert.get("claims");
    if (cl.t != J::OBJ) { errs.push_back("claims must be an object"); return; }
    auto fmt = cl.get("format");
    if (fmt && fmt->t == J::STR && fmt->s != "polyline3d-v1")
        errs.push_back("unsupported claims format");
    auto dig = cl.get("output_digest");
    if (!dig || dig->t != J::STR || dig->s.size() != 64) {
        errs.push_back("output_digest must be 64 lowercase hex characters");
        return;
    }
    for (char c : dig->s)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
            errs.push_back("output_digest must be 64 lowercase hex characters");
            return;
        }
    auto poly = cert.get("occt_polyline");
    if (!poly) { errs.push_back("claims present but no occt_polyline"); return; }
    std::string actual = sha256_hex(canon(*poly));
    if (actual != dig->s) errs.push_back("output digest mismatch");
}

static bool verify_cert(const J& cert, Stats& st, std::vector<std::string>& errs) {
    check_claims(cert, errs);
    const J& patches = *cert.get("patches");
    Patch A = load_patch(patches.a[0]);
    Patch B = load_patch(patches.a[1]);
    Patch dAu = deriv(A, 0), dAv = deriv(A, 1), dBu = deriv(B, 0), dBv = deriv(B, 1);
    if (!positive_weights(A) || !positive_weights(B))
        errs.push_back("patch weights must be positive");
    if (cert.has("nurbs") && cert.get("nurbs")->t == J::ARR) {
        const J& nurbs = *cert.get("nurbs");
        for (size_t k = 0; k < nurbs.a.size(); ++k) {
            Blocks blk = p7_convert(nurbs.a[k]);
            const Patch& claimed = (k == 0) ? A : B;
            bool match = false;
            for (auto& b : blk.blocks) if (patch_equals(claimed, b, blk.du, blk.dv)) { match = true; break; }
            if (match) st.p7_ok++;
            else { st.p7_fail++; errs.push_back("P7 mismatch"); }
        }
    }
    if (cert.has("tree")) {
        const J& root = *cert.get("tree");
        if (root.get("box") && !unit_box(parse_box(*root.get("box"))))
            errs.push_back("tree root is not [0,1]^4");
        Ctx c{A, B, dAu, dAv, dBu, dBv, &st, &errs};
        verify_tree(c, root, A, B);
    }
    if (cert.has("certified")) {
        for (auto& item : cert.get("certified")->a) {
            auto box = parse_box(*item.get("box"));
            I U(box[0][0], box[0][1]), V(box[1][0], box[1][1]);
            I S(box[2][0], box[2][1]), T(box[3][0], box[3][1]);
            if (krawczyk_slice(A, B, dAu, dAv, dBu, dBv, U, V, S, T)) st.cert_ok++;
            else { st.cert_fail++; errs.push_back("certified box fails"); }
        }
    }
    if (cert.has("coincidence")) check_coincidence(*cert.get("coincidence"), A, B, st, errs);
    if (cert.has("tangent_points"))
        check_tangency(*cert.get("tangent_points"), A, B, dAu, dAv, dBu, dBv, st, errs);
    if (cert.has("occt_polyline")) {
        if (cert.has("solution_points")) {
            if (cert.has("delta2")) check_miss(cert, A, B, st, errs);
            else check_cover(cert, A, B, st, errs);
        } else if (cert.has("kind") && cert.get("kind")->s == "branch_coverage") {
            check_branch(cert, A, st, errs);
        }
    }
    if (cert.has("chunks")) errs.push_back("chunked certificate: pass --chunked");
    return errs.empty();
}

// ------------------------------------------------------------ chunked mode
static std::vector<std::pair<int, int>> parse_path(const J& p) {
    std::vector<std::pair<int, int>> out;
    for (auto& e : p.a) out.emplace_back((int)e.a[0].num, (int)e.a[1].num);
    return out;
}
static std::string pkey(const std::vector<std::pair<int, int>>& p) {
    std::string s;
    for (auto& pr : p) { s += std::to_string(pr.first); s += ":"; s += std::to_string(pr.second); s += ","; }
    return s;
}
struct ChunkTask {
    const J* ch;
    std::array<std::array<Q, 2>, 4> box;
    Patch A, B;
};
static void verify_chunked_cert(const J& cert, Stats& st, std::vector<std::string>& errs, int jobs) {
    check_claims(cert, errs);
    Patch A = load_patch((*cert.get("patches")).a[0]);
    Patch B = load_patch((*cert.get("patches")).a[1]);
    Patch dAu = deriv(A, 0), dAv = deriv(A, 1), dBu = deriv(B, 0), dBv = deriv(B, 1);
    const J& chunks = *cert.get("chunks");
    std::map<std::string, const J*> cmap;
    std::map<std::string, std::vector<std::pair<int, int>>> pmap;
    for (auto& ch : chunks.a) {
        auto path = parse_path(*ch.get("path"));
        for (auto& pr : path)
            if (pr.first < 0 || pr.first > 3 || pr.second < 0 || pr.second > 1)
                errs.push_back("invalid chunk path");
        std::string key = pkey(path);
        cmap[key] = &ch;
        pmap[key] = path;
    }
    if (cmap.size() != chunks.a.size()) errs.push_back("duplicate chunk path");
    std::vector<ChunkTask> tasks;
    std::vector<std::string> perr;
    std::function<void(std::vector<std::pair<int, int>>, std::array<std::array<Q, 2>, 4>, Patch, Patch)> walk;
    walk = [&](std::vector<std::pair<int, int>> path, std::array<std::array<Q, 2>, 4> box, Patch Pa, Patch Pb) {
        std::string key = pkey(path);
        auto it = cmap.find(key);
        if (it != cmap.end()) {
            auto rb = parse_box(*it->second->get("root"));
            for (int k = 0; k < 4; ++k)
                if (rb[k][0] != box[k][0] || rb[k][1] != box[k][1])
                    perr.push_back("chunk root does not match frontier box");
            tasks.push_back({it->second, box, Pa, Pb});
            return;
        }
        int axis = -1;
        for (auto& kv : pmap)
            if (kv.second.size() > path.size() &&
                std::equal(path.begin(), path.end(), kv.second.begin())) { axis = kv.second[path.size()].first; break; }
        if (axis < 0 || axis > 3) { perr.push_back("frontier hole or invalid axis"); return; }
        for (int side = 0; side < 2; ++side) {
            auto nb = box;
            Q mid = (nb[axis][0] + nb[axis][1]) / 2;
            if (side == 0) nb[axis][1] = mid; else nb[axis][0] = mid;
            auto np = path;
            np.emplace_back(axis, side);
            Patch na = Pa, nb2 = Pb;
            if (axis == 0) na = split_u(Pa, side == 1);
            else if (axis == 1) na = split_v(Pa, side == 1);
            else if (axis == 2) nb2 = split_u(Pb, side == 1);
            else nb2 = split_v(Pb, side == 1);
            walk(np, nb, na, nb2);
        }
    };
    std::array<std::array<Q, 2>, 4> root;
    for (int k = 0; k < 4; ++k) { root[k][0] = 0; root[k][1] = 1; }
    walk({}, root, A, B);
    if (tasks.size() != chunks.a.size()) perr.push_back("unused chunks");
    errs.insert(errs.end(), perr.begin(), perr.end());

    auto run_task = [&](const ChunkTask& t, Stats& out, std::vector<std::string>& e) {
        std::string h = sha256_hex(canon(*t.ch->get("tree")));
        if (h != t.ch->get("sha256")->s) e.push_back("chunk sha256 mismatch");
        Ctx c{A, B, dAu, dAv, dBu, dBv, &out, &e};
        verify_tree(c, *t.ch->get("tree"), t.A, t.B);
    };
    if (jobs <= 1) {
        for (auto& t : tasks) run_task(t, st, errs);
    } else {
        std::vector<Stats> sts(tasks.size());
        std::vector<std::vector<std::string>> errss(tasks.size());
        std::atomic<size_t> idx{0};
        std::vector<std::thread> ths;
        for (int j = 0; j < jobs; ++j)
            ths.emplace_back([&]() {
                while (true) {
                    size_t i = idx.fetch_add(1);
                    if (i >= tasks.size()) break;
                    run_task(tasks[i], sts[i], errss[i]);
                }
            });
        for (auto& th : ths) th.join();
        for (size_t i = 0; i < tasks.size(); ++i) {
            st.nodes += sts[i].nodes; st.excluded += sts[i].excluded; st.output += sts[i].output;
            st.cert_ok += sts[i].cert_ok; st.cert_fail += sts[i].cert_fail;
            st.cert_ok += sts[i].cert_ok;
            st.p7_ok += sts[i].p7_ok; st.p7_fail += sts[i].p7_fail;
            st.k4_ok += sts[i].k4_ok; st.k4_fail += sts[i].k4_fail;
            st.coinc_ok += sts[i].coinc_ok; st.tan_ok += sts[i].tan_ok;
            st.miss_ok += sts[i].miss_ok; st.cover_ok += sts[i].cover_ok;
            st.branch_ok += sts[i].branch_ok; st.branch_fail += sts[i].branch_fail;
            errs.insert(errs.end(), errss[i].begin(), errss[i].end());
        }
    }
    if (cert.has("certified")) {
        for (auto& item : cert.get("certified")->a) {
            auto box = parse_box(*item.get("box"));
            I U(box[0][0], box[0][1]), V(box[1][0], box[1][1]);
            I S(box[2][0], box[2][1]), T(box[3][0], box[3][1]);
            if (krawczyk_slice(A, B, dAu, dAv, dBu, dBv, U, V, S, T)) st.cert_ok++;
            else { st.cert_fail++; errs.push_back("certified box fails"); }
        }
    }
    if (cert.has("nurbs") && cert.get("nurbs")->t == J::ARR) {
        for (size_t k = 0; k < cert.get("nurbs")->a.size(); ++k) {
            Blocks blk = p7_convert(cert.get("nurbs")->a[k]);
            const Patch& claimed = (k == 0) ? A : B;
            bool match = false;
            for (auto& b : blk.blocks) if (patch_equals(claimed, b, blk.du, blk.dv)) { match = true; break; }
            if (match) st.p7_ok++;
            else { st.p7_fail++; errs.push_back("P7 mismatch (chunked)"); }
        }
    }
}
// -------------------------------------------------------------------- main
static void print_stats(const Stats& s, const std::string& tag) {
    std::printf("%s nodes=%ld excluded=%ld output=%ld cert_ok=%ld cert_fail=%ld p7_ok=%ld p7_fail=%ld "
                "k4_ok=%ld k4_fail=%ld coinc_ok=%ld tan_ok=%ld miss_ok=%ld cover_ok=%ld "
                "branch_ok=%ld branch_fail=%ld\n",
                tag.c_str(), s.nodes, s.excluded, s.output, s.cert_ok, s.cert_fail, s.p7_ok, s.p7_fail,
                s.k4_ok, s.k4_fail, s.coinc_ok, s.tan_ok, s.miss_ok, s.cover_ok, s.branch_ok,
                s.branch_fail);
}
static double now_ms() {
    using namespace std::chrono;
    return duration<double, std::milli>(steady_clock::now().time_since_epoch()).count();
}
static bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}
static void list_json(const std::string& dir, std::vector<std::string>& out) {
    std::string cmd = "ls -1 " + dir + "/*.json 2>/dev/null";
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) return;
    char buf[4096];
    while (std::fgets(buf, sizeof buf, p)) {
        std::string s(buf);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        if (ends_with(s, ".json")) out.push_back(s);
    }
    pclose(p);
}
int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: %s <cert.json> | --batch <dir>\n", argv[0]);
        return 2;
    }
    std::string mode = argv[1];
    if (mode == "--batch" && argc >= 3) {
        std::vector<std::string> files;
        list_json(argv[2], files);
        long total = 0, accept = 0;
        double t0 = now_ms();
        Stats agg;
        for (auto& f : files) {
            Stats st;
            std::vector<std::string> errs;
            try {
                J cert = parse_file(f);
                if (cert.has("chunks")) verify_chunked_cert(cert, st, errs, 1);
                else verify_cert(cert, st, errs);
            } catch (std::exception& e) {
                errs.push_back(std::string("exception: ") + e.what());
            }
            ++total;
            if (errs.empty()) ++accept;
            agg.nodes += st.nodes; agg.excluded += st.excluded; agg.output += st.output;
            agg.cert_ok += st.cert_ok; agg.cert_fail += st.cert_fail; agg.p7_ok += st.p7_ok;
            agg.k4_ok += st.k4_ok; agg.k4_fail += st.k4_fail; agg.coinc_ok += st.coinc_ok;
            agg.tan_ok += st.tan_ok; agg.miss_ok += st.miss_ok; agg.cover_ok += st.cover_ok;
            agg.branch_ok += st.branch_ok; agg.branch_fail += st.branch_fail;
        }
        double dt = now_ms() - t0;
        std::printf("BATCH files=%ld accept=%ld reject=%ld time=%.1f ms (%.3f ms/file)\n", total, accept,
                    total - accept, dt, dt / std::max(1L, total));
        print_stats(agg, "TOTAL");
        return accept == total ? 0 : 1;
    }
    const char* path = argv[1];
    bool chunked = false;
    int jobs = 1;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--chunked") chunked = true;
        else if (a == "--jobs" && i + 1 < argc) jobs = std::atoi(argv[++i]);
    }
    double t0 = now_ms();
    Stats st;
    std::vector<std::string> errs;
    try {
        J cert = parse_file(path);
        if (cert.has("chunks")) chunked = true;
        if (chunked) verify_chunked_cert(cert, st, errs, jobs);
        else verify_cert(cert, st, errs);
    } catch (std::exception& e) {
        errs.push_back(std::string("exception: ") + e.what());
    }
    double dt = now_ms() - t0;
    print_stats(st, "RESULT");
    for (size_t i = 0; i < errs.size() && i < 5; ++i) std::printf("  ERROR: %s\n", errs[i].c_str());
    std::printf("time=%.2f ms\n", dt);
    std::printf("FINAL: %s\n", errs.empty() ? "ACCEPT" : "REJECT");
    return errs.empty() ? 0 : 1;
}








