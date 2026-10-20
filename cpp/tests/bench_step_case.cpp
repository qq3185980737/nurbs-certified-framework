// ============================================================================
// bench_step_case.cpp — 真实 STEP 零件上的认证曲线-曲面求交
//   读取 step_case/case_data.txt：Rational Bezier 曲线段 + Rational Bezier 曲面片
//   + OCCT 参考交点；对每个 (边,面) 对枚举参数盒并做 3x3 区间 Krawczyk 认证。
//   自包含（仅标准库），不依赖主库。
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <numeric>

struct IV {
    double lo, hi;
    IV(double l = 0, double h = 0) : lo(l), hi(h) {}
    double mid() const { return 0.5 * (lo + hi); }
    double wid() const { return hi - lo; }
};
static inline IV I(double l, double h) { return IV(l, h); }
static inline IV add(const IV& a, const IV& b) { return IV(std::nextafter(a.lo + b.lo, -1e300), std::nextafter(a.hi + b.hi, 1e300)); }
static inline IV sub(const IV& a, const IV& b) { return IV(std::nextafter(a.lo - b.hi, -1e300), std::nextafter(a.hi - b.lo, 1e300)); }
static inline IV mul(const IV& a, const IV& b) {
    double p[4] = { a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi };
    double mn = p[0], mx = p[0];
    for (int i = 1; i < 4; ++i) { mn = std::min(mn, p[i]); mx = std::max(mx, p[i]); }
    return IV(std::nextafter(mn, -1e300), std::nextafter(mx, 1e300));
}
static inline IV divi(const IV& a, const IV& b, bool* ok) {
    if (b.lo <= 0.0 && b.hi >= 0.0) { if (ok) *ok = false; return IV(0, 0); }
    double p[4] = { a.lo / b.lo, a.lo / b.hi, a.hi / b.lo, a.hi / b.hi };
    double mn = p[0], mx = p[0];
    for (int i = 1; i < 4; ++i) { mn = std::min(mn, p[i]); mx = std::max(mx, p[i]); }
    return IV(std::nextafter(mn, -1e300), std::nextafter(mx, 1e300));
}
static inline bool contains0(const IV& a) { return a.lo <= 0.0 && a.hi >= 0.0; }

using H4 = std::array<double, 4>;   // homogeneous (wx, wy, wz, w)

// ---- Rational Bezier 曲线 ----
struct RCurve {
    int deg = 0;
    std::vector<H4> P;
    H4 heval_iv(const IV& t) const {
        std::vector<std::array<IV, 4>> c(P.size());
        for (size_t i = 0; i < P.size(); ++i)
            for (int k = 0; k < 4; ++k) c[i][k] = I(P[i][k], P[i][k]);
        for (size_t r = 1; r < c.size(); ++r)
            for (size_t i = 0; i + r < c.size(); ++i)
                for (int k = 0; k < 4; ++k)
                    c[i][k] = add(mul(c[i][k], I(1.0 - t.hi, 1.0 - t.lo)), mul(c[i + 1][k], t));
        H4 out{};
        for (int k = 0; k < 4; ++k) { out[k] = c[0][k].lo; }
        // 需要区间：返回中点用于直接计算，区间用别处实现
        return out;
    }
    // 区间求值：返回区间四元组
    std::array<IV, 4> heval(const IV& t) const {
        std::vector<std::array<IV, 4>> c(P.size());
        for (size_t i = 0; i < P.size(); ++i)
            for (int k = 0; k < 4; ++k) c[i][k] = I(P[i][k], P[i][k]);
        IV one(1.0, 1.0);
        for (size_t r = 1; r < c.size(); ++r)
            for (size_t i = 0; i + r < c.size(); ++i)
                for (int k = 0; k < 4; ++k)
                    c[i][k] = add(mul(c[i][k], sub(one, t)), mul(c[i + 1][k], t));
        return c[0];
    }
    H4 heval_pt(double t) const {
        if (P.size() < 2) return P.empty() ? H4{0,0,0,1} : P[0];
        std::vector<H4> c = P;
        for (size_t r = 1; r < c.size(); ++r)
            for (size_t i = 0; i + r < c.size(); ++i)
                for (int k = 0; k < 4; ++k) c[i][k] = (1 - t) * c[i][k] + t * c[i + 1][k];
        return c[0];
    }
    std::array<IV, 4> hderiv(const IV& t) const {   // 齐次导数控制网
        int p = deg;
        std::vector<std::array<IV, 4>> c(p);
        for (int i = 0; i < p; ++i)
            for (int k = 0; k < 4; ++k)
                c[i][k] = mul(I((double)p, (double)p), sub(I(P[i + 1][k], P[i + 1][k]), I(P[i][k], P[i][k])));
        IV one(1.0, 1.0);
        for (size_t r = 1; r < c.size(); ++r)
            for (size_t i = 0; i + r < c.size(); ++i)
                for (int k = 0; k < 4; ++k)
                    c[i][k] = add(mul(c[i][k], sub(one, t)), mul(c[i + 1][k], t));
        return c[0];
    }
};

// ---- Rational Bezier 曲面（行优先 [i*nv+j]）----
struct RSurf {
    int du = 0, dv = 0, nu = 0, nv = 0;
    std::vector<H4> P;
    std::array<IV, 4> heval(const IV& s, const IV& t) const {
        // 先沿 v 方向 de Casteljau，再沿 u
        std::vector<std::array<IV, 4>> col(nu);
        IV one(1.0, 1.0);
        for (int i = 0; i < nu; ++i) {
            std::vector<std::array<IV, 4>> row(nv);
            for (int j = 0; j < nv; ++j) {
                const H4& p = P[i * nv + j];
                for (int k = 0; k < 4; ++k) row[j][k] = I(p[k], p[k]);
            }
            for (int r2 = 1; r2 < nv; ++r2)
                for (int j = 0; j + r2 < nv; ++j)
                    for (int k = 0; k < 4; ++k)
                        row[j][k] = add(mul(row[j][k], sub(one, t)), mul(row[j + 1][k], t));
            col[i] = row[0];
        }
        for (int r2 = 1; r2 < nu; ++r2)
            for (int i = 0; i + r2 < nu; ++i)
                for (int k = 0; k < 4; ++k)
                    col[i][k] = add(mul(col[i][k], sub(one, s)), mul(col[i + 1][k], s));
        return col[0];
    }
    H4 heval_pt(double s, double t) const {
        if (nu < 2 || nv < 2) return H4{0,0,0,1};
        std::vector<H4> col(nu);
        for (int i = 0; i < nu; ++i) {
            std::vector<H4> row(nv);
            for (int j = 0; j < nv; ++j) row[j] = P[i * nv + j];
            for (int r2 = 1; r2 < nv; ++r2)
                for (int j = 0; j + r2 < nv; ++j)
                    for (int k = 0; k < 4; ++k) row[j][k] = (1 - t) * row[j][k] + t * row[j + 1][k];
            col[i] = row[0];
        }
        for (int r2 = 1; r2 < nu; ++r2)
            for (int i = 0; i + r2 < nu; ++i)
                for (int k = 0; k < 4; ++k) col[i][k] = (1 - s) * col[i][k] + s * col[i + 1][k];
        return col[0];
    }
    std::array<IV, 4> hderiv_s(const IV& s, const IV& t) const {
        if (nu <= 1) { std::array<IV, 4> z{}; for (int k = 0; k < 4; ++k) z[k] = I(0, 0); return z; }
        RSurf d = *this;
        d.nu = nu - 1; d.du = du - 1;
        d.P.clear();
        for (int i = 0; i < nu - 1; ++i)
            for (int j = 0; j < nv; ++j) {
                H4 h{};
                for (int k = 0; k < 4; ++k) h[k] = (double)du * (P[(i + 1) * nv + j][k] - P[i * nv + j][k]);
                d.P.push_back(h);
            }
        return d.heval(s, t);
    }
    std::array<IV, 4> hderiv_t(const IV& s, const IV& t) const {
        if (nv <= 1) { std::array<IV, 4> z{}; for (int k = 0; k < 4; ++k) z[k] = I(0, 0); return z; }
        RSurf d = *this;
        d.nv = nv - 1; d.dv = dv - 1;
        d.P.clear();
        for (int i = 0; i < nu; ++i)
            for (int j = 0; j < nv - 1; ++j) {
                H4 h{};
                for (int k = 0; k < 4; ++k) h[k] = (double)dv * (P[i * nv + j + 1][k] - P[i * nv + j][k]);
                d.P.push_back(h);
            }
        return d.heval(s, t);
    }
};

// ---- 有理除法：C=H/W, C'=(H' - C W')/W ----
static inline void rat_point(const std::array<IV, 4>& h, double out[3], bool* ok) {
    IV w = h[3];
    if (contains0(w)) { if (ok) *ok = false; out[0] = out[1] = out[2] = 0; return; }
    for (int i = 0; i < 3; ++i) out[i] = divi(h[i], w, ok).mid();
}
static inline void rat_iv(const std::array<IV, 4>& h, IV out[3], bool* ok) {
    IV w = h[3];
    for (int i = 0; i < 3; ++i) out[i] = divi(h[i], w, ok);
}
static inline void rat_deriv_iv(const std::array<IV, 4>& h, const std::array<IV, 4>& dh,
                                const IV C[3], IV out[3], bool* ok) {
    for (int i = 0; i < 3; ++i) {
        IV num = sub(dh[i], mul(C[i], dh[3]));
        out[i] = divi(num, h[3], ok);
    }
}

struct Box3 { IV u, s, t; };

// 3x3 Krawczyk
static bool krawczyk(const RCurve& C, const RSurf& S, const Box3& X, double* outPt) {
    IV Cc[3]; bool ok = true;
    auto hc = C.heval(X.u);
    rat_iv(hc, Cc, &ok); if (!ok) return false;
    IV Sc[3]; auto hs = S.heval(X.s, X.t);
    rat_iv(hs, Sc, &ok); if (!ok) return false;
    IV F[3];
    for (int i = 0; i < 3; ++i) F[i] = sub(Cc[i], Sc[i]);
    for (int i = 0; i < 3; ++i) if (!contains0(F[i])) return false;   // 排除

    // 区间 Jacobian
    IV dC[3]; auto dhc = C.hderiv(X.u);
    rat_deriv_iv(hc, dhc, Cc, dC, &ok); if (!ok) return false;
    IV dS_s[3], dS_t[3];
    auto hs_s = S.hderiv_s(X.s, X.t);
    auto hs_t = S.hderiv_t(X.s, X.t);
    rat_deriv_iv(hs, hs_s, Sc, dS_s, &ok); if (!ok) return false;
    rat_deriv_iv(hs, hs_t, Sc, dS_t, &ok); if (!ok) return false;
    IV J[3][3];
    for (int i = 0; i < 3; ++i) { J[i][0] = dC[i]; J[i][1] = IV(-dS_s[i].hi, -dS_s[i].lo); J[i][2] = IV(-dS_t[i].hi, -dS_t[i].lo); }

    // 中点浮点求值
    double mu = X.u.mid(), ms = X.s.mid(), mt = X.t.mid();
    H4 hc0 = C.heval_pt(mu); H4 hs0 = S.heval_pt(ms, mt);
    if (std::fabs(hc0[3]) < 1e-14 || std::fabs(hs0[3]) < 1e-14) return false;
    double Cp[3], Sp[3];
    for (int i = 0; i < 3; ++i) { Cp[i] = hc0[i] / hc0[3]; Sp[i] = hs0[i] / hs0[3]; }
    double Fm[3] = { Cp[0] - Sp[0], Cp[1] - Sp[1], Cp[2] - Sp[2] };
    // 浮点导数
    double dCp[3], dSp_s[3], dSp_t[3];
    {
        H4 d = C.heval_pt(mu); // 已用
        // 曲线导数：解析控制点差分
        int p = C.deg;
        std::vector<H4> Q(p);
        for (int i = 0; i < p; ++i) for (int k = 0; k < 4; ++k) Q[i][k] = (double)p * (C.P[i + 1][k] - C.P[i][k]);
        H4 hd{};
        { std::vector<H4> c = Q; for (size_t r = 1; r < c.size(); ++r) for (size_t i = 0; i + r < c.size(); ++i) for (int k = 0; k < 4; ++k) c[i][k] = (1 - mu) * c[i][k] + mu * c[i + 1][k]; hd = c.empty() ? H4{0,0,0,0} : c[0]; }
        double w = hc0[3];
        for (int i = 0; i < 3; ++i) dCp[i] = (hd[i] - Cp[i] * hd[3]) / w;
        (void)d;
    }
    {
        // 曲面导数：数值中心差分（足够用于 Y 中点的近似逆）
        double h = 1e-7;
        H4 a = S.heval_pt(std::min(1.0, ms + h), mt), b = S.heval_pt(std::max(0.0, ms - h), mt);
        double A[3], B[3];
        for (int i = 0; i < 3; ++i) { A[i] = a[i] / a[3]; B[i] = b[i] / b[3]; }
        for (int i = 0; i < 3; ++i) dSp_s[i] = (A[i] - B[i]) / (std::min(1.0, ms + h) - std::max(0.0, ms - h));
        H4 c = S.heval_pt(ms, std::min(1.0, mt + h)), d2 = S.heval_pt(ms, std::max(0.0, mt - h));
        double Cc2[3], D2[3];
        for (int i = 0; i < 3; ++i) { Cc2[i] = c[i] / c[3]; D2[i] = d2[i] / d2[3]; }
        for (int i = 0; i < 3; ++i) dSp_t[i] = (Cc2[i] - D2[i]) / (std::min(1.0, mt + h) - std::max(0.0, mt - h));
    }
    double Jm[3][3] = { { dCp[0], -dSp_s[0], -dSp_t[0] },
                        { dCp[1], -dSp_s[1], -dSp_t[1] },
                        { dCp[2], -dSp_s[2], -dSp_t[2] } };
    // 3x3 逆（Gauss-Jordan）
    double M[3][6];
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) M[i][j] = Jm[i][j]; for (int j = 0; j < 3; ++j) M[i][3 + j] = (i == j) ? 1.0 : 0.0; }
    for (int col = 0; col < 3; ++col) {
        int piv = col; double best = std::fabs(M[col][col]);
        for (int r2 = col + 1; r2 < 3; ++r2) if (std::fabs(M[r2][col]) > best) { best = std::fabs(M[r2][col]); piv = r2; }
        if (best < 1e-14) return false;
        if (piv != col) for (int j = 0; j < 6; ++j) std::swap(M[piv][j], M[col][j]);
        double d = M[col][col];
        for (int j = 0; j < 6; ++j) M[col][j] /= d;
        for (int r2 = 0; r2 < 3; ++r2) if (r2 != col) { double f = M[r2][col]; for (int j = 0; j < 6; ++j) M[r2][j] -= f * M[col][j]; }
    }
    double Y[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Y[i][j] = M[i][3 + j];

    // K = m - Y*Fm + (I - Y*J(X))*(X - m)
    Box3 mbox = { I(mu, mu), I(ms, ms), I(mt, mt) };
    IV v[3] = { sub(X.u, mbox.u), sub(X.s, mbox.s), sub(X.t, mbox.t) };
    IV K[3];
    for (int i = 0; i < 3; ++i) {
        double yf = 0.0;
        for (int j = 0; j < 3; ++j) yf += Y[i][j] * Fm[j];
        IV acc = sub(mbox.u, I(yf, yf));   // 后面按坐标加
        IV sum = IV(0, 0);
        for (int j = 0; j < 3; ++j) {
            // (I - Y*J)[i][j]
            IV yj = IV(0, 0);
            for (int l = 0; l < 3; ++l) yj = add(yj, mul(I(Y[i][l], Y[i][l]), J[l][j]));
            IV mij = sub(I((i == j) ? 1.0 : 0.0, (i == j) ? 1.0 : 0.0), yj);
            sum = add(sum, mul(mij, v[j]));
        }
        K[i] = add(acc, sum);
    }
    // 接受判定：K ⊂ int(X)，三坐标分别
    IV Xc[3] = { X.u, X.s, X.t };
    for (int i = 0; i < 3; ++i) {
        if (!(Xc[i].lo < K[i].lo && K[i].hi < Xc[i].hi)) return false;
    }
    // 输出点：C(u) 浮点
    H4 hp = C.heval_pt(mu);
    for (int i = 0; i < 3; ++i) outPt[i] = hp[i] / hp[3];
    return true;
}

// 枚举
static void enumerate(const RCurve& C, const RSurf& S, Box3 X, int depth,
                      std::vector<std::array<double, 3>>& roots, long& undecided) {
    double pts[3];
    if (depth >= 26) { ++undecided; return; }
    // 若 Krawczyk 接受
    if (krawczyk(C, S, X, pts)) {
        roots.push_back({ pts[0], pts[1], pts[2] });
        return;
    }
    double wmax = std::max(X.u.wid(), std::max(X.s.wid(), X.t.wid()));
    if (wmax <= 1e-4) { ++undecided; return; }
    // 二分最长边
    if (X.u.wid() >= X.s.wid() && X.u.wid() >= X.t.wid()) {
        double m = X.u.mid();
        enumerate(C, S, Box3{ I(X.u.lo, m), X.s, X.t }, depth + 1, roots, undecided);
        enumerate(C, S, Box3{ I(m, X.u.hi), X.s, X.t }, depth + 1, roots, undecided);
    } else if (X.s.wid() >= X.t.wid()) {
        double m = X.s.mid();
        enumerate(C, S, Box3{ X.u, I(X.s.lo, m), X.t }, depth + 1, roots, undecided);
        enumerate(C, S, Box3{ X.u, I(m, X.s.hi), X.t }, depth + 1, roots, undecided);
    } else {
        double m = X.t.mid();
        enumerate(C, S, Box3{ X.u, X.s, I(X.t.lo, m) }, depth + 1, roots, undecided);
        enumerate(C, S, Box3{ X.u, X.s, I(m, X.t.hi) }, depth + 1, roots, undecided);
    }
}

int main(int argc, char** argv) {
    const char* path = argc > 1 ? argv[1] : "step_case/case_data.txt";
    std::ifstream f(path);
    if (!f) { std::printf("无法打开 %s\n", path); return 1; }
    std::vector<std::vector<RCurve>> curves;
    std::vector<std::vector<RSurf>> surfs;
    struct Pair { int e, s; std::vector<std::array<double, 3>> pts; };
    std::vector<Pair> pairs;
    std::string tok;
    long ntok = 0;
    while (f >> tok) {
        if (++ntok % 200 == 0) { std::printf("  parse %ld tok=%s\n", ntok, tok.c_str()); std::fflush(stdout); }
        if (tok == "NSURF") { int n; f >> n; surfs.resize(n); }
        else if (tok == "NCURVE") { int n; f >> n; curves.resize(n); }
        else if (tok == "SURF") {
            int fid, pi, du, dv, nu, nv; f >> fid >> pi >> du >> dv >> nu >> nv;
            RSurf s; s.du = du; s.dv = dv; s.nu = nu; s.nv = nv;
            s.P.resize((size_t)nu * nv);
            for (int i = 0; i < nu; ++i) for (int j = 0; j < nv; ++j) {
                double x, y, z, w; f >> x >> y >> z >> w;
                s.P[(size_t)i * nv + j] = H4{ x * w, y * w, z * w, w };
            }
            surfs[fid].push_back(s);
            std::string end; f >> end;
        }
        else if (tok == "CURVE") {
            int eid, si, deg, n; f >> eid >> si >> deg >> n;
            RCurve c; c.deg = deg; c.P.resize(n);
            for (int i = 0; i < n; ++i) {
                double x, y, z, w; f >> x >> y >> z >> w;
                c.P[i] = H4{ x * w, y * w, z * w, w };
            }
            curves[eid].push_back(c);
            std::string end; f >> end;
        }
        else if (tok == "NPAIR") { int n; f >> n; pairs.resize(n); }
        else if (tok == "PAIR") {
            static int pi = 0;
            int e, s, m; f >> e >> s >> m;
            pairs[pi].e = e; pairs[pi].s = s;
            for (int i = 0; i < m; ++i) { double x, y, z; f >> x >> y >> z; pairs[pi].pts.push_back({ x, y, z }); }
            std::string end; f >> end;
            ++pi;
        }
    }
    std::printf("读取: 曲线组=%zu 曲面组=%zu 边-面对=%zu 期望交点=%zu\n",
                curves.size(), surfs.size(), pairs.size(),
                (size_t)std::accumulate(pairs.begin(), pairs.end(), 0,
                    [](int a, const Pair& p) { return a + (int)p.pts.size(); }));
    auto t0 = std::chrono::high_resolution_clock::now();
    int matched = 0, extra = 0, missing = 0;
    long undecided = 0;
    int pidx = 0;
    for (auto& p : pairs) {
        std::printf("[%d/%zu] pair edge=%d face=%d segs=%zu patches=%zu\n", ++pidx, pairs.size(), p.e, p.s, (size_t)curves[p.e].size(), (size_t)surfs[p.s].size());
        std::fflush(stdout);
        std::vector<std::array<double, 3>> roots;
        for (auto& c : curves[p.e]) for (auto& s : surfs[p.s])
            enumerate(c, s, Box3{ I(0, 1), I(0, 1), I(0, 1) }, 0, roots, undecided);
        // 去重
        std::vector<std::array<double, 3>> uniq;
        for (auto& r : roots) {
            bool dup = false;
            for (auto& q : uniq) {
                double d = std::sqrt((r[0] - q[0]) * (r[0] - q[0]) + (r[1] - q[1]) * (r[1] - q[1]) + (r[2] - q[2]) * (r[2] - q[2]));
                if (d < 1e-6) { dup = true; break; }
            }
            if (!dup) uniq.push_back(r);
        }
        // 匹配
        std::vector<bool> usedExp(p.pts.size(), false), usedGot(uniq.size(), false);
        for (size_t i = 0; i < p.pts.size(); ++i) {
            for (size_t j = 0; j < uniq.size(); ++j) {
                if (usedGot[j]) continue;
                double d = std::sqrt((p.pts[i][0] - uniq[j][0]) * (p.pts[i][0] - uniq[j][0]) +
                                     (p.pts[i][1] - uniq[j][1]) * (p.pts[i][1] - uniq[j][1]) +
                                     (p.pts[i][2] - uniq[j][2]) * (p.pts[i][2] - uniq[j][2]));
                if (d < 1e-4) { usedExp[i] = true; usedGot[j] = true; ++matched; break; }
            }
        }
        for (size_t i = 0; i < p.pts.size(); ++i) if (!usedExp[i]) ++missing;
        for (size_t j = 0; j < uniq.size(); ++j) if (!usedGot[j]) ++extra;
        std::printf("  pair (edge %d, face %d): 期望 %zu, 认证 %zu\n", p.e, p.s, p.pts.size(), uniq.size());
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::printf("\n汇总: 边-面对 %zu | 期望交点 %d | 认证匹配 %d | 漏检 %d | 多检 %d | 未定盒 %ld | 耗时 %.1f ms\n",
                pairs.size(), matched + missing, matched, missing, extra, undecided, ms);
    return (missing == 0 && extra == 0) ? 0 : 2;
}
