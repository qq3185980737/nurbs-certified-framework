// Certified column for the shared SSI fixtures A-E.
//
// This is the *same* interval-arithmetic exclusion + subdivision cover as
// cpp/tests/bench_ssi.cpp (the routine is copied verbatim from there, only the
// fixture polynomials are the ones shared with the SISL/OCCT baselines), so the
// numbers are directly comparable with the in-repo certified SSI experiment.
// Output per fixture: excluded boxes / cover boxes / subdivisions / wall time.
// tol is the parameter-domain box size (bench_ssi uses 1e-3).
#include <cstdio>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <chrono>

struct Ival { double lo, hi; Ival(double l = 0, double h = 0) : lo(l), hi(h) {} };
inline Ival iadd(Ival a, Ival b) {
  return Ival(std::nextafter(a.lo + b.lo, -std::numeric_limits<double>::infinity()),
              std::nextafter(a.hi + b.hi,  std::numeric_limits<double>::infinity()));
}
inline Ival imul(Ival a, Ival b) {
  double x1 = a.lo * b.lo, x2 = a.lo * b.hi, x3 = a.hi * b.lo, x4 = a.hi * b.hi;
  return Ival(std::nextafter(std::min(std::min(x1, x2), std::min(x3, x4)),
                             -std::numeric_limits<double>::infinity()),
              std::nextafter(std::max(std::max(x1, x2), std::max(x3, x4)),
                              std::numeric_limits<double>::infinity()));
}
inline Ival ipow(Ival x, int n) { Ival r(1, 1); for (int i = 0; i < n; ++i) r = imul(r, x); return r; }

struct Poly2 {
  std::vector<std::vector<double>> a;  // a[i][j] coefficient of s^i t^j
  int maxI() const { return (int)a.size() - 1; }
  double eval(double s, double t) const {
    double r = 0;
    for (int i = 0; i <= maxI(); ++i) {
      double ri = 0;
      for (int j = 0; j < (int)a[i].size(); ++j) ri += a[i][j] * std::pow(t, j);
      r += ri * std::pow(s, i);
    }
    return r;
  }
  Ival ieval(Ival s, Ival t) const {  // rigorous interval evaluation (Horner)
    Ival r(0, 0), si(1, 1);
    for (int i = 0; i <= maxI(); ++i) {
      Ival ri(0, 0), tj(1, 1);
      for (int j = 0; j < (int)a[i].size(); ++j) {
        ri = iadd(ri, imul(Ival(a[i][j], a[i][j]), tj));
        tj = imul(tj, t);
      }
      r = iadd(r, imul(ri, si));
      si = imul(si, s);
    }
    return r;
  }
};

struct Stats { long excl = 0, cover = 0, subdiv = 0; };

static void cover(const Poly2& f, double sl, double sh, double tl, double th,
                  double tol, Stats& st) {
  Ival F = f.ieval(Ival(sl, sh), Ival(tl, th));
  if (!(F.lo <= 0.0 && F.hi >= 0.0)) { ++st.excl; return; }
  if ((sh - sl) <= tol && (th - tl) <= tol) { ++st.cover; return; }
  ++st.subdiv;
  if ((sh - sl) >= (th - tl)) {
    const double m = 0.5 * (sl + sh);
    cover(f, sl, m, tl, th, tol, st); cover(f, m, sh, tl, th, tol, st);
  } else {
    const double m = 0.5 * (tl + th);
    cover(f, sl, sh, tl, m, tol, st); cover(f, sl, sh, m, th, tol, st);
  }
}

int main() {
  const double tol = 1e-3;
  struct Case { const char* name; Poly2 f; const char* truth; };
  std::vector<Case> cases = {
      {"A", {{{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},}, "1 curve: diagonal s=t"},
      {"B", {{{0.3}},}, "0 intersections"},
      {"C", {{{0.0, 0.5}, {1.0, 0.0}},}, "1 point: corner (0,0) tangency"},
      {"D", {{{0, 0, 1}, {0, -2, 0}, {1, 0, 0}},}, "1 tangential curve s=t"},
      {"E", {{{-0.1, 0, 0, 1}, {0}, {0}, {1}},}, "1 interior curve: u^3+v^3 = 0.1"},
  };
  printf("# certified interval cover for the shared SSI fixtures, parameter tol=%.0e\n",
         tol);
  printf("# method: cpp/tests/bench_ssi.cpp (interval exclusion + subdivision)\n");
  printf("fixture,excluded_boxes,cover_boxes,subdivisions,time_ms,truth\n");
  for (auto& c : cases) {
    // E: g = u^3 + v^3 - 0.1?  The baseline fixture is u^3 = 0.1 - v^3, i.e.
    // g = u^3 + v^3 - 0.1, with a[0][0] = -0.1.
    Poly2 f = c.f;
    Stats st;
    auto t0 = std::chrono::steady_clock::now();
    cover(f, 0.0, 1.0, 0.0, 1.0, tol, st);
    double dt = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - t0).count();
    printf("%s,%ld,%ld,%ld,%.2f,\"%s\"\n", c.name, st.excl, st.cover, st.subdiv,
           dt, c.truth);
    fflush(stdout);
  }
  return 0;
}
