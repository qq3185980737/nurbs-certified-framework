// SISL failure-mode probes (companion to the A-F baseline table):
//  (1) multi-root curve/curve: quintic Bezier curve C(t) = (t, 30*prod(t-k/6))
//      with 5 transverse crossings of the x-axis, intersected with the x-axis
//      via s1857.  Reports the number of roots found, the parameter values and
//      the residual |C(t).y|.
//  (2) surface self-intersection: SISL exposes no surface self-intersection
//      operator; s1859(S,S) on the same surface is probed to document what the
//      coincidence code path returns.
//  (3) tolerance behaviour on the degenerate tangential fixture D:
//      s1859/s1310 is called with epsge = 1e-6, 1e-9 and 1e-12; the traced 3D
//      curve is densely sampled and the max distance between the curve and the
//      two surfaces is reported.  If that distance stays at the epsge=1e-6
//      level while a stricter tolerance is requested, the solver does not honour
//      its own tolerance (i.e. it returns geometry without a tolerance
//      certificate).
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "sisl.h"
using std::vector;

static vector<double> bezier_knots(int k) {
  vector<double> et(2 * k);
  for (int i = 0; i < k; ++i) { et[i] = 0.0; et[k + i] = 1.0; }
  return et;
}

static double binom(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  double r = 1.0;
  for (int i = 0; i < k; ++i) r = r * (n - i) / (i + 1);
  return r;
}

static vector<double> mono1d_to_bezier(const vector<double>& a, int n) {
  vector<double> c(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = 0; k <= i && k < (int)a.size(); ++k)
      s += a[k] * binom(i, k) / binom(n - 1, k);
    c[i] = s;
  }
  return c;
}

static SISLCurve* make_bezier3d(const vector<double>& ax,
                                const vector<double>& ay, int n) {
  vector<double> et = bezier_knots(n);
  vector<double> cx = mono1d_to_bezier(ax, n);
  vector<double> cy = mono1d_to_bezier(ay, n);
  vector<double> coef((size_t)n * 3, 0.0);
  for (int i = 0; i < n; ++i) {
    coef[3 * i] = cx[i];
    coef[3 * i + 1] = cy[i];
    coef[3 * i + 2] = 0.0;
  }
  return newCurve(n, n, et.data(), coef.data(), 3, 3, 1);
}

// monomial (u,v) -> Bernstein net (n x n)
static vector<vector<double>> mono2d_to_bezier(const vector<vector<double>>& a,
                                               int n) {
  vector<vector<double>> z(n, vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double s = 0.0;
      for (int k = 0; k <= i; ++k)
        for (int l = 0; l <= j; ++l) {
          double ak = (k < (int)a.size() && l < (int)a[k].size()) ? a[k][l] : 0.0;
          s += ak * binom(i, k) / binom(n - 1, k) * binom(j, l) / binom(n - 1, l);
        }
      z[i][j] = s;
    }
  return z;
}

static SISLSurf* make_graph(const vector<vector<double>>& mono, int n) {
  vector<vector<double>> z = mono2d_to_bezier(mono, n);
  vector<double> et = bezier_knots(n);
  vector<double> coef((size_t)n * n * 3, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      size_t idx = ((size_t)j * n + i) * 3;
      coef[idx] = double(i) / (n - 1);
      coef[idx + 1] = double(j) / (n - 1);
      coef[idx + 2] = z[i][j];
    }
  return newSurf(n, n, n, n, et.data(), et.data(), coef.data(), 3, 3, 1);
}

static void eval_curve(SISLCurve* c, double t, double o[3]) {
  int lk = 0, st = 0; double d[3];
  s1221(c, 0, t, &lk, d, &st);
  o[0] = d[0]; o[1] = d[1]; o[2] = d[2];
}
static void eval_surf(SISLSurf* s, double u, double v, double o[3]) {
  int l1 = 0, l2 = 0, st = 0; double par[2] = {u, v}, d[3];
  s1421(s, 0, par, &l1, &l2, d, nullptr, &st);
  o[0] = d[0]; o[1] = d[1]; o[2] = d[2];
}

// ---------------------------------------------------------------- (1) multi-root
static void probe_multirroot() {
  const double roots[5] = {1.0 / 6, 2.0 / 6, 3.0 / 6, 4.0 / 6, 5.0 / 6};
  // g(t) = 30 * prod_k (t - roots[k]) ; build monomial coefficients by convolution
  vector<double> g(1, 30.0);
  for (int k = 0; k < 5; ++k) {
    vector<double> next(g.size() + 1, 0.0);
    for (size_t i = 0; i < g.size(); ++i) {
      next[i] += -roots[k] * g[i];
      next[i + 1] += g[i];
    }
    g = next;
  }
  vector<double> ax(6, 0.0);  // x(t) = t
  ax[1] = 1.0;
  SISLCurve* c = make_bezier3d(ax, g, 6);
  SISLCurve* line = make_bezier3d({0.0, 1.0}, {0.0, 0.0}, 2);  // x-axis segment

  int npt = 0, ncr = 0, st = 0;
  double *p1 = nullptr, *p2 = nullptr;
  SISLIntcurve** ic = nullptr;
  const double epsge = 1e-6;
  s1857(c, line, 1e-15, epsge, &npt, &p1, &p2, &ncr, &ic, &st);
  printf("[multi-root] s1857: quintic C(t)=(t, 30*prod_{k=1..5}(t-k/6)) vs x-axis\n");
  printf("  epsge=1e-6  stat=%d  num_points=%d  num_curves=%d   (ground truth: 5 simple roots)\n",
         st, npt, ncr);
  double maxres = 0.0;
  for (int i = 0; i < npt; ++i) {
    double t = p1[i];  // curve/curve: one parameter value per intersection point
    double P[3];
    eval_curve(c, t, P);
    double res = std::fabs(P[1]);
    if (res > maxres) maxres = res;
    printf("  root[%2d]  t=%.9f  C(t).y=%+.3e  P=(%.6f,%.6f,%.6f)\n", i, t, P[1],
           P[0], P[1], P[2]);
  }
  printf("  max |C(t).y| over reported roots = %.3e\n", maxres);
  printf("  -> SISL returns the geometric intersections but nothing about root\n"
         "     multiplicity/existence: no certified count.\n");
  if (p1) free(p1);
  if (p2) free(p2);
  if (ncr) freeIntcrvlist(ic, ncr);
  freeCurve(c); freeCurve(line);
}

// ------------------------------------------------------- (2) surface self-intersect
static void probe_selfint() {
  vector<vector<double>> a = {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}};  // u^2 - v^2
  SISLSurf* s = make_graph(a, 3);
  int npt = 0, ncr = 0, st = 0;
  double *p1 = nullptr, *p2 = nullptr;
  SISLIntcurve** ic = nullptr;
  s1859(s, s, 1e-15, 1e-6, &npt, &p1, &p2, &ncr, &ic, &st);
  printf("[self-intersection] s1859(S,S) on S=(u,v,u^2-v^2): stat=%d npt=%d ncr=%d\n",
         st, npt, ncr);
  printf("  SISL has no surface self-intersection operator; the coincidence code\n"
         "  path returns the boundary of the parameter square as 'curves':\n");
  for (int i = 0; i < ncr; ++i)
    printf("    curve[%d] ipoint=%d itype=%d start=(%.3f,%.3f) end=(%.3f,%.3f)\n", i,
           ic[i]->ipoint, ic[i]->itype, ic[i]->epar1[0], ic[i]->epar1[1],
           ic[i]->epar1[2 * (ic[i]->ipoint - 1)],
           ic[i]->epar1[2 * (ic[i]->ipoint - 1) + 1]);
  if (p1) free(p1);
  if (p2) free(p2);
  if (ncr) freeIntcrvlist(ic, ncr);
  freeSurf(s);
}

// --------------------------------------------------- (3) D tolerance behaviour
static void probe_D_tolerance(double epsge) {
  vector<vector<double>> a = {{0, 0, 1}, {0, -2, 0}, {1, 0, 0}};  // (u-v)^2
  SISLSurf* s1 = make_graph(a, 3);
  vector<vector<double>> z0 = {{0}};
  SISLSurf* s2 = make_graph(z0, 2);
  int npt = 0, ncr = 0, st = 0;
  double *p1 = nullptr, *p2 = nullptr;
  SISLIntcurve** ic = nullptr;
  s1859(s1, s2, 1e-15, epsge, &npt, &p1, &p2, &ncr, &ic, &st);
  printf("[D tolerance] epsge=%.0e : s1859 stat=%d npt=%d ncr=%d\n", epsge, st, npt, ncr);
  int traced = 0, tstat_max = 0;
  double max_gap = 0.0, max_param_dev = 0.0, max_f = 0.0;
  double wgP[3] = {0,0,0}, wgA[3] = {0,0,0}, wgB[3] = {0,0,0};
  double wdP[3] = {0,0,0};
  long nsamp = 0;
  for (int i = 0; i < ncr; ++i) {
    int tst = 0;
    s1310(s1, s2, ic[i], epsge, 0.0, 2, 0, &tst);
    if (tst > tstat_max) tstat_max = tst;
    if (!ic[i]->pgeom) continue;
    ++traced;
    SISLCurve* pg = ic[i]->pgeom;
    double lo = pg->et[pg->ik - 1], hi = pg->et[pg->in];
    const int N = 400;
    for (int k = 0; k <= N; ++k) {
      double t = lo + (hi - lo) * double(k) / N;
      double P[3], A[3], B[3];
      eval_curve(pg, t, P);
      // graph patches: (x,y) of P are the (u,v) parameters of both patches
      eval_surf(s1, P[0], P[1], A);
      eval_surf(s2, P[0], P[1], B);
      double f = std::fabs((P[0] - P[1]) * (P[0] - P[1]));  // |f1-f2| at (x,y)
      if (f > max_f) max_f = f;
      double d1 = std::sqrt((P[0]-A[0])*(P[0]-A[0]) + (P[1]-A[1])*(P[1]-A[1]) +
                            (P[2]-A[2])*(P[2]-A[2]));
      double d2 = std::sqrt((P[0]-B[0])*(P[0]-B[0]) + (P[1]-B[1])*(P[1]-B[1]) +
                            (P[2]-B[2])*(P[2]-B[2]));
      double gap = d1 > d2 ? d1 : d2;
      if (gap > max_gap) {
        max_gap = gap;
        wgP[0] = P[0]; wgP[1] = P[1]; wgP[2] = P[2];
        wgA[0] = A[0]; wgA[1] = A[1]; wgA[2] = A[2];
        wgB[0] = B[0]; wgB[1] = B[1]; wgB[2] = B[2];
      }
      double dev = std::fabs(P[0] - P[1]);
      if (dev > max_param_dev) {
        max_param_dev = dev;
        wdP[0] = P[0]; wdP[1] = P[1]; wdP[2] = P[2];
      }
      ++nsamp;
    }
  }
  printf("  traced_curves=%d trace_stat_max=%d samples=%ld\n", traced, tstat_max, nsamp);
  printf("  max dist(curve, S1 or S2) = %.3e   (requested epsge = %.0e) -> %s\n",
         max_gap, epsge, max_gap <= epsge ? "WITHIN tolerance" : "EXCEEDS tolerance");
  printf("  max |f1-f2| at curve points = %.3e   max |u-v| = %.3e\n", max_f, max_param_dev);
  printf("  worst-gap sample: P=(%.6f,%.6f,%.6f) S1=(%.6f,%.6f,%.6f) S2=(%.6f,%.6f,%.6f)\n",
         wgP[0], wgP[1], wgP[2], wgA[0], wgA[1], wgA[2], wgB[0], wgB[1], wgB[2]);
  printf("  worst-dev sample: P=(%.6f,%.6f,%.6f)\n", wdP[0], wdP[1], wdP[2]);
  if (p1) free(p1);
  if (p2) free(p2);
  if (ncr) freeIntcrvlist(ic, ncr);
  freeSurf(s1); freeSurf(s2);
}

int main(int argc, char** argv) {
  std::string mode = (argc > 1) ? argv[1] : "all";
  if (mode == "multi") {
    probe_multirroot();
  } else if (mode == "self") {
    probe_selfint();
  } else if (mode == "D") {
    double epsge = (argc > 2) ? atof(argv[2]) : 1e-6;
    probe_D_tolerance(epsge);
  } else {
    probe_multirroot();
    printf("\n");
    probe_selfint();
    printf("\n");
    probe_D_tolerance(1e-6);
    probe_D_tolerance(1e-4);
  }
  return 0;
}
