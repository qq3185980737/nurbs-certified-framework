// SISL 4.7 surface-surface / surface-curve intersection baseline harness
// Fixtures A-F, same as the certified framework experiments (graph patches
// S1=(s,t,f1), S2=(s,t,f2); intersection is the zero set of f1-f2).
//
// API used (verified against SISL 4.7 doc/manual/func/*.tex & app/surf_surf_int.C):
//   newSurf / newCurve     - construct Bezier patch / curve
//   s1421                  - evaluate surface (position)
//   s1221                  - evaluate curve (position)
//   s1859                  - surface/surface intersection topology
//   s1310                  - march intersection curve (3D + param-plane curves)
//   s1858                  - surface/curve intersection topology
//
// Residual definitions:
//   res_pt    : max over isolated intersection points of |S1(p1)-S2(p2)|
//   res_guide : max over guide points of intersection curves   |S1(p1)-S2(p2)|
//   res_trace : max over traced-curve control points of the distance between
//               the 3D control point and S1/S2 evaluated at the corresponding
//               parameter-plane control point (consistency of s1310 output)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>

#include "sisl.h"

using std::vector;
using std::string;

static double now_ms() {
  using namespace std::chrono;
  return duration_cast<duration<double, std::milli>>(
             steady_clock::now().time_since_epoch()).count();
}

static double binom(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  double r = 1.0;
  for (int i = 0; i < k; ++i) r = r * (n - i) / (i + 1);
  return r;
}

static vector<vector<double>> mono_to_bezier(
    const vector<vector<double>>& a, int n1, int n2) {
  vector<vector<double>> z(n1, vector<double>(n2, 0.0));
  for (int i = 0; i < n1; ++i)
    for (int j = 0; j < n2; ++j) {
      double s = 0.0;
      for (int k = 0; k <= i; ++k)
        for (int l = 0; l <= j; ++l) {
          double ak = (k < (int)a.size() && l < (int)a[k].size()) ? a[k][l] : 0.0;
          s += ak * binom(i, k) / binom(n1 - 1, k) * binom(j, l) / binom(n2 - 1, l);
        }
      z[i][j] = s;
    }
  return z;
}

static vector<double> bezier_knots(int k) {
  vector<double> et(2 * k);
  for (int i = 0; i < k; ++i) { et[i] = 0.0; et[k + i] = 1.0; }
  return et;
}

static SISLSurf* make_graph_surface(const vector<vector<double>>& z, int n1,
                                   int n2) {
  vector<double> et1 = bezier_knots(n1), et2 = bezier_knots(n2);
  vector<double> coef((size_t)n1 * n2 * 3, 0.0);
  for (int i1 = 0; i1 < n1; ++i1)
    for (int i2 = 0; i2 < n2; ++i2) {
      size_t idx = ((size_t)i2 * n1 + i1) * 3;  // SISL layout (i2*in1 + i1)
      coef[idx + 0] = (n1 > 1) ? double(i1) / (n1 - 1) : 0.0;
      coef[idx + 1] = (n2 > 1) ? double(i2) / (n2 - 1) : 0.0;
      coef[idx + 2] = z[i1][i2];
    }
  return newSurf(n1, n2, n1, n2, et1.data(), et2.data(), coef.data(), 3, 3, 1);
}

static SISLSurf* make_graph_mono(const vector<vector<double>>& a, int n1,
                                 int n2) {
  return make_graph_surface(mono_to_bezier(a, n1, n2), n1, n2);
}

// exact rational quarter circle/cylinder: x=cos(theta), y=sin(theta), z=t
static SISLSurf* make_quarter_cylinder() {
  double wx[3] = {1.0, std::sqrt(2.0) / 2.0, 1.0};
  double px[3] = {1.0, 1.0, 0.0};
  double py[3] = {0.0, 1.0, 1.0};
  vector<double> et1 = bezier_knots(3), et2 = bezier_knots(2);
  vector<double> coef((size_t)3 * 2 * 4, 0.0);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 2; ++j) {
      size_t idx = ((size_t)j * 3 + i) * 4;
      coef[idx + 0] = wx[i] * px[i];
      coef[idx + 1] = wx[i] * py[i];
      coef[idx + 2] = wx[i] * double(j);
      coef[idx + 3] = wx[i];
    }
  return newSurf(3, 2, 3, 2, et1.data(), et2.data(), coef.data(), 4, 3, 1);
}

static SISLCurve* make_line(const double p0[3], const double p1[3]) {
  vector<double> et = {0.0, 0.0, 1.0, 1.0};
  vector<double> coef(6);
  for (int k = 0; k < 3; ++k) { coef[k] = p0[k]; coef[3 + k] = p1[k]; }
  return newCurve(2, 2, et.data(), coef.data(), 1, 3, 1);
}

static void eval_surf(SISLSurf* s, double u, double v, double out[3]) {
  int l1 = 0, l2 = 0, stat = 0;
  double par[2] = {u, v}, der[3];
  s1421(s, 0, par, &l1, &l2, der, nullptr, &stat);
  out[0] = der[0]; out[1] = der[1]; out[2] = der[2];
}

static void eval_curve3(SISLCurve* c, double t, double out[3]) {
  int lk = 0, stat = 0;
  double der[3];
  s1221(c, 0, t, &lk, der, &stat);
  out[0] = der[0]; out[1] = der[1]; out[2] = der[2];
}

static double dist3(const double a[3], const double b[3]) {
  return std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) +
                   (a[2]-b[2])*(a[2]-b[2]));
}

// valid parameter interval of a SISL curve: [et[ik-1], et[in]]
static double curve_param_lo(SISLCurve* c) { return c->et[c->ik - 1]; }
static double curve_param_hi(SISLCurve* c) { return c->et[c->in]; }

struct SSResult {
  int num_int_pt = 0, num_int_cr = 0, topo_stat = 0;
  long guide_pts = 0;
  int traced_curves = 0, traced_pts = 0, trace_err = 0, trace_sing = 0;
  int trace_samples = 0;
  double res_pt = 0, res_guide = 0, res_trace = 0;
  double t_topo = 0, t_trace = 0;
  string detail;
};

static SSResult run_ssi(SISLSurf* s1, SISLSurf* s2, double epsge, bool trace) {
  SSResult r;
  const double epsco = 1.0e-15;
  int numintpt = 0, numintcr = 0, stat = 0;
  double *pp1 = nullptr, *pp2 = nullptr;
  SISLIntcurve** ic = nullptr;
  char buf[512];

  double t0 = now_ms();
  s1859(s1, s2, epsco, epsge, &numintpt, &pp1, &pp2, &numintcr, &ic, &stat);
  r.t_topo = now_ms() - t0;
  r.topo_stat = stat;
  r.num_int_pt = numintpt;
  r.num_int_cr = numintcr;

  for (int i = 0; i < numintpt; ++i) {
    double a[3], b[3];
    eval_surf(s1, pp1[2 * i], pp1[2 * i + 1], a);
    eval_surf(s2, pp2[2 * i], pp2[2 * i + 1], b);
    double d = dist3(a, b);
    if (d > r.res_pt) r.res_pt = d;
    snprintf(buf, sizeof buf, "  point[%d] par1=(%.6f,%.6f) par2=(%.6f,%.6f) P=(%.6f,%.6f,%.6f) res=%.3e\n",
             i, pp1[2*i], pp1[2*i+1], pp2[2*i], pp2[2*i+1], a[0], a[1], a[2], d);
    r.detail += buf;
  }

  for (int i = 0; i < numintcr; ++i) {
    SISLIntcurve* c = ic[i];
    r.guide_pts += c->ipoint;
    for (int j = 0; j < c->ipoint; ++j) {
      double a[3], b[3];
      eval_surf(s1, c->epar1[2*j], c->epar1[2*j+1], a);
      eval_surf(s2, c->epar2[2*j], c->epar2[2*j+1], b);
      double d = dist3(a, b);
      if (d > r.res_guide) r.res_guide = d;
    }
    snprintf(buf, sizeof buf, "  curve[%d] ipoint=%d itype=%d start1=(%.4f,%.4f) end1=(%.4f,%.4f) start2=(%.4f,%.4f) end2=(%.4f,%.4f)\n",
             i, c->ipoint, c->itype,
             c->epar1[0], c->epar1[1],
             c->epar1[2*(c->ipoint-1)], c->epar1[2*(c->ipoint-1)+1],
             c->epar2 ? c->epar2[0] : 0.0, c->epar2 ? c->epar2[1] : 0.0,
             c->epar2 ? c->epar2[2*(c->ipoint-1)] : 0.0,
             c->epar2 ? c->epar2[2*(c->ipoint-1)+1] : 0.0);
    r.detail += buf;
  }

  if (trace) {
    double t1 = now_ms();
    for (int i = 0; i < numintcr; ++i) {
      SISLIntcurve* c = ic[i];
      int tstat = 0;
      s1310(s1, s2, c, epsge, 0.0, 2, 0, &tstat);
      if (tstat < 0) { ++r.trace_err; continue; }
      if (tstat == 3) ++r.trace_sing;
      if (!c->pgeom) continue;
      ++r.traced_curves;
      r.traced_pts += c->pgeom->in;
      SISLCurve* pg = c->pgeom;
      int np = pg->in;
      // endpoints of the 3D curve
      snprintf(buf, sizeof buf, "  traced[%d] stat=%d ctrl=%d P_start=(%.4f,%.4f,%.4f) P_end=(%.4f,%.4f,%.4f)\n",
               i, tstat, np,
               pg->ecoef[0], pg->ecoef[1], pg->ecoef[2],
               pg->ecoef[3*(np-1)], pg->ecoef[3*(np-1)+1], pg->ecoef[3*(np-1)+2]);
      r.detail += buf;
      // Residual of the traced curve.  For graph patches S=(u,v,f) the x,y
      // coordinates of a 3D intersection point are exactly the (u,v) parameters
      // of BOTH patches, so the traced 3D curve can be checked directly against
      // both surfaces without relying on the parameter-plane curves:
      //   residual = max_k max(|P_k - S1(Px,Py)|, |P_k - S2(Px,Py)|)
      {
        double lo = curve_param_lo(pg), hi = curve_param_hi(pg);
        for (int k = 0; k <= 40; ++k) {
          double tt = lo + (hi - lo) * double(k) / 40.0;
          double P[3], a[3], b[3];
          eval_curve3(pg, tt, P);
          eval_surf(s1, P[0], P[1], a);
          eval_surf(s2, P[0], P[1], b);
          double d1 = dist3(a, P), d2 = dist3(b, P);
          double d = d1 > d2 ? d1 : d2;
          if (d > r.res_trace) r.res_trace = d;
          ++r.trace_samples;
        }
      }
    }
    r.t_trace = now_ms() - t1;
  }

  if (pp1) free(pp1);
  if (pp2) free(pp2);
  if (numintcr > 0) freeIntcrvlist(ic, numintcr);
  return r;
}

struct SCResult {
  int num_int_pt = 0, num_int_cr = 0, stat = 0;
  double maxres = 0.0, time_ms = 0.0;
  string detail;
};

static SCResult run_sci(SISLSurf* s, SISLCurve* c, double epsge) {
  SCResult r;
  const double epsco = 1.0e-15;
  int numintpt = 0, numintcr = 0, stat = 0;
  double *pp1 = nullptr, *pp2 = nullptr;
  SISLIntcurve** ic = nullptr;
  char buf[512];
  double t0 = now_ms();
  s1858(s, c, epsco, epsge, &numintpt, &pp1, &pp2, &numintcr, &ic, &stat);
  r.time_ms = now_ms() - t0;
  r.stat = stat;
  r.num_int_pt = numintpt;
  r.num_int_cr = numintcr;
  for (int i = 0; i < numintpt; ++i) {
    double a[3], b[3];
    eval_surf(s, pp1[2*i], pp1[2*i+1], a);
    eval_curve3(c, pp2[i], b);
    double d = dist3(a, b);
    if (d > r.maxres) r.maxres = d;
    snprintf(buf, sizeof buf, "  point[%d] par_surf=(%.6f,%.6f) par_crv=%.6f P=(%.6f,%.6f,%.6f) res=%.3e\n",
             i, pp1[2*i], pp1[2*i+1], pp2[i], a[0], a[1], a[2], d);
    r.detail += buf;
  }
  if (pp1) free(pp1);
  if (pp2) free(pp2);
  if (numintcr > 0) freeIntcrvlist(ic, numintcr);
  return r;
}

// ---------- fixtures ----------
static const char* kTruth[6] = {
    "1 curve: diagonal s=t (transverse except at origin)",
    "0 intersections",
    "1 point: corner (0,0) tangency",
    "1 tangential/degenerate curve s=t",
    "1 interior curve: u^3 + v^3 = 0.1 (bicubic x bicubic)",
    "1 point: (0.866025,0.5,0.5)"};

static SSResult run_fixture(int f, double epsge, bool trace) {
  SSResult r;
  if (f == 0) {  // A: f = u^2 - v^2
    vector<vector<double>> a = {{0,0,-1},{0,0,0},{1,0,0}};
    SISLSurf* s1 = make_graph_mono(a, 3, 3);
    vector<vector<double>> z = {{0}};
    SISLSurf* s2 = make_graph_mono(z, 2, 2);
    r = run_ssi(s1, s2, epsge, trace);
    freeSurf(s1); freeSurf(s2);
  } else if (f == 1) {  // B: f = 0.3
    vector<vector<double>> a = {{0.3}};
    SISLSurf* s1 = make_graph_mono(a, 2, 2);
    vector<vector<double>> z = {{0}};
    SISLSurf* s2 = make_graph_mono(z, 2, 2);
    r = run_ssi(s1, s2, epsge, trace);
    freeSurf(s1); freeSurf(s2);
  } else if (f == 2) {  // C: f = u + 0.5 v
    vector<vector<double>> a = {{0.0, 0.5}, {1.0, 0.0}};
    SISLSurf* s1 = make_graph_mono(a, 2, 2);
    vector<vector<double>> z = {{0}};
    SISLSurf* s2 = make_graph_mono(z, 2, 2);
    r = run_ssi(s1, s2, epsge, trace);
    freeSurf(s1); freeSurf(s2);
  } else if (f == 3) {  // D: f = (u-v)^2 = u^2 - 2uv + v^2
    vector<vector<double>> a = {{0,0,1},{0,-2,0},{1,0,0}};
    SISLSurf* s1 = make_graph_mono(a, 3, 3);
    vector<vector<double>> z = {{0}};
    SISLSurf* s2 = make_graph_mono(z, 2, 2);
    r = run_ssi(s1, s2, epsge, trace);
    freeSurf(s1); freeSurf(s2);
  } else if (f == 4) {  // E: u^3 - v^3 + 0.1
    vector<vector<double>> a1(4, vector<double>(4, 0.0));
    a1[3][0] = 1.0;
    SISLSurf* s1 = make_graph_mono(a1, 4, 4);
    vector<vector<double>> a2(4, vector<double>(4, 0.0));
    a2[0][0] = 0.1; a2[0][3] = -1.0;
    SISLSurf* s2 = make_graph_mono(a2, 4, 4);
    r = run_ssi(s1, s2, epsge, trace);
    freeSurf(s1); freeSurf(s2);
  } else {  // F handled separately (surface/curve)
    r.topo_stat = -999;
  }
  return r;
}

static int selftest() {
  int fail = 0;
  double p[3];
  vector<vector<double>> au = {{0.0}, {1.0}};
  SISLSurf* su = make_graph_mono(au, 2, 2);
  eval_surf(su, 0.5, 0.5, p);
  printf("selftest f=u   at (0.5,0.5): (%.6f %.6f %.6f) expect (0.5,0.5,0.5)\n", p[0], p[1], p[2]);
  if (std::fabs(p[0]-0.5) > 1e-12 || std::fabs(p[1]-0.5) > 1e-12 || std::fabs(p[2]-0.5) > 1e-12) ++fail;

  vector<vector<double>> a2 = {{0,0,0},{0,0,0},{1,0,0}};
  SISLSurf* s2 = make_graph_mono(a2, 3, 3);
  eval_surf(s2, 0.5, 0.25, p);
  printf("selftest f=u^2 at (0.5,0.25): (%.6f %.6f %.6f) expect (0.5,0.25,0.25)\n", p[0], p[1], p[2]);
  if (std::fabs(p[0]-0.5) > 1e-12 || std::fabs(p[1]-0.25) > 1e-12 || std::fabs(p[2]-0.25) > 1e-12) ++fail;

  SISLSurf* cyl = make_quarter_cylinder();
  eval_surf(cyl, 0.5, 0.5, p);  // rational param u=0.5 is the 45 deg point
  printf("selftest cyl at (0.5,0.5): (%.6f %.6f %.6f) expect (0.707107,0.707107,0.500000)\n", p[0], p[1], p[2]);
  if (std::fabs(p[0]-std::sqrt(2.0)/2) > 1e-6 || std::fabs(p[1]-std::sqrt(2.0)/2) > 1e-6 ||
      std::fabs(p[2]-0.5) > 1e-12) ++fail;

  double q0[3] = {0.5, 0.5, 0.5}, q1[3] = {1.1, 0.5, 0.5};
  SISLCurve* ln = make_line(q0, q1);
  eval_curve3(ln, 0.5, p);
  printf("selftest line at 0.5: (%.6f %.6f %.6f) expect (0.800000,0.500000,0.500000)\n", p[0], p[1], p[2]);
  if (std::fabs(p[0]-0.8) > 1e-12) ++fail;

  freeSurf(su); freeSurf(s2); freeSurf(cyl); freeCurve(ln);
  printf("selftest: %s\n", fail ? "FAIL" : "PASS");
  return fail;
}

int main(int argc, char** argv) {
  bool detail = (argc > 1 && string(argv[1]) == "--detail");
  if (argc > 1 && string(argv[1]) == "--selftest") return selftest();

  const double epsges[2] = {1.0e-6, 1.0e-4};
  const char* names[6] = {"A", "B", "C", "D", "E", "F"};

  // warm-up so that the reported timings exclude one-time library initialisation
  {
    SSResult w = run_fixture(1, 1.0e-6, false);  // B, cheapest fixture
    (void)w;
    SISLSurf* cyl = make_quarter_cylinder();
    double q0[3] = {0.5, 0.5, 0.5}, q1[3] = {1.1, 0.5, 0.5};
    SISLCurve* ln = make_line(q0, q1);
    SCResult wc = run_sci(cyl, ln, 1.0e-6);
    (void)wc;
    freeSurf(cyl); freeCurve(ln);
  }

  if (!detail) {
    printf("fixture,epsge,topology_stat,num_int_points,num_int_curves,"
           "guide_points,traced_curves,traced_ctrl_pts,trace_errors,"
           "res_point,res_guide,res_trace,time_topology_ms,time_trace_ms,truth\n");
    for (int f = 0; f < 6; ++f) {
      for (int ei = 0; ei < 2; ++ei) {
        double epsge = epsges[ei];
        if (f == 5) {
          SISLSurf* cyl = make_quarter_cylinder();
          double q0[3] = {0.5, 0.5, 0.5}, q1[3] = {1.1, 0.5, 0.5};
          SISLCurve* ln = make_line(q0, q1);
          SCResult sc = run_sci(cyl, ln, epsge);
          printf("F,%.1e,%d,%d,%d,0,0,0,0,%.3e,0,0,%.3f,0,\"%s\"\n",
                 epsge, sc.stat, sc.num_int_pt, sc.num_int_cr, sc.maxres,
                 sc.time_ms, kTruth[f]);
          freeSurf(cyl); freeCurve(ln);
          continue;
        }
        SSResult r = run_fixture(f, epsge, true);
        printf("%s,%.1e,%d,%d,%d,%ld,%d,%d,%d,%.3e,%.3e,%.3e,%.3f,%.3f,\"%s\"\n",
               names[f], epsge, r.topo_stat, r.num_int_pt, r.num_int_cr,
               r.guide_pts, r.traced_curves, r.traced_pts,
               r.trace_err + r.trace_sing, r.res_pt, r.res_guide, r.res_trace,
               r.t_topo, r.t_trace, kTruth[f]);
        fflush(stdout);
      }
    }
  } else {
    for (int f = 0; f < 6; ++f) {
      printf("### fixture %s  (truth: %s)\n", names[f], kTruth[f]);
      if (f == 5) {
        SISLSurf* cyl = make_quarter_cylinder();
        double q0[3] = {0.5, 0.5, 0.5}, q1[3] = {1.1, 0.5, 0.5};
        SISLCurve* ln = make_line(q0, q1);
        SCResult sc = run_sci(cyl, ln, 1.0e-6);
        printf("s1858 stat=%d num_pt=%d num_cr=%d res=%.3e t=%.3fms\n",
               sc.stat, sc.num_int_pt, sc.num_int_cr, sc.maxres, sc.time_ms);
        fputs(sc.detail.c_str(), stdout);
        freeSurf(cyl); freeCurve(ln);
      } else {
        SSResult r = run_fixture(f, 1.0e-6, true);
        printf("s1859 stat=%d num_pt=%d num_cr=%d guide=%ld traced=%d err=%d sing=%d\n",
               r.topo_stat, r.num_int_pt, r.num_int_cr, r.guide_pts,
               r.traced_curves, r.trace_err, r.trace_sing);
        printf("res_point=%.3e res_guide=%.3e res_trace=%.3e t_topo=%.3fms t_trace=%.3fms\n",
               r.res_pt, r.res_guide, r.res_trace, r.t_topo, r.t_trace);
        fputs(r.detail.c_str(), stdout);
      }
      fflush(stdout);
    }
  }
  return 0;
}
