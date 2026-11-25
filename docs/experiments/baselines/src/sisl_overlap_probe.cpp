// Supplementary SISL probes: 2D-overlap (coplanar / identical patch) cases.
// These correspond to the certified framework's "coplanar overlap" case, which
// is not a curve intersection.  Used to document SISL's failure mode on the
// paper's overlap fixture (see comparison.md).
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include "sisl.h"
using std::vector;

static vector<double> bezier_knots(int k) {
  vector<double> et(2 * k);
  for (int i = 0; i < k; ++i) { et[i] = 0.0; et[k + i] = 1.0; }
  return et;
}

// graph surface over [0,1]^2: S(u,v) = (u0 + su*u, v0 + sv*v, z(u,v))
// z given by monomial coefficients a[i][j] of u^i v^j
static SISLSurf* graph(const vector<vector<double>>& a, int n,
                       double u0 = 0.0, double su = 1.0,
                       double v0 = 0.0, double sv = 1.0) {
  auto bin = [](int N, int K) {
    double r = 1.0;
    for (int i = 0; i < K; ++i) r = r * (N - i) / (i + 1);
    return r;
  };
  vector<vector<double>> z(n, vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double s = 0;
      for (int k = 0; k <= i; ++k)
        for (int l = 0; l <= j; ++l) {
          double ak = (k < (int)a.size() && l < (int)a[k].size()) ? a[k][l] : 0.0;
          s += ak * bin(i, k) / bin(n - 1, k) * bin(j, l) / bin(n - 1, l);
        }
      z[i][j] = s;
    }
  vector<double> et = bezier_knots(n);
  vector<double> coef((size_t)n * n * 3);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      size_t idx = ((size_t)j * n + i) * 3;
      coef[idx] = u0 + su * double(i) / (n - 1);
      coef[idx + 1] = v0 + sv * double(j) / (n - 1);
      coef[idx + 2] = z[i][j];
    }
  return newSurf(n, n, n, n, et.data(), et.data(), coef.data(), 3, 3, 1);
}

static void probe(const char* name, SISLSurf* s1, SISLSurf* s2, double epsge) {
  int npt = 0, ncr = 0, st = 0;
  double *p1 = nullptr, *p2 = nullptr;
  SISLIntcurve** ic = nullptr;
  s1859(s1, s2, 1e-15, epsge, &npt, &p1, &p2, &ncr, &ic, &st);
  printf("%-52s epsge=%.0e stat=%d npt=%d ncr=%d", name, epsge, st, npt, ncr);
  for (int i = 0; i < ncr; ++i) {
    printf(" | crv%d: ipoint=%d itype=%d ep1=(%.3f,%.3f)->(%.3f,%.3f)", i,
           ic[i]->ipoint, ic[i]->itype, ic[i]->epar1[0], ic[i]->epar1[1],
           ic[i]->epar1[2 * (ic[i]->ipoint - 1)],
           ic[i]->epar1[2 * (ic[i]->ipoint - 1) + 1]);
  }
  printf("\n");
  if (p1) free(p1);
  if (p2) free(p2);
  if (ncr) freeIntcrvlist(ic, ncr);
}

int main() {
  const double eps = 1e-6;
  // G: identical graph patches (same parametrization) -> full 2D overlap
  {
    vector<vector<double>> a = {{0,0,-1},{0,0,0},{1,0,0}};
    SISLSurf* s1 = graph(a, 3);
    SISLSurf* s2 = graph(a, 3);
    probe("G  identical patches z=u^2-v^2 (full 2D overlap)", s1, s2, eps);
    freeSurf(s1); freeSurf(s2);
  }
  // G2: identical planes z=0 -> full 2D overlap
  {
    vector<vector<double>> a = {{0}};
    SISLSurf* s1 = graph(a, 2);
    SISLSurf* s2 = graph(a, 2);
    probe("G2 identical planes z=0 (full 2D overlap)", s1, s2, eps);
    freeSurf(s1); freeSurf(s2);
  }
  // H: two coplanar patches with a partial 2D overlap (shifted in x)
  {
    vector<vector<double>> a = {{0}};
    SISLSurf* s1 = graph(a, 2);
    SISLSurf* s2 = graph(a, 2, 0.5, 1.0);  // x in [0.5,1.5]
    probe("H  coplanar planes z=0, x in [0,1] vs [0.5,1.5]", s1, s2, eps);
    freeSurf(s1); freeSurf(s2);
  }
  // H2: coplanar planes meeting only along a boundary edge (x=1)
  {
    vector<vector<double>> a = {{0}};
    SISLSurf* s1 = graph(a, 2);
    SISLSurf* s2 = graph(a, 2, 1.0, 1.0);  // x in [1,2]
    probe("H2 coplanar planes z=0, x in [0,1] vs [1,2] (edge contact)", s1, s2, eps);
    freeSurf(s1); freeSurf(s2);
  }
  // I: same plane, transposed parametrization (same point set)
  {
    vector<vector<double>> a = {{0}};
    SISLSurf* s1 = graph(a, 2);
    // transposed: swap u and v -> build manually with x=v, y=u
    vector<double> et = bezier_knots(2);
    vector<double> coef = {0,0,0, 0,1,0, 1,0,0, 1,1,0};
    SISLSurf* s2 = newSurf(2, 2, 2, 2, et.data(), et.data(), coef.data(), 3, 3, 1);
    probe("I  plane z=0 with transposed (u,v) parametrization", s1, s2, eps);
    freeSurf(s1); freeSurf(s2);
  }
  // J: near-coplanar (z = 1e-7 * u), below/at epsge
  {
    vector<vector<double>> a0 = {{0}};
    vector<vector<double>> a1 = {{0}, {1e-7}};
    SISLSurf* s1 = graph(a0, 2);
    SISLSurf* s2 = graph(a1, 2);
    probe("J  near-coplanar z=0 vs z=1e-7*u", s1, s2, eps);
    freeSurf(s1); freeSurf(s2);
  }
  return 0;
}
