#include <dde/dde.h>
#include <stddef.h>
#include <math.h>

// Parameters (passed in order):
//   p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau
void bc_model(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars  = (double*) data;
  double p0     = pars[0];
  double q0     = pars[1];
  double v0     = pars[2];
  double d0     = pars[3];
  double p1     = pars[4];
  double q1     = pars[5];
  double v1     = pars[6];
  double d1     = pars[7];
  double d2     = pars[8];
  double beta0  = pars[9];
  double beta1  = pars[10];
  double tau    = pars[11];

  // Lag: y[2] (third state, 0-indexed) at t - tau
  static const int idx[1] = {2};
  double hist3[1];
  ylag_vec_int(t - tau, idx, 1, hist3);

  double h3sq   = hist3[0] * hist3[0];
  double frac0  = v0 / (1.0 + beta0 * h3sq);
  double frac1  = v1 / (1.0 + beta1 * h3sq);

  dydt[0] = frac0 * (p0 - q0) * y[0] - d0 * y[0];
  dydt[1] = frac0 * (1.0 - p0 + q0) * y[0]
          + frac1 * (p1 - q1) * y[1]
          - d1 * y[1];
  dydt[2] = frac1 * (1.0 - p1 + q1) * y[1] - d2 * y[2];
}

// This needs to be included exactly once per shared library.
#include <dde/dde.c>