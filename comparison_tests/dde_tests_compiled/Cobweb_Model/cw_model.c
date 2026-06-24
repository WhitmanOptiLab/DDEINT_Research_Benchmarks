#include <dde/dde.h>
#include <stddef.h>
#include <math.h>

// Parameters (passed in order):
//   a, b, c, d, tau, beta, speed
void cw_model(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double a     = pars[0];
  double b     = pars[1];
  double c     = pars[2];
  double d     = pars[3];
  double tau   = pars[4];
  double beta  = pars[5];
  double speed = pars[6];

  // Lag: y[1] (p_expected, 0-indexed) at t - tau
  static const int idx[1] = {1};
  double p_exp_lag[1];
  ylag_vec_int(t - tau, idx, 1, p_exp_lag);

  double demand = a * y[0] + b;
  double supply = c * p_exp_lag[0] + d;

  dydt[0] = speed * (demand - supply);
  dydt[1] = beta  * (y[0] - y[1]);
}

#include <dde/dde.c>