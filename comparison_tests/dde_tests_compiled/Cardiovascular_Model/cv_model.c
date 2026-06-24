#include <dde/dde.h>
#include <stddef.h>
#include <math.h>

// Parameters (passed in order):
//   ca, cv, R, r, Vstr, alpha0, alphas, alphap, alphaH,
//   beta0, betas, betap, betaH, gammaH, tau
void cv_model(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars   = (double*) data;
  double ca      = pars[0];
  double cv_     = pars[1];
  // pars[2] is R0 (unused; R is computed dynamically below)
  double r       = pars[3];
  double Vstr    = pars[4];
  double alphas  = pars[6];
  double alphap  = pars[7];
  double alphaH  = pars[8];
  double betas   = pars[10];
  double betap   = pars[11];
  double betaH   = pars[12];
  double gammaH  = pars[13];
  double tau     = pars[14];

  // Time-varying resistance (matches R model exactly)
  double R = (t <= 600.0) ? 1.05 : 0.21 * exp(600.0 - t) + 0.84;

  // Lag: y[0] (Pa, 0-indexed) at t - tau
  static const int idx[1] = {0};
  double Patau[1];
  ylag_vec_int(t - tau, idx, 1, Patau);

  double Pa = y[0];
  double Pv = y[1];
  double H  = y[2];

  double dydt0 = -(1.0 / (ca * R)) * Pa
               +  (1.0 / (ca * R)) * Pv
               +  (1.0 / ca) * Vstr * H;

  double dydt1 =  (1.0 / (cv_ * R)) * Pa
               - ((1.0 / (cv_ * R)) + (1.0 / (cv_ * r))) * Pv;

  double Ts   = 1.0 / (1.0 + pow(Patau[0] / alphas, betas));
  double Tp   = 1.0 / (1.0 + pow(alphap / Pa,       betap));
  double dydt2 = (alphaH * Ts) / (1.0 + gammaH * Tp) - betaH * Tp;

  dydt[0] = dydt0;
  dydt[1] = dydt1;
  dydt[2] = dydt2;
}

#include <dde/dde.c>