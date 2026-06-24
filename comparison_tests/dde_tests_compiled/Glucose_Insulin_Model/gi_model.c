#include <dde/dde.h>
#include <stddef.h>
#include <math.h>

// Parameters (passed in order):
//   V_G, V_I, S_I, G_b, I_b, n, gamma, tau
void gi_model(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars  = (double*) data;
  double V_G    = pars[0];
  double V_I    = pars[1];
  double S_I    = pars[2];
  double G_b    = pars[3];
  double I_b    = pars[4];
  double n_par  = pars[5];
  double gamma  = pars[6];
  double tau    = pars[7];

  // Lag: y[0] (G) and y[1] (I) at t - tau
  static const int idx[2] = {0, 1};
  double ylags[2];
  ylag_vec_int(t - tau, idx, 2, ylags);
  double G_tau = ylags[0];
  double I_tau = ylags[1];

  dydt[0] = (y[0] - G_b) / V_G
           - (S_I * I_tau * y[0]) / V_G
           + 10.0;

  dydt[1] = -n_par * (y[1] - I_b)
           + (gamma * G_tau) / V_I;
}

#include <dde/dde.c>