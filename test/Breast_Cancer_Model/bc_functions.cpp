#include "bc_functions.hpp"
#include <cmath>

BCParams bc_p = {0.2, 0.3, 1.0, 5.0, 0.2, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

void bc_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history)
{
    double hist3 = history.at_time(t - bc_p.tau, 2);
    du[0] = (bc_p.v0 / (1 + bc_p.beta0 * std::pow(hist3, 2))) * (bc_p.p0 - bc_p.q0) * u[0] - bc_p.d0 * u[0];
    du[1] = (bc_p.v0 / (1 + bc_p.beta0 * std::pow(hist3, 2))) * (1 - bc_p.p0 + bc_p.q0) * u[0] +
            (bc_p.v1 / (1 + bc_p.beta1 * std::pow(hist3, 2))) * (bc_p.p1 - bc_p.q1) * u[1] - bc_p.d1 * u[1];
    du[2] = (bc_p.v1 / (1 + bc_p.beta1 * std::pow(hist3, 2))) * (1 - bc_p.p1 + bc_p.q1) * u[1] - bc_p.d2 * u[2];
}

double history_bc(double t)
{
    return 1.0;
}
