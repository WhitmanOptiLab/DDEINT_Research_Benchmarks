#include "gi_functions.hpp"

GIParams gi_p = {10.0, 5.0, 1.0, 90.0, 10.0, 0.1, 0.05, 5.0};

void gi_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history)
{
    double G_Tau = history.at_time(t - gi_p.tau, 0);
    double I_Tau = history.at_time(t - gi_p.tau, 1);

    du[0] = (u[0] - gi_p.G_b) / gi_p.V_G - (gi_p.S_I * I_Tau * u[0]) / gi_p.V_G + 10.0;
    du[1] = - gi_p.n * (u[1] - gi_p.I_b) + (gi_p.gamma * G_Tau) / gi_p.V_I;
}

double history_gi(double t)
{
    return 1.0;
}
