#include "cv_functions.hpp"
#include <cmath>

CVParams cv_p = {1.55, 519, 1.05, 0.068, 67.9, 93, 93, 93, 0.84, 7, 7, 7, 1.17, 0, 4.0};

void cv_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history) 
{
    double R = t <= 600 ? 1.05 : 0.21 * std::exp(600 - t) + 0.84;
    double Patau = history.at_time(t - cv_p.tau, 0);
    double Paoft = u[0];
    double Pvoft = u[1];
    double Hoft = u[2];

    du[0] = - (1 / (cv_p.ca * R)) * Paoft + (1 / (cv_p.ca * R)) * Pvoft + (1 / cv_p.ca) * cv_p.Vstr * Hoft;
    du[1] = (1 / (cv_p.cv * R)) * Paoft - (1 / (cv_p.cv * R) + 1 / (cv_p.cv * cv_p.r)) * Pvoft;
    double Ts = 1 / (1 + std::pow((Patau / cv_p.alphas), cv_p.betas));
    double Tp = 1 / (1 + std::pow((cv_p.alphap / Paoft), cv_p.betap));
    du[2] = (cv_p.alphaH * Ts) / (1 + cv_p.gammaH * Tp) - cv_p.betaH * Tp;
}

double history_Pa(double t) 
{
    return 93;
}

double history_Pv(double t) 
{
    return (1 / (1 + cv_p.R / cv_p.r)) * 93;
}

double history_H(double t) 
{
    return (1 / (cv_p.R * cv_p.Vstr)) * (1 / (1 + cv_p.r / cv_p.R)) * 93;
}
