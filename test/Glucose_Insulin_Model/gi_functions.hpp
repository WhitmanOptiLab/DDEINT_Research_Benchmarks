#pragma once

#include "../../DDEINT/history/history.hpp"

#include <vector>

struct GIParams
{
    double V_G;
    double V_I;
    double S_I;
    double G_b;
    double I_b;
    double n;
    double gamma;
    double tau;
};

extern GIParams gi_p;

void gi_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history);

double history_gi(double t);
