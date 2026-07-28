#pragma once

#include "../../DDEINT/PolicyInterface/policy.hpp"
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

void gi_dde(double t, std::vector<double>& u, std::vector<double>& du, SimContext& ctx);

double history_gi(double t);
