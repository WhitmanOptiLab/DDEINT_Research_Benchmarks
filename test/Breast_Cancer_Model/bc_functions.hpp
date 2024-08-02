#pragma once

#include "../../DDEINT/history/history.hpp"

#include <vector>

struct BCParams
{
    double p0;
    double q0;
    double v0;
    double d0;
    double p1;
    double q1;
    double v1;
    double d1;
    double d2;
    double beta0;
    double beta1;
    double tau;
};

extern BCParams bc_p;

void bc_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history);

double history_bc(double t);
