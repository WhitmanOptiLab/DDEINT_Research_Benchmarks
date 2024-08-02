#pragma once 

#include "../../DDEINT/history/history.hpp"
#include <vector>

struct CVParams {
    double ca;
    double cv;
    double R;
    double r;
    double Vstr;
    double alpha0;
    double alphas;
    double alphap;
    double alphaH;
    double beta0;
    double betas;
    double betap;
    double betaH;
    double gammaH;
    double tau;
};

extern CVParams cv_p;

void cv_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history);

double history_Pa(double t);

double history_Pv(double t);

double history_H(double t);  
