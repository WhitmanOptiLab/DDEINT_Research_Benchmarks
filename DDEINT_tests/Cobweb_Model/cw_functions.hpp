#pragma once 

#include "../../DDEINT/history/history.hpp"
#include <vector>

struct CWParams
{   
    // Configured with Thiemer's baseline parameters
    // Demand slope: a, Base demand: b, Supply slope: c, Base supply: d
    double a, b;  // Demand Parameters: x_D = a*p + b 
    double c, d;  // Supply parameters: x_S = c*p + d
    double tau;   // Production Delay (1.0 cycle)
    double beta; 
    double speed; // Price adjustment speed (How fast the market is updating)

};

extern CWParams cw_p;

void cw_dde(double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history);

double history_cw(double t);
