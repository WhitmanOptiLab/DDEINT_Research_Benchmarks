#include "cw_functions.hpp"

// cw params 
CWParams cw_p = {-0.5, 1.5, 0.8, 0.0, 1.5, 0.6, 1.2};

void cw_dde (double t, std::vector<double>& u, std::vector<double>& du, History<double,double>& history) 
{ 
    double p_now = u[0];
    
    // Past expected price: p^E(t - tau) used to calculate delayed market supply
    double p_expected_delayed = HIST_AT_TIME(history, t - cw_p.tau, 1);


    double demand = cw_p.a * p_now + cw_p.b;
    double supply =  cw_p.c * p_expected_delayed + cw_p.d;

    // dp/dt = speed * (Demand - Supply)
    du[0] = cw_p.speed * (demand - supply);

    // dp^E/dt = beta (real price now - thought of current price)
    du[1] = cw_p.beta * (p_now - u[1]);
}

double history_cw(double t)
{   
    return 0.4;
}