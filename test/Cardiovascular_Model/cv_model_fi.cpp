#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

#include <iostream>
#include <ctime>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

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

int main()
{
    std::cout << "Running DDE solver..." << std::endl;
    std::cout << "Running the Cardiovascular Model..." << std::endl;
    
    // Initial conditions and time span
    std::vector<double> u0 = {93, (1 / (1 + cv_p.R / cv_p.r)) * 93, (1 / (cv_p.R * cv_p.Vstr)) * (1 / (1 + cv_p.r / cv_p.R)) * 93};
    double t_initial = 0.0;
    double t_final = 1000.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_Pa, history_Pv, history_H};
    std::vector<double> max_delays = {4.0, 4.0, 4.0}; // Ensure the size matches the number of equations

    DDEint_dopri_5<cv_dde> dde_solver(3, max_delays, prehistory);

    // Measure the time taken to solve the problem
    std::clock_t start = std::clock();
    std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 20000, ABS_TOL, REL_TOL);
    // Measure the time taken to solve the problem
    std::clock_t end = std::clock();

    double elapsed_time = (end - start) / (double) CLOCKS_PER_SEC;
    // conver to milliseconds
    elapsed_time = elapsed_time * 1000;
    std::cout << "Cardiovascular Model took " << elapsed_time << " milliseconds to solve." << std::endl;

    return 0;
}
