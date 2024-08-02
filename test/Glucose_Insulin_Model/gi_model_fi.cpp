#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

#include <iostream>
#include <ctime>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

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

int main()
{
    std::cout << "Running DDE solver..." << std::endl;
    std::cout << "Running the Glucose-Insulin Regulation Model..." << std::endl;

    // Initial conditions and time span
    std::vector<double> u0 = {gi_p.G_b, gi_p.I_b};
    double t_initial = 0.0;
    double t_final = 200.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_gi, history_gi};
    std::vector<double> max_delays = {gi_p.tau, gi_p.tau}; // Ensure the size matches the number of equations

    DDEint_dopri_5<gi_dde> dde_solver(2, max_delays, prehistory);

    // Measure the time taken to solve the problem
    std::clock_t start = std::clock();
    std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 10000, ABS_TOL, REL_TOL);
    // Measure the time taken to solve the problem
    std::clock_t end = std::clock();

    double elapsed_time = (end - start) / (double) CLOCKS_PER_SEC * 1000;  // Convert to milliseconds
    std::cout << "Glucose-Insulnn Regulatiom Model took " << elapsed_time << " milliseconds to solve." << std::endl;

    return 0;
}
