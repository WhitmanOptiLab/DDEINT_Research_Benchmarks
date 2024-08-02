#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

#include <iostream>
#include <ctime>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

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

BCParams bc_p = {0.2, 0.3, 1.0, 5.0, 0.2, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

void bc_dde(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history)
{
    double hist3 = history.at_time(t - bc_p.tau, 2);
    du[0] = (bc_p.v0 / (1 + bc_p.beta0 * std::pow(hist3, 2))) * (bc_p.p0 - bc_p.q0) * u[0] - bc_p.d0 * u[0];
    du[1] = (bc_p.v0 / (1 + bc_p.beta0 * std::pow(hist3, 2))) * (1 - bc_p.p0 + bc_p.q0) * u[0] +
            (bc_p.v1 / (1 + bc_p.beta1 * std::pow(hist3, 2))) * (bc_p.p1 - bc_p.q1) * u[1] - bc_p.d1 * u[1];
    du[2] = (bc_p.v1 / (1 + bc_p.beta1 * std::pow(hist3, 2))) * (1 - bc_p.p1 + bc_p.q1) * u[1] - bc_p.d2 * u[2];
}

double history_bc(double t)
{
    return 1.0;
}

int main()
{
    std::cout << "Running DDE solver..." << std::endl;
    std::cout << "Running the Breast Cancer Model..." << std::endl;
    
    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.0};
    double t_i = 0.0;
    double t_f = 10.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_bc, history_bc, history_bc};
    std::vector<double> max_delays = {1.0, 1.0, 1.0}; // Ensure the size matches the number of equations

    DDEint_dopri_5<bc_dde> dde_solver(3, max_delays, prehistory);

    // Measure the time taken to solve the problem
    std::clock_t start = std::clock();
    std::vector<std::vector<double>> solution = dde_solver.run(t_i, t_f, u0, 0.1, 1e-5, 10000, ABS_TOL, REL_TOL);
    // Measure the time taken to solve the problem
    std::clock_t end = std::clock();

    double elapsed_time = (end - start) / (double) CLOCKS_PER_SEC;
    // convert to milliseconds
    elapsed_time = elapsed_time * 1000;
    std::cout << "Breast Cancer Model took " << elapsed_time << " milliseconds to solve." << std::endl;

    return 0;
}
