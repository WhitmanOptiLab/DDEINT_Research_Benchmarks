#include "Methods/Dormand_Prince/DoPri_5.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>

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

void bc_dde(double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history)
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
    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.0};
    double t_i = 0.0;
    double t_f = 10.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_bc, history_bc, history_bc};
    std::vector<double> max_delays = {bc_p.tau, bc_p.tau, bc_p.tau}; // Ensure the size matches the number of equations

    DoPri_5<bc_dde> dde_solver(3, max_delays, prehistory);
    Results results;
    dde_solver.initialize(0, 0.1, 1e-5, u0, ABS_TOL, REL_TOL);
    // Start clock
    const auto start{std::chrono::high_resolution_clock::now()};
    results = dde_solver.solve(t_i, t_f, 200, 100000);
    // Stop timer
    const auto end{std::chrono::high_resolution_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    // Print time(diff between two times)
    std::cout << "Time Taken: " << elapsed_seconds.count() << std::endl;
    

    return 0;
}
