#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <ctime>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

// Breast cancer model

// Define the delay
const double tau = 1.0;
// const std::vector<double> lags = {tau};

// Parameters for the model
const double p0 = 0.2;
const double q0 = 0.3;
const double v0 = 1.0;
const double d0 = 5.0;
const double p1 = 0.2;
const double q1 = 0.3;
const double v1 = 1.0;
const double d1 = 1.0;
const double d2 = 1.0;
const double beta0 = 1.0;
const double beta1 = 1.0;

// Define the DDE function
void bc_ddefun(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history) 
{
    double hist3 = history.at_time(t - tau, 2);
    du[0] = (v0 / (1 + beta0 * std::pow(hist3, 2))) * (p0 - q0) * u[0] - d0 * u[0];
    du[1] = (v0 / (1 + beta0 * std::pow(hist3, 2))) * (1 - p0 + q0) * u[0] +
            (v1 / (1 + beta1 * std::pow(hist3, 2))) * (p1 - q1) * u[1] - d1 * u[1];
    du[2] = (v1 / (1 + beta1 * std::pow(hist3, 2))) * (1 - p1 + q1) * u[1] - d2 * u[2];
}

// Define the history function
double history_function(double t) 
{
    return 1.0;
}

double BC_Model() 
{
    std::cout << "Running DDE solver..." << std::endl;
    std::cout << "Running the Breast Cancer Model..." << std::endl;
    
    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.0};
    double t_i = 0.0;
    double t_f = 10.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_function, history_function, history_function};
    std::vector<double> max_delays = {tau, tau, tau}; // Ensure the size matches the number of equations

    DDEint_dopri_5<bc_ddefun> dde_solver(3, max_delays, prehistory);

    // Measure the time taken to solve the problem
    std::clock_t start = std::clock();
    std::vector<std::vector<double>> solution = dde_solver.run(t_i, t_f, u0, 0.1, 1e-5, 10000, ABS_TOL, REL_TOL);
    // Measure the time taken to solve the problem
    std::clock_t end = std::clock();

    double elapsed_time = (end - start) / (double) CLOCKS_PER_SEC;
    // convert to milliseconds
    elapsed_time = elapsed_time * 1000;
    std::cout << "Breast Cancer Model took " << elapsed_time << " milliseconds to solve." << std::endl;

    // Prepare data for plotting
    std::vector<double> time_points;
    std::vector<double> u1_values;
    std::vector<double> u2_values;
    std::vector<double> u3_values;

    for (const auto& row : solution) 
    {
        time_points.push_back(row[0]);
        u1_values.push_back(row[1]);
        u2_values.push_back(row[2]);
        u3_values.push_back(row[3]);
    }

    // save the data to a csv file
    std::ofstream file("data/bc_model_cpp.csv");
    file << "Time,u1,u2,u3\n";
    for (size_t i = 0; i < time_points.size(); ++i) 
    {
        file << time_points[i] << "," << u1_values[i] << "," << u2_values[i] << "," << u3_values[i] << "\n";
    }
    file.close();

    return elapsed_time;
}