#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <ctime>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

// Repressilator model

// Repressilator parameters
const double tau3 = 0.1;
const double beta = 50;
const double n = 2;
const double k = 1;
const double gamma_r = 1;

// Define the DDE function for the repressilator with delay
void rs_ddefun(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history) 
{
    double x3_tau = history.at_time(t - tau3, 2);
    double x1_tau = history.at_time(t - tau3, 0);
    double x2_tau = history.at_time(t - tau3, 1);

    du[0] = beta / (1 + std::pow(x3_tau / k, n)) - gamma_r * u[0];
    du[1] = beta / (1 + std::pow(x1_tau / k, n)) - gamma_r * u[1];
    du[2] = beta / (1 + std::pow(x2_tau / k, n)) - gamma_r * u[2];
}

// Define the history function
double history_func(double t) 
{
    return 1.0;
}

double RS_Model() 
{
    std::cout << "Running DDE solver..." << std::endl;

    std::cout << "Running the Repressilator Model..." << std::endl;

    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.2};
    double t_initial = 0.0;
    double t_final = 50.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_func, history_func, history_func};
    std::vector<double> max_delays = {tau3, tau3, tau3}; // Ensure the size matches the number of equations

    DDEint_dopri_5<rs_ddefun> dde_solver(3, max_delays, prehistory);

    // Measure the time taken to solve the problem
    std::clock_t start = std::clock();
    std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 10000, ABS_TOL, REL_TOL);
    // Measure the time taken to solve the problem
    std::clock_t end = std::clock();

    double elapsed_time = (end - start) / (double) CLOCKS_PER_SEC;
    // conver to milliseconds
    elapsed_time = elapsed_time * 1000;
    std::cout << "Repressilator Model took " << elapsed_time << " milliseconds to solve." << std::endl;

    // Prepare data for plotting
    std::vector<double> time_points;
    std::vector<double> x1_values;
    std::vector<double> x2_values;
    std::vector<double> x3_values;

    for (const auto& row : solution) 
    {
        time_points.push_back(row[0]);
        x1_values.push_back(row[1]);
        x2_values.push_back(row[2]);
        x3_values.push_back(row[3]);
    }

    // Save the data to a CSV file
    std::ofstream file("data/repressilator_model_cpp.csv");
    file << "Time,x1,x2,x3\n";
    for (size_t i = 0; i < time_points.size(); ++i) 
    {
        file << time_points[i] << "," << x1_values[i] << "," << x2_values[i] << "," << x3_values[i] << "\n";
    }
    file.close();
    
    return elapsed_time;
}