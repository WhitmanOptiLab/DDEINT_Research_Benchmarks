#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <ctime>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

// Cardiovascular model

// Define physical parameters
struct Params {
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
};

Params p = {1.55, 519, 1.05, 0.068, 67.9, 93, 93, 93, 0.84, 7, 7, 7, 1.17, 0};

// Define delay
const double tau2 = 4.0;
// const std::vector<double> Lags = {tau2};

// Define the DDE function
void cv_ddefun(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history) 
{
    double R = t <= 600 ? 1.05 : 0.21 * std::exp(600 - t) + 0.84;
    double Patau = history.at_time(t - tau2, 0);
    double Paoft = u[0];
    double Pvoft = u[1];
    double Hoft = u[2];

    du[0] = - (1 / (p.ca * R)) * Paoft + (1 / (p.ca * R)) * Pvoft + (1 / p.ca) * p.Vstr * Hoft;
    du[1] = (1 / (p.cv * R)) * Paoft - (1 / (p.cv * R) + 1 / (p.cv * p.r)) * Pvoft;
    double Ts = 1 / (1 + std::pow((Patau / p.alphas), p.betas));
    double Tp = 1 / (1 + std::pow((p.alphap / Paoft), p.betap));
    du[2] = (p.alphaH * Ts) / (1 + p.gammaH * Tp) - p.betaH * Tp;
}

// Define the history functions
double history_Pa(double t) 
{
    return 93;
}

double history_Pv(double t) 
{
    return (1 / (1 + p.R / p.r)) * 93;
}

double history_H(double t) 
{
    return (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * 93;
}

void CV_Model()
{
    std::cout << "Running DDE solver..." << std::endl;
  
    std::cout << "Running the Cardiovascular Model..." << std::endl;
    
        // Initial conditions and time span
        std::vector<double> u0 = {93, (1 / (1 + p.R / p.r)) * 93, (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * 93};
        double t_initial = 0.0;
        double t_final = 1000.0;

        // Create the DDE problem and solve it
        std::vector<std::function<double(double)>> prehistory = {history_Pa, history_Pv, history_H};
        std::vector<double> max_delays = {tau2, tau2, tau2}; // Ensure the size matches the number of equations

        DDEint_dopri_5<cv_ddefun> dde_solver(3, max_delays, prehistory);

        // Measure the time taken to solve the problem
        std::clock_t start = std::clock();
        std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 20000, ABS_TOL, REL_TOL);
        // Measure the time taken to solve the problem
        std::clock_t end = std::clock();

        double elapsed_time = (end - start) / (double) CLOCKS_PER_SEC;
        // conver to milliseconds
        elapsed_time = elapsed_time * 1000;
        std::cout << "Cardiovascular Model took " << elapsed_time << " milliseconds to solve." << std::endl;


        // Prepare data for plotting
        std::vector<double> time_points;
        std::vector<double> Pa_values;
        std::vector<double> Pv_values;
        std::vector<double> heart_rate;

        for (const auto& row : solution) 
        {
            time_points.push_back(row[0]);
            Pa_values.push_back(row[1]);
            Pv_values.push_back(row[2]);
            heart_rate.push_back(row[3]); // H(t) is the third element
        }

        // save the data to a csv file
        std::ofstream file("data/cardiovascular_model_cpp.csv");
        file << "Time,Pa,Pv,HR\n";
        for (size_t i = 0; i < time_points.size(); ++i) 
        {
            file << time_points[i] << "," << Pa_values[i] << "," << Pv_values[i] << "," << heart_rate[i] << "\n";
        }
        file.close();
}