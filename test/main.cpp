#include "DDEINT/dopri/ddeint_dopri_5.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "cpp_utils/plot.hpp"


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
void bc_model(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history) 
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
const double Tau = 4.0;
// const std::vector<double> Lags = {Tau};

// Define the DDE function
void ddefun(size_t num_eq, double t, std::vector<double>& u, std::vector<double>& du, History<double, double>& history) {
    double R = t <= 600 ? 1.05 : 0.21 * std::exp(600 - t) + 0.84;
    double Patau = history.at_time(t - Tau, 0);
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
double history_Pa(double t) {
    return 93;
}

double history_Pv(double t) {
    return (1 / (1 + p.R / p.r)) * 93;
}

double history_H(double t) {
    return (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * 93;
}

int main() 
{
    std::cout << "Running DDE solvers..." << std::endl;
    std::cout << "Running the Breast Cancer Model..." << std::endl;
    {
        // Initial conditions and time span
        std::vector<double> u0 = {1.0, 1.0, 1.0};
        double t_initial = 0.0;
        double t_final = 10.0;

        // Create the DDE problem and solve it
        std::vector<std::function<double(double)>> prehistory = {history_function, history_function, history_function};
        std::vector<double> max_delays = {tau, tau, tau}; // Ensure the size matches the number of equations

        DDEint_dopri_5<bc_model> dde_solver(3, max_delays, prehistory);

        std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 10000, 1e-6, 1e-6);

        // Prepare data for plotting
        std::vector<double> time_points;
        std::vector<double> u1_values;
        std::vector<double> u2_values;
        std::vector<double> u3_values;

        for (const auto& row : solution) {
            time_points.push_back(row[0]);
            u1_values.push_back(row[1]);
            u2_values.push_back(row[2]);
            u3_values.push_back(row[3]);
        }

        // Plot the solution
        // Plot plot("bc_model_cpp");
        // plot.setXlabel("Time t");
        // plot.setYlabel("Solution u");
        // plot.plot(time_points, u1_values, "data_u1.dat", "lines", "blue", "u1");
        // plot.plot(time_points, u2_values, "data_u2.dat", "lines", "red", "u2");
        // plot.plot(time_points, u3_values, "data_u3.dat", "lines", "green", "u3");
        // plot.show();

        // save the data to a csv file
        std::ofstream file("data/bc_model_cpp.csv");
        file << "Time,u1,u2,u3\n";
        for (size_t i = 0; i < time_points.size(); ++i) {
            file << time_points[i] << "," << u1_values[i] << "," << u2_values[i] << "," << u3_values[i] << "\n";
        }
        file.close();

    }

    std::cout << "Running the Cardiovascular Model..." << std::endl;
    {
        // Initial conditions and time span
        std::vector<double> u0 = {93, (1 / (1 + p.R / p.r)) * 93, (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * 93};
        double t_initial = 0.0;
        double t_final = 1000.0;

        // Create the DDE problem and solve it
        std::vector<std::function<double(double)>> prehistory = {history_Pa, history_Pv, history_H};
        std::vector<double> max_delays = {Tau, Tau, Tau}; // Ensure the size matches the number of equations

        DDEint_dopri_5<ddefun> dde_solver(3, max_delays, prehistory);

        std::vector<std::vector<double>> solution = dde_solver.run(t_initial, t_final, u0, 0.1, 1e-5, 10000, 1e-6, 1e-6);

        // Prepare data for plotting
        std::vector<double> time_points;
        std::vector<double> heart_rate;

        for (const auto& row : solution) {
            time_points.push_back(row[0]);
            heart_rate.push_back(row[3]); // H(t) is the third element
        }

        // Plot the solution
        // Plot plot("cardiovascular_model_cpp");
        // plot.setXlabel("Time t");
        // plot.setYlabel("Heart Rate H(t)");
        // plot.plot(time_points, heart_rate, "data_heart_rate.dat", "lines", "blue", "Heart Rate");
        // plot.show();

        // save the data to a csv file
        std::ofstream file("data/cardiovascular_model_cpp.csv");
        file << "Time,Heart Rate\n";
        for (size_t i = 0; i < time_points.size(); ++i) {
            file << time_points[i] << "," << heart_rate[i] << "\n";
        }
        file.close();
    }
    return 0;
}