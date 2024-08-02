#include "gi_functions.hpp"
#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

#include <iostream>
#include <ctime>
#include <fstream>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

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

    // Prepare data for plotting
    std::vector<double> time_points;
    std::vector<double> glucose_values;
    std::vector<double> insulin_values;

    for (const auto& row : solution) 
    {
        time_points.push_back(row[0]);
        glucose_values.push_back(row[1]);
        insulin_values.push_back(row[2]); // G(t) and I(t)
    }

    // Save the data to a CSV file
    std::ofstream file("../data/gi_model_cpp.csv");
    file << "Time,Glucose,Insulin\n";
    for (size_t i = 0; i < time_points.size(); ++i) 
    {
        file << std::setprecision(10) <<  time_points[i] << "," << glucose_values[i] << "," << insulin_values[i] << "\n";
    }
    file.close();

    return 0;
}
