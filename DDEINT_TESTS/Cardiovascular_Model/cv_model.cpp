#include "cv_functions.hpp"
#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"

#include <iostream>
#include <ctime>
#include <fstream>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

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
    std::ofstream file("../data/cv_model_cpp.csv");
    file << "Time,Pa,Pv,HR\n";
    for (size_t i = 0; i < time_points.size(); ++i) 
    {
        // make the precision 10 decimal places
        file << std::setprecision(10) << time_points[i] << "," << Pa_values[i] << "," << Pv_values[i] << "," << heart_rate[i] << "\n";
    }
    file.close();

    return 0;
}
