#include "cv_functions.hpp"
#include "../../DDEINT/Methods/Dormand_Prince/DoPri_5.hpp"

#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <filesystem>

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

    DoPri_5<cv_dde> dde_solver(3,max_delays, prehistory);
    dde_solver.initialize(0, 0.1, 1e-5, u0, ABS_TOL, REL_TOL, false, true);
    Results results = dde_solver.solve(t_initial, t_final, 500, 20000);


    // Prepare data for plotting
    std::vector<double> time_points;
    std::vector<double> Pa_values;
    std::vector<double> Pv_values;
    std::vector<double> heart_rate;

    for (size_t i = 0; i < results.times.size(); i++)
    {
        time_points.push_back(results.times[i]);
        Pa_values.push_back(results.solutions[i][0]);
        Pv_values.push_back(results.solutions[i][1]);
        heart_rate.push_back(results.solutions[i][2]); // H(t) is the third element
    }

    // save the data to a csv file
    std::filesystem::path data_path = std::filesystem::canonical("/proc/self/exe").parent_path().parent_path() / "data" / "csv_files" / "cv_model.csv";
    std::filesystem::create_directories(data_path.parent_path());
    std::ofstream file(data_path);
    file << "Time,Pa,Pv,HR\n";
    for (size_t i = 0; i < time_points.size(); ++i) 
    {
        // make the precision 10 decimal places
        file << std::setprecision(10) << time_points[i] << "," << Pa_values[i] << "," << Pv_values[i] << "," << heart_rate[i] << "\n";
    }
    file.close();

    return 0;
}
