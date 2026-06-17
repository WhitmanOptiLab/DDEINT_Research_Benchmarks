#include "../../DDEINT/Methods/Dormand_Prince/DoPri_5.hpp"
#include "bc_functions.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <filesystem>

#define ABS_TOL 1e-9
#define REL_TOL 1e-9

int main()
{
    std::cout << "Running DDE solver..." << std::endl;
    std::cout << "Running the Breast Cancer Model..." << std::endl;
    
    // Initial conditions and time span
    std::vector<double> u0 = {1.0, 1.0, 1.0};
    double t_initial = 0.0;
    double t_final = 10.0;

    // Create the DDE problem and solve it
    std::vector<std::function<double(double)>> prehistory = {history_bc, history_bc, history_bc};
    std::vector<double> max_delays = {1.0, 1.0, 1.0}; // Ensure the size matches the number of equations

    DoPri_5<bc_dde> dde_solver(3, max_delays, prehistory);
    dde_solver.initialize(0, 0.1, 1e-5, u0, ABS_TOL, REL_TOL, false, true);
    Results results = dde_solver.solve(t_initial, t_final, 500, 10000);

    // Prepare data for plotting
    std::vector<double> time_points;
    std::vector<double> u1_values;
    std::vector<double> u2_values;
    std::vector<double> u3_values;

    for (size_t i; i < results.times.size(); i++)
    {
        time_points.push_back(results.times[i]);
        u1_values.push_back(results.solutions[i][0]);
        u2_values.push_back(results.solutions[i][1]);
        u3_values.push_back(results.solutions[i][2]);
    }

    // save the data to a csv file
    std::filesystem::path data_path = std::filesystem::canonical("/proc/self/exe").parent_path().parent_path() / "data" / "csv_files" / "bc_model_cpp.csv";
    std::filesystem::create_directories(data_path.parent_path());
    std::ofstream file(data_path);
    file << "Time,u1,u2,u3\n";
    for (size_t i = 0; i < time_points.size(); ++i) 
    {
        file << std::setprecision(10) << time_points[i] << "," << u1_values[i] << "," << u2_values[i] << "," << u3_values[i] << "\n";
    }
    file.close();

    return 0;
}
