#include "heat_equation.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>

#include <fstream>

// Function to initialize the state vector y based on the initial condition
void initialize_y(size_t nx_loc, size_t ny_loc, double dx, double dy, std::vector<double> &y) {
    // Resize y to fit the grid
    y.resize(nx_loc * ny_loc);

    // Loop over each grid point and initialize y
    for (size_t j = 0; j < ny_loc; ++j) {
        for (size_t i = 0; i < nx_loc; ++i) {
            double x = i * dx;
            double y_val = j * dy;
            size_t idx = i + j * nx_loc;
            y[idx] = std::sin(M_PI * x) * std::sin(M_PI * x) * std::sin(M_PI * y_val) * std::sin(M_PI * y_val);
        }
    }
}


void write_to_csv(const std::vector<double>& time ,const std::vector<double>& data, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return; 

    }

    file << "time,||u||_rms" << std::endl;
    for (size_t i = 0; i < time.size(); ++i) {
        file << time[i] << "," << data[i] << std::endl;
    }

    
    
    file << std::endl;

    file.close();
}

double calculate_rms(const std::vector<double>& values) {
    double sum_of_squares = 0.0;
    for (double value : values) {
        sum_of_squares += value * value;
    }
    double mean = sum_of_squares / values.size();
    return sqrt(mean);
}

//update DDE library so that it does  not return a vector of all the snapshots but only the last snapshot of the time it finished. Outside of the ring buffer we are not saving any other snapshots.
//optimize m_output array. check to see where its being updated. 
// change test case to take snapshots 1/20th for a 1sec. basically call run 20 times and get the last snapshot of each run

int main() {
    // Dummy prehistory function (returns 0.0 for all inputs)
    double (*dummy_prehist)(double) = [](double) { return 0.0; };
    std::vector<std::function<double(double)>> no_prehist(100, dummy_prehist);

    // Initial conditions for the Laplacian
    std::vector<double> init_cond_laplacian;
    size_t nx_loc = 32;
    size_t ny_loc = 32;
    const double dx = 1.0 / (nx_loc - 1);
    const double dy = 1.0 / (ny_loc - 1);
    int num_eq = nx_loc * ny_loc;

    // Initialize the initial condition vector
    std::vector<double> hist_init_cond(num_eq, 0.0);
    DDEint_dopri_5<laplacian> test_laplacian(num_eq, hist_init_cond, no_prehist);

    // Initialize the state
    initialize_y(nx_loc, ny_loc, dx, dy, init_cond_laplacian);
    
    std::vector<double> output, time;
    std::vector<std::vector<double>> snapshots;
    snapshots.reserve(20);

    auto start_time = std::chrono::high_resolution_clock::now();

    double dt = 0.05; // Time step
    double tf = 1.0; // Final time
    int num_steps = tf / dt; // Number of time steps

    for (int step = 0; step < num_steps; ++step) {
        double t_start = step * dt;
        double t_end = (step + 1) * dt;

        if (step == 0) {
            // Initialize the integrator at the first step
            output = test_laplacian.run(t_start, t_end, init_cond_laplacian, 0.005, 0.0005, 50000, 1e-10, 1e-5, false);
            snapshots.push_back(output);
            time.push_back(test_laplacian.get_t());
        } else {
            // Continue the integration for subsequent steps
           output = test_laplacian.continue_integration(t_end, 0.0005, 50000, 1e-10, 1e-5, false);
           snapshots.push_back(output);
           time.push_back(test_laplacian.get_t());
        }
    }


    auto end_time = std::chrono::high_resolution_clock::now();
    double execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << "Execution time: " << execution_time << " milliseconds\n";
    //output.insert(output.begin(), tf);
    while (!output.empty()) {
        output.pop_back();
    }
    // Calculate the rms of each snapshot
    int col_width = 25; 
    std::cout << std::left;
    std::cout << std::setw(col_width) << "time" << std::setw(col_width) << "||u||_rms " << std::endl;
    for (size_t i = 0; i < snapshots.size(); ++i) {
        double rms = calculate_rms(snapshots[i]);
        output.push_back(rms);
        std::cout << std::scientific;
        std::cout << std::setprecision(16);
        std::cout << std::setw(col_width) << time[i] << std::setw(col_width) << rms << std::endl;

        
    }
     

    
    write_to_csv(time, output, "data/heat_equation_output.csv");


    return 0;
}













