#ifndef HEAT_EQUATION_HPP
#define HEAT_EQUATION_HPP

#include <vector>
#include <cmath>
#include "../../DDEINT/dopri/DOPRI_5.hpp"  




// Forward declaration of the History class
template <typename S, typename T> class History;

// // Struct to hold problem parameters
// struct LaplacianParams {
//     double xl = 0.0;
//     double xu = 1.0;
//     double yl = 0.0;
//     double yu = 1.0;
//     double kx = 1.0; 
//     double ky = 1.0; 
//     size_t nx_loc = 32; 
//     size_t ny_loc = 32; 
// };

// Function declarations
void laplacian(size_t n, double t, std::vector<double>& y, std::vector<double>& dydt, History<double, double>& hist);
void initialize_y(size_t nx_loc, size_t ny_loc, double dx, double dy, std::vector<double>& y);
void print_statistics(const std::vector<std::vector<double>>& results, double execution_time);

#endif // HEAT_EQUATION_HPP
