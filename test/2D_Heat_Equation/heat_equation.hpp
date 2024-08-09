#ifndef HEAT_EQUATION_HPP
#define HEAT_EQUATION_HPP

#include <vector>
#include <cmath>
#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"  
#include "laplacian.hpp"



// Function declarations
void laplacian(size_t n, double t, std::vector<double>& y, std::vector<double>& dydt, History<double, double>& hist);
void initialize_y(size_t nx_loc, size_t ny_loc, double dx, double dy, std::vector<double>& y);
void print_statistics(const std::vector<std::vector<double>>& results, double execution_time);

#endif // HEAT_EQUATION_HPP
