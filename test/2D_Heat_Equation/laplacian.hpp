#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include <vector>
#include <cmath>
#include <iostream>

#include "../../DDEINT/dopri/ddeint_dopri_5.hpp"  

template <typename S, typename T> class History;

void laplacian(size_t n, double t, std::vector<double>& y, std::vector<double>& dydt, History<double, double>& hist);

#endif // LAPLACIAN_HPP
