
#include "laplacian.hpp"
#include <iomanip>

// 2d heat eqaution
void laplacian(size_t n, double t, std::vector<double> &y, std::vector<double> &dydt, History<double, double>& hist) {
    
    // Check for size mismatch
    if (n != y.size() || n != dydt.size()) {
        throw std::runtime_error("Size mismatch in laplacian function");
    }
    //Boundary conditions
    double xl = 0.0;
    double xu = 1.0;
    double yl = 0.0;
    double yu = 1.0;

    // Diffusion coefficients in the x and y directions
    double kx = 1.0; 
    double ky = 1.0; 

    //grid size
    size_t nx_loc = 32; 
    size_t ny_loc = 32; 


    // Grid spacing in the x and y directions
    const double dx = xu / (nx_loc - 1); 
    const double dy = yu / (ny_loc - 1); 


    // diffusion coefficients in the x and y directions
    const double cx = kx / (dx * dx); 
    const double cy = ky / (dy * dy);
    const double cc = -2.0 * (cx + cy);

    // Initialize dydt to zero
    std::fill(dydt.begin(), dydt.end(), 0.0);

    // precompute time-dependent terms
    double sin_t_cos_t = sin(M_PI * t) * cos(M_PI * t);
    double cos_sqr_t = cos(M_PI * t) * cos(M_PI * t);

    for (size_t j = 1; j < ny_loc - 1; j++) {
        for (size_t i = 1; i < nx_loc - 1; i++) {
            size_t idx = i + j * nx_loc; // Index of current grid point

            //coordinates of current grid point
            double x = i * dx;
            double y_val = j * dy; // 'y' renamed to 'y_val' to avoid name clash with input vector

            double sin_sqr_x = sin(M_PI * x) * sin(M_PI * x);
            double sin_sqr_y = sin(M_PI * y_val) * sin(M_PI * y_val);
            double cos_sqr_x = cos(M_PI * x) * cos(M_PI * x);
            double cos_sqr_y = cos(M_PI * y_val) * cos(M_PI * y_val);

            // Compute the forcing term
            dydt[idx] = -2.0 * M_PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t -
                        kx * 2.0 * M_PI * M_PI * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t -
                        ky * 2.0 * M_PI * M_PI * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;

            // Add Diffusion term
            dydt[idx] += cc * y[idx] +
                         cx * (y[idx - 1] + y[idx + 1]) +
                         cy * (y[idx - nx_loc] + y[idx + nx_loc]);
        }
    }
    // Set boundary conditions
    for (size_t j = 0; j < ny_loc; j++) {
        dydt[j * nx_loc] = 0.0;               // Left boundary (x = 0)
        dydt[(j + 1) * nx_loc - 1] = 0.0;     // Right boundary (x = x_max)
    }
    for (size_t i = 0; i < nx_loc; i++) {
        dydt[i] = 0.0;                        // Bottom boundary (y = 0)
        dydt[(ny_loc - 1) * nx_loc + i] = 0.0; // Top boundary (y = y_max)
    }
}