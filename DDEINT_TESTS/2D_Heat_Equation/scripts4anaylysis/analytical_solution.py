import numpy as np
import matplotlib.pyplot as plt
import csv

# Boundary conditions
xl = 0.0
xu = 1.0
yl = 0.0
yu = 1.0

# Diffusion coefficients in the x and y directions
kx = 1.0
ky = 1.0

# Grid size
nx_loc = 32
ny_loc = 32

# Grid spacing in the x and y directions
dx = xu / (nx_loc - 1)
dy = yu / (ny_loc - 1)

# Diffusion coefficients in the x and y directions
cx = kx / (dx * dx)
cy = ky / (dy * dy)
cc = -2.0 * (cx + cy)

# Define the time at which to evaluate the solution
t = 1

# Generate the grid points
x = np.linspace(xl, xu, nx_loc)
y = np.linspace(yl, yu, ny_loc)
X, Y = np.meshgrid(x, y)

# Compute the analytical solution
PI = np.pi
cos_sqr_t = np.cos(PI * t) ** 2
sin_sqr_x = np.sin(PI * X) ** 2
sin_sqr_y = np.sin(PI * Y) ** 2
u = sin_sqr_x * sin_sqr_y * cos_sqr_t





# Plot the analytical solution
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, u, levels=50, cmap='viridis')
plt.colorbar(label='u(x, y, t)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Analytical Solution of 2D Heat Equation at t={t}')
plt.savefig('plots/analytical_solution.png')
plt.show()


# Save the analytical solution to a CSV file
with open('data/analytical_solution.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y', 'u'])
    for i in range(nx_loc):
        for j in range(ny_loc):
            writer.writerow([X[i, j], Y[i, j], u[i, j]])

print("Analytical solution saved to analytical_solution.csv")