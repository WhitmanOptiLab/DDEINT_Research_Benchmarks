import numpy as np
import matplotlib.pyplot as plt
import csv

# Read the numerical solution from the provided CSV file
def read_numerical_solution(filename):
    x = []
    y = []
    u = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
            u.append(float(row[2]))
    return np.array(x), np.array(y), np.array(u)

# Read numerical solution from CSV file
x_numerical, y_numerical, u_numerical = read_numerical_solution('data/laplacian_out.csv')

# Reshape the data to match the grid size
nx_loc = 32
ny_loc = 32
x_numerical = x_numerical.reshape((nx_loc, ny_loc))
y_numerical = y_numerical.reshape((nx_loc, ny_loc))
u_numerical = u_numerical.reshape((nx_loc, ny_loc))

# Generate the grid points
xl = 0.0
xu = 1.0
yl = 0.0
yu = 1.0
x = np.linspace(xl, xu, nx_loc)
y = np.linspace(yl, yu, ny_loc)
X, Y = np.meshgrid(x, y)

# Plot the numerical solution
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, u_numerical, levels=50, cmap='viridis')
plt.colorbar(label='u(x, y, t)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Numerical Solution of 2D Heat Equation')
plt.savefig('plots/numerical_solution.png')
plt.show()


