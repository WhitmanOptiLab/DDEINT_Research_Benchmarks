import csv 
import numpy as np

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
_, _, u_numerical = read_numerical_solution('data/laplacian_out.csv')

def read_analytical_solution(filename):
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

# Read analytical solution from CSV file
_, _, u_analytical = read_analytical_solution('data/analytical_solution.csv')

# Compute error
error = np.mean(np.abs(u_analytical - u_numerical))
print(f'Abs Error: {error}')
