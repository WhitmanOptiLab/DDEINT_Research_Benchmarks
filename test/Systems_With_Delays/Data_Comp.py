import pandas as pd
import matplotlib.pyplot as plt

bc_cpp_sol = pd.read_csv('data/bc_model_cpp.csv')
bc_jl_sol = pd.read_csv('data/bc_model_jl.csv')

cv_cpp_sol = pd.read_csv('data/cardiovascular_model_cpp.csv')
cv_jl_sol = pd.read_csv('data/cardiovascular_model_jl.csv')

repressilator_cpp_sol = pd.read_csv('data/repressilator_model_cpp.csv')
repressilator_jl_sol = pd.read_csv('data/repressilator_model_jl.csv')

# sir_cpp_sol = pd.read_csv('data/sir_model_cpp.csv')
# sir_jl_sol = pd.read_csv('data/sir_model_jl.csv')

# create two plots with two subplots each
bc_fig, bc_axs = plt.subplots(1, 2, figsize=(20, 10))
cv_fig, cv_axs = plt.subplots(1, 2, figsize=(20, 10))
rm_fig, rm_axs = plt.subplots(1, 2, figsize=(20, 10))
# sir_fig, sir_axs = plt.subplots(1, 2, figsize=(20, 10))

# remove side margins
bc_fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
cv_fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
rm_fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
# sir_fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

# plot the results
# bc_model_cpp.csv has Time, u1, u2, u3 
# bc_model_jl.csv has timestap, value1, value2, value3
bc_axs[0].plot(bc_cpp_sol['Time'], bc_cpp_sol['u1'], label='u1')
bc_axs[0].plot(bc_cpp_sol['Time'], bc_cpp_sol['u2'], label='u2')
bc_axs[0].plot(bc_cpp_sol['Time'], bc_cpp_sol['u3'], label='u3')
bc_axs[0].set_title('DDEInt Breast Cancer Model Solution')
bc_axs[0].set_xlabel('Time')
bc_axs[0].set_ylabel('Values')
bc_axs[0].legend()

bc_axs[1].plot(bc_jl_sol['time'], bc_jl_sol['u1'], label='u1')
bc_axs[1].plot(bc_jl_sol['time'], bc_jl_sol['u2'], label='u2')
bc_axs[1].plot(bc_jl_sol['time'], bc_jl_sol['u3'], label='u3')
bc_axs[1].set_title('DifferenitalEquations.jl Breast Cancer Model Solution')
bc_axs[1].set_xlabel('Time')
bc_axs[1].set_ylabel('Values')
bc_axs[1].legend()

# cv_model_cpp.csv has Time, Heart Rate
# cv_model_jl.csv has timestamp, value1, value2, value3
cv_axs[0].plot(cv_cpp_sol['Time'], cv_cpp_sol['HR'], label='Heart Rate')
cv_axs[0].set_title('DDEInt Cardiovascular Model Solution')
cv_axs[0].set_xlabel('Time')
cv_axs[0].set_ylabel('Heart Rate')
cv_axs[0].legend()

cv_axs[1].plot(cv_jl_sol['time'], cv_jl_sol['HR'], label='Heart Rate')
cv_axs[1].set_title('DifferenitalEquations.jl Cardiovascular Model Solution')
cv_axs[1].set_xlabel('Time')
cv_axs[1].set_ylabel('Heart Rate')
cv_axs[1].legend()

# repressilator_model_cpp.csv has Time, x1, x2, x3
# repressilator_model_jl.csv has timestamp, value1, value2, value3
rm_axs[0].plot(repressilator_cpp_sol['Time'], repressilator_cpp_sol['x1'], label='x1')
rm_axs[0].plot(repressilator_cpp_sol['Time'], repressilator_cpp_sol['x2'], label='x2')
rm_axs[0].plot(repressilator_cpp_sol['Time'], repressilator_cpp_sol['x3'], label='x3')
rm_axs[0].set_title('DDEInt Repressilator Model Solution')
rm_axs[0].set_xlabel('Time')
rm_axs[0].set_ylabel('Values')
rm_axs[0].legend()

rm_axs[1].plot(repressilator_jl_sol['time'], repressilator_jl_sol['x1'], label='x1')
rm_axs[1].plot(repressilator_jl_sol['time'], repressilator_jl_sol['x2'], label='x2')
rm_axs[1].plot(repressilator_jl_sol['time'], repressilator_jl_sol['x3'], label='x3')
rm_axs[1].set_title('DifferenitalEquations.jl Repressilator Model Solution')
rm_axs[1].set_xlabel('Time')
rm_axs[1].set_ylabel('Values')
rm_axs[1].legend()

# sir_model_cpp.csv has Time, S, I, R
# sir_model_jl.csv has timestamp, S, I, R
# sir_axs[0].plot(sir_cpp_sol['Time'], sir_cpp_sol['S'], label='S')
# sir_axs[0].plot(sir_cpp_sol['Time'], sir_cpp_sol['I'], label='I')
# sir_axs[0].plot(sir_cpp_sol['Time'], sir_cpp_sol['R'], label='R')
# sir_axs[0].set_title('DDEInt SIR Model Solution')
# sir_axs[0].set_xlabel('Time')
# sir_axs[0].set_ylabel('Values')
# sir_axs[0].legend()

# sir_axs[1].plot(sir_jl_sol['time'], sir_jl_sol['S'], label='S')
# sir_axs[1].plot(sir_jl_sol['time'], sir_jl_sol['I'], label='I')
# sir_axs[1].plot(sir_jl_sol['time'], sir_jl_sol['R'], label='R')
# sir_axs[1].set_title('DifferenitalEquations.jl SIR Model Solution')
# sir_axs[1].set_xlabel('Time')
# sir_axs[1].set_ylabel('Values')
# sir_axs[1].legend()

# save the plots
bc_fig.savefig('plots/bc_model_comp.png')
cv_fig.savefig('plots/cv_model_comp.png')
rm_fig.savefig('plots/repressilator_model_comp.png')
# sir_fig.savefig('plots/sir_model_comp.png')