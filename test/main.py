import pandas as pd
import matplotlib.pyplot as plt

bc_cpp_sol = pd.read_csv('data/bc_model_cpp.csv')
bc_jl_sol = pd.read_csv('data/bc_model_jl.csv')

cv_cpp_sol = pd.read_csv('data/cardiovascular_model_cpp.csv')
cv_jl_sol = pd.read_csv('data/cardiovascular_model_jl.csv')

# create two plots with two subplots each
bc_fig, bc_axs = plt.subplots(1, 2, figsize=(20, 10))
cv_fig, cv_axs = plt.subplots(1, 2, figsize=(20, 10))

# remove side margins
bc_fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
cv_fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

# plot the results
# bc_model_cpp.csv has Time, u1, u2, u3 
# bc_model_jl.csv has timestap, value1, value2, value3
bc_axs[0].plot(bc_cpp_sol['Time'], bc_cpp_sol['u1'], label='u1')
bc_axs[0].plot(bc_cpp_sol['Time'], bc_cpp_sol['u2'], label='u2')
bc_axs[0].plot(bc_cpp_sol['Time'], bc_cpp_sol['u3'], label='u3')
bc_axs[0].set_title('DDEInt BC Model Solution')
bc_axs[0].set_xlabel('Time')
bc_axs[0].set_ylabel('Values')
bc_axs[0].legend()

bc_axs[1].plot(bc_jl_sol['timestamp'], bc_jl_sol['value1'], label='u1')
bc_axs[1].plot(bc_jl_sol['timestamp'], bc_jl_sol['value2'], label='u2')
bc_axs[1].plot(bc_jl_sol['timestamp'], bc_jl_sol['value3'], label='u3')
bc_axs[1].set_title('DifferenitalEquations.jl BC Model Solution')
bc_axs[1].set_xlabel('Time')
bc_axs[1].set_ylabel('Values')
bc_axs[1].legend()

# cv_model_cpp.csv has Time, Heart Rate
# cv_model_jl.csv has timestamp, value1, value2, value3
cv_axs[0].plot(cv_cpp_sol['Time'], cv_cpp_sol['Heart Rate'], label='Heart Rate')
cv_axs[0].set_title('DDEInt Cardiovascular Model Solution')
cv_axs[0].set_xlabel('Time')
cv_axs[0].set_ylabel('Heart Rate')
cv_axs[0].legend()

cv_axs[1].plot(cv_jl_sol['timestamp'], cv_jl_sol['value3'], label='Heart Rate')
cv_axs[1].set_title('DifferenitalEquations.jl Cardiovascular Model Solution')
cv_axs[1].set_xlabel('Time')
cv_axs[1].set_ylabel('Heart Rate')
cv_axs[1].legend()

# save the plots
bc_fig.savefig('plots/bc_model_comp.png')
cv_fig.savefig('plots/cv_model_comp.png')
