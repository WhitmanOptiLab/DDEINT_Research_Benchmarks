import numpy as np
import matplotlib.pyplot as plt

# Data, csv files
bc_model_lsoda = np.loadtxt('data/bc_model_lsoda.csv', delimiter='\t', skiprows=25, usecols=(0, 1, 2, 3)) # Breast Cancer Model
cv_model_lsoda = np.loadtxt('data/cv_model_lsoda.csv', delimiter='\t', skiprows=25, usecols=(0, 1, 2, 3)) # Cardiovascular Model
cw_model_lsoda = np.loadtxt('data/cw_model_lsoda.csv', delimiter='\t', skiprows=25, usecols=(0, 1))       # Cobweb Model

bc_model_ddeint = np.loadtxt('data/bc_model_ddeint.csv', delimiter='\t', skiprows=25, usecols=(0, 1, 2, 3)) # Breast Cancer Model
cv_model_ddeint = np.loadtxt('data/cv_model_ddeint.csv', delimiter='\t', skiprows=25, usecols=(0, 1, 2, 3)) # Cardiovascular Model
cw_model_ddeint = np.loadtxt('data/cw_model_ddeint.csv', delimiter='\t', skiprows=25, usecols=(0, 1))       # Cobweb Model

# Retrieving time, derivative(s)
# Breast Cancer Model -- data --
bc_lsoda_t = bc_model_lsoda[:, 0]
bc_lsoda_u = bc_model_lsoda[:, 1]
bc_lsoda_u1 = bc_model_lsoda[:, 2]
bc_lsoda_u2 = bc_model_lsoda[:, 3]
bc_ddeint_t = bc_model_ddeint[:, 0]
bc_ddeint_u = bc_model_ddeint[:, 1]
bc_ddeint_u1 = bc_model_ddeint[:, 2]
bc_ddeint_u2 = bc_model_ddeint[:, 3]
# Cardiovascular Model -- data --
cv_lsoda_t = cv_model_lsoda[:, 0]
cv_lsoda_u = cv_model_lsoda[:, 1]
cv_lsoda_u1 = cv_model_lsoda[:, 2]
cv_lsoda_u2 = cv_model_lsoda[:, 3]
cv_ddeint_t = cv_model_ddeint[:, 0]
cv_ddeint_u = cv_model_ddeint[:, 1]
cv_ddeint_u1 = cv_model_ddeint[:, 2]
cv_ddeint_u2 = cv_model_ddeint[:, 3]
# Cobweb Model --data--
cw_lsoda_t = cw_model_lsoda[:, 0]
cw_lsoda_u = cw_model_lsoda[:, 1]
cw_ddeint_t = cw_model_ddeint[:, 0]
cw_ddeint_u = cw_model_ddeint[:, 1]

# Interpolation DDEINT_Method onto LSODA_METHOD time grid
# BC --interpolation--
bc_u_interp = np.interp(bc_ddeint_t, bc_lsoda_t, bc_lsoda_u)
bc_u1_interp = np.interp(bc_ddeint_t, bc_lsoda_t, bc_lsoda_u1)
bc_u2_interp = np.interp(bc_ddeint_t, bc_lsoda_t, bc_lsoda_u2)
# CV --interpolation--
cv_u_interp = np.interp(cv_ddeint_t, cv_lsoda_t, cv_lsoda_u)
cv_u1_interp = np.interp(cv_ddeint_t, cv_lsoda_t, cv_lsoda_u1)
cv_u2_interp = np.interp(cv_ddeint_t, cv_lsoda_t, cv_lsoda_u2)
# CW --interpolation--
cw_u_interp = np.interp(cw_ddeint_t, cw_lsoda_t, cw_lsoda_u)

# RMSE (root mean square error) and RRMSE (relative root mean square error)
# BC --rmse-- and --rrmse--
bc_rmse_u = np.sqrt(np.mean((bc_u_interp - bc_ddeint_u)**2))
bc_rmse_u1 = np.sqrt(np.mean((bc_u1_interp - bc_ddeint_u1)**2))
bc_rmse_u2 = np.sqrt(np.mean((bc_u2_interp - bc_ddeint_u2)**2))
bc_rrmse_u = bc_rmse_u / np.mean(np.abs(bc_ddeint_u)) * 100 # percentage
bc_rrmse_u1 = bc_rmse_u1 / np.mean(np.abs(bc_ddeint_u1)) * 100 
bc_rrmse_u2 = bc_rmse_u2 / np.mean(np.abs(bc_ddeint_u2)) * 100 
# CV --rmse-- and --rrmse--
cv_rmse_u = np.sqrt(np.mean((cv_u_interp - cv_ddeint_u)**2))
cv_rmse_u1 = np.sqrt(np.mean((cv_u1_interp - cv_ddeint_u1)**2))
cv_rmse_u2 = np.sqrt(np.mean((cv_u2_interp - cv_ddeint_u2)**2))
cv_rrmse_u = cv_rmse_u / np.mean(np.abs(cv_ddeint_u)) * 100 # percentage
cv_rrmse_u1 = cv_rmse_u1 / np.mean(np.abs(cv_ddeint_u1)) * 100 
cv_rrmse_u2 = cv_rmse_u2 / np.mean(np.abs(cv_ddeint_u2)) * 100 
# CW --rmse-- and --rrmse--
cw_rmse_u = np.sqrt(np.mean((cw_u_interp  - cw_ddeint_u)**2))
cw_rrmse_u = cw_rmse_u / np.mean(np.abs(cw_ddeint_u)) * 100 # percentage

# Test Printing 
print("------Breast Cancer Model------")
print(f"RMSE u: {bc_rmse_u:.6f}" )
print(f"RMSE u1: {bc_rmse_u1:.6f}")
print(f"RMSE u2: {bc_rmse_u2:.6f}")
print(f"Relative RMSE u : {bc_rrmse_u:.4f}%")
print(f"Relative RMSE u1: {bc_rrmse_u1:.4f}%")
print(f"Relative RMSE u2:  {bc_rrmse_u2:.4f}%")
print("------Cardiovascular Model------")
print(f"RMSE u: {cv_rmse_u:.6f}" )
print(f"RMSE u1: {cv_rmse_u1:.6f}")
print(f"RMSE u2: {cv_rmse_u2:.6f}")
print(f"Relative RMSE u : {cv_rrmse_u:.4f}%")
print(f"Relative RMSE u1: {cv_rrmse_u1:.4f}%")
print(f"Relative RMSE u2:  {cv_rrmse_u2:.4f}%")
print("------Cobweb Model------")
print(f"RMSE u: {cw_rmse_u:.6f}" )
print(f"Relative RMSE u : {cw_rrmse_u:.4f}%")

# Plotting
# Plotting Breast Cancer Model
fig_bc, axes_bc = plt.subplots(3,1, figsize=(8,6), sharex=True)

axes_bc[0].plot(bc_ddeint_t, bc_ddeint_u, label='DDEINT (u)', linewidth=1)
axes_bc[0].plot(bc_ddeint_t, bc_u_interp, label='COPASI (u)', linewidth=1, linestyle='-.')
axes_bc[0].set_ylabel('u')
axes_bc[0].legend()
axes_bc[0].set_title(f'RMSE: {bc_rmse_u:.6f} RRMSE: {bc_rrmse_u:.4f}%')

axes_bc[1].plot(bc_ddeint_t, bc_ddeint_u1, label='DDEINT (u1)', linewidth=1)
axes_bc[1].plot(bc_ddeint_t, bc_u1_interp, label='COPASI (u1)', linewidth=1, linestyle='-.')
axes_bc[1].set_ylabel('u1')
axes_bc[1].legend()
axes_bc[1].set_title(f'RMSE: {bc_rmse_u1:.6f} RRMSE: {bc_rrmse_u1:.4f}%')

axes_bc[2].plot(bc_ddeint_t, bc_ddeint_u2, label='DDEINT (u2)', linewidth=1)
axes_bc[2].plot(bc_ddeint_t, bc_u2_interp, label='COPASI (u2)', linewidth=1, linestyle='-.')
axes_bc[2].set_ylabel('u2')
axes_bc[2].legend()
axes_bc[2].set_title(f'RMSE: {bc_rmse_u2:.6f} RRMSE: {bc_rrmse_u2:.4f}%')

fig_bc.tight_layout()
fig_bc.savefig('validation/breast_cancer_model_comparison.png', dpi=150)
fig_bc.show()

# Plotting Cardiovascular Model
fig_cv, axes_cv = plt.subplots(3,1, figsize=(8,6), sharex=True)

axes_cv[0].plot(cv_ddeint_t, cv_ddeint_u, label='DDEINT (u)', linewidth=1)
axes_cv[0].plot(cv_ddeint_t, cv_u_interp, label='COPASI (u)', linewidth=1, linestyle='-.')
axes_cv[0].set_ylabel('u')
axes_cv[0].legend()
axes_cv[0].set_title(f'RMSE: {cv_rmse_u:.6f} RRMSE: {cv_rrmse_u:.4f}%')

axes_cv[1].plot(cv_ddeint_t, cv_ddeint_u1, label='DDEINT (u1)', linewidth=1)
axes_cv[1].plot(cv_ddeint_t, cv_u1_interp, label='COPASI (u1)', linewidth=1, linestyle='-.')
axes_cv[1].set_ylabel('u1')
axes_cv[1].legend()
axes_cv[1].set_title(f'RMSE: {cv_rmse_u1:.6f} RRMSE: {cv_rrmse_u1:.4f}%')

axes_cv[2].plot(cv_ddeint_t, cv_ddeint_u2, label='DDEINT (u2)', linewidth=1)
axes_cv[2].plot(cv_ddeint_t, cv_u2_interp, label='COPASI (u2)', linewidth=1, linestyle='-.')
axes_cv[2].set_ylabel('u2')
axes_cv[2].legend()
axes_cv[2].set_title(f'RMSE: {cv_rmse_u2:.6f} RRMSE: {cv_rrmse_u2:.4f}%')

fig_cv.tight_layout()
fig_cv.savefig('validation/cardiovascular_model_comparison.png', dpi=150)

# Plotting Cobweb Model
fig_cw, axes_cw = plt.subplots(1,1, figsize=(8,6), sharex=True)

axes_cw.plot(cw_ddeint_t, cw_ddeint_u, label='DDEINT (u)', linewidth=1)
axes_cw.plot(cw_ddeint_t, cw_u_interp, label='COPASI (u)', linewidth=1, linestyle='-.')
axes_cw.set_ylabel('u')
axes_cw.legend()
axes_cw.set_title(f'RMSE: {cw_rmse_u:.6f} RRMSE: {cw_rrmse_u:.4f}%')

fig_cw.tight_layout()
fig_cw.savefig('validation/cobweb_model_comparison.png', dpi=150)


plt.show()