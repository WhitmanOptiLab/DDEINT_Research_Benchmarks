import numpy as np
import matplotlib.pyplot as plt

# Load data
copasi = np.loadtxt('data/copasi_cv_model.csv', delimiter='\t', skiprows=25, usecols=(0, 1, 2, 3))
dde = np.loadtxt('data/ddeint_cv_model.csv', delimiter=',', skiprows=1)

copasi_t  = copasi[:, 0]
copasi_u  = copasi[:, 1]
copasi_u1 = copasi[:, 2]
copasi_u2 = copasi[:, 3]

ddeint_t  = dde[:, 0]
ddeint_y1 = dde[:, 1]
ddeint_y2 = dde[:, 2]
ddeint_y3 = dde[:, 3]

# Interpolate COPASI onto DDE time grid
u_interp  = np.interp(ddeint_t, copasi_t, copasi_u)
u1_interp = np.interp(ddeint_t, copasi_t, copasi_u1)
u2_interp = np.interp(ddeint_t, copasi_t, copasi_u2)

# RMSE (ROOT MEAN SQUARE ERROR)
rmse_u  = np.sqrt(np.mean((u_interp  - ddeint_y1)**2))
rmse_u1 = np.sqrt(np.mean((u1_interp - ddeint_y2)**2))
rmse_u2 = np.sqrt(np.mean((u2_interp - ddeint_y3)**2))

print(f"RMSE u: {rmse_u:.3f}" )
print(f"RMSE u1: {rmse_u1:.3f}")
print(f"RMSE u2: {rmse_u2:.3f}")

# Relative RMSE
rrmse_u  = rmse_u  / np.mean(np.abs(ddeint_y1)) * 100
rrmse_u1 = rmse_u1 / np.mean(np.abs(ddeint_y2)) * 100
rrmse_u2 = rmse_u2 / np.mean(np.abs(ddeint_y3)) * 100

print(f"Relative RMSE u : {rrmse_u:.6f}%")
print(f"Relative RMSE u1: {rrmse_u1:.6f}%")
print(f"Relative RMSE u2:  {rrmse_u2:.6f}%")



# Plot
fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

axes[0].plot(ddeint_t, ddeint_y1, label='DDE (y1)', linewidth=1)
axes[0].plot(ddeint_t, u_interp, label='COPASI (u)', linewidth=1, linestyle='-.')
axes[0].set_ylabel('Pa')
axes[0].legend()
axes[0].set_title(f'RMSE: {rmse_u:.3f} RRMSE:{rrmse_u:.4f}%')

axes[1].plot(ddeint_t, ddeint_y2, label='DDE (y2)', linewidth=1)
axes[1].plot(ddeint_t, u1_interp, label='COPASI (u1)', linewidth=1, linestyle='-.')
axes[1].set_ylabel('Pv')
axes[1].legend()
axes[1].set_title(f'RMSE: {rmse_u1:.3f} RRMSE:{rrmse_u1:.4f}%')

axes[2].plot(ddeint_t, ddeint_y3, label='DDE (y2)', linewidth=1)
axes[2].plot(ddeint_t, u2_interp, label='COPASI (u2)', linewidth=1, linestyle='-.')
axes[2].set_ylabel('H')
axes[2].legend()
axes[2].set_title(f'RMSE: {rmse_u2:.3f} RRMSE:{rrmse_u2:.4f}%')

plt.tight_layout()
plt.savefig('validation/cv_model_comparison.png', dpi=150)
plt.show()

#mask_copasi = copasi_t <= 200
#mask_dde = ddeint_t <= 200

#fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

#axes[0].plot(ddeint_t[mask_dde], ddeint_y1[mask_dde], label='DDE (y1)', linewidth=1)
#axes[0].plot(copasi_t[mask_copasi], copasi_u[mask_copasi], label='COPASI (u)', linewidth=1, linestyle='-.')
#axes[0].set_ylabel('u0')
#axes[0].legend()

#axes[1].plot(ddeint_t[mask_dde], ddeint_y2[mask_dde], label='DDE (y2)', linewidth=1)
#axes[1].plot(copasi_t[mask_copasi], copasi_u1[mask_copasi], label='COPASI (u1)', linewidth=1, linestyle='-.')
#axes[1].set_ylabel('u1')
#axes[1].legend()

#axes[2].plot(ddeint_t[mask_dde], ddeint_y3[mask_dde], label='DDE (y3)', linewidth=1)
#axes[2].plot(copasi_t[mask_copasi], copasi_u2[mask_copasi], label='COPASI (u2)', linewidth=1, linestyle='-.')
#axes[2].set_ylabel('u2')
#axes[2].set_xlabel('Time')
#axes[2].legend()

#plt.tight_layout()
#plt.savefig('validation/bc_model_comparison.png', dpi=150)
#plt.show()
