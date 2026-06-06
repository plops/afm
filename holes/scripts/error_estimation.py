"""
Purpose: Perform bootstrap error estimation on calibrated hole coordinates
and CMM geometric parameters. Uses residual bootstrapping to calculate
parameter uncertainties and analyzes multivariate interdependencies.
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import lsqr
import os

# Set style for professional look
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

# Load processed dataset
ds = xr.open_dataset("holes/processed_drill_data.nc")
block_names = ds['block'].values
n_blocks = len(block_names)
u_vals = ds['u'].values
v_vals = ds['v'].values

# Extract coordinates
v_grid, u_grid = np.meshgrid(ds['v'].values, ds['u'].values, indexing='ij')
u_flat = u_grid.flatten()
v_flat = v_grid.flatten()
N = len(u_flat) # 943

# Build global sparse system
n_block_params = 2 * N + 10
total_params = n_blocks * n_block_params + 3
total_eqs = n_blocks * 4 * N + n_blocks * 2 * N

A_mat = lil_matrix((total_eqs, total_params))
y_val = np.zeros(total_eqs)

def get_param_idx(b_idx, p_offset):
    return b_idx * n_block_params + p_offset

idx_sx = n_blocks * n_block_params
idx_sy = n_blocks * n_block_params + 1
idx_alpha = n_blocks * n_block_params + 2

print("Setting up global sparse system...")
for b_idx, b in enumerate(block_names):
    dx_u = ds['dx_unrot'].sel(block=b).values.flatten()
    dy_u = ds['dy_unrot'].sel(block=b).values.flatten()
    dx_r = ds['dx_rot'].sel(block=b).values.flatten()
    dy_r = ds['dy_rot'].sel(block=b).values.flatten()

    t_u_raw = ds['time_unrot'].sel(block=b).values.flatten()
    t_r_raw = ds['time_rot'].sel(block=b).values.flatten()
    t_u = (t_u_raw - t_u_raw.min()) / np.timedelta64(1, 's')
    t_r = (t_r_raw - t_r_raw.min()) / np.timedelta64(1, 's')

    for i in range(N):
        u = u_flat[i]
        v = v_flat[i]
        tu = t_u[i]
        tr = t_r[i]
        
        # eq 1: dx_u[i]
        A_mat[4*i + b_idx*4*N, get_param_idx(b_idx, i)] = 1.0
        A_mat[4*i + b_idx*4*N, idx_sx] = u
        A_mat[4*i + b_idx*4*N, get_param_idx(b_idx, 2*N)] = -v
        A_mat[4*i + b_idx*4*N, get_param_idx(b_idx, 2*N+1)] = 1.0
        A_mat[4*i + b_idx*4*N, get_param_idx(b_idx, 2*N+6)] = tu
        y_val[4*i + b_idx*4*N] = dx_u[i]
        
        # eq 2: dy_u[i]
        A_mat[4*i+1 + b_idx*4*N, get_param_idx(b_idx, N + i)] = 1.0
        A_mat[4*i+1 + b_idx*4*N, idx_sy] = v
        A_mat[4*i+1 + b_idx*4*N, idx_alpha] = u
        A_mat[4*i+1 + b_idx*4*N, get_param_idx(b_idx, 2*N)] = u
        A_mat[4*i+1 + b_idx*4*N, get_param_idx(b_idx, 2*N+2)] = 1.0
        A_mat[4*i+1 + b_idx*4*N, get_param_idx(b_idx, 2*N+7)] = tu
        y_val[4*i+1 + b_idx*4*N] = dy_u[i]
        
        # eq 3: dx_r[i]
        A_mat[4*i+2 + b_idx*4*N, get_param_idx(b_idx, N + i)] = 1.0
        A_mat[4*i+2 + b_idx*4*N, idx_sx] = v
        A_mat[4*i+2 + b_idx*4*N, get_param_idx(b_idx, 2*N+3)] = u
        A_mat[4*i+2 + b_idx*4*N, get_param_idx(b_idx, 2*N+4)] = 1.0
        A_mat[4*i+2 + b_idx*4*N, get_param_idx(b_idx, 2*N+8)] = tr
        y_val[4*i+2 + b_idx*4*N] = dx_r[i]
        
        # eq 4: dy_r[i]
        A_mat[4*i+3 + b_idx*4*N, get_param_idx(b_idx, i)] = -1.0
        A_mat[4*i+3 + b_idx*4*N, idx_sy] = -u
        A_mat[4*i+3 + b_idx*4*N, idx_alpha] = v
        A_mat[4*i+3 + b_idx*4*N, get_param_idx(b_idx, 2*N+3)] = v
        A_mat[4*i+3 + b_idx*4*N, get_param_idx(b_idx, 2*N+5)] = 1.0
        A_mat[4*i+3 + b_idx*4*N, get_param_idx(b_idx, 2*N+9)] = tr
        y_val[4*i+3 + b_idx*4*N] = dy_r[i]

# Add regularization equations
eq_idx = n_blocks * 4 * N
lambda_reg = 1e-6
for b_idx in range(n_blocks):
    for i in range(2 * N):
        A_mat[eq_idx, get_param_idx(b_idx, i)] = lambda_reg
        y_val[eq_idx] = 0.0
        eq_idx += 1

A_mat = A_mat.tocsr()

# Solve original system
print("Solving original system...")
res_lsqr = lsqr(A_mat, y_val, damp=0.0)
x_sol = res_lsqr[0]

# Compute fitted values and residuals
y_fit = A_mat.dot(x_sol)
residuals = y_val - y_fit
# We only want residuals of the measurement equations (not the regularization equations)
meas_residuals = residuals[:n_blocks * 4 * N]
print(f"Residual standard deviation of measurements: {meas_residuals.std() * 1000:.4f} um")

# Run Residual Bootstrapping
B = 50
print(f"Running {B} bootstrap iterations...")
bootstrap_params = []

# Set random seed for reproducibility
np.random.seed(42)

for b_run in range(B):
    # Resample residuals with replacement
    boot_res = np.random.choice(meas_residuals, size=len(meas_residuals), replace=True)
    # The regularization equations still target 0.0
    boot_y = np.concatenate([y_fit[:n_blocks * 4 * N] + boot_res, np.zeros(n_blocks * 2 * N)])
    
    # Solve bootstrap system
    res_boot = lsqr(A_mat, boot_y, damp=0.0)
    bootstrap_params.append(res_boot[0])
    
    if (b_run + 1) % 10 == 0:
        print(f"  Completed {b_run + 1}/{B} iterations")

bootstrap_params = np.array(bootstrap_params) # Shape: (B, total_params)

# Extract bootstrap parameters for statistical analysis
s_x_boot = bootstrap_params[:, idx_sx] * 1e6 # ppm
s_y_boot = bootstrap_params[:, idx_sy] * 1e6 # ppm
alpha_boot = bootstrap_params[:, idx_alpha] * 1e6 # urad

# Standard errors
se_sx = s_x_boot.std()
se_sy = s_y_boot.std()
se_alpha = alpha_boot.std()

print("\n--- Calibration Parameter Uncertainty (1-Sigma Standard Error) ---")
print(f"CMM X scale error (s_x):      {x_sol[idx_sx]*1e6:8.4f} +/- {se_sx:6.4f} ppm")
print(f"CMM Y scale error (s_y):      {x_sol[idx_sy]*1e6:8.4f} +/- {se_sy:6.4f} ppm")
print(f"CMM squareness error (alpha): {x_sol[idx_alpha]*1e6:8.4f} +/- {se_alpha:6.4f} urad")

# Calculate correlation matrix to inspect interdependencies
df_boot = pd.DataFrame({
    's_x': s_x_boot,
    's_y': s_y_boot,
    'alpha': alpha_boot
})
corr_matrix = df_boot.corr()
print("\n--- Correlation Matrix of Global Parameters (Interdependency) ---")
print(corr_matrix.to_string())

# Compute spatial uncertainty maps for physical deviations (standard error)
se_du = bootstrap_params[:, :N].std(axis=0) * 1000 # um
se_dv = bootstrap_params[:, N:2*N].std(axis=0) * 1000 # um

se_du_grid = se_du.reshape((23, 41))
se_dv_grid = se_dv.reshape((23, 41))

# Save outputs to CSV
df_se = pd.DataFrame({
    'u_nominal_mm': u_flat,
    'v_nominal_mm': v_flat,
    'se_du_um': se_du,
    'se_dv_um': se_dv
})
df_se.to_csv("holes/calibrated_deviations_uncertainty.csv", index=False)
print("\nSaved coordinate uncertainties to holes/calibrated_deviations_uncertainty.csv")

# Plot 1: Bootstrap distributions of global CMM parameters
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].hist(s_x_boot, bins=15, color='#2b5c8f', alpha=0.7, edgecolor='black')
axes[0].axvline(x_sol[idx_sx]*1e6, color='red', linestyle='--', linewidth=2, label='Fitted Value')
axes[0].set_title(f"$s_x$ Distribution\n{x_sol[idx_sx]*1e6:.2f} $\pm$ {se_sx:.2f} ppm")
axes[0].set_xlabel("X Scale Error (ppm)")
axes[0].set_ylabel("Frequency")
axes[0].legend()

axes[1].hist(s_y_boot, bins=15, color='#d95f02', alpha=0.7, edgecolor='black')
axes[1].axvline(x_sol[idx_sy]*1e6, color='red', linestyle='--', linewidth=2, label='Fitted Value')
axes[1].set_title(f"$s_y$ Distribution\n{x_sol[idx_sy]*1e6:.2f} $\pm$ {se_sy:.2f} ppm")
axes[1].set_xlabel("Y Scale Error (ppm)")
axes[1].legend()

axes[2].hist(alpha_boot, bins=15, color='#7570b3', alpha=0.7, edgecolor='black')
axes[2].axvline(x_sol[idx_alpha]*1e6, color='red', linestyle='--', linewidth=2, label='Fitted Value')
axes[2].set_title(f"$\\alpha$ Distribution\n{x_sol[idx_alpha]*1e6:.2f} $\pm$ {se_alpha:.2f} $\\mu$rad")
axes[2].set_xlabel("Squareness Error ($\\mu$rad)")
axes[2].legend()

plt.suptitle("Bootstrap Parameter Distributions & Standard Errors (B=50)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig("holes/cmm_bootstrap_distributions.png", dpi=150)
plt.close()

# Plot 2: 2D Spatial uncertainty (standard error) map on the plate
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

im0 = axes[0].imshow(se_du_grid, origin='lower', extent=[u_vals.min(), u_vals.max(), v_vals.min(), v_vals.max()],
                     cmap='plasma', aspect='auto')
axes[0].set_title("Standard Error in $u$ Coordinate ($\\mu$m)")
axes[0].set_xlabel("U Coordinate (mm)")
axes[0].set_ylabel("V Coordinate (mm)")
fig.colorbar(im0, ax=axes[0], label="Uncertainty ($\\mu$m)")

im1 = axes[1].imshow(se_dv_grid, origin='lower', extent=[u_vals.min(), u_vals.max(), v_vals.min(), v_vals.max()],
                     cmap='plasma', aspect='auto')
axes[1].set_title("Standard Error in $v$ Coordinate ($\\mu$m)")
axes[1].set_xlabel("U Coordinate (mm)")
axes[1].set_ylabel("V Coordinate (mm)")
fig.colorbar(im1, ax=axes[1], label="Uncertainty ($\\mu$m)")

plt.suptitle("2D Spatial Coordinate Uncertainty Map (1-Sigma)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig("holes/coordinate_uncertainty_map.png", dpi=150)
plt.close()

print("Plots saved to holes/cmm_bootstrap_distributions.png and holes/coordinate_uncertainty_map.png")
