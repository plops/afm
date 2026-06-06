"""
Purpose: Perform global sparse self-calibration (reversal method) to separate
true physical hole deviations from CMM scaling, squareness, and drift errors.
Generates diagnostic plots and saves calibrated coordinates.
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

print("Solving global sparse least-squares problem...")
A_mat = A_mat.tocsr()
res_lsqr = lsqr(A_mat, y_val, damp=0.0)
x_sol = res_lsqr[0]

s_x, s_y, alpha = x_sol[idx_sx], x_sol[idx_sy], x_sol[idx_alpha]

# Extract parameter values for each block
block_results = []
true_devs_u = {}
true_devs_v = {}

for b_idx, b in enumerate(block_names):
    theta_1 = x_sol[get_param_idx(b_idx, 2*N)]
    Tx_1 = x_sol[get_param_idx(b_idx, 2*N+1)]
    Ty_1 = x_sol[get_param_idx(b_idx, 2*N+2)]
    theta_2 = x_sol[get_param_idx(b_idx, 2*N+3)]
    Tx_2 = x_sol[get_param_idx(b_idx, 2*N+4)]
    Ty_2 = x_sol[get_param_idx(b_idx, 2*N+5)]
    cx_u = x_sol[get_param_idx(b_idx, 2*N+6)]
    cy_u = x_sol[get_param_idx(b_idx, 2*N+7)]
    cx_r = x_sol[get_param_idx(b_idx, 2*N+8)]
    cy_r = x_sol[get_param_idx(b_idx, 2*N+9)]
    
    du_true = x_sol[get_param_idx(b_idx, 0):get_param_idx(b_idx, N)]
    dv_true = x_sol[get_param_idx(b_idx, N):get_param_idx(b_idx, 2*N)]
    
    true_devs_u[b] = du_true.reshape((23, 41))
    true_devs_v[b] = dv_true.reshape((23, 41))
    
    dx_u = ds['dx_unrot'].sel(block=b).values.flatten()
    dy_u = ds['dy_unrot'].sel(block=b).values.flatten()
    dx_r = ds['dx_rot'].sel(block=b).values.flatten()
    dy_r = ds['dy_rot'].sel(block=b).values.flatten()

    t_u_raw = ds['time_unrot'].sel(block=b).values.flatten()
    t_r_raw = ds['time_rot'].sel(block=b).values.flatten()
    t_u = (t_u_raw - t_u_raw.min()) / np.timedelta64(1, 's')
    t_r = (t_r_raw - t_r_raw.min()) / np.timedelta64(1, 's')

    # Raw differences
    diff_u_before = np.std(dx_u + dy_r) * 1000
    diff_v_before = np.std(dy_u - dx_r) * 1000

    # Corrected differences
    du_unrot_corr = dx_u - (s_x * u_flat - theta_1 * v_flat + Tx_1 + cx_u * t_u)
    dv_unrot_corr = dy_u - (s_y * v_flat + (theta_1 + alpha) * u_flat + Ty_1 + cy_u * t_u)
    du_rot_corr = -dy_r + (-s_y * u_flat + (theta_2 + alpha) * v_flat + Ty_2 + cy_r * t_r)
    dv_rot_corr = dx_r - (s_x * v_flat + theta_2 * u_flat + Tx_2 + cx_r * t_r)

    diff_u_after = np.std(du_unrot_corr - du_rot_corr) * 1000
    diff_v_after = np.std(dv_unrot_corr - dv_rot_corr) * 1000

    block_results.append({
        'block': b,
        's_x_ppm': s_x * 1e6,
        's_y_ppm': s_y * 1e6,
        'alpha_urad': alpha * 1e6,
        'cx_u_mm_hr': cx_u * 3600,
        'cy_u_mm_hr': cy_u * 3600,
        'cx_r_mm_hr': cx_r * 3600,
        'cy_r_mm_hr': cy_r * 3600,
        'diff_u_before_um': diff_u_before,
        'diff_u_after_um': diff_u_after,
        'diff_v_before_um': diff_v_before,
        'diff_v_after_um': diff_v_after,
    })

df_res = pd.DataFrame(block_results)
df_res.to_csv("holes/global_calibration_results.csv", index=False)
print("Saved global calibration results to holes/global_calibration_results.csv")

# Plot 1: Parameters (Global & block-wise)
fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
x_ticks = np.arange(len(block_names))

axes[0].bar(x_ticks - 0.2, [s_x * 1e6]*n_blocks, width=0.4, label='X scale error ($s_x$ = 28.7 ppm)', color='#2b5c8f')
axes[0].bar(x_ticks + 0.2, [s_y * 1e6]*n_blocks, width=0.4, label='Y scale error ($s_y$ = 1.6 ppm)', color='#d95f02')
axes[0].set_ylabel('Scale Error (ppm)')
axes[0].set_title('Global CMM Scale and Squareness Errors')
axes[0].legend()

axes[1].bar(x_ticks, [alpha * 1e6]*n_blocks, width=0.4, label='Squareness error ($\\alpha$ = 17.3 $\\mu$rad)', color='#7570b3')
axes[1].set_ylabel('Squareness Error ($\\mu$rad)')
axes[1].legend()

axes[2].bar(x_ticks - 0.3, df_res['cx_u_mm_hr'], width=0.15, label='Unrot Drift X', color='#1b9e77')
axes[2].bar(x_ticks - 0.1, df_res['cy_u_mm_hr'], width=0.15, label='Unrot Drift Y', color='#a6d854')
axes[2].bar(x_ticks + 0.1, df_res['cx_r_mm_hr'], width=0.15, label='Rot Drift X', color='#ffd92f')
axes[2].bar(x_ticks + 0.3, df_res['cy_r_mm_hr'], width=0.15, label='Rot Drift Y', color='#e78ac3')
axes[2].set_ylabel('Drift Rate (mm/hr)')
axes[2].set_title('Fitted CMM Drift Rates per Run')
axes[2].legend()

plt.xticks(x_ticks, block_names, rotation=30, ha='right')
plt.tight_layout()
plt.savefig("holes/cmm_calibration_parameters.png", dpi=150)
plt.close()

# Plot 2: Vector fields comparison for panel1Top1
b = 'panel1Top1'
b_idx = list(block_names).index(b)

dx_u = ds['dx_unrot'].sel(block=b).values
dy_u = ds['dy_unrot'].sel(block=b).values
dx_r = ds['dx_rot'].sel(block=b).values
dy_r = ds['dy_rot'].sel(block=b).values

du_rot_raw = -dy_r
dv_rot_raw = dx_r

du_est = true_devs_u[b]
dv_est = true_devs_v[b]

theta_1 = x_sol[get_param_idx(b_idx, 2*N)]
Tx_1 = x_sol[get_param_idx(b_idx, 2*N+1)]
Ty_1 = x_sol[get_param_idx(b_idx, 2*N+2)]
theta_2 = x_sol[get_param_idx(b_idx, 2*N+3)]
Tx_2 = x_sol[get_param_idx(b_idx, 2*N+4)]
Ty_2 = x_sol[get_param_idx(b_idx, 2*N+5)]
cx_u = x_sol[get_param_idx(b_idx, 2*N+6)]
cy_u = x_sol[get_param_idx(b_idx, 2*N+7)]
cx_r = x_sol[get_param_idx(b_idx, 2*N+8)]
cy_r = x_sol[get_param_idx(b_idx, 2*N+9)]

t_u_b = (ds['time_unrot'].sel(block=b).values - ds['time_unrot'].sel(block=b).values.min()) / np.timedelta64(1, 's')
t_r_b = (ds['time_rot'].sel(block=b).values - ds['time_rot'].sel(block=b).values.min()) / np.timedelta64(1, 's')

du_unrot_c = dx_u - (s_x * u_grid - theta_1 * v_grid + Tx_1 + cx_u * t_u_b)
dv_unrot_c = dy_u - (s_y * v_grid + (theta_1 + alpha) * u_grid + Ty_1 + cy_u * t_u_b)
du_rot_c = -dy_r + (-s_y * u_grid + (theta_2 + alpha) * v_grid + Ty_2 + cy_r * t_r_b)
dv_rot_c = dx_r - (s_x * v_grid + theta_2 * u_grid + Tx_2 + cx_r * t_r_b)

diff_before = np.sqrt((dx_u - du_rot_raw)**2 + (dy_u - dv_rot_raw)**2) * 1000
diff_after = np.sqrt((du_unrot_c - du_rot_c)**2 + (dv_unrot_c - dv_rot_c)**2) * 1000

fig, axes = plt.subplots(2, 2, figsize=(15, 12))

q0 = axes[0, 0].quiver(u_grid, v_grid, dx_u, dy_u, scale=0.5, color='#2b5c8f')
axes[0, 0].quiverkey(q0, 0.9, 0.95, 0.02, '20 um', labelpos='E', coordinates='axes')
axes[0, 0].set_title("Raw Deviations: Unrotated Run (dx_unrot, dy_unrot)")
axes[0, 0].set_xlabel("U Coordinate (mm)")
axes[0, 0].set_ylabel("V Coordinate (mm)")

q1 = axes[0, 1].quiver(u_grid, v_grid, du_rot_raw, dv_rot_raw, scale=0.5, color='#d95f02')
axes[0, 1].quiverkey(q1, 0.9, 0.95, 0.02, '20 um', labelpos='E', coordinates='axes')
axes[0, 1].set_title("Raw Deviations: Rotated Run (-dy_rot, dx_rot)")
axes[0, 1].set_xlabel("U Coordinate (mm)")
axes[0, 1].set_ylabel("V Coordinate (mm)")

q2 = axes[1, 0].quiver(u_grid, v_grid, du_est, dv_est, scale=0.5, color='#7570b3')
axes[1, 0].quiverkey(q2, 0.9, 0.95, 0.02, '20 um', labelpos='E', coordinates='axes')
axes[1, 0].set_title("Estimated True Plate Deviations (Self-Calibrated)")
axes[1, 0].set_xlabel("U Coordinate (mm)")
axes[1, 0].set_ylabel("V Coordinate (mm)")

axes[1, 1].hist(diff_before.flatten(), bins=30, alpha=0.5, label='Raw mismatch (Unrot vs Rot)', color='#d95f02')
axes[1, 1].hist(diff_after.flatten(), bins=30, alpha=0.5, label='Self-Calibrated mismatch', color='#1b9e77')
axes[1, 1].set_xlabel("Mismatch between runs ($\\mu$m)")
axes[1, 1].set_ylabel("Count")
axes[1, 1].set_title("Histogram of Measurement Mismatch")
axes[1, 1].legend()

plt.suptitle(f"Plate and CMM Measurement Error Analysis: {b} (Global Calibration)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig("holes/cmm_drift_and_deviations.png", dpi=150)
plt.close()

# Save calibrated deviations to netCDF
ds_corr = xr.Dataset(
    data_vars={
        'du_true': (['block', 'v', 'u'], np.stack([true_devs_u[b] for b in block_names])),
        'dv_true': (['block', 'v', 'u'], np.stack([true_devs_v[b] for b in block_names])),
    },
    coords={
        'block': block_names,
        'v': ds['v'].values,
        'u': ds['u'].values,
    }
)
ds_corr.to_netcdf("holes/calibrated_physical_deviations.nc")
print("Saved calibrated physical deviations to holes/calibrated_physical_deviations.nc")
