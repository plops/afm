# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "lmfit>=1.3.0",
#     "numpy>=1.26.0",
#     "matplotlib>=3.8.0",
# ]
# ///

import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameters

# --- 1. SIMULATE LOADING YOUR DATA ---
def load_mock_data(timestep):
    grid_x, grid_y = np.meshgrid(np.linspace(0, 10, 20), np.linspace(0, 10, 20))
    x = grid_x.flatten()
    y = grid_y.flatten()
    
    # Structural error + global thermal drift per timestep
    drift_x = timestep * 0.4
    drift_y = timestep * -0.15 # added y-drift for a realistic vector chart
    
    # X and Y distortion components
    dx = 0.05 * x + 0.01 * x**2 + 0.02 * y + drift_x + np.random.normal(0, 0.02, len(x))
    dy = -0.02 * x + 0.04 * y + drift_y + np.random.normal(0, 0.02, len(y))
    return x, y, dx, dy

x1, y1, dx1, dy1 = load_mock_data(timestep=1)
x2, y2, dx2, dy2 = load_mock_data(timestep=2)
x3, y3, dx3, dy3 = load_mock_data(timestep=3)

# --- 2. DECOUPLE GLOBAL THERMAL DRIFT ---
# Tracking spatial means for both axes
mean_dx = [np.mean(dx1), np.mean(dx2), np.mean(dx3)]
mean_dy = [np.mean(dy1), np.mean(dy2), np.mean(dy3)]

dx1_corr = dx1 - mean_dx[0]
dx2_corr = dx2 - mean_dx[1]
dx3_corr = dx3 - mean_dx[2]

dy1_corr = dy1 - mean_dy[0]
dy2_corr = dy2 - mean_dy[1]
dy3_corr = dy3 - mean_dy[2]

# --- 3. AVERAGE TIME STEPS TO ISOLATE GEOMETRIC ERROR ---
dx_geometric = (dx1_corr + dx2_corr + dx3_corr) / 3.0
dy_geometric = (dy1_corr + dy2_corr + dy3_corr) / 3.0

# --- 4. FIT SCANNER ARTIFACTS USING LMFIT ---
def scanner_error_model(x, y, scale_x, non_lin_x, cross_talk_y):
    return (scale_x * x) + (non_lin_x * x**2) + (cross_talk_y * y)

model = Model(scanner_error_model, independent_vars=['x', 'y'])
params = Parameters()
params.add('scale_x', value=0.01, min=-1.0, max=1.0)
params.add('non_lin_x', value=0.001, min=-0.1, max=0.1)
params.add('cross_talk_y', value=0.001, min=-0.5, max=0.5)

fit_result = model.fit(dx_geometric, params, x=x1, y=y1)
print(fit_result.fit_report())

# Calculate noise floor residuals
residuals_t1 = dx1_corr - dx_geometric
all_residuals = np.concatenate([residuals_t1, dx2_corr - dx_geometric, dx3_corr - dx_geometric])
print(f"\nInstrument Precision Limit (Noise Floor RMS): {np.sqrt(np.mean(all_residuals**2)):.4f} units")

# --- 5. GENERATE DIAGNOSTIC DIAGRAMS ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Plot 1: Quiver Plot (Vector field of pure geometric scanner distortion)
# This shows exactly how and where the piezo grid pushes features out of bounds.
q = axes[0].quiver(x1, y1, dx_geometric, dy_geometric, color='teal', scale_units='xy', angles='xy')
axes[0].quiverkey(q, X=0.1, Y=0.92, U=1, label='Error Scale (1 Unit)', labelpos='E')
axes[0].set_title("Geometric Scanner Distortion Field")
axes[0].set_xlabel("X Position")
axes[0].set_ylabel("Y Position")
axes[0].grid(True, linestyle='--', alpha=0.5)

# Plot 2: Heatmap of Fit Residuals (Targeted to your X-axis fit)
# If this look completely random (like salt and pepper), your polynomial model successfully decoupled the error.
res_grid = fit_result.residual.reshape(20, 20)
im = axes[1].imshow(res_grid, extent=[0, 10, 0, 10], origin='lower', cmap='bwr', aspect='auto')
fig.colorbar(im, ax=axes[1], label="Residual Error")
axes[1].set_title("X-Fit Spatial Residuals (Noise Profile)")
axes[1].set_xlabel("X Position")
axes[1].set_ylabel("Y Position")

# Plot 3: Thermal Drift Profile over Time
# Quantifies frame movement across the 2-hour waiting periods.
timesteps = [0, 2, 4] # 2 hours wait between measurements
axes[2].plot(timesteps, mean_dx, 'o-', label='X Drift Baseline', color='crimson')
axes[2].plot(timesteps, mean_dy, 's-', label='Y Drift Baseline', color='navy')
axes[2].set_title("Global Thermal Drift Track")
axes[2].set_xlabel("Time (Hours)")
axes[2].set_ylabel("Absolute Sensor Offset")
axes[2].legend()
axes[2].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig("afm_diagnostic_plots.png", dpi=200)
print("\n[Success] Diagrams successfully saved to 'afm_diagnostic_plots.png'")
