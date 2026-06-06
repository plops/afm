# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "lmfit>=1.3.0",
#     "numpy>=1.26.0",
# ]
# ///

# run with: uv run cal.py

import numpy as np
from lmfit import Model, Parameters

# --- 1. SIMULATE LOADING YOUR DATA ---
def load_mock_data(timestep):
    grid_x, grid_y = np.meshgrid(np.linspace(0, 10, 20), np.linspace(0, 10, 20))
    x = grid_x.flatten()
    y = grid_y.flatten()
    
    # Simulating structural error + global thermal drift per timestep
    drift_x = timestep * 0.4
    dx = 0.05 * x + 0.01 * x**2 + 0.02 * y + drift_x + np.random.normal(0, 0.02, len(x))
    return x, y, dx

x1, y1, dx1 = load_mock_data(timestep=1)
x2, y2, dx2 = load_mock_data(timestep=2)
x3, y3, dx3 = load_mock_data(timestep=3)

# --- 2. DECOUPLE GLOBAL THERMAL DRIFT ---
mean_dx1 = np.mean(dx1)
mean_dx2 = np.mean(dx2)
mean_dx3 = np.mean(dx3)

dx1_drift_corrected = dx1 - mean_dx1
dx2_drift_corrected = dx2 - mean_dx2
dx3_drift_corrected = dx3 - mean_dx3

# --- 3. AVERAGE TIME STEPS TO ISOLATE GEOMETRIC ERROR ---
dx_geometric_pure = (dx1_drift_corrected + dx2_drift_corrected + dx3_drift_corrected) / 3.0

# --- 4. FIT SCANNER ARTIFACTS USING LMFIT ---
def scanner_error_model(x, y, scale_x, non_lin_x, cross_talk_y):
    return (scale_x * x) + (non_lin_x * x**2) + (cross_talk_y * y)

model = Model(scanner_error_model, independent_vars=['x', 'y'])

params = Parameters()
params.add('scale_x', value=0.01, min=-1.0, max=1.0)
params.add('non_lin_x', value=0.001, min=-0.1, max=0.1)
params.add('cross_talk_y', value=0.001, min=-0.5, max=0.5)

fit_result = model.fit(dx_geometric_pure, params, x=x1, y=y1)

# --- 5. PRINT DECOUPLED CALIBRATION RESULTS ---
print(fit_result.fit_report())

# --- 6. QUANTIFY RANDOM NOISE FLOOR ---
residuals_t1 = dx1_drift_corrected - dx_geometric_pure
residuals_t2 = dx2_drift_corrected - dx_geometric_pure
residuals_t3 = dx3_drift_corrected - dx_geometric_pure

all_residuals = np.concatenate([residuals_t1, residuals_t2, residuals_t3])
noise_floor_rms = np.sqrt(np.mean(all_residuals**2))
print(f"\nInstrument Precision Limit (Noise Floor RMS): {noise_floor_rms:.4f} units")
