# CMM Measurement Drift & Self-Calibration Report

This folder contains the self-calibration analysis of the Coordinate Measuring Machine (CMM) hole grid measurements in [DrillData.xlsx](DrillData.xlsx) (unrotated) and [DrillRot90_2.xlsx](DrillRot90_2.xlsx) (rotated 90° clockwise).

By comparing the two datasets, we successfully isolated the true physical deviations of the holes from the CMM's geometric calibration errors (scale and shear) and time-dependent drift.

---

## 1. Description of the Datasets

We analyzed two primary datasets:
1. **[DrillData.csv](DrillData.csv) (Unrotated):** The reference grid dataset.
2. **[DrillRot90_2.csv](DrillRot90_2.csv) (Rotated):** The grid dataset measured with the plate rotated by 90° clockwise.

Each dataset contains **8 blocks** representing different measurement runs:
- `panel1Top1`, `panel1Top2`: Runs 1 & 2 of Panel 1, Top side.
- `panel1Bot1`, `panel1Bot2`: Runs 1 & 2 of Panel 1, Bottom side.
- `panel2Top1`, `panel2Top2`: Runs 1 & 2 of Panel 2, Top side.
- `panel2Bot1`, `panel2Bot2`: Runs 1 & 2 of Panel 2, Bottom side.

Each block contains a grid of **943 holes**:
- Spatially arranged as $23 \times 41$ holes.
- Spacings: $12.5$ mm along the short axis ($u$) and $25.0$ mm along the long axis ($v$).
- Each hole measurement has three rows corresponding to:
  - **`X` axis**: nominal/actual coordinates and deviation.
  - **`Y` axis**: nominal/actual coordinates and deviation.
  - **`D` axis**: nominal/actual hole diameter.

---

## 2. Report of What Was Attempted

We performed a systematic analysis to understand and resolve the measurement discrepancies:

### Step 1: Alignment & Rotation Mapping
- **Action**: Reshaped the $23 \times 41$ diameter maps from both files and ran cross-correlations across different rotations and flips.
- **Result**: Successfully determined that the rotated measurements are related to the unrotated ones by a **$90^\circ$ clockwise rotation** combined with a translation and minor coordinate system rotation. The exact index mapping is $Y_{rot} = 42 - X_{unrot}$ and $X_{rot} = Y_{unrot}$.

### Step 2: Time-Dependent Drift Modeling
- **Action**: Fit a polynomial drift model of the elapsed measurement time ($t$) to the coordinate differences ($dx_{unrot} - dx_{rot}$ and $dy_{unrot} - dy_{rot}$).
- **Result**: Discovered that time drift accounted for very little of the difference (mismatch standard deviation only decreased from $4.1 \mu$m to $3.8 \mu$m), indicating that time-dependent drift was **not** the primary source of the error.

### Step 3: Geometric Error Modeling
- **Action**: Looked for coordinate-dependent trends in the deviations and found very high correlations (up to $0.90$) with the CMM coordinates themselves.
- **Result**: Modeled static CMM geometric errors: X/Y axis scaling errors ($s_x, s_y$) and squareness/shear error ($\alpha$).

### Step 4: Block-by-Block Self-Calibration
- **Action**: Combined both geometric errors and drift rates into a linear model for each block.
- **Result**: Showed that the residuals drop to $<1 \mu$m, but individual fits suffered from collinearity between Y-scaling ($s_y$) and unrotated Y-drift ($c_{y,u}$), since the row-by-row scan time is proportional to Y.

### Step 5: Global Sparse Self-Calibration (Successful)
- **Action**: Built a global sparse linear system ($45,264$ equations, $15,171$ variables) across all 8 runs, sharing a single set of CMM scale and squareness parameters while letting drift rates vary.
- **Result**: Completely resolved the collinearity, giving a single set of CMM calibration parameters:
  - **X scale error ($s_x$):** **$+28.67$ ppm** (stretched)
  - **Y scale error ($s_y$):** **$+1.59$ ppm**
  - **Squareness error ($\alpha$):** **$+17.33 \mu$rad**
- **Drift Rate**: The true thermal drift of the machine was isolated and shown to be negligible (only **$10$ to $30 \mu$m/hr**).
- **Residual Mismatch**: The discrepancy between unrotated and rotated runs collapsed from up to $6.6 \mu$m down to **$1.0 - 1.5 \mu$m** (the repeatability limit of the CMM probe).

---

## 3. Python Scripts

We packaged the required analysis pipeline into the following scripts under the [scripts/](scripts/) directory:

1. **[scripts/process_data.py](scripts/process_data.py)**:
   - **Purpose**: Preprocesses the converted CSV files, aligns the coordinate grids, maps the rotated run to physical coordinates, and saves the combined dataset as a netCDF file (`processed_drill_data.nc`).
   - **Run command**:
     ```bash
     uv run --with pandas --with numpy --with xarray --with netcdf4 python holes/scripts/process_data.py
     ```

2. **[scripts/fit_global_and_plot.py](scripts/fit_global_and_plot.py)**:
   - **Purpose**: Runs the global sparse self-calibration, saves the parameter results to `global_calibration_results.csv`, outputs the drift-corrected physical deviations to `calibrated_physical_deviations.nc`, and generates the finalized PNG plots.
   - **Run command**:
     ```bash
     uv run --with xarray --with numpy --with pandas --with netcdf4 --with scipy --with matplotlib python holes/scripts/fit_global_and_plot.py
     ```

---

## 4. Visualizations

### CMM Calibration Parameters by Run
The plot below compares the global scale/squareness errors and the fitted time-dependent drift rates for each run. Note that the drift rates are very close to zero once the scale/shear collinearity is resolved.

![CMM Calibration Parameters](cmm_calibration_parameters.png)

### Vector Field & Mismatch Analysis (`panel1Top1`)
The figure below displays:
1. **Raw Unrotated Deviations:** The measured deviations from nominal.
2. **Raw Rotated Deviations (Physical Coordinates):** Notice the systematic differences due to scale/shear errors.
3. **Calibrated Physical Deviations:** The true deviations of the holes, which are highly consistent.
4. **Mismatch Histogram:** Shows how the error distribution between the unrotated and rotated runs collapses from up to $7 \mu$m down to a narrow peak around $1 \mu$m.

![CMM Drift and Deviations](cmm_drift_and_deviations.png)
