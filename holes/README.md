# CMM Measurement Calibration & Self-Calibration Report

This folder contains the calibration analysis of Coordinate Measuring Machine (CMM) hole grid measurements. It details how we compared measurements of the same plate in two different orientations to isolate and correct CMM errors down to the micrometer level.

---

## 1. Description of the Datasets

We analyzed two primary datasets:
1. **[DrillData.csv](DrillData.csv) (Unrotated):** Measured with the plate aligned with the CMM axes.
2. **[DrillRot90_2.csv](DrillRot90_2.csv) (Rotated):** Measured with the plate rotated by 90° clockwise.

Each dataset contains **8 blocks** representing different measurement runs:
- `panel1Top1`, `panel1Top2`: Runs 1 & 2 of Panel 1, Top side.
- `panel1Bot1`, `panel1Bot2`: Runs 1 & 2 of Panel 1, Bottom side.
- `panel2Top1`, `panel2Top2`: Runs 1 & 2 of Panel 2, Top side.
- `panel2Bot1`, `panel2Bot2`: Runs 1 & 2 of Panel 2, Bottom side.

Each run measures a grid of **943 holes**:
- Spatially arranged as $23 \times 41$ holes.
- Nominal spacing: $12.5$ mm along the short axis ($u$) and $25.0$ mm along the long axis ($v$).
- For each hole, the CMM records coordinates along the **`X` axis**, **`Y` axis**, and the hole diameter **`D`**.

---

## 2. Coordinate Mapping & Reversal Physics

To separate CMM errors from actual plate manufacturing deviations, we use the **reversal method** (self-calibration). 

Let us define:
- $(u, v)$ as the nominal, plate-fixed coordinate system of the holes, where:
  - $u$ ranges from $0$ to $500$ mm ($41$ holes, $12.5$ mm spacing).
  - $v$ ranges from $0$ to $550$ mm ($23$ holes, $25$ mm spacing).
- $(\Delta u, \Delta v)$ as the **true physical deviations** of the holes (manufacturing errors). These are properties of the plate and do not change.

### 2.1 Unrotated Setup
In the unrotated run, the plate is aligned directly with the CMM:
- CMM X-axis is along plate $u$ ($X_{nom} = u$).
- CMM Y-axis is along plate $v$ ($Y_{nom} = v$).
- The physical deviations map directly:
  - $\Delta X_{phys} = \Delta u$
  - $\Delta Y_{phys} = \Delta v$

### 2.2 Rotated Setup (90° Clockwise)
In the rotated run, the plate is rotated 90° clockwise on the CMM bed. Comparing hole diameters confirmed the physical mapping:
- CMM X-axis is along plate $v$ ($X_{nom} = v$).
- CMM Y-axis is along plate $-u$ with an offset ($Y_{nom} = 500 - u$).
- Because the plate is rotated, the physical deviations rotate too:
  - A deviation in the plate's $v$ direction now points along CMM X: $\Delta X_{phys} = \Delta v$
  - A deviation in the plate's $u$ direction now points along negative CMM Y: $\Delta Y_{phys} = -\Delta u$

---

## 3. CMM Error Sources and Projections

The coordinate values reported by the CMM are corrupted by three sources of error:
1. **Geometric Axis Errors (Scale & Shear):**
   - **X scale error ($s_x$):** Adds a coordinate-dependent error $s_x \cdot X_{nom}$.
   - **Y scale error ($s_y$):** Adds a coordinate-dependent error $s_y \cdot Y_{nom}$.
   - **Squareness error ($\alpha$):** Non-perpendicularity between X and Y axes. Adds a shear error $\alpha \cdot X_{nom}$ to the Y axis.
2. **Alignment Errors (Translation & Rotation):**
   - Small fixturing rotations ($\theta$) and translations ($T_x, T_y$).
3. **Time-Dependent Drift ($c_x, c_y$):**
   - Thermal expansion or servo drift over time. Adds $c_x \cdot t$ and $c_y \cdot t$, where $t$ is the elapsed time since the start of the run.

---

## 4. Mathematical Equations of the Self-Calibration Model

By combining the physical deviations, geometric errors, alignments, and drift, we obtain the full mathematical model for the coordinate deviations.

### 4.1 Unrotated Model Equations
For each hole $i$:
- **X deviation ($dx_{unrot,i}$):**
  - *LaTeX*: $\Delta x_{unrot,i} = \Delta u_i + s_x \cdot u_i - \theta_1 \cdot v_i + T_{x1} + c_{x,u} \cdot t_{u,i}$
  - *Plain text*: `dx_unrot = du + sx * u - theta_1 * v + Tx_1 + cx_u * t_u`
- **Y deviation ($dy_{unrot,i}$):**
  - *LaTeX*: $\Delta y_{unrot,i} = \Delta v_i + s_y \cdot v_i + (\theta_1 + \alpha) \cdot u_i + T_{y1} + c_{y,u} \cdot t_{u,i}$
  - *Plain text*: `dy_unrot = dv + sy * v + (theta_1 + alpha) * u + Ty_1 + cy_u * t_u`

### 4.2 Rotated Model Equations
For each hole $i$ (with $X_{nom} = v_i$ and $Y_{nom} = 500 - u_i$):
- **X deviation ($dx_{rot,i}$):**
  - *LaTeX*: $\Delta x_{rot,i} = \Delta v_i + s_x \cdot v_i + \theta_2 \cdot u_i + T_{x2} + c_{x,r} \cdot t_{r,i}$
  - *Plain text*: `dx_rot = dv + sx * v + theta_2 * u + Tx_2 + cx_r * t_r`
- **Y deviation ($dy_{rot,i}$):**
  - *LaTeX*: $\Delta y_{rot,i} = -\Delta u_i - s_y \cdot u_i + (\theta_2 + \alpha) \cdot v_i + T'_{y2} + c_{y,r} \cdot t_{r,i}$
  - *Plain text*: `dy_rot = -du - sy * u + (theta_2 + alpha) * v + Ty_2 + cy_r * t_r`
  *(Note: The constant term $500 \cdot s_y$ has been absorbed into the translation offset parameter $T'_{y2}$.)*

---

## 5. Solving the System (Separating Scale from Drift)

If we only analyze the unrotated run, the measurement time is highly correlated with the $v$ coordinate because the CMM scans row-by-row. Consequently, Y-drift ($c_{y,u} \cdot t$) and Y-scaling ($s_y \cdot v$) are collinear and cannot be separated.

The rotated run breaks this degeneracy:
- In the unrotated run: Time is collinear with physical $v$ (long axis).
- In the rotated run: Time is collinear with physical $u$ (short axis).

By solving the equations for all 8 runs simultaneously using a global sparse linear solver (`scipy.sparse.linalg.lsqr`), we isolated the true CMM parameters.

---

## 6. Results & Visualizations

### Global CMM Calibration Parameters:
- **X scale error ($s_x$):** **$+28.67$ ppm** (stretched by $28.67 \mu$m per meter)
- **Y scale error ($s_y$):** **$+1.59$ ppm** (virtually zero error)
- **Squareness error ($\alpha$):** **$+17.33 \mu$rad**

### Dynamic Drift Rates:
The true time-dependent drift rates are extremely small: **$10$ to $40 \mu$m/hr**. The apparent large "drift" was a signature of the static $28.67$ ppm scale error on CMM X, which rotated relative to the plate.

### Mismatch Reduction:
Applying this calibration reduced the coordinate mismatch between the unrotated and rotated runs from up to **$6.6 \mu$m** down to **$1.0 - 1.5 \mu$m** (the repeatability limit of the CMM probe).

### Parameter Plots
![CMM Calibration Parameters](cmm_calibration_parameters.png)

### Vector Field & Mismatch Analysis (`panel1Top1`)
![CMM Drift and Deviations](cmm_drift_and_deviations.png)

---

## 7. Error Estimation & Parameter Interdependencies

To evaluate the precision of our calibrated physical coordinates and CMM parameters, we performed **Residual Bootstrapping** over $B=50$ runs. In each iteration, we resampled the residuals of our global fit (measurement noise $\sigma \approx 0.63\ \mu$m) with replacement and refit the global linear system.

### 7.1 Parameter Standard Errors (1-Sigma Confidence)
- **X scale error ($s_x$):** $28.6695 \pm 0.1256$ ppm
- **Y scale error ($s_y$):** $1.5939 \pm 0.1946$ ppm
- **Squareness error ($\alpha$):** $17.3279 \pm 0.8229\ \mu$rad

These values prove that the CMM geometric errors are determined with exceptionally high precision (errors of $<0.2$ ppm for scale and $<0.85\ \mu$rad for squareness).

### 7.2 Parameter Interdependency (Correlation Matrix)
The bootstrap run revealed the following correlation coefficients:

| Parameter | $s_x$ | $s_y$ | $\alpha$ |
| :--- | :---: | :---: | :---: |
| **$s_x$** | 1.000 | 0.583 | 0.182 |
| **$s_y$** | 0.583 | 1.000 | **0.797** |
| **$\alpha$** | 0.182 | **0.797** | 1.000 |

> [!IMPORTANT]
> The **high correlation of $0.80$** between Y-scale ($s_y$) and squareness ($\alpha$) highlights a strong multivariate interdependency. In the rotated coordinate equations, $s_y$ and $\alpha$ are coupled through the alignment rotation angles ($\theta_1, \theta_2$). Despite this coupling, the global combination of unrotated and rotated measurements is strong enough to cleanly isolate both parameters with very small uncertainties.

### 7.3 Coordinate Uncertainty (2D Spatial Map)
The standard error of the calibrated true physical coordinate deviations ($\Delta u, \Delta v$) is extremely low and uniform across the plate:
- **Average standard error in $u$**: **$0.12\ \mu$m**
- **Average standard error in $v$**: **$0.13\ \mu$m**

The uncertainty is uniform across the entire surface of the plate, indicating that the self-calibration remains highly stable even near the edges and corners.

### Bootstrap Parameter Distributions
![Bootstrap Parameter Distributions](cmm_bootstrap_distributions.png)

### 2D Spatial Coordinate Uncertainty Map
![2D Coordinate Uncertainty Map](coordinate_uncertainty_map.png)

---

## 8. Python Analysis Scripts

All data analysis was executed using the scripts in the [scripts/](scripts/) directory:

1. **[scripts/process_data.py](scripts/process_data.py)**: Aligns and structures the unrotated and rotated CSV files into a single netCDF dataset (`processed_drill_data.nc`).
2. **[scripts/fit_global_and_plot.py](scripts/fit_global_and_plot.py)**: Performs the global self-calibration, outputs the calibration results and corrected coordinates, and generates the vector plots.
3. **[scripts/error_estimation.py](scripts/error_estimation.py)**: Performs the residual bootstrapping, computes standard errors and parameter correlations, and plots the uncertainty distributions.
