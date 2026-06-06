[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/plops/afm)

# AFM Calibration Project README

## Overview

This project provides a computational approach to decoupling systemic errors in piezoceramic scanner tubes used in Atomic Force Microscopes (AFM). It implements a five-stage pipeline to isolate geometric artifacts from temporal drift and random noise, ensuring high-precision positioning and instrument characterization.

The software targets three primary error sources:
- **Geometric Distortion**: Non-linear scaling, hysteresis, and cross-talk between axes
- **Thermal Drift**: Uniform translation of the scanner frame due to temperature changes
- **Instrument Noise Floor**: Random measurement uncertainty limit 

## Installation

This project uses the `uv` package manager with PEP 723 inline metadata for dependency management, allowing for reproducible execution without manual environment setup.

### Prerequisites
- Python 3.10 or higher
- uv package manager

### Dependencies
- `numpy>=1.26.0` - Numerical computing
- `lmfit>=1.3.0` - Physics-oriented fitting with parameter bounds
- `matplotlib>=3.8.0` - Diagnostic visualizations (cal2.py only) 

## Usage

### Core Calibration
Run the core calibration script to process displacement data and generate fit statistics:
```bash
uv run cal.py
``` 

### Visualization Suite
Run the extended script with Y-axis tracking and diagnostic visualizations:
```bash
uv run cal2.py
``` 

## Project Structure

| File | Role | Key Features |
|------|------|--------------|
| `cal.py` | Core calibration script | Numerical results, lmfit reports, drift correction, geometric modeling |
| `cal2.py` | Extended script | Y-axis tracking, diagnostic visualizations (quiver plots, heatmaps, drift profiles) |
| `prompt.txt` | Documentation | Background on AFM physics, error sources, and methodology |
| `output.txt` | Sample output | Levenberg-Marquardt fit statistics example |

## Calibration Pipeline

The calibration logic follows a sequential pipeline:

1. **Data Ingestion**: Load multi-timestep displacement data via `load_mock_data()`
2. **Drift Correction**: Compute and subtract spatial means (`np.mean(dx)`) to remove thermal drift
3. **Geometric Isolation**: Average drift-corrected timesteps to isolate geometric error from random noise
4. **Model Fitting**: Use `lmfit` to fit a 2D polynomial model extracting scaling errors and cross-talk
5. **Noise Analysis**: Calculate residual RMS to quantify instrument precision limit

### Key Functions

- `load_mock_data(timestep)`: Simulates loading AFM scan data with structural error and thermal drift 
- `scanner_error_model(x, y, scale_x, non_lin_x, cross_talk_y)`: 2D polynomial model for spatial distortions

## Output Interpretation

The calibration output includes:
- **Fit Statistics**: Chi-square, reduced chi-square, AIC, BIC, R-squared
- **Parameter Estimates**: scale_x, non_lin_x, cross_talk_y with confidence intervals
- **Correlations**: Parameter correlation matrix
- **Noise Floor RMS**: Instrument precision limit in units 

## Domain Context

Understanding the output requires knowledge of AFM-specific phenomena:
- Difference between open-loop and closed-loop scanners
- Why thermal equilibrium is necessary before calibration
- Physical meaning of piezo scaling factors and cross-talk coefficients 

## Notes

This README is based on the project's Overview wiki page and code analysis. The project uses mock data for demonstration - replace `load_mock_data()` with actual file loading logic for real calibration data. 

