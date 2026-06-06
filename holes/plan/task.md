# Worklog & Task List

- `[x]` Implement bootstrap error estimation script [error_estimation.py](file:///home/kiel/stage/afm/holes/scripts/error_estimation.py)
  - `[x]` Load aligned dataset from `processed_drill_data.nc`
  - `[x]` Implement block bootstrapping (resample blocks with replacement)
  - `[x]` Solve global sparse system for each bootstrap sample
  - `[x]` Calculate standard errors (uncertainties) of physical deviations and CMM parameters
  - `[x]` Calculate correlation matrix of parameters to analyze interdependency
- `[x]` Generate diagnostic visualizations for error estimation
  - `[x]` Plot bootstrap distributions of global CMM parameters
  - `[x]` Plot 2D spatial map of coordinate uncertainty (standard error) on the plate
- `[x]` Update repository documentation
  - `[x]` Add error estimation section and findings to the holes README
  - `[x]` Commit new script and updated files to the repository and push
