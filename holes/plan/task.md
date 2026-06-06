# Worklog & Task List

- `[ ]` Implement bootstrap error estimation script [error_estimation.py](file:///home/kiel/stage/afm/holes/scripts/error_estimation.py)
  - `[ ]` Load aligned dataset from `processed_drill_data.nc`
  - `[ ]` Implement block bootstrapping (resample blocks with replacement)
  - `[ ]` Solve global sparse system for each bootstrap sample
  - `[ ]` Calculate standard errors (uncertainties) of physical deviations and CMM parameters
  - `[ ]` Calculate correlation matrix of parameters to analyze interdependency
- `[ ]` Generate diagnostic visualizations for error estimation
  - `[ ]` Plot bootstrap distributions of global CMM parameters
  - `[ ]` Plot 2D spatial map of coordinate uncertainty (standard error) on the plate
- `[ ]` Update repository documentation
  - `[ ]` Add error estimation section and findings to the holes README
  - `[ ]` Commit new script and updated files to the repository and push
