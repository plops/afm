"""
Purpose: Preprocess CMM grid hole measurements from raw CSV format, align
coordinate grids, map rotated runs to physical coordinates, and save as netCDF.
"""

import pandas as pd
import numpy as np
import xarray as xr

# Load CSVs
df_data = pd.read_csv("holes/DrillData.csv")
df_rot = pd.read_csv("holes/DrillRot90_2.csv")

# Constants
N_BLOCKS = 8
N_HOLES = 943 # 23 * 41

# Blocks metadata
block_names = [
    'panel1Top1', 'panel1Top2', 'panel1Bot1', 'panel1Bot2',
    'panel2Top1', 'panel2Top2', 'panel2Bot1', 'panel2Bot2'
]

# We will define a grid with coordinates:
# u (short axis, 41 values: 0, 12.5, ..., 500)
# v (long axis, 23 values: 0, 25, ..., 550)
u_vals = np.arange(41) * 12.5
v_vals = np.arange(23) * 25.0

# Initialize arrays for unrotated data (shape: 8, 23, 41)
dx_unrot = np.zeros((N_BLOCKS, 23, 41))
dy_unrot = np.zeros((N_BLOCKS, 23, 41))
d_unrot = np.zeros((N_BLOCKS, 23, 41))
time_unrot = np.zeros((N_BLOCKS, 23, 41), dtype='datetime64[ns]')

# Initialize arrays for rotated data (shape: 8, 23, 41)
dx_rot = np.zeros((N_BLOCKS, 23, 41))
dy_rot = np.zeros((N_BLOCKS, 23, 41))
d_rot = np.zeros((N_BLOCKS, 23, 41))
time_rot = np.zeros((N_BLOCKS, 23, 41), dtype='datetime64[ns]')

# Fill unrotated data
for b in range(N_BLOCKS):
    block = df_data.iloc[b*2829:(b+1)*2829]
    for (y_idx, x_idx), group in block.groupby(['Unnamed: 2', 'Unnamed: 3']):
        v_idx = int(y_idx) - 1
        u_idx = int(x_idx) - 1
        
        row_x = group[group['Unnamed: 4'] == 'X'].iloc[0]
        row_y = group[group['Unnamed: 4'] == 'Y'].iloc[0]
        row_d = group[group['Unnamed: 4'] == 'D'].iloc[0]
        
        dx_unrot[b, v_idx, u_idx] = row_x['Unnamed: 9']
        dy_unrot[b, v_idx, u_idx] = row_y['Unnamed: 9']
        d_unrot[b, v_idx, u_idx] = row_d['Unnamed: 6']
        
        dt_str = f"{row_x['Unnamed: 15']} {row_x['Unnamed: 14']}"
        time_unrot[b, v_idx, u_idx] = np.datetime64(dt_str)

# Fill rotated data
for b in range(N_BLOCKS):
    block = df_rot.iloc[b*2829:(b+1)*2829]
    for (y_coord, x_coord), group in block.groupby(['Y Coordinate', 'X Coordinate']):
        v_idx = int(x_coord) - 1
        u_idx = 42 - int(y_coord) - 1
        
        row_x = group[group['axis'] == 'X'].iloc[0]
        row_y = group[group['axis'] == 'Y'].iloc[0]
        row_d = group[group['axis'] == 'D'].iloc[0]
        
        dx_rot[b, v_idx, u_idx] = row_x['diff(mm)']
        dy_rot[b, v_idx, u_idx] = row_y['diff(mm)']
        d_rot[b, v_idx, u_idx] = row_d['realValue']
        
        dt_str = f"{row_x['Unnamed: 14'].split()[0]} {row_x['Unnamed: 13']}"
        time_rot[b, v_idx, u_idx] = np.datetime64(dt_str)

# Create xarray dataset
ds = xr.Dataset(
    data_vars={
        'dx_unrot': (['block', 'v', 'u'], dx_unrot),
        'dy_unrot': (['block', 'v', 'u'], dy_unrot),
        'd_unrot': (['block', 'v', 'u'], d_unrot),
        'time_unrot': (['block', 'v', 'u'], time_unrot),
        'dx_rot': (['block', 'v', 'u'], dx_rot),
        'dy_rot': (['block', 'v', 'u'], dy_rot),
        'd_rot': (['block', 'v', 'u'], d_rot),
        'time_rot': (['block', 'v', 'u'], time_rot),
    },
    coords={
        'block': block_names,
        'v': v_vals,
        'u': u_vals,
    }
)

ds.to_netcdf("holes/processed_drill_data.nc")
print("Saved netCDF dataset to holes/processed_drill_data.nc")
