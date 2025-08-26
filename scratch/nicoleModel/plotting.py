from glob import glob
import os
import pandas as pd
from Python_H2model_gridsearch import get_H2model, grid_from_radius, vrad_to_vobs, prog_14_datavel, prog_14_dataflux, prog_14_dataerr
import matplotlib.pyplot as plt
import math
import numpy as np

folder_path = "run_outputs"

# Get all CSV files in the folder
csv_files = glob(os.path.join(folder_path, "*.txt"))

# Column names (as in your example)
columns = ["chi^2", "z", "gamma", "T", "q", "rchar", "MH2", "LyA17flux"]

# Read and concatenate all CSV files
# df_list = []
# for file in csv_files:
#     data_table = pd.read_csv(file, sep="/", header=None)
#     data_table.columns = columns
#     df_list.append(data_table)

# combined_df = pd.concat(df_list, ignore_index=True)

# # Save the combined dataframe to the same folder
# output_file = os.path.join(folder_path, "combined_data.csv")
# combined_df.to_csv(output_file, index=False)

data_table = pd.read_csv(os.path.join(folder_path, 'combined_data.csv'))
# Plot best fit result
sorted_data_table = data_table.assign(x=np.abs(np.log(data_table['chi^2'])).replace(-np.inf, np.nan)).sort_values("x", ascending=True, na_position='last').drop('x',axis=1)
index_sorted = sorted_data_table.index
best_fit_idx = index_sorted[0]
print(sorted_data_table.loc[best_fit_idx])

#minimum target information
targ_inclination = 19.
targ_d = 154
theta_keys = ['z', 'gamma', 'T', 'q', 'rchar', 'MH2', 'LyA17flux']
param_dict = {"Inclination": targ_inclination * math.pi / 180.,
              "flux_level": 0.,
              "lnf": 0.,
              "Molecule": "H2"
              }

LyA_keys = ["LyA01flux", "LyA02flux", "LyA14flux", "LyA17flux"]
LyA_pumping = [1217.205, 1217.643, 1216.070, 1215.726]
LyA_file = "V4046Sgr_LyAprof_France2014.csv"
LyA_df = pd.read_csv(LyA_file)

for pumpwave_idx, pump_wave in enumerate(LyA_pumping):
    LyA_tokeep = np.array(LyA_df['ry_out'].loc[(LyA_df['lambda'] >= pump_wave - 0.01)
                                               & (LyA_df['lambda'] <= pump_wave + 0.01)])[0]
    param_dict[LyA_keys[pumpwave_idx]] = LyA_tokeep

grid_dict = grid_from_radius(targ_d, param_dict)
vobs = vrad_to_vobs(grid_dict["rgrid"], grid_dict["phigrid"])

vobs_flat = vobs.flatten() #same as np.ravel with order 'C'
bins, edges = np.histogram(vobs_flat, bins=599, range=(-300., 300.))
inds = np.digitize(vobs_flat, edges)

theta = [sorted_data_table.loc[best_fit_idx]["z"], sorted_data_table.loc[best_fit_idx]["gamma"], sorted_data_table.loc[best_fit_idx]["T"], sorted_data_table.loc[best_fit_idx]["q"], sorted_data_table.loc[best_fit_idx]["rchar"], sorted_data_table.loc[best_fit_idx]["MH2"], sorted_data_table.loc[best_fit_idx]["LyA17flux"]]

full_interp_model = get_H2model(theta, theta_keys, param_dict, inds)

plt.clf()
plt.close()
plt.title("Ru Lupi: Best-Fit H2 Model Emission Lines for [1,7]")
plt.plot(np.arange(len(prog_14_dataflux)),
         prog_14_dataflux,
         color="k", lw=2.)
plt.plot(np.arange(len(prog_14_dataflux)),
         full_interp_model,
         color="teal", lw=1.)
plt.gca().minorticks_on()
plt.show(block=True)