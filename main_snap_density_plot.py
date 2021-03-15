'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_ics_density_plot.py: produce a large number of plots of density and velocities in the ICs
'''

import numpy as np
import read_files as rf
import plot_utils as pu
import config as cfg
import pandas as pd
import tools as t
import os

# Full dataset
file_path = '/home/edoardo/CLUES/DATA/LGF/512/05_14/snapshot_054'
#file_path = '/home/edoardo/CLUES/DATA/LGF/1024/SNAPS/00_06/snapshot_054'

# Plot properties
center = np.array([50.0e+3] * 3)
side_size = 100.0e+3
thickness = 8.0e+3
n_files = 1
frac = 1.0
units = 'kpc'
part_type = [0,1,2.0,3.0,4,5]
grid_size = 600
rand_seed = 1
fig_size = 10
augment = False
legend = True
hex_plot = True
show_plot = False
velocity = False

part_df_1 = rf.read_snap(file_name=file_path, part_types=[1], n_files=1)
part_df_2 = rf.read_snap(file_name=file_path, part_types=[2], n_files=1)

n_resample = int(len(part_df_1) / 8.0)

part_df_1r = part_df_1.sample(n_resample)
part_df = pd.concat([part_df_1r, part_df_2], axis=0)
slab_part_df = t.find_slab(part_df=part_df, side=side_size, thick=thickness, center=center, reduction_factor=frac, z_axis=2, rand_seed=rand_seed)

for col in ['X', 'Y']:
    slab_part_df[col] = slab_part_df[col].values - 50.0

#print(part_df.head())

this_fout = 'output/dens_plot.png'
# Feed the previously chosen dataframe and plot its 2D density projection
pu.plot_density(data=slab_part_df, axes_plot=[0, 1], file_name=this_fout, show_plot=show_plot, legend=legend, hex_plot=hex_plot,
    grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size, velocity=velocity)






