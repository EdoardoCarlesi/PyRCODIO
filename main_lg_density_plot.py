'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_lg_density_plot.py: produce a large number of plots centered around Local Groups
'''

import read_files as rf
import halo_utils as hu
import plot_utils as pu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

# Configure the LG model and subpaths
[model_run, dict_model] = cfg.lg_models()
code_run = cfg.simu_runs()
sub_run = cfg.sub_runs()

# Local data path
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'

# Full dataset
#base_path = '/media/edoardo/data1/DATA/'
base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'
snap_path = 'snapshot_054'

# Plot properties
side_size = 2.0e+3
thickness = 1.0e+3
n_files = 1
frac = 0.95
units = 'kpc'
part_type = 1
grid_size = 600
rand_seed = 1
fig_size = 12
show_plot = False
do_velocity = False
augment = False

# Input file containing all the main properties needed to do the plots
input_all_csv = 'output/lg_pairs_2048.csv'
data_all = pd.read_csv(input_all_csv)

code_run = data_all['simu_code'].unique()
sub_run = data_all['sub_code'].unique()

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
        com = ['Xc_LG', 'Yc_LG', 'Zc_LG']
        this_com = data_all[(data_all['sub_code'] == sub) & (data_all['simu_code'] == code)][com].values

        sub = '%02d' % sub
        this_path = base_path + code + '/' + sub + '/'
        this_fout = 'output/lg_' + code + '_' + sub + '_'
        this_snap = this_path + snap_path

        # Check that file exists
        if os.path.isfile(this_snap):

            for z_axis in range(0, 3):
                # Select particles for the plot, do a selection first. z_axis is orthogonal to the plane of the projection
                ax0 = (z_axis + 1) % 3
                ax1 = (z_axis + 2) % 3
                center = this_com[0]

                print('Found: ', this_snap, ', plotting around: ', center)
                print('Plot axes, z=', z_axis, ' ax0=', ax0, ' ax1=', ax1)

                # Select a slab around a given axis, this function returns a dataframe
                part_df = pu.find_slab(file_name=this_snap, side=side_size, thick=thickness, center=center, reduction_factor=frac, 
                            z_axis=z_axis, velocity=do_velocity, rand_seed=rand_seed)

                #pu.plot_density(data=part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, 
                #            grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size)
                pu.plot_density(data=part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, 
                            grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size, )


