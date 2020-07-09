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
frac = 1.0
units = 'kpc'
part_type = [1, 2]
grid_size = 600
rand_seed = 1
fig_size = 15
show_plot = False
velocity = False
augment = False
legend = True
vel_components = ['Vx', 'Vy', 'Vz']

# Input file containing all the main properties needed to do the plots
input_all_csv = 'output/lg_pairs_2048.csv'
data_all = pd.read_csv(input_all_csv)

code_run = data_all['simu_code'].unique()
sub_run = data_all['sub_code'].unique()

# Now loop on all the simulations and gather data
for code in code_run[:1]:

    for sub in sub_run[:2]:
        com = ['Xc_LG', 'Yc_LG', 'Zc_LG']
        this_com = data_all[(data_all['sub_code'] == sub) & (data_all['simu_code'] == code)][com].values

        sub = '%02d' % sub
        this_path = base_path + code + '/' + sub + '/'
        this_fout = 'output/lg_' + code + '_' + sub + '_'
        this_snap = this_path + snap_path

        # Check that file exists
        if os.path.isfile(this_snap):
            center = this_com[0]
            print('Found: ', this_snap, ', plotting around: ', center)

            # Select a slab around a given axis, this function returns a dataframe
            part_df = rf.read_snap(file_name=this_snap, velocity=velocity, part_type=1, n_files=1)

            # Select particles for the plot, do a selection first. z_axis is orthogonal to the plane of the projection
            for z_axis in range(0, 3):
                ax0 = (z_axis + 1) % 3
                ax1 = (z_axis + 2) % 3

                print('Plot axes, z=', z_axis, ' ax0=', ax0, ' ax1=', ax1)

                slab_part_df = t.find_slab(part_df=part_df, side=side_size, thick=thickness, center=center, reduction_factor=frac, z_axis=z_axis, rand_seed=rand_seed)

                # Feed the previously chosen dataframe and plot its 2D density projection
                pu.plot_density(data=slab_part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, legend=legend,
                            grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size, velocity=velocity, vel=vel_components)


