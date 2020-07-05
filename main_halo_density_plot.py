'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_halo_density_plot.py: produce a large number of plots centered around individual objects (halos, clusters...)
'''

import read_files as rf
import halo_utils as hu
import plot_utils as pu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

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
frac = 1.00
units = 'kpc'
part_type = 1
grid_size = 400
rand_seed = 1
fig_size = 1.2
show_plot = False
do_velocity = False
augment = True
shift = False

# Configure the LG model and subpaths
code_run = cfg.simu_runs()
sub_run = cfg.sub_runs()

# Now loop on all the simulations and gather data
#for code in code_run:
for code in ['09_18']:

    for sub in sub_run:
        input_all_csv = 'output/clusters_' + code + '_' + sub + '.csv' 
        data_all = pd.read_csv(input_all_csv)
        com = ['Xc(6)', 'Yc(7)', 'Zc(8)']
        these_com = data_all[com].values

        this_path = base_path + code + '/' + sub + '/'
        this_snap = this_path + snap_path

        # Check that file exists
        if os.path.isfile(this_snap):

            for i, this_com in enumerate(these_com):

                for z_axis in range(0, 3):

                    # Select particles for the plot, do a selection first. z_axis is orthogonal to the plane of the projection
                    ax0 = (z_axis + 1) % 3
                    ax1 = (z_axis + 2) % 3

                    if shift == True:
                        'Implement this function'
                        center = t.shift(this_com, side_size * 0.5)
                        this_fout = 'output/cluster_shift_' + code + '_' + sub + '.' + str(i) + '.'
                    else:
                        center = this_com
                        this_fout = 'output/cluster_' + code + '_' + sub + '.' + str(i) + '.'

                    print(i, 'Found: ', this_snap, ', plotting around: ', center)
                    print('Plot axes, z=', z_axis, ' ax0=', ax0, ' ax1=', ax1)

                    # Select a slab around a given axis, this function returns a dataframe
                    part_df = pu.find_slab(file_name=this_snap, side=side_size, thick=thickness, center=center, reduction_factor=frac, 
                            z_axis=z_axis, velocity=do_velocity, rand_seed=rand_seed, part_type=4)

                    pu.plot_density(data=part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, 
                            grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size)


