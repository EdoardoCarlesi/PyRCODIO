'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_density_plot.py: produce a large number of plots centered around individual objects (LGs, halos, clusters...)
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
frac = 0.3
units = 'kpc'
part_type = 1
doVelocity = True

# Production mode: read LG from several snapshots OR single Halos from one snapshots
LG_mode = True; Halo_mode = False
#LG_mode = False; Halo_mode = True

if LG_mode == True:
    
    # Input file containing all the main properties needed to do the plots
    input_all_csv = 'output/lg_pairs_2048.csv'
    data_all = pd.read_csv(input_all_csv)

    code_run = data_all['simu_code'].unique()
    sub_run = data_all['sub_code'].unique()

    # Now loop on all the simulations and gather data
    for code in code_run[:2]:

        for sub in sub_run:
            com = ['Xc_LG', 'Yc_LG', 'Zc_LG']
            this_com = data_all[(data_all['sub_code'] == sub) & (data_all['simu_code'] == code)][com].values

            sub = '%02d' % sub
            this_path = base_path + code + '/' + sub + '/'
            this_fout = 'output/lg_' + code + '_' + sub + '_'
            this_snap = this_path + snap_path

            # Check that file exists
            if os.path.isfile(this_snap):

                # TODO loop on the axis
                # Select particles for the plot, do a selection first. z_axis is orthogonal to the plane of the projection
                z_axis = 2
                center = this_com[0]
                print('Found: ', this_snap, ', plotting around: ', center)

                # Select a slab around a given axis, this function returns a dataframe
                part_df = pu.find_slab(file_name=this_snap, side=side_size, thick=thickness, center=center, reduction_factor=frac, z_axis=z_axis, velocity=doVelocity, rand_seed=1)

                pu.plot_density(data=part_df, axes_plot=[0, 1], file_name=this_fout)

# In this kind of production mode we select individual halos from a single snapshot
elif Halo_mode == True:

    input_halo_csv = 'specify_a_file'
    data_all = pd.read_csv(input_halo_csv)
    





