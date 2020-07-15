'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_ics_density_plot.py: produce a large number of plots of density and velocities in the ICs
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
code_run = cfg.gen_runs(0, 1)
sub_run = cfg.gen_runs(0, 1)

# Full dataset
#base_path = '/media/edoardo/data1/DATA/'
#base_path = '/media/edoardo/Elements/CLUES/DATA/ICs/'
#base_path = '/z/carlesi/CLUES/ginnungagap/ginnungagap/ICs/'
base_path = '/z/carlesi/HestiaNoam/RE_SIMS/512/DM_ONLY/'
ic_path = 'zoom_cf2_512_100.000_'
ic_extension = '.1'

# Plot properties
side_size = 10.0e+3
thickness = 5.0e+3
n_files = 1
frac = 1.0
units = 'kpc'
part_type = 1
grid_size = 50
log_plot = False
rand_seed = 1
fig_size = 1.2
show_plot = False
velocity = False
augment = False
legend = False
hex_plot = False
vel_components = ['Vx', 'Vy', 'Vz']

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
        com = ['Xc_LG', 'Yc_LG', 'Zc_LG']
        #this_ic = base_path + ic_path + code + '_' + sub + ic_extension
        this_ic = base_path + code + '_' + sub + '/' + ic_path + code + '_' + sub + ic_extension
        this_fout = 'output/ic_' + ic_path + code + '_' + sub 

        if os.path.isfile(this_ic):
            part_df = rf.read_snap(file_name=this_ic, velocity=velocity, part_types=part_type, n_files=1)
            center = t.particles_com(part_df)

            print('Found: ', this_ic, ', plotting around: ', center)

            for z_axis in range(0, 3):
                # Select particles for the plot, do a selection first. z_axis is orthogonal to the plane of the projection
                ax0 = (z_axis + 1) % 3
                ax1 = (z_axis + 2) % 3

                print('Plot axes, z=', z_axis, ' ax0=', ax0, ' ax1=', ax1)

                # Select a slab around a given axis, this function returns a dataframe
                slab_part_df = t.find_slab(part_df=part_df, side=side_size, thick=thickness, center=center, reduction_factor=frac, z_axis=z_axis, rand_seed=rand_seed)

                try:
                    # Feed the previously chosen dataframe and plot its 2D density projection
                    pu.plot_density(data=slab_part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, legend=legend, hex_plot=hex_plot,
                        grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size, velocity=velocity, vel=vel_components)
                except:
                    print('Error. Could not plot density slabs for file: ', this_ic)

