'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_particles_extract.py: extract particle data within a high-res region in the ICs or in a snapshot
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
code_run = cfg.gen_runs(0, 80)
sub_run = cfg.gen_runs(0, 40)

# IC data
base_path = '/z/carlesi/CLUES/ginnungagap/ginnungagap/ICs/'
ic_path = 'zoom_cf2_512_100.000_'
ic_extension = '.1'

# Snapshot data


# Plot properties
n_files = 1
frac = 1.0
units = 'kpc'
part_type = 1
grid_size = 50
log_plot = False
rand_seed = 1
fig_size = 1.2
show_plot = False
velocity = True
augment = False
legend = False
hex_plot = False
vel_components = ['Vx', 'Vy', 'Vz']

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
        com = ['Xc_LG', 'Yc_LG', 'Zc_LG']
        this_ic = base_path + ic_path + code + '_' + sub + ic_extension
        this_fout = 'output/ic_' + ic_path + code + '_' + sub 

        if os.path.isfile(this_ic):
            part_df = rf.read_snap(file_name=this_ic, velocity=velocity, part_types=part_type, n_files=1)
            center = t.particles_com(part_df)


