'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_particles_features.py = extract features (triaxialities, densities, overdensities etc.) from a set of particles
'''

import os
import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t

# Configure the LG model and subpaths
code_run = cfg.gen_runs(0, 2)
sub_run = cfg.gen_runs(0, 10)

# Data (we assume it's all pkl / csv files, extracted in the same folder)
base_path = '/home/edoardo/CLUES/DATA/Particles/512/'

# Plot properties
velocity = True
rand_state = 1

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:

        this_ic = base_path + '/ic_' + code + '_' + sub + '.pkl'
        this_snap = base_path + '/snap_' + code + '_' + sub + '.pkl'
        this_lg = base_path + '/lg_' + code + '_' + sub + '.pkl'
    
        # Check if files exist
        if os.path.isfile(this_ic) and os.path.isfile(this_snap) and os.path.isfile(this_lg):
            part_ic = pd.read_pickle(this_ic)
            part_snap = pd.read_pickle(this_snap)
            lg = pd.read_pickle(this_lg)

            lg.info(dump=True)
            print('ICs head:')
            print(part_ic.head())
            print('snap head:')
            print(part_snap.head())
#            halo_ahf = pd.
            #print(this_ic)
            #print(part_ic.head())
            #f_ic.close()
