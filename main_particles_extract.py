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
#base_ic_path = '/home/edoardo/CLUES/DATA/ICs/'
base_ic_path = '/z/carlesi/CLUES/ginnungagap/ginnungagap/ICs/'
n_ic_files = 2; ic_root = 'zoom_cf2_512_100.000_'

# Snapshot data
base_snap_path = '/z/carlesi/CLUES/DATA/512/'
n_snap_files = 1; snap_root = 'snapshot_054'
#n_snap_files = 8; snap_root = 'snapshot_127'

# Save files in PKL format
out_extension = '.pkl'

# Look at snapshots OR ICs?
#snapshot = False
snapshot = True

# Plot properties
velocity = True
n_files_ics = 2
n_files_snap = 1
part_type = [1]
rand_state = 1
reduce_factor = 0.33

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:

        if snapshot == False:
            this_file = base_ic_path + ic_root + code + '_' + sub
            this_fout = 'output/ic_' + code + '_' + sub + out_extension
            n_files = n_ic_files
    
        # If this is not a snapshot, try reading ICs
        else:
            this_file = base_snap_path + code + '_' + sub + '/' + snap_root
            this_fout = 'output/snap_' + code + '_' + sub + out_extension
            n_files = n_snap_files
    
        if n_files > 1:
            this_file_test = this_file + '.0'
        else:
            this_file_test = this_file

            # First check if file exists
            if os.path.isfile(this_file_test):
                part_df = rf.read_snap(file_name=this_file, velocity=velocity, part_types=part_type, n_files=n_files)
                #print(part_df.head())

                # Then compress the data and save only a subset of the total particles
                if reduce_factor < 1.0 and len(part_df) > 1000:
                    part_df = part_df.sample(frac=reduce_factor, random_state=rand_state)
            
                if len(part_df) > 1000:
                    print('Saving file to: ', this_fout)
                    part_df.to_pickle(this_fout)
    
