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
import particles as ps

# Configure the LG model and subpaths
code_run = cfg.gen_runs(0, 1)
sub_run = cfg.gen_runs(0, 5)

# Data (we assume it's all pkl / csv files, extracted in the same folder)
base_path = '/home/edoardo/CLUES/DATA/Particles/512/'

# Plot properties
velocity = True
radius = 1.0e+3

rand_state = 1
resize = 0.05

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:

        this_ic = base_path + '/ic_' + code + '_' + sub + '.pkl'
        this_snap = base_path + '/snap_' + code + '_' + sub + '.pkl'
        this_lg = base_path + '/lg_' + code + '_' + sub + '.pkl'
        this_ahf = base_path + '/ahf_' + code + '_' + sub + '.csv'
    
        # Check if files exist
        if os.path.isfile(this_ic) and os.path.isfile(this_snap) and os.path.isfile(this_lg):
            print('Found: ', this_ic, this_snap, this_lg)

            # Read the ICs, the snapshots, the local group and AHF catalogs
            part_ic = pd.read_pickle(this_ic)
            part_snap = pd.read_pickle(this_snap)
            lg = pd.read_pickle(this_lg)
            print('Files have been read in. Extracting particles from a set of ', len(part_snap))    
        
            if resize < 1.0:
                part_snap = part_snap.sample(frac=resize, random_state=rand_state)

            m31 = lg.LG1
            mw = lg.LG2

            # Select a few particles around the very center of the LG
            print('Looking for M31 particles, resampled to ', len(part_snap), m31.r() * 0.5)
            part_m31 = ps.find_particles(part_snap, m31.pos(), m31.r())
            print('Done.\n', part_m31.head(), '\n N sampled particles :', len(part_m31))

            # Try to find the selected particles in the ICs
            m31_ics = ps.match_particle_ids(data=part_ic, ids=part_m31['ID'].values)

            #print(part_ic['ID'])

            print(m31_ics)
            

            '''
            # Sanity check
            print(m31.r())
            print(m31.pos())
            '''

            '''
            # Find the LG candidates in the AHF catalog - not really needed but who knows
            ahf = pd.read_csv(this_ahf)
            lg_halos = hu.find_halos(ahf, lg.geo_com(), radius)
            lg_halos.sort_values(by=['Mvir(4)'], inplace=True, ascending=False)
            m31 = hu.Halo(lg_halos.iloc[0])
            print(m31.info())
            '''

            '''
            # Dump some information
            lg.info(dump=True)
            print('ICs head:')
            print(part_ic.head())
            print('snap head:')
            print(part_snap.head())
            print('AHF head:')
            print(ahf.head())
            print('LG head:')
            print(lg_halos.head(2))
            '''

            




