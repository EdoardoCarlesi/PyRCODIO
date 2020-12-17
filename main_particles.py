'''
    Python Routines for COsmology and Data I/O (PyRCODIO) v0.2
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_particles.py = extract features (triaxialities, densities, overdensities etc.) from a set of particles and also analyze them 
'''

import os
import read_files as rf
import halo_utils as hu
import plot_utils as pu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import particles as ps


def particles_extract():
    '''
    Find a series of particles in gadget format files and extract their properties
    '''

    # IC data
    #base_ic_path = '/home/edoardo/CLUES/DATA/ICs/'
    base_ic_path = '/z/carlesi/CLUES/ginnungagap/ginnungagap/ICs/'
    #base_ic_path = '/z/carlesi/HestiaNoam/RE_SIMS/512/DM_ONLY/'
    #00_10/ICs'
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
    reduce_factor = 1.0

    # Now loop on all the simulations and gather data
    for code in code_run:

        for sub in sub_run:

            if snapshot == False:
                #this_file = base_ic_path + code + '_' + sub + '/ICs/' + ic_root + code + '_' + sub
                #this_fout = 'output/ic_hestia_' + code + '_' + sub + out_extension
                this_file = base_ic_path + ic_root + code + '_' + sub
                this_fout = 'output/ic_' + code + '_' + sub + out_extension
                n_files = n_ic_files
        
            # If this is not ICs, try reading a snapshot
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

                # Then compress the data and save only a subset of the total particles
                if reduce_factor < 1.0 and len(part_df) > 1000:
                    part_df = part_df.sample(frac=reduce_factor, random_state=rand_state)
                
                if len(part_df) > 1000:
                    print('Saving file to: ', this_fout)
                    part_df.to_pickle(this_fout)

     return None


def particles_features():
    '''
    Read a bunch of files that contain processed information about particles, then do 
    '''

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
    return None


if __name__ == '__main__':
    ''' Wrapper '''
    
    global code_run
    global sub_run

    # Configure the LG model and subpaths
    code_run = cfg.gen_runs(0, 1)
    sub_run = cfg.gen_runs(0, 5)

    particles_extract()
    particles_features()




