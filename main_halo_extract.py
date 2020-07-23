'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_halo_pca_extract.py: extract halo properties and dump them to a separate file
'''

import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

# Configure the LG model and subpaths
code_run = cfg.simu_runs()
sub_run = cfg.sub_runs()

# Initialize this as an empty dataframe
all_halo_df = pd.DataFrame()

# Output file base path, save all relevant halo MAHs here
out_base_pkl = 'saved/lg_pair_'

# Output file containing all the LG properties 
out_base = 'output/clusters_'

# File and path specifications
base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'
file_ahf = 'snapshot_full_054.z0.000.AHF_halos'
#file_ahf = 'snapshot_full_054.z0.001.AHF_halos'

# Cluster mass threshold
m_thresh = 1.e+14

# Loop on all the simulations and gather data
#for code in code_run:
for code in ['01_12']:
#for code in ['45_17']:

    for sub in sub_run[:]:
        this_path = base_path + code + '/' + sub + '/'
        this_ahf = this_path + file_ahf
        #print(this_ahf)

        # Check that file exists
        if os.path.isfile(this_ahf):
            print('Reading AHF file: ', this_ahf)
            halos, halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
            halo_select_df = halo_df[halo_df['Mvir(4)'] > m_thresh]
            out_all_halo_csv = out_base + code + '_' + sub + '.csv'
            halo_select_df.to_csv(out_all_halo_csv)


