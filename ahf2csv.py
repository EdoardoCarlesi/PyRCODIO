'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    ahf2csv.py: convert (and compress) AHF halo catalogs to csv files
'''

import read_files as rf
import halo_utils as hu
import seaborn as sns
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

# Configure the LG model and subpaths
#code_run = cfg.simu_runs()
#sub_run = cfg.sub_runs()

code_run = cfg.gen_runs(0, 80)
sub_run = cfg.gen_runs(0, 40)

# Local data path, file names and file format
file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'
#file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'
base_ahf = 'HESTIA_100Mpc_512_'
format_ahf = '.127.z0.000.AHF_halos'

# Full dataset
#base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'
base_path = '/z/carlesi/CLUES/DATA/512/'

#kpcFac = 1.0
kpcFac = 1.e+3

# Select a subsample from the full catalog to look for local groups
cat_radius = 10.0 * kpcFac
cat_center = [50.0 * kpcFac] * 3

# Output file base path, save all relevant halo MAHs here
out_base = 'output/ahf_'

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
#        this_path = base_path + code + '/' + sub + '/'
        this_path = base_path + code + '_' + sub + '/'
        #this_ahf = this_path + file_ahf
        this_ahf = this_path + base_ahf + code + '_' + sub + format_ahf
        out_file = out_base + code + '_' + sub + '.csv'
#        print(this_ahf)

        # Check that file exists
        if os.path.isfile(this_ahf):
            print('Reading AHF file: ', this_ahf)
            halos, halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
#            print('AHF file: ', halo_df.head())
            select_halo_df = hu.find_halos(halo_df, cat_center, cat_radius)
            print('Selected ', len(select_halo_df), ' halos, saving to', out_file)           
#            print('AHF file: ', select_halo_df.head())
            select_halo_df.to_csv(out_file)




