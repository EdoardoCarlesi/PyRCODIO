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
code_run = cfg.simu_runs()
sub_run = cfg.sub_runs()

#code_run = cfg.gen_runs(0, 10)
#sub_run = cfg.gen_runs(0, 40)

# Local data path, file names and file format
file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'

# Full dataset
base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'

# Select a subsample from the full catalog to look for local groups
cat_radius = 10.0e+3
cat_center = [50.e+3] * 3

# Output file base path, save all relevant halo MAHs here
out_base = 'output/ahf_'

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
        this_path = base_path + code + '/' + sub + '/'
        this_ahf = this_path + file_ahf
        out_file = out_base + code + '_' + sub + '.csv'

        # Check that file exists
        if os.path.isfile(this_ahf):
            print('Reading AHF file: ', this_ahf)
            halos, halo_df = rf.read_ahf_halo(this_ahf)
            select_halo_df = hu.find_halos(halo_df, cat_center, cat_radius)
            print('Selected ', len(select_halo_df), ' halos, saving to', out_file)           
            select_halo_df.to_csv(out_file)




