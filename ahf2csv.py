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

# Basic properties of the files to analyze
resolution = '1024'
#resolution = '2048'
#resolution = '512'
file_format = 'AHF'

# Configure the LG model and subpaths
if resolution == '2048':
    code_run = cfg.simu_runs()
    sub_run = cfg.sub_runs()
else:
    code_run = cfg.gen_runs(0, 80)
    sub_run = cfg.gen_runs(0, 40)

# Default value
mpi = False

# Local data path, file names and file format for both standard and NIL format
if resolution == '1024':
    file_ahf = 'snapshot_054.z0.000.AHF_halos'
    mpi = True
elif resolution == '512':
    file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'

base_ahf = 'AHF_output/HESTIA_100Mpc_512_'
format_ahf = '.127.z0.000.AHF_halos'

# Full dataset path
#base_path = '/z/carlesi/CLUES/DATA/512/'
#base_path = '/z/carlesi/HestiaNoam/RE_SIMS/512/DM_ONLY/'
base_path = '/media/edoardo/Elements/CLUES/DATA/LGF/' + resolution + '/'
out_path = '/home/edoardo/CLUES/DATA/LGF/1024/CSV/'

#kpcFac = 1.0
kpcFac = 1.e+3

# Select a subsample from the full catalog to look for local groups
cat_radius = 10.0 * kpcFac
cat_center = [50.0 * kpcFac] * 3

# Output file base path, save all relevant halo MAHs here
out_base = out_path + 'ahf_'

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
   
        # Which path
        if resolution == '2048':
            this_path = base_path + code + '/' + sub + '/'
        else:
            this_path = base_path + code + '_' + sub + '/'

        # Format of the AHF file
        if file_format == 'AHF':
            this_ahf = this_path + file_ahf
        elif file_format == 'NIL':
            this_ahf = this_path + base_ahf + code + '_' + sub + format_ahf

        out_file = out_base + code + '_' + sub + '.csv'

        # Check that file exists
        if os.path.isfile(this_ahf):
            print('Reading AHF file: ', this_ahf)
            halo_df = rf.read_ahf_halo(this_ahf, file_mpi=mpi)
            select_halo_df = hu.find_halos(halo_df, cat_center, cat_radius)
            print('Selected ', len(select_halo_df), ' halos, saving to', out_file)           
            select_halo_df.to_csv(out_file)
        #else:
        #    print('AHF file: ', this_ahf, ' not found.')




