'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com
    
    main_extract_vweb.py: extract v-web data around a given point
'''

import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

# Configure subpaths
code_run = cfg.gen_runs(0, 80)
sub_run = cfg.gen_runs(0, 40)

# Local data path, file names and file format
base_vweb = 'vweb_'
format_vweb = '.000128.Vweb-csv'

# Full dataset
base_path = '/z/carlesi/CLUES/DATA/512/Vweb/'

kpcFac = 1.0
#kpcFac = 1.e+3

# Select a subsample of nodes of the web around a given point
radius = 10.0 * kpcFac
center = [50.0 * kpcFac] * 3

# Output file base path
out_base = base_path + 'vweb_center_0128_'

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
        this_path = base_path + code + '_' + sub + '/'
        this_vweb = this_path + base_vweb + code + '_' + sub + format_vweb
        out_file = out_base + code + '_' + sub + '.pkl'

        # Check that file exists
        if os.path.isfile(this_vweb):
            select_vweb = rf.extract_vweb(file_name=this_web, center=center, radius=radius)
            select_vweb.to_pickle(out_file)

#            print('Reading AHF file: ', this_ahf)
#            halos, halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
#            print('AHF file: ', halo_df.head())
#            select_halo_df = hu.find_halos(halo_df, cat_center, cat_radius)
#            print('Selected ', len(select_halo_df), ' halos, saving to', out_file)           
#            print('AHF file: ', select_halo_df.head())
#            select_halo_df.to_csv(out_file)




