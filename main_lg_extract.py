'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_lg_extract.py: find LGs and extract their properties 
'''

import read_files as rf
import halo_utils as hu
import seaborn as sns
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

# Use AHF / csv catalogs
csvAhf = False

# Configure the LG model and subpaths
if csvAhf == True:
    code_run = cfg.gen_runs(0, 80)
    sub_run = cfg.gen_runs(0, 40)
    [model_run, dict_model] = cfg.lg_models()

else:
    [model_run, dict_model] = cfg.lg_models()
    code_run = cfg.simu_runs()
    sub_run = cfg.sub_runs()

# Local data path, file names and file format
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'
file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'

# Full dataset base path
if csvAhf == True:
    base_path = '/media/edoardo/Elements/CLUES/DATA/Particles/512/'
else:
    base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'

# Select a subsample from the full catalog to look for local groups
cat_radius = 10.0e+3
cat_center = [50.e+3] * 3

# Read the | Gyr / z / a | time conversion table
time = rf.read_time(data_path)

all_halo_mah = []

# Output files 
if csvAhf == True:
    out_base_pkl = base_path + 'lg_'
    out_all_lgs_csv = 'output/lg_pairs_512.csv'
else:
    out_base_pkl = 'saved/lg_pair_'
    out_all_lgs_csv = 'output/lg_pairs_2048.csv'


# Write the file header
out_all_lgs = open(out_all_lgs_csv, 'w')
h = hu.Halo('void')
out_all_lgs.write(hu.LocalGroup(h, h).header()+'\n')
out_all_lgs.close()

# Re-open the file in append mode
out_all_lgs = open(out_all_lgs_csv, 'a')

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:
        if csvAhf == True:
            this_ahf = base_path + 'ahf_' + code + '_' + sub + '.csv'
        else:
            this_path = base_path + code + '/' + sub + '/'
            this_ahf = this_path + file_ahf

        # Check that file exists
        if os.path.isfile(this_ahf):
            out_file_pkl = out_base_pkl + code + '_' + sub + '.pkl'

            if os.path.isfile(out_file_pkl):
                print('Loading MAHs from .pkl file...')
                out_f_pkl = open(out_file_pkl, 'rb')
                lg = pkl.load(out_f_pkl)
                lg.code_simu = code
                lg.code_sub = sub
                out_f_pkl.close()
                print('Done.')
                print('Writing LG properties... ')
                out_all_lgs.write(lg.info()+'\n')
            else:
                print('Reading AHF file: ', this_ahf)

                if csvAhf == True:
                    halo_df = pd.read_csv(this_ahf)
                    this_model = model_run[dict_model['GENERIC']]

                else:
                    halos, halo_df = rf.read_ahf_halo(this_ahf)

                    print('Looking for Local Group candidates...')
                    this_model = model_run[dict_model[code]]

                if len(halo_df) > 0:
                    # The local group model might eventually find more than one pair in the given volume
                    these_lgs = hu.find_lg(halo_df, this_model, cat_center, cat_radius)
                else:
                    these_lgs = []

                # Check if there is at least one LG in the selection
                if len(these_lgs) > 0:
                    print('Found a LG candidate with properties: ')

                    # Choose only one LG 
                    these_lgs[0].info()
                    lg = these_lgs[0]
                    lg.code_simu = code
                    lg.code_sub = sub

                    print('Writing LG properties... ')
                    out_all_lgs.write(lg.info(dump=False)+'\n')

                    print('Saving LG output to file: ', out_file_pkl)
                    out_f_pkl = open(out_file_pkl, 'wb')
                    pkl.dump(lg, out_f_pkl)
                    out_f_pkl.close()

# Close the CSV file containing all the LG infos
out_all_lgs.close()




