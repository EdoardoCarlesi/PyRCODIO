'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_sub_mass_function.py: wrapper for utilites to compute and plot subhalo mass functions
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
[model_run, dict_model] = cfg.lg_models()
code_run = cfg.simu_runs()
sub_run = cfg.sub_runs()

# Local data path, file names and file format
file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'
file_mah = 'halo_'
form_mah = '.allinfo'

# Local test files
#base_path = '/media/edoardo/data1/DATA/'
#base_path_mah = '/media/edoardo/data1/DATA/'

# Full dataset
base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'
base_path_mah = '/media/edoardo/Elements/CLUES/DATA/trees/2048/'
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'

# Select a subsample from the full catalog to look for local groups
cat_radius = 10.0e+3
cat_center = [50.e+3] * 3

# Read the | Gyr / z / a | time conversion table
time = rf.read_time(data_path)

all_halo_mah = []

# Output file base path, save all relevant halo MAHs here
out_base_pkl = 'saved/lg_mahs_'

# Now loop on all the simulations and gather data
for code in code_run[1:2]:

    for sub in sub_run:
        this_path = base_path + code + '/' + sub + '/'
        this_path_mah = base_path_mah + code + '/' + sub + '/'
        this_ahf = this_path + file_ahf
        this_mah = this_path_mah + file_mah

        # Check that file exists
        if os.path.isfile(this_ahf):
            out_file_pkl = out_base_pkl + code + '_' + sub + '.pkl'

            if os.path.isfile(out_file_pkl):
                print('Loading MAHs from .pkl file...')
                out_f_pkl = open(out_file_pkl, 'rb')
                mahs, lg = pkl.load(out_f_pkl)
                out_f_pkl.close()
                print('Done.')
                print('LG properties: ')
                lg.info()
            else:
                print('Reading AHF file: ', this_ahf)
                halos, halo_df = rf.read_ahf_halo(this_ahf)

                print('Looking for Local Group candidates...')
                this_model = model_run[dict_model[code]]

                # The local group model might eventually find more than one pair in the given volume
                these_lgs = hu.find_lg(halo_df, this_model, cat_center, cat_radius)

                # SELECT ONE LG
                print('Found a LG candidate with properties: ')
                these_lgs[0].info()
                lg = these_lgs[0]

                # Define a gather radius for halos around the LG
                lg_radius = lg.LG1.r() + lg.LG2.r() + lg.r_halos()
                lg_radius = 1.2 * lg_radius

                # TODO there might be more LGs in the area but here I am selecting only one
                #for this_lg in these_lgs:
                lg_center = lg.geo_com()
                id_list = hu.halo_ids_around_center(halos, lg_center, lg_radius)

                print('Reading MAHs for ', len(id_list), ' halos at ', lg_radius, ' kpc/h around the LG center of mass... ')
                mahs = rf.read_mah_halo(id_list, this_mah, time)
                print('Done.')
            
                print('Saving MAHs output to file: ', out_file_pkl)
                out_f_pkl = open(out_file_pkl, 'wb')
                pkl.dump([mahs, lg], out_f_pkl)
                out_f_pkl.close()

                '''
                print(mahs[0].m_t())
                print(mahs[0].m_max())
                print(mahs[0].m_t_norm())
                print(mahs[0].formation_time())
                print(mahs[0].last_major_merger())
                '''
