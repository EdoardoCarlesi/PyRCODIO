'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_density_plot.py: produce a large number of plots centered around individual objects (LGs, halos, clusters...)
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
[model_run, dict_model] = cfg.lg_models()
code_run = cfg.simu_runs()
sub_run = cfg.sub_runs()

# Local data path
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'

# Full dataset
#base_path = '/media/edoardo/data1/DATA/'
base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'
snap_path = 'snapshot_054'

# Plot properties
side_size = 2.0e+3
thickness = 1.0e+3
n_files = 1
reduce_factor = 1.0
units = 'kpc'
part_type = 1

# Production mode: read LG from several snapshots OR single Halos from one snapshots
LG_mode = True; Halo_mode = False
#LG_mode = False; Halo_mode = True

if LG_mode == True:
    
    # Input file containing all the main properties needed to do the plots
    input_all_csv = 'output/lg_pairs_2048.csv'
    data_all = pd.read_csv(input_all_csv)

    code_run = data_all['simu_code'].unique()
    sub_run = data_all['sub_code'].unique()

    # Now loop on all the simulations and gather data
    for code in code_run[:2]:

        for sub in sub_run:
            com = ['Xc_LG', 'Yc_LG', 'Zc_LG']
            this_com = data_all[(data_all['sub_code'] == sub) & (data_all['simu_code'] == code)][com].values

            sub = '%02d' % sub
            this_path = base_path + code + '/' + sub + '/'
            this_snap = this_path + snap_path

            # Check that file exists
            if os.path.isfile(this_snap):

                # Do the plot
                center = this_com[0]
                print('Found: ', this_snap, ', plotting around: ', center)

                pu.find_slab(file_name=this_snap, side=side_size, thick=thickness, center=center)

#                x_tmp, y_tmp = pu.find_slab(this_snap, )

'''
            out_file_pkl = out_base_pkl + code + '_' + sub + '.pkl'
            if os.path.isfile(out_file_pkl):
                print('Loading MAHs from .pkl file...')
                out_f_pkl = open(out_file_pkl, 'rb')
                lg = pkl.load(out_f_pkl)
                out_f_pkl.close()
                print('Done.')
                print('Writing LG properties... ')
                out_all_lgs.write(lg.info()+'\n')
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

                print('Writing LG properties... ')
                out_all_lgs.write(lg.info()+'\n')

                print('Saving MAHs output to file: ', out_file_pkl)
                out_f_pkl = open(out_file_pkl, 'wb')
                pkl.dump(lg, out_f_pkl)
                out_f_pkl.close()

# Close the CSV file containing all the LG infos
out_all_lgs.close()
'''



