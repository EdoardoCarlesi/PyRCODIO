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
csvAhf = True; hiResAhf = False
#csvAhf = False; hiResAhf = True

# Configure the LG model and subpaths
if csvAhf == True:
    code_run = cfg.gen_runs(0, 100)
    sub_run = cfg.gen_runs(0, 100)
    [model_run, dict_model] = cfg.lg_models()
elif hiResAhf == True:
    [model_run, dict_model] = cfg.lg_models()
    code_run = cfg.simu_runs()
    sub_run = cfg.sub_runs()

# Local data path, file names and file format
data_path = '/home/edoardo/CLUES/PyRCODIO/data/'
file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'

#resolution = '512'
resolution = '1024'

# Full dataset base path
if csvAhf == True:
    #base_path = '/home/edoardo/Elements/CLUES/DATA/CSV/' + resolution + '/'
    base_path = '/home/edoardo/CLUES/DATA/LGF/' + resolution + '/CSV/'
elif hiResAhf == True:
    base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'

# Select a subsample from the full catalog to look for local groups
facKpc = 1.0e+3
radius = 7.0 * facKpc
center = [50.0 * facKpc] * 3

# Read the | Gyr / z / a | time conversion table
time = rf.read_time(data_path)

all_halo_mah = []

# Output files 
if csvAhf == True:
    out_base_pkl = base_path + 'lg_'
    out_all_lgs_csv = 'output/lg_pairs_' + resolution + '.csv'
if hiResAhf == True:
    out_base_pkl = 'saved/lg_pair_'
    out_all_lgs_csv = 'output/lg_pairs_2048.csv'

# All LG dataframe (to be dumped to csv later on), use a dummy halo to properly generate the file header, use a dummy halo to properly generate the file header
h = hu.Halo('void')
cols = hu.LocalGroup(h, h).header(dump=False)

all_lgs_df = pd.DataFrame(columns=cols)

# Now loop on all the simulations and gather data
for code in code_run:

    for sub in sub_run:

        out_file_pkl = 'output/lg_' + resolution + '_' + code + '_' + sub + '.pkl'
        these_lgs = []

        if csvAhf == True:
            this_ahf = base_path + 'ahf_' + code + '_' + sub + '.csv'
        elif hiResAhf == True:
            this_path = base_path + code + '/' + sub + '/'
            this_ahf = this_path + file_ahf

        # Check that file exists
        if os.path.isfile(this_ahf):
            print('Reading AHF file: ', this_ahf)

            if os.path.isfile(out_file_pkl):
                with open(out_file_pkl, 'rb') as f_pkl:
                    these_lgs = pkl.load(f_pkl)
            else:
                if csvAhf == True:
                    halo_df = pd.read_csv(this_ahf)
                    this_model = model_run[dict_model['GENERIC']]

                elif hiResAhf == True:
                    halos, halo_df = rf.read_ahf_halo(this_ahf)
                    this_model = model_run[dict_model[code]]

                print('Looking for Local Group candidates...')

                if len(halo_df) > 0:
                    these_lgs = hu.find_lg(halo_df, this_model, center, radius, center_cut=True, search='Box')
                else:
                    print('Problems with file: ', this_ahf, '. Reading the file results in zero length.')

            # Check if there is at least one LG in the selection
            if len(these_lgs) > 0:

                # Save simu info in LG object
                for i in range(0, len(these_lgs)):
                    these_lgs[i].code_simu = code
                    these_lgs[i].code_sub = sub

                print('Saving LG output to file: ', out_file_pkl)
                out_f_pkl = open(out_file_pkl, 'wb')
                pkl.dump(these_lgs, out_f_pkl)
                out_f_pkl.close()

                for lg in these_lgs:
                    try:
                        print(code, '_', sub, ' --> found a LG candidate with properties: ')
                        lg.info()
    
                        # Append the LG to the dataframe
                        this_row = lg.info(dump=False)
                        this_series = pd.Series(this_row, index=cols)
                        all_lgs_df = all_lgs_df.append(this_series, ignore_index=True)
                    except:
                        print('Error on simu ', code, '_', sub, ', skipping LG')

# Save CSV output file
all_lgs_df.to_csv(out_all_lgs_csv)




