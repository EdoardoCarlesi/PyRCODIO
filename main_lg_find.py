import numpy as np
import config as cfg
import read_files as rf
import time
import pickle

# Simulation & catalog
file_single='snapshot_054.z0.000.AHF_halos'
box_size = 100000.0
#base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'
base_path = '/home/edoardo/CLUES/DATA/FullBox/'

sub_runs = cfg.gen_runs(0, 4)

#refineSelection = False
refineSelection = True
n_models = 5

m_thresh = 4.0e+11
models = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6']

lg_models, index = cfg.lg_models()

model = lg_models[index['M1']]

kpcFac = 1.0e+3
cat_center = [50.0 * kpcFac] * 3
cat_radius = 200.0 * kpcFac

for run in sub_runs:

    this_ahf = base_path + run + '/' + file_single
    print('Reading file: ', this_ahf)
    halos, halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
    print('Found: ', len(halos), ' objects.')
    halo_select_df = halo_df[halo_df['Mvir(4)'] > m_thresh]
    print('Cut to ', len(halo_select_df), ' objects.')
    these_lgs = hu.find_lg(halo_select_df, this_model, cat_center, cat_radius)
    print('Found a total of', len(these_lgs), ' LGs.')


    # TODO: speedup the algorithm by taking subvolumes!
    # TODO: overlap all the spheres, then use the unique() function of DataFrames to remove duplicates


#    lg_model = lg_models[index[model]]
#    print(index[model], model)
#    this_ahf_file = base_path + sub_path + '/' + file_single
