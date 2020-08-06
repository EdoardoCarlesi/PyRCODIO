import pandas as pd
import numpy as np
import halo_utils as hu
import config as cfg
import read_files as rf
import time
import pickle as pkl

# Simulation & catalog
#file_single='snapshot_054.z0.000.AHF_halos'
file_single='snapshot_full_054.z0.000.AHF_halos'
box_size = 100000.0
#base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'
#base_path = '/home/edoardo/CLUES/DATA/FullBox/'
base_path='/media/edoardo/Elements/CLUES/DATA/2048/00_06/'

sub_runs = cfg.gen_runs(0, 5)

#refineSelection = False
refineSelection = True
n_models = 5

m_thresh = 4.0e+11
models = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6']

lg_models, index = cfg.lg_models()

#this_model = lg_models[index['M1']]
this_model = lg_models[index['GENERIC']]

kpcFac = 1.0e+3
radius = 25.0 * kpcFac
side_buffer = 2.0 * kpcFac

n_sub = int(box_size / radius)
print('Subdivision in ', n_sub, ' subcubes per axis, radius: ', radius, ' and side_buffer: ', side_buffer)

for run in sub_runs:
    all_lgs = []
    this_ahf = base_path + run + '/' + file_single
    print('Reading file: ', this_ahf)
    halos, halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
    print('Found: ', len(halos), ' objects.')

    for ix in range(0, n_sub):
        for iy in range(0, n_sub):
            for iz in range(0, n_sub):
                this_center = np.array([radius * (0.5 + ix), radius * (0.5 + iy), radius * (0.5 + iz)])
                this_radius = radius * 0.5 + side_buffer
                print('Subbox around center: ', this_center, ' rad: ', this_radius)
                these_lgs = hu.find_lg(halo_df, this_model, this_center, this_radius, center_cut=True, search='Box')

                for this_lg in these_lgs:
                    this_lg.code_simu = 'FB'
                    this_lg.code_sub = run

                    all_lgs.append(this_lg)

    this_lg_df = pd.DataFrame(columns = this_lg.header(dump=False))

    for lg in all_lgs:
        this_row = lg.info(dump=False)
        this_series = pd.Series(this_row, index = this_lg_df.columns)
        this_lg_df = this_lg_df.append(this_series, ignore_index=True)

    print(this_lg_df.head())
    this_csv = 'output/lg_fullbox_' + run + '.csv'
    this_lg_df.drop_duplicates(inplace = True)
    this_lg_df.to_csv(this_csv)

    this_pkl = 'output/lg_fullbox_' + run + '.pkl'
    f_pkl = open(this_pkl, 'wb')
    pkl.dump(all_lgs, f_pkl)


