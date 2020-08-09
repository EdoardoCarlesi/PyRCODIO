'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_lg_fb_find: find local groups in full box simulations
'''

import dask.dataframe as dd
import pandas as pd
import numpy as np
import halo_utils as hu
import config as cfg
import read_files as rf
import tools as t
import time
import pickle as pkl

# Choose catalog type
#use_ahf = True; use_rs = False
use_ahf = False; use_rs = True

# Simulation & catalog
if use_ahf == True:
    #file_single='snapshot_054.z0.000.AHF_halos'
    file_single='snapshot_full_054.z0.000.AHF_halos'
    base_file_out = 'output/lg_fullbox_'    
    box_size = 100000.0
    #base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'
    #base_path = '/home/edoardo/CLUES/DATA/FullBox/'
    base_path='/media/edoardo/Elements/CLUES/DATA/2048/00_06/'
    sub_runs = cfg.gen_runs(0, 5)

elif use_rs == True:
    box_size = 2500000.0
    file_single = '.part'
    base_file_out = 'output/lg_fullbox_rs_'    
    base_path = '/home/edoardo/CLUES/DATA/RS/out_79_csv/'
    #base_path = '/z/carlesi/STORE/MultiDark/RockStarCSV/BigMD_3840_Planck1/out_79_csv/'
    sub_runs = []

    '''
    n_parts = 10
    for i in range(0, n_parts):
        sub = '%04d' % i
        sub_runs.append(sub)
    '''

    sub_runs = ['0011']

lg_models, index = cfg.lg_models()
this_model = lg_models[index['GENERIC']]

kpcFac = 1.0e+3
radius = 15.0 * kpcFac
side_buffer = 2.0 * kpcFac

n_sub_x = np.ceil(box_size / radius)
n_sub_y = n_sub_x
n_sub_z = n_sub_y

print('Subdivision in ', n_sub_x, ' subcubes per axis, radius: ', radius, ' and side_buffer: ', side_buffer)

for run in sub_runs:

    all_lgs = []
    if use_ahf:
        this_ahf = base_path + run + '/' + file_single
        print('Reading file: ', this_ahf)
        halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
        print('Found: ', len(halo_df), ' objects.')

        x_min = 0.0
        y_min = 0.0
        z_min = 0.0

    elif use_rs:
        this_rs = base_path + run + file_single
        
        print('Reading file: ', this_rs)
        halo_df = pd.read_csv(this_rs)
        #halo_df = rf.read_rs_halo(read_file=this_rs, with_dask=True)
        halo_df.columns = t.header_rs2ahf(halo_df.columns)
        #print('Found: ', len(halo_df), ' objects.')

        # The file is divided in sub files, determine the extension
        x_max = halo_df['Xc(6)'].max() * kpcFac
        y_max = halo_df['Yc(7)'].max() * kpcFac
        z_max = halo_df['Zc(8)'].max() * kpcFac
        x_min = halo_df['Xc(6)'].min() * kpcFac
        y_min = halo_df['Yc(7)'].min() * kpcFac
        z_min = halo_df['Zc(8)'].min() * kpcFac

        print(x_min, y_min, z_min)
        print(x_max, y_max, z_max)

        n_sub_x = int(np.ceil((x_max - x_min) / radius))
        n_sub_y = int(np.ceil((y_max - y_min) / radius))
        n_sub_z = int(np.ceil((z_max - z_min) / radius))

        print('N subcubes ', n_sub_x, n_sub_y, n_sub_z)

    for ix in range(0, n_sub_x):
        for iy in range(0, n_sub_y):
            for iz in range(0, n_sub_z):
                this_center = np.array([radius * (0.5 + ix)+x_min, radius * (0.5 + iy)+y_min, radius * (0.5 + iz)+z_min])
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
    this_csv = base_file_out + run + '.csv'
    this_lg_df.drop_duplicates(inplace = True)
    this_lg_df.to_csv(this_csv, float_format='%.3f')

    this_pkl = base_file_out + run + '.pkl'
    f_pkl = open(this_pkl, 'wb')
    pkl.dump(all_lgs, f_pkl)


