'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_halo_pca_extract.py: extract halo properties and dump them to a separate file
'''

import numpy as np
import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

# Resimulation mode
#mode='2048'
mode='1024'

base_path = '/media/edoardo/Elements/CLUES/DATA/LGF/' + mode + '/'

# Configure the LG model and subpaths
if mode == '2048':
    code_run = cfg.simu_runs()
    sub_run = cfg.sub_runs()

elif mode == '1024':
    num_run = cfg.gen_runs(0, 80)
    sub_run = cfg.gen_runs(0, 30)

# Read csv list of LGs
lgs_csv = 'output/lg_pairs_' + mode + '.csv'
df_lgs = pd.read_csv(lgs_csv)
print('TotLen: ', len(df_lgs))

x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
x_cols_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
box_center = np.array([5e+4] * 3)

r_max = 10000.0

'''
# Do some LG filtering
v_max = - 0.0
d_max = 7000.0

df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
print('TotLen: ', len(df_lgs))

# Select LGs depending on their distance from the box center
'''

#file_ahf = 'snapshot_full_054.z0.000.AHF_halos'
#file_ahf = 'snapshot_full_054.z0.001.AHF_halos'
file_ahf = 'snapshot_054.z0.000.AHF_halos'

# Loop on all the simulations and gather data
if mode == '2048':
    for code in code_run:

        for sub in sub_run:
            this_path = base_path + code + '/' + sub + '/'
            this_ahf = this_path + file_ahf

            # Check that file exists
            if os.path.isfile(this_ahf):
                print('Reading AHF file: ', this_ahf)
                halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
                halo_select_df = halo_df[halo_df['Mvir(4)'] > m_thresh]
                out_all_halo_csv = out_base + code + '_' + sub + '.csv'
                halo_select_df.to_csv(out_all_halo_csv)

elif mode == '1024':
    
    for ilg, row in df_lgs.iterrows():
        num = str('%02d' % int(row['simu_code']))
        sub = str('%02d' % int(row['sub_code']))

        this_ahf = base_path + num + '_'  + sub + '/' + file_ahf
        this_x = row[x_cols].values
        #print(this_x)

        if os.path.isfile(this_ahf):
            
            halo_df = rf.read_ahf_halo(this_ahf, file_mpi=True)
            halo_df['Dist'] = halo_df[x_cols_ahf].T.apply(lambda x: t.distance(x, this_x))
            halo_df = halo_df[halo_df['Dist'] < r_max]
            #print(halo_df.head())

            f_out = 'output/LG_' + mode + '/lg_center_' + num + '_' + sub + '.' + str(ilg) + '.csv'

            print('Saving to: ', f_out)
            halo_df.to_csv(f_out)


