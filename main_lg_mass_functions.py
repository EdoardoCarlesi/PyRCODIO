import read_files as rf
import halo_utils as hu
import plot_utils as pu
import config as cfg
import pickle as pkl
import pandas as pd
import numpy as np
import tools as t
import swifter
import os

# Base path for halo catalogs and identified LGs
ahf_base = '/media/edoardo/Elements/CLUES/DATA/FullBox/'
ahf_file = 'snapshot_054.z0.000.AHF_halos'
ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.csv'
lgs_base = 'output/lg_fullbox_'

# Loop over these catalogs
i_ini = 0
i_end = 1

# Read AHF from original and export to CSV or not
export_csv = False

# Max and min radius to look for stuff
R_max = 8000.0

# Column coordinates
x_col = ['Xc_LG', 'Yc_LG', 'Zc_LG']
x_col_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
m_col = ['Mvir(4)']

# Loop over halo catalogs and lg lists
for i_cat in range(i_ini, i_end):
    i_str = '%02d' % i_cat

    this_ahf = ahf_base + i_str + '/' + ahf_file
    this_ahf_csv = 'output/' + ahf_file_csv
    this_lgs = lgs_base + i_str + '.csv'  

    if export_csv == True: 
        print('AHF FILE: ', this_ahf)
        df_ahf = rf.read_ahf_halo(this_ahf)
        print('AHF FILE TO CSV: ', this_ahf_csv)
        df_ahf.to_csv(this_ahf_csv)
    else:
        print('AHF CSV FILE: ', this_ahf_csv)
        df_ahf = pd.read_csv(this_ahf_csv)
        print(df_ahf[x_col_ahf])

    print('LGS FILE: ', this_lgs)
    df_lgs = pd.read_csv(this_lgs)
    print(df_lgs[x_col])
    
    for i_lg, row in df_lgs.iterrows():

        x_lg = row[x_col].values
        print(i_lg, x_lg)

        lg_csv = 'output/halos_simu_' + i_str + '_lg_' + str(i_lg) + '.csv'  
        #df_ahf['Dist'] = df_ahf[x_col_ahf].T.apply(lambda x: t.distance(x, x_lg))  
        df_ahf['Dist'] = df_ahf[x_col_ahf].swifter.apply(lambda x: t.distance(x, x_lg), axis=1)  
        df_tmp = pd.DataFrame()
        
        print(df_ahf['Dist'])

    '''
        df_tmp['R'] = df_ahf[df_ahf['Dist'] < R_max]
        print(df_tmp.head())

        print('Exporting: ', lg_csv)
        df_tmp.to_csv(lg_csv)
        
        # Find objects within r_min
        for rad in r_max:
            r_str = 'R_' + str(rad)

            

    '''



