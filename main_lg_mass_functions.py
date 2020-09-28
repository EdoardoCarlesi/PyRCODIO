import seaborn as sns
import matplotlib.pyplot as plt
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
ahf_base = '/home/edoardo/CLUES/DATA/FullBox/'
#ahf_file = 'snapshot_054.z0.000.AHF_halos'

#ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.csv'; lg_csv_file = 'output/LG_HALOS/halos_simu_'
ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.csv'; lg_csv_file = 'output/LG_HALOS/halos_simu_periodic_'
#ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.periodic.csv'; lg_csv_file = 'output/LG_HALOS/halos_simu_periodic_'

# Loop over these catalogs
i_ini = 0
i_end = 5

# Max and min radius to look for stuff
r_max = 5000.0
box_size = 1.0e+5
radii = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], dtype='float') * 1.e+3 

# Column coordinates
x_col = ['Xc_LG', 'Yc_LG', 'Zc_LG']
x_col_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
m_col = ['Mvir(4)']
d_col = 'Dist'
m_min = 4.0e+11
v_max = -0.0
d_max = 5000.0

n_bin = 100
m_max = 0.75e+14
m_min = 1.59e+9

mbins = t.gen_bins(nbins=n_bin, binmax=m_max, binmin=m_min) 

all_bins = [[] for i in range(0, n_bin-1)]

#mode='fullbox'
mode='lgf'
#mode='plots' 

if mode == 'fullbox':

    # Loop over halo catalogs and lg lists
    for i_cat in range(i_ini, i_end):

        i_str = '%02d' % i_cat
        print('Analyzing fullbox LG candidates mass functions ...')
        this_lgs = lgs_base + i_str + '.csv'  
 
        print('LGS FILE: ', this_lgs)
        df_lgs_orig = pd.read_csv(this_lgs)
        n_lgs_pre = len(df_lgs_orig)

        print('N pre: ', n_lgs_pre)
        df_lgs = df_lgs_orig.drop_duplicates(subset=['M_M31'])
        df_lgs = df_lgs[df_lgs['M_M31'] > m_min]
        df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
        print('N post: ', len(df_lgs))

        # Loop over the mass functions already extracted 
        for i_lg, row in df_lgs.iterrows():
            lg_csv = lg_csv_file + i_str + '_lg_' + str(i_lg) + '.csv'  

            if os.path.isfile(lg_csv):
                df_tmp = pd.read_csv(lg_csv)

                # Find the halo list
                df_tmp = df_tmp[df_tmp[d_col] < r_max]

                x_bin, y_bin = t.check_df_bins(data=df_tmp, bins=mbins, binmode='log')

                for ib in range(0, n_bin-1):

                    v = y_bin[ib]

                    if v > 0:
                        all_bins[ib].append(v)

    f_out = 'output/all_mass_functions_fullbox.pkl'
    print('Saving output to: ', f_out)
    pkl.dump(all_bins, open(f_out, 'wb'))

elif mode == 'lgf':

    lgs_lgf_base = '/home/edoardo/CLUES/PyRCODIO/output/LG_1024/lg_center_'

    # Read the full list of halos
    df_lgs = pd.read_csv('output/lg_pairs_1024.csv')
    print('TotLen: ', len(df_lgs))

    # Select a halo subsample
    df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
    print('TotLen: ', len(df_lgs))

    x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
    df_lgs['Dist'] = df_lgs[x_cols].T.apply(lambda x: t.distance(x, box_center))
    df_lgs = df_lgs[df_lgs['Dist'] < d_max]
    print('TotLen: ', len(df_lgs))

    print('Analyzing fullbox LG candidates mass functions ...')
    
    # Loop over halo catalogs and lg lists
    for ilg, row in df_lgs.iterrows():
        num = str('%02d' % int(row['simu_code']))
        sub = str('%02d' % int(row['sub_code']))
        
        this_lg = lgs_lgf_base + num + '_' + sub + '.' + str(ilg) + '.csv'

        if os.path.isfile(this_lg):
            df_tmp = pd.read_csv(this_lg)
            #df_tmp = pkl.load(open(this_lg, 'rb'))
            print(df_tmp)
            #print(df_tmp.head())

            # Find the halo list
            #df_tmp = df_tmp[df_tmp[d_col] < r_max]

            '''
                x_bin, y_bin = t.check_df_bins(data=df_tmp, bins=mbins, binmode='log')

                for ib in range(0, n_bin-1):

                    v = y_bin[ib]

                    if v > 0:
                        all_bins[ib].append(v)

    f_out = 'output/all_mass_functions_lgf.pkl'
    print('Saving output to: ', f_out)
    pkl.dump(all_bins, open(f_out, 'wb'))

            '''
