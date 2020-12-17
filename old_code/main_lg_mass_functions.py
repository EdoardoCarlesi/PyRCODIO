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
import pickle
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
v_max = 10.0
d_max = 7000.0

n_bin = 100
m_max = 0.75e+14
m_min = 1.59e+9

masscut_pca = 5.0e+10

mbins = t.gen_bins(nbins=n_bin, binmax=m_max, binmin=m_min) 

all_bins = [[] for i in range(0, n_bin-1)]

#mode='fullbox'
#mode='lgf'
#mode='plots_mf'
mode='plots_pca' 

if mode == 'fullbox':
    
    # FIXME
    lgs_base = '/home/edoardo/CLUES/PyRCODIO/output/lg_fb_new_'
    #lgs_base = '/home/edoardo/CLUES/PyRCODIO/output/LG_HALOS/halos_simu_periodic_'

    all_pca = []

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
            this_x = row[x_col].values

            if os.path.isfile(lg_csv):
                df_tmp = pd.read_csv(lg_csv)

                # Find the halo list
                df_tmp = df_tmp[df_tmp[d_col] < r_max]

                df_tmp[d_col] = df_tmp[x_col_ahf].T.apply(lambda x: t.distance(x, this_x))
                df_tmp = df_tmp[df_tmp[d_col] < r_max]

                pca = t.spatial_pca(data=df_tmp, cols=x_col_ahf)
                print('PCA result: ', pca)
            
                df_tmp = df_tmp[df_tmp['Mvir(4)'] > masscut_pca] 
                pca = t.spatial_pca(data=df_tmp, cols=x_col_ahf)
                print('PCA result mass: ', pca)
        
                all_pca.append(pca)

    f_pca = 'output/pca_fullbox.pkl'
    pickle.dump(all_pca, open(f_pca, 'wb'))

    print('Output PCA written to: ', f_pca)

    '''
                x_bin, y_bin = t.check_df_bins(data=df_tmp, bins=mbins, binmode='log')

                for ib in range(0, n_bin-1):

                    v = y_bin[ib]

                    if v > 0:
                        all_bins[ib].append(v)
                        break

    f_out = 'output/all_mass_functions_fullbox.pkl'
    print('Saving output to: ', f_out)
    pkl.dump(all_bins, open(f_out, 'wb'))

    f_out = 'output/all_mass_functions_x_bin_fullbox.pkl'
    pkl.dump(all_bins, open(f_out, 'wb'))
    '''

elif mode == 'lgf':

    lgs_lgf_base = '/home/edoardo/CLUES/PyRCODIO/output/LG_1024/lg_center_'

    # Read the full list of halos
    df_lgs = pd.read_csv('output/lg_pairs_1024.csv')
    print('TotLen: ', len(df_lgs))

    # Select a halo subsample
    df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
    print('TotLen: ', len(df_lgs))
    box_center = [5.0e+4]*3

    df_lgs[d_col] = df_lgs[x_col].T.apply(lambda x: t.distance(x, box_center))
    df_lgs = df_lgs[df_lgs[d_col] < d_max]
    print('TotLen: ', len(df_lgs))

    n_lgs = len(df_lgs)
    print('Analyzing fullbox LG candidates mass functions ...')
    
    all_pca = []

    # Loop over halo catalogs and lg lists
    for ilg, row in df_lgs.iterrows():
        num = str('%02d' % int(row['simu_code']))
        sub = str('%02d' % int(row['sub_code']))
        this_x = row[x_col].values
        
        this_lg = lgs_lgf_base + num + '_' + sub + '.' + str(ilg) + '.csv'

        if os.path.isfile(this_lg):
            df_tmp = pd.read_csv(this_lg)

            print('Reading ', this_lg)
            df_tmp[d_col] = df_tmp[x_col_ahf].T.apply(lambda x: t.distance(x, this_x))
            df_tmp = df_tmp[df_tmp[d_col] < r_max]

            pca = t.spatial_pca(data=df_tmp, cols=x_col_ahf)
            print('PCA result: ', pca)
            
            df_tmp = df_tmp[df_tmp['Mvir(4)'] > masscut_pca]
            pca = t.spatial_pca(data=df_tmp, cols=x_col_ahf)
            print('PCA result mass: ', pca)
        
            all_pca.append(pca)

    f_pca = 'output/pca_lgf.pkl'
    pickle.dump(all_pca, open(f_pca, 'wb'))

    print('Output PCA written to: ', f_pca)

    '''
        # FIXME this is commented just for the moment when computing spatial distribution of halos

            x_bin, y_bin = t.check_df_bins(data=df_tmp, bins=mbins, binmode='log')

            for ib in range(0, n_bin-1):
                v = y_bin[ib]

                if v > 0:
                    all_bins[ib].append(v)

    #print(all_bins[0])
    #f_out = 'output/all_mass_functions_lgf.pkl'
    f_out = 'output/all_mass_functions_lgf.' + str(n_lgs) + '.pkl'
    print('Saving output to: ', f_out)
    #pkl.dump(all_bins, open(f_out, 'wb'))
    
    # Save the x bin values 
    f_out = 'output/all_mass_functions_x_bin_lgf.pkl'
    pkl.dump(x_bin, open(f_out, 'wb'))
    '''

elif mode == 'plots_mf':
  
    y_fb = 'output/all_mass_functions_fullbox.pkl'
    x_lg = 'output/all_mass_functions_x_bin_lgf.pkl'
    y_lg = 'output/all_mass_functions_lgf.160.pkl'
    #y_lg = 'output/all_mass_functions_lgf.338.pkl'
    #y_fb = 'output/all_mass_functions_lgf.223.pkl'

#    x_bin_fb = pkl.load(open(x_fb, 'rb'))
#    print(x_bin_fb)
    #print(y_bin_fb)
    x_bin_lg = pkl.load(open(x_lg, 'rb'))
    #print(x_bin_lg)
    y_bin_lg = pkl.load(open(y_lg, 'rb'))
    y_bin_fb = pkl.load(open(y_fb, 'rb'))
    #print(y_bin_lg)

    n_bins = len(x_bin_lg)
    perc_bins_y_lg = np.zeros((n_bins, 3))
    perc_bins_y_fb = np.zeros((n_bins, 3))
    percs = [10, 50, 90]

    radius = 5.0
    vol = 4.0 / 3.0 * np.pi * (radius ** 3.0)

    for i in range(0, n_bins):
        this_bin_lg = y_bin_lg[i]
        this_bin_fb = y_bin_fb[i]

        if len(this_bin_lg) > 2:
            pbin_lg = np.percentile(this_bin_lg, percs)
            perc_bins_y_lg[i] = pbin_lg / vol

        if len(this_bin_fb) > 2:
            pbin_fb = np.percentile(this_bin_fb, percs)
            perc_bins_y_fb[i] = pbin_fb / vol

    #plt.axis([9.5, 13.0, -0.1/vol, 3.25/vol])
    plt.axis([9.5, 13.0, -2.6, 0.65])
    plt.xlabel(r'$\log_{10} M \quad [M_{\odot} h^{-1}]$')
    plt.ylabel(r'$\log_{10} n \quad [h^3 Mpc^{-3}]$')
    plt.plot(x_bin_lg, np.log10(perc_bins_y_lg[:,1]), color='black')
    plt.plot(x_bin_lg, np.log10(perc_bins_y_fb[:,1]), color='red')
    plt.fill_between(x_bin_lg, np.log10(perc_bins_y_fb[:, 0]), np.log10(perc_bins_y_fb[:, 2]), facecolor='orange', alpha=0.3)
    plt.fill_between(x_bin_lg, np.log10(perc_bins_y_lg[:, 0]), np.log10(perc_bins_y_lg[:, 2]), facecolor='grey', alpha=0.2)
    plt.tight_layout()

    out_fig = 'output/mass_functions_5mpc_cs_vs_fb.png'
    plt.savefig(out_fig)
    plt.show(block=False)
    plt.pause(3)
    plt.close()
    
elif mode == 'plots_pca':

    f_pca_lgf = 'output/pca_lgf.pkl'
    f_pca_fullbox = 'output/pca_fullbox.pkl'
    pca_lg = pickle.load(open(f_pca_lgf, 'rb'))
    pca_fb = pickle.load(open(f_pca_fullbox, 'rb'))

    pca_lg = np.array(pca_lg)
    pca_fb = np.array(pca_fb)

    sns.distplot(pca_lg[:,1])
    sns.distplot(pca_fb[:,1])
    sns.distplot(pca_lg[:,2])
    sns.distplot(pca_fb[:,2])
    plt.show()



