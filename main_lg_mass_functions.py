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

def apply_dist(data=None, center=None, cols=None):
    dist = np.zeros((len(data)))
    for i, d in enumerate((data[cols].values - center)):
        dist[i] = np.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
        #print(i, dist[i])
    return dist


# Base path for halo catalogs and identified LGs
#ahf_base = '/media/edoardo/Elements/CLUES/DATA/FullBox/'
ahf_base = '/home/edoardo/CLUES/DATA/FullBox/'
ahf_file = 'snapshot_054.z0.000.AHF_halos'
ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.csv'
#lgs_base = 'output/lg_fullbox_'
lgs_base = 'output/lg_fb_new_'

# Loop over these catalogs
i_ini = 0
i_end = 5

# Read AHF from original and export to CSV or not
export_csv = False
analyze_only = True

# Max and min radius to look for stuff
r_max = 8000.0
radii = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0], dtype='float') * 1.e+3 

# Column coordinates
x_col = ['Xc_LG', 'Yc_LG', 'Zc_LG']
x_col_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
m_col = ['Mvir(4)']
d_col = 'Dist'
m_min = 5.0e+11
v_max = -1.0

# Loop over halo catalogs and lg lists
for i_cat in range(i_ini, i_end):
    i_str = '%02d' % i_cat

    this_ahf = ahf_base + i_str + '/' + ahf_file
    this_ahf_csv = 'output/full.' + i_str + '.' + ahf_file_csv
    this_lgs = lgs_base + i_str + '.csv'  
    this_lgs_rad = lgs_base + i_str + '_radii.csv'  

    if export_csv == True: 
        print('AHF FILE: ', this_ahf)
        df_ahf = rf.read_ahf_halo(this_ahf)
        print('AHF FILE TO CSV: ', this_ahf_csv)
        df_ahf.to_csv(this_ahf_csv)
    elif analyze_only == False:
        print('AHF CSV FILE: ', this_ahf_csv)
        df_ahf = pd.read_csv(this_ahf_csv)
        print(df_ahf[x_col_ahf])
    
    print('LGS FILE: ', this_lgs)
    df_lgs_orig = pd.read_csv(this_lgs)
    n_lgs_pre = len(df_lgs_orig)

    print('N pre: ', n_lgs_pre)
    df_lgs = df_lgs_orig.drop_duplicates(subset=['M_M31'])
    df_lgs = df_lgs[df_lgs['M_M31'] > m_min]
    df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
    print('N post: ', len(df_lgs))

    # Skip the AHF file reading and read the individual mass functions only, which have already been stored
    if analyze_only == False:
    
        # initialize a string containing all the radii - this will be the dataframe column with the max halo mass
        radii_str = []
        for rad in radii:
                r_str = 'R_' + str(rad)
                radii_str.append(r_str)

        m_rad = np.zeros((len(radii), len(df_lgs)))

        for i_lg, row in df_lgs.iterrows():

            x_lg = row[x_col].values
            print(i_lg, x_lg)

            df_ahf[d_col] = apply_dist(data=df_ahf, center=x_lg, cols=x_col_ahf)
            lg_csv = 'output/halos_simu_' + i_str + '_lg_' + str(i_lg) + '.csv'  
        
            df_tmp = pd.DataFrame()
            df_tmp = df_ahf[df_ahf[d_col] < r_max]

            print('Exporting: ', lg_csv)
            df_tmp.to_csv(lg_csv)

    # If only doing the analysis we won't need to dump all the stuff, just read it
    else:
        
        print('Running analysis only...')
        
        r_str = []
        for r in radii:
            r_str.append('R_' + str(r))

        m_rad = np.zeros((n_lgs_pre, len(radii)))
        
        for i_lg, row in df_lgs.iterrows():
            lg_csv = 'output/halos_simu_' + i_str + '_lg_' + str(i_lg) + '.csv'  

            if os.path.isfile(lg_csv):
                df_tmp = pd.read_csv(lg_csv)
                #print('Reading: ', lg_csv)
                    
                for i_r, rad in enumerate(radii):
                    df_tmp = df_tmp[df_tmp[d_col] > 2000.0]
                    m_max = df_tmp[df_tmp[d_col] < rad][m_col].max()
                    m_rad[i_lg, i_r] = m_max

                #print('LG(', i_lg, ') ', m_rad[i_lg, :])
                
        print(m_rad[0:10, :])
        print(m_rad.shape)

        for i_r, r_col in enumerate(r_str):
            df_lgs_orig[r_col] = m_rad[:, i_r]
            #print(m_rad[:, i_r])

        print('Saving to: ', this_lgs_rad)
        df_lgs_orig.to_csv(this_lgs_rad)



