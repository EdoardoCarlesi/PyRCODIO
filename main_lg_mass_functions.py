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


def apply_dist(data=None, center=None, cols=None, box=100.0):
    dist = np.zeros((len(data)))

    dd = data[cols].values - center
    dp = (data[cols].values + box - center) % box 
    
    # FIXME check this implementation of periodic boundary conditions
    for i, d in enumerate(dd):
        d0 = np.sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2])
        d1 = np.sqrt(dp[i][0]**2 + dp[i][1]**2 + dp[i][2]**2)

        dist[i] = np.min([d0, d1])
    
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
#mode = 'export_csv'
#mode = 'read_csv'
#mode = 'analysis'
mode = 'plots'

# Max and min radius to look for stuff
r_max = 15000.0
radii = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0], dtype='float') * 1.e+3 
 
r_str = []
for r in radii:
    r_str.append('R_' + str(r))

df_rm_all = []

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
    this_lgs = lgs_base + i_str + '.csv'  
 
    print('LGS FILE: ', this_lgs)
    df_lgs_orig = pd.read_csv(this_lgs)
    n_lgs_pre = len(df_lgs_orig)

    print('N pre: ', n_lgs_pre)
    df_lgs = df_lgs_orig.drop_duplicates(subset=['M_M31'])
    df_lgs = df_lgs[df_lgs['M_M31'] > m_min]
    df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
    print('N post: ', len(df_lgs))

    # Read AHF and export it to CSV 
    if mode == 'export_csv': 
        print('AHF FILE: ', this_ahf)
        df_ahf = rf.read_ahf_halo(this_ahf)
        print('AHF FILE TO CSV: ', this_ahf_csv)
        df_ahf.to_csv(this_ahf_csv)

    # Read CSV ahf files and extract mass functions around LG candidates
    elif mode == 'read_csv':
        this_ahf_csv = 'output/full.' + i_str + '.' + ahf_file_csv

        print('AHF CSV FILE: ', this_ahf_csv)
        df_ahf = pd.read_csv(this_ahf_csv)
        print(df_ahf[x_col_ahf])
    
        # Initialize a string containing all the radii - this will be the dataframe column with the max halo mass
        radii_str = []
        for rad in radii:
                r_str = 'R_' + str(rad)
                radii_str.append(r_str)

        m_rad = np.zeros((len(radii), len(df_lgs)))

        # Loop over LGs
        for i_lg, row in df_lgs.iterrows():

            x_lg = row[x_col].values
            df_ahf[d_col] = apply_dist(data=df_ahf, center=x_lg, cols=x_col_ahf)
            lg_csv = 'output/LG_HALOS/halos_simu_' + i_str + '_lg_' + str(i_lg) + '.csv'  
        
            df_tmp = pd.DataFrame()
            df_tmp = df_ahf[df_ahf[d_col] < r_max]

            print('Exporting: ', lg_csv)
            df_tmp.to_csv(lg_csv)

    # If only doing the analysis we won't need to dump all the stuff, just read it
    elif mode == 'analysis':
        
        print('Running analysis only...')
        m_rad = np.zeros((n_lgs_pre, len(radii)))
        
        # Loop over the mass functions already extracted 
        for i_lg, row in df_lgs.iterrows():
            lg_csv = 'output/LG_HALOS/halos_simu_' + i_str + '_lg_' + str(i_lg) + '.csv'  

            if os.path.isfile(lg_csv):
                df_tmp = pd.read_csv(lg_csv)
                print('Reading: ', lg_csv)
                    
                for i_r, rad in enumerate(radii):
                    df_tmp = df_tmp[df_tmp[d_col] > 2000.0]
                    m_max = df_tmp[df_tmp[d_col] < rad][m_col].max()
                    m_rad[i_lg, i_r] = m_max

        for i_r, r_col in enumerate(r_str):
            df_lgs_orig[r_col] = m_rad[:, i_r]

        this_lgs_rad = lgs_base + i_str + '_radii.csv'  
        print('Saving to: ', this_lgs_rad)
        df_lgs_orig.to_csv(this_lgs_rad)

    elif mode == 'plots':

        this_lgs_rad = lgs_base + i_str + '_radii.csv'  
        print('Reading from : ', this_lgs_rad)
        df_rad = pd.read_csv(this_lgs_rad)
        df_rad = df_rad[df_rad[r_str[0]] > 0.0]

        df_rm_all.append(df_rad)

        #print(df_rad.head())


if mode == 'plots':
    
    print('Merge all DFs...')
    df_rm = pd.concat(df_rm_all)
    print(df_rm.head())

    slopes = []
    masses = []
    vtans = []
    for rs in r_str[0:6]:
        f_out = 'output/masses_lg_r' + rs + '.png'
        sns.distplot(np.log10(df_rm[rs]))
        title = 'LG masses at R = ' + rs 
        plt.title(title)
        plt.tight_layout()
        print('Saving file to: ', f_out)
        plt.savefig(f_out)
        plt.cla()
        plt.clf()

        med = np.percentile(np.log10(df_rm[rs]), [20, 50, 80])
        medians = '%.3f_{%.3f}^{+%.3f}' % (med[1], med[0]-med[1], med[2]-med[1]) 
        print(medians)
    
        m_cluster = 1.0e+14
        n_cluster = len(df_rm[df_rm[rs] > m_cluster])

        print('Fractions of LG at ', rs, ' from a cluster: ',  n_cluster/len(df_rm), ' (', n_cluster, '/', len(df_rm), ')')
        
        x = np.log10(df_rm[rs]) #/1.0e+12) 
        y = np.log10(df_rm['Vtan']/100.0) 
        #y = np.log10(-df_rm['Vrad']/100.0) 
        #y = np.log10((df_rm['M_MW'] + df_rm['M_M31'])/1.0e+12) 

        f_out = 'output/masses_kde_lg_' + rs + '.png'
        sns.scatterplot(x, y)
        sns.kdeplot(x, y)

        #plt.plot()
        
        slope = np.polyfit(x, y, 1)
        print('params: ', slope)
        #print('slope: ', slope[0])
        slopes.append(slope[0])
        masses.append(med[0])

        plt.xlabel(r'$\log_{10} M_{max}$')
        plt.ylabel(r'$\log_{10} (V_{tan}/100)$')
        plt.title(rs)
        plt.tight_layout()
        '''
        plt.show(block=False)
        plt.pause(3)
        plt.close()
        '''
        plt.savefig(f_out)
        plt.cla()
        plt.clf()

print(masses)
print(slopes)

plt.plot(masses, slopes, color='black')
plt.xlabel(r'$\log_{10} M_{max}$')
plt.ylabel(r'b')
plt.tight_layout()
plt.savefig('output/masses_vs_vtancorr.png')
plt.pause(3)
plt.close()
plt.cla()
plt.clf()

