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


def apply_dist(data=None, center=None, cols=None):
    dist = np.zeros((len(data)))

    dd = data[cols].values - center
    
    for i, d in enumerate(dd):
        dist[i] = np.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
    
    return dist

#simu = 'fullbox'
simu = 'lgf'

print('SimuType: ', simu)

# Base path for halo catalogs and identified LGs
if simu == 'fullbox':
    #ahf_base = '/media/edoardo/Elements/CLUES/DATA/FullBox/'
    ahf_base = '/home/edoardo/CLUES/DATA/FullBox/'
    ahf_file = 'snapshot_054.z0.000.AHF_halos'
    #ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.csv'; lg_csv_file = 'output/LG_HALOS/halos_simu_'
    ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.periodic.csv'; lg_csv_file = 'output/LG_HALOS/halos_simu_periodic_'
    #ahf_file_csv = 'snapshot_054.z0.000.AHF_halos.periodic.csv'; lg_csv_file = 'output/LG_HALOS/halos_simu_periodic_25mpc_'
    #lgs_base = 'output/lg_fullbox_'
    lgs_base = 'output/lg_fb_new_'

    # Loop over these catalogs
    i_ini = 0
    i_end = 1
    r_max = 25000.0
    box_size = 1.0e+5
    #radii = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 
    #               21.0, 22.0, 23.0, 24.0, 25.0], dtype='float') * 1.e+3 
    radii = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], dtype='float') * 1.e+3 

elif simu == 'lgf':
    i_ini = 0
    i_end = 1

    ahf_base = '/home/edoardo/CLUES/DATA/FullBox/'
    lg_csv_file = 'output/LG_1024/lg_center_'
    lgs_base = 'output/lg_pairs_1024'
    box_size = 1.0e+5
    r_max = 10000.0
    d_max = 6000.0
    radii = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], dtype='float') * 1.e+3 

# Read AHF from original and export to CSV or not
#mode = 'export_csv'
#mode = 'read_csv'
#mode = 'analysis'
mode = 'plots'

r_str = []
for r in radii:
    r_str.append('R_' + str(r))

df_rm_all = []

# These data will be used for the 'plots'
lgs_base_fb = 'output/lg_fb_new_'
i_ini_tot = 0
i_end_tot = 5

# Column coordinates
x_col = ['Xc_LG', 'Yc_LG', 'Zc_LG']
x_col_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
m_col = ['Mvir(4)']
d_col = 'Dist'
m_min = 4.0e+11
v_max = -0.0

# Loop over halo catalogs and lg lists
for i_cat in range(i_ini, i_end):
    i_str = '%02d' % i_cat

    if simu == 'fullbox':
        this_ahf = ahf_base + i_str + '/' + ahf_file
        this_lgs = lgs_base + i_str + '.csv'  
        this_ahf_csv = 'output/full.' + i_str + '.' + ahf_file_csv

    elif simu == 'lgf':
        this_lgs = lgs_base + '.csv'

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
        df_ahf = t.periodic_boundaries(data=df_ahf, slab_size=r_max, box=box_size)
        print('AHF FILE TO CSV: ', this_ahf_csv)
        df_ahf.to_csv(this_ahf_csv)

    # Read CSV ahf files and extract mass functions around LG candidates
    elif mode == 'read_csv':
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
            lg_csv = lg_csv_file + i_str + '_lg_' + str(i_lg) + '.csv'  
        
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

            if simu == 'fullbox':
                lg_csv = lg_csv_file + i_str + '_lg_' + str(i_lg) + '.csv'  
                this_lgs_rad = lgs_base + i_str + '_radii.csv'  
                
            elif simu == 'lgf':
                num = str('%02d' % int(row['simu_code']))
                sub = str('%02d' % int(row['sub_code']))
                lg_x = row[x_col].values
                lg_csv = lg_csv_file + num + '_' + sub + '.' + str(i_lg) + '.csv'
                this_lgs_rad = lgs_base + '_radii.csv'  

            if os.path.isfile(lg_csv):
                df_tmp = pd.read_csv(lg_csv)
                print('Reading: ', lg_csv)
                    
                for i_r, rad in enumerate(radii):
                    df_tmp = df_tmp[df_tmp[d_col] > 2000.0]
                    m_max = df_tmp[df_tmp[d_col] < rad][m_col].max()
                    m_rad[i_lg, i_r] = m_max

        for i_r, r_col in enumerate(r_str):
            df_lgs_orig[r_col] = m_rad[:, i_r]

        print('Saving to: ', this_lgs_rad)
        df_lgs_orig.to_csv(this_lgs_rad)

if mode == 'plots':

    label_fb = 'RAND'
    label_lg = 'CS'

    # Load FB data
    for i_cat in range(i_ini_tot, i_end_tot):
        i_str = '%02d' % i_cat
        this_lgs_rad = lgs_base_fb + i_str + '_radii.csv'  

        print('Reading from : ', this_lgs_rad)
        df_rad = pd.read_csv(this_lgs_rad)
        df_rad = df_rad[df_rad[r_str[0]] > 0.0]
        df_rm_all.append(df_rad)

    print('Merge all DFs...')
    df_fb = pd.concat(df_rm_all)
    print(df_fb.head())

    # Load CS data
    lgs_rad = lgs_base + '_radii.csv'  

    print('Reading from : ', lgs_rad, ' n_tot: ', len(lgs_rad))
    df_rad = pd.read_csv(lgs_rad)
    df_lg = df_rad[df_rad[r_str[0]] > 0.0]
    box_center = [5.0e+4] * 3
    df_lg['Dist'] = df_lg[x_col].T.apply(lambda x: t.distance(x, box_center))
    df_lg = df_lg[df_lg['Dist'] < d_max]
    print('New number: ', len(df_lg))

    slopes_fb = []; masses_fb = []; fracs_fb = []
    slopes_lg = []; masses_lg = []; fracs_lg = []
    meds_lg = np.zeros((len(r_str), 3))
    meds_fb = np.zeros((len(r_str), 3)) 

    for ir, rs in enumerate(r_str):
        f_out = 'output/masses_lg_cs_fb_' + rs + '.png'
        sns.distplot(np.log10(df_fb[rs]), color='orange', label=label_fb)
        sns.distplot(np.log10(df_lg[rs]), color='grey', label=label_lg)
        plt.xlabel(r'$\log_{10}M\quad [M_{\odot} h^{-1}]$')
        plt.legend()
        title = r'Max mass at R: ' + str(radii[ir]/1.e+3) + ' $h^{-1}$ Mpc'
        plt.title(title)
        plt.tight_layout()
        print('Saving file to: ', f_out)
        plt.savefig(f_out)
        plt.cla()
        plt.clf()

        percs = [10, 50, 90]

        med = np.percentile(np.log10(df_fb[rs]), percs)
        meds_fb[ir] = med
        masses_fb.append(med[1])
        medians_fb = '%.3f_{%.3f}^{+%.3f}' % (med[1], med[0]-med[1], med[2]-med[1]) 

        med = np.percentile(np.log10(df_lg[rs]), percs)
        meds_lg[ir] = med
        masses_lg.append(med[1])
        medians_lg = '%.3f_{%.3f}^{+%.3f}' % (med[1], med[0]-med[1], med[2]-med[1]) 
        print('LG:', medians_lg)
        print('FB:', medians_fb)
    
        m_cluster = 1.0e+14
        n_cluster = len(df_fb[df_fb[rs] > m_cluster])
        f_cluster = n_cluster / len(df_fb)
        print('Fractions of FB at ', rs, ' from a cluster: ',  f_cluster) 
        fracs_fb.append(f_cluster)

        n_cluster = len(df_lg[df_lg[rs] > m_cluster])
        f_cluster = n_cluster / len(df_lg)
        print('Fractions of LG at ', rs, ' from a cluster: ',  f_cluster) 
        fracs_lg.append(f_cluster)

        # FIXME: these fits are fishy and inconsistent, do not use them atm

        '''
        x_fb = np.log10(df_fb[rs]) 
        y_fb = np.log10(df_fb['Vtan']/100.0) 
        x_lg = np.log10(df_lg[rs])
        y_lg = np.log10(df_lg['Vtan']/100.0) 

        slope = np.polyfit(x_fb, y_fb, 1)
        slopes_fb.append(slope[0])

        f_out = 'output/masses_kde_FB_' + rs + '.png'
        plt.xlabel(r'$\log_{10} M_{max}$')
        plt.ylabel(r'$\log_{10} (V_{tan}/100)$')
        plt.tight_layout()
        sns.scatterplot(x_fb, y_fb)
        sns.kdeplot(x_fb, y_fb)
        plt.savefig(f_out)
        plt.cla()
        plt.clf()

        slope = np.polyfit(x_lg, y_lg, 1)
        slopes_lg.append(slope[0])

        f_out = 'output/masses_kde_LG_' + rs + '.png'
        plt.xlabel(r'$\log_{10} M_{max}$')
        plt.ylabel(r'$\log_{10} (V_{tan}/100)$')
        plt.tight_layout()
        sns.scatterplot(x_lg, y_lg)
        sns.kdeplot(x_lg, y_lg)
        plt.savefig(f_out)
        plt.cla()
        plt.clf()
        '''

    # Still in 'plots' mode
    print('FB: ')
    print(masses_fb)
    print(slopes_fb)
    print(fracs_fb)

    print('LG: ')
    print(masses_lg)
    print(slopes_lg)
    print(fracs_lg)

    '''
    plt.plot(masses, slopes, color='black')
    plt.xlabel(r'$\log_{10} M_{max}$')
    plt.ylabel(r'b')
    plt.tight_layout()
    plt.savefig('output/masses_vs_vtancorr.png')
    #plt.pause(3)
    plt.close()
    plt.cla()
    plt.clf()
    '''

    plt.plot(radii, masses_fb, color='red', label = label_fb)
    plt.fill_between(radii, meds_fb[:, 0], meds_fb[:, 2], color='orange', alpha=0.3)
    plt.plot(radii, masses_lg, color='black', label = label_lg)
    plt.fill_between(radii, meds_lg[:, 0], meds_lg[:, 2], color='grey', alpha=0.2)
    plt.xlabel(r'$R \quad [Mpc h^{-1}]$')
    plt.ylabel(r'$\log_{10}M \quad [M_{\odot} h^{-1}]$')
    plt.tight_layout()
    f_out = 'output/m_max_R_cs_vs_fb.png'
    plt.savefig(f_out) 

    delta_masses = (np.array(masses_fb) - np.array(masses_lg)) / np.array(masses_lg)
    delta_m = (np.array(masses_fb) - np.array(masses_lg)) 
    ratio_m = (np.array(masses_fb) / np.array(masses_lg)) 
    print(delta_masses)
    print(ratio_m)
    print(10 ** delta_m)
    #plt.savefig('output/frac_cl_vs_r.png')
    #plt.pause(3)
    #plt.close()
    #plt.cla()
    #plt.clf()




