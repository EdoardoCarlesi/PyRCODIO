'''
    Python Routines for COsmology and Data I/O (PyRCODIO) v0.2
    Edoardo Carlesi (2020)
    ecarlesi83@gmail.com

    main_lg: find LGs, extract and analyze their properties 
'''

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_mass_max():
    """ Plot the results of the analysis """

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


def plots_pca():
    '''
    Once the PCA has been applied to the halos' coordinates to look for flatness of the distribution, do a plot
    '''
    # TODO this does the full pca thing now
    '''
            MOVE THE PCA STUFF SOMEWHERE ELSE
        '''

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


def density_plots():
    '''
    Do a series of plots centered around LG-like pairs
    '''

    # Choose catalog type
    use_ahf = True; use_rs = False
    #use_ahf = False; use_rs = True

    # Simulation & catalog
    if use_ahf == True:
        file_single='snapshot_054.z0.000.AHF_halos'
        #file_single='snapshot_full_054.z0.000.AHF_halos'
        base_file_out = 'output/lg_fb_new_'    
        box_size = 100000.0
        #base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'
        base_path = '/home/edoardo/CLUES/DATA/FullBox/'
        #base_path='/media/edoardo/Elements/CLUES/DATA/2048/00_06/'
        sub_runs = cfg.gen_runs(0, 5)

    elif use_rs == True:
        box_size = 2500000.0
        file_single = '.part'
        base_file_out = 'output/lg_fullbox_rs_'    
        base_path = '/home/edoardo/CLUES/DATA/RS/out_79_csv/'
        #base_path = '/z/carlesi/STORE/MultiDark/RockStarCSV/BigMD_3840_Planck1/out_79_csv/'
        sub_runs = []

        n_start = 10
        n_parts = 20
        for i in range(n_start, n_parts):
            sub = '%04d' % i
            sub_runs.append(sub)
        '''
        sub_runs = ['0011']
        '''

    lg_models, index = cfg.lg_models()
    this_model = lg_models[index['GENERIC']]

    kpcFac = 1.0e+3
    radius = 7.0 * kpcFac
    side_buffer = 1.0 * kpcFac

    n_sub_x = int(np.ceil(box_size / radius))
    n_sub_y = int(n_sub_x)
    n_sub_z = int(n_sub_y)
    n_tot = np.power(n_sub_x, 3)

    print('Subdivision in ', n_sub_x, ' subcubes per axis, radius: ', radius, ' and side_buffer: ', side_buffer)

    for run in sub_runs:

        all_lgs = []
        if use_ahf:
            this_ahf = base_path + run + '/' + file_single
            print('Reading file: ', this_ahf)
            #halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
            halo_df = rf.read_ahf_halo(this_ahf, file_mpi=True)
            print('Found: ', len(halo_df), ' objects.')

            x_min = 0.0;         y_min = 0.0;         z_min = 0.0

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
            n_tot = n_sub_x * n_sub_y * n_sub_z

            print('N subcubes ', n_sub_x, n_sub_y, n_sub_z, ' Ntot: ', n_tot)

        #print(halo_df.head())
        n_count = 0
        old_time = time.time()
        for ix in range(0, int(n_sub_x)):
            for iy in range(0, n_sub_y):
                new_time = time.time()
                dif_time = '%.3f' % (new_time - old_time)
                percent = '%.3f' % (100.0 * n_count/n_tot) 
                print('Done: ', percent, '% in ', dif_time, ' seconds. Tot LGs: ', len(all_lgs), flush = True)
                #old_time = time.time()

                for iz in range(0, n_sub_z):
                    n_count += 1
                    this_center = np.array([radius * (0.5 + ix)+x_min, radius * (0.5 + iy)+y_min, radius * (0.5 + iz)+z_min])
                    this_radius = radius * 0.5 + side_buffer
                    #print('Subbox around center: ', this_center, ' rad: ', this_radius, flush=True)
                    these_lgs = hu.find_lg(halo_df, this_model, this_center, this_radius, center_cut=True, search='Box', verbose=False)
                    #print('FindLG: ', len(these_lgs))

                    for this_lg in these_lgs:
                        this_lg.code_simu = 'FB'
                        this_lg.code_sub = run

                        all_lgs.append(this_lg)

            #for ii in range(0, 10):
            #    print(ii, all_lgs[len(all_lgs) - ii - 1].info())

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


def plot_mass_functions():
    """ Take the binned mass functions and do some plots """ 

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


def plot_mass_max():
    """ Plot the results of the analysis """

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


