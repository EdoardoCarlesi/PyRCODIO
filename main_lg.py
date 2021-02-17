'''
    Python Routines for COsmology and Data I/O (PyRCODIO) v0.2
    Edoardo Carlesi (2020)
    ecarlesi83@gmail.com

    main_lg: find LGs, extract and analyze their properties
'''

import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import seaborn as sns
import read_files as rf
import halo_utils as hu
import numpy as np
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os
from tqdm import tqdm
plt.rcParams.update({'font.size': 15})

# Set some important global variables
global simu, ahf_base, ahf_file, ahf_file_csv, lgs_base
global i_ini_tot, i_end_tot, m_min, v_max, masscut_pca
global x_col, x_col_ahf, m_col, dens_col, dist_col, mw_col, m31_col
global radii_cols, lg_models

# Define them here and use them through all the program
x_col = ['Xc_LG', 'Yc_LG', 'Zc_LG']
x_col_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
lg_models = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6']
dens_col = 'Dens'
dist_col = 'Distance'
m31_col = 'M_M31'
mw_col = 'M_MW'
m_col = 'Mvir(4)'
radii_cols = ['R', 'Dens', 'Mtot', 'Mmax', 'I_a', 'I_b', 'I_c', 'Iw_a',  'Iw_b', 'Iw_c', 'PCA_a', 'PCA_b', 'PCA_c', 'ID', 'simu']


def lg_density(simu='512', dist=5.0e+3):
    """ Compute the density of the LGs """

    if simu == '512':
        f_lg_csv = 'output/lg_pairs_512.csv'
        lgf=True

    elif simu == 'FB':
        f_lg_csv = 'output/lg_pairs_FB.csv'
        lgf=True

    elif simu == '1024':
        f_lg_csv = 'output/lg_pairs_1024.csv'
        lgf=False

    else:
        print(f'Simulation type {simu} does not exist.')
        exit()

    all_lg_models, all_lg_index = cfg.lg_models()
    data = pd.read_csv(f_lg_csv)
    simu_codes = data['simu_code'].values
    sub_codes = data['sub_code'].values
    all_codes = []

    for sim, sub in zip(simu_codes, sub_codes):
        str_sim = '%02d' % sim
        str_sub = '%02d' % sub
        all_codes.append(str_sim + '_' + str_sub)

    all_codes_set = set(all_codes)
    n_sims = len(all_codes_set)
    #print(len(all_codes), len(all_codes_set))

    for model in lg_models:
        this_index = all_lg_index[model]
        this_model = all_lg_models[this_index]

        these_lgs = hu.select_lgs(data=data, lg_model=this_model, lgf=lgf, dist=dist)
        n_lgs = len(these_lgs)

        print(f'{model} {n_lgs/n_sims} {n_lgs}')


def halos_around_lg(verbose=False, simu='fullbox', res='512', run_min=0, run_max=5):
    """
    Extract the mass functions around LG candidates in FB simulations, the most massive member at a given radius and also the density

    Output:

    - one list of dataframes with properties as a function of radius:
        RADIUS, MaxMass, MatterDensity, PCA params (1, 2, 3), Triaxiality (1, 2, 3)
    - one list of dataframes including the distance from the LG center:
        AHF HALO DATA, distance from the center

    """

    if simu == 'fullbox':
        print('halos_around_lg() running for fullbox simulation...')

        lg_file_base = '/home/edoardo/CLUES/PyRCODIO/output/LG_HALOS/halos_around_lg_'
        halos_suffix = '.csv'
        lg_suffix = '.pkl'

        runs = cfg.gen_runs(run_min, run_max)
        n_lg_max = 1000

        # Output file prefixes
        df_list_out = 'output/lg_fb_df.pkl'
        mf_list_out = 'output/lg_fb_mf.pkl'
        mf_x_list_out = 'output/lg_fb_mf_x' + res + '.pkl'

    elif simu == 'lgf':
        print(f'halos_around_lg() running for costrained simulation {res}')

        lg_file_base = '/home/edoardo/CLUES/PyRCODIO/output/LG_' +  res + '/halos_around_lg_'
        halos_suffix = '.csv'
        lg_suffix = '.pkl'

        runs = cfg.gen_all_runs(run_min, run_max, 0, 40)
        n_lg_max = 7

        # Output file prefixes
        df_list_out = 'output/lg_lgf_df' + res + '.pkl'
        mf_list_out = 'output/lg_lgf_mf' + res + '.pkl'
        mf_x_list_out = 'output/lg_lgf_mf_x' + res + '.pkl'

    # We will append all dataframes to this list
    df_list = []
    mf_list = []
    lg_list = []

    # This variable is used to keep track of all the indexes and rescale them
    all_n_tot = []

    # These are the radii at which we want to take different shells and check for the most massive halos therein
    r_min, r_max = 2, 13
    radii = [1000.0 * j for j in range(r_min, r_max)]

    this_lg_df = pd.DataFrame(columns = hu.LocalGroup(None, None).header(dump=False))
    print('Read in all files and removing duplicates...')

    # We first read in all LGs and remove duplicates
    for run in tqdm(runs):
        all_n_tot.append(len(lg_list))

        # Use tqdm to produce a nice progress bar
        for i in range(0, n_lg_max):

            lg_file_pkl = lg_file_base + run + '_' + str(i) + lg_suffix

            # Check that the file exists
            if os.path.isfile(lg_file_pkl):
                this_lg = pkl.load(open(lg_file_pkl,'rb'))
                this_lg.code_simu = run
                lg_list.append(this_lg)
                this_row = this_lg.info(dump=False)
                this_series = pd.Series(this_row, index = this_lg_df.columns)
                this_lg_df = this_lg_df.append(this_series, ignore_index=True)

    all_n_tot = np.sort(np.array(list(set(all_n_tot)), dtype=int))

    # Compute the raw number of LG candidates
    n_lg = len(lg_list)
    print(f'LG Total found: {n_lg}')

    # Clean the list
    this_lg_df = this_lg_df.drop_duplicates()
    n_lg_new = len(this_lg_df)
    print(f'LG Total found after cleaning: {n_lg_new}')

    # Sanity check
    if verbose:
        print(this_lg_df.columns)
        print(this_lg_df.index)
        print(this_lg_df.info())
        print(this_lg_df.head())

    # How to bin the mass function - settings
    n_bins = 25
    bin_min = 1.5e+9
    bin_max = 5.0e+13
    rad_mf = 5000.0

    # Other parameters
    box_size = 20.0e+3  # this is the smaller box, which is twixe the radius in L
    rad_max = radii[-1]

    # These mass bins will be used for all LGs
    x_bins = t.gen_bins(nbins=n_bins, binmax=bin_max, binmin=bin_min, binmode='log')
    pkl.dump(x_bins, open(mf_x_list_out, 'wb'))

    codes = this_lg_df['simu_code'].unique()

    # Loop on all the rows of the dataframe
    for i, row in tqdm(this_lg_df.iterrows()):

        run = row['simu_code']
        masses_r = []

        ind_run = np.where(codes == run)

        # We need to rescale the index of the dataframe
        i = i - all_n_tot[ind_run[0][0]]
        lg_list_csv = lg_file_base + run + '_' + str(i) + halos_suffix

        this_halos = pd.read_csv(lg_list_csv)

        if simu == 'lgf':
            dens0 = np.sum(this_halos[m_col].values) / (box_size ** 3.0)
        elif simu == 'fullbox':
            dens0 = 41.0

        this_x = np.reshape(row[x_col].values, (3, 1))

        # First, select only a subset of halos
        this_halos = t.select_sphere(data=this_halos, radius=rad_max, col=dist_col, x_col=x_col_ahf, center=this_x)
        this_df = pd.DataFrame(columns = radii_cols)

        # Now for each halo pair we look for the properties around it, at different radii
        for rad in radii[::-1]:
            new_row = []
            this_halos = this_halos[this_halos[dist_col] < rad]

            #print(dist_col, rad, len(this_halos))

            # Keep track of mass bins only at a given radius
            if rad == rad_mf:
                mass = this_halos[m_col].values
                y_bins = t.bin_df(data=this_halos, x_bins=x_bins, binmode='log')
                masses_r.append(y_bins)

            # Once the halos within the radius have been selected, compute some quantities
            vol = (np.pi * (rad ** 3.0)) * 4.0 / 3.0
            mmax = this_halos[m_col].max()
            mtot = np.sum(this_halos[m_col].values)
            dens = mtot / vol 
            fac = 1.0e+14
            #print(f'r={rad}, dens:{dens}, dens0:{dens0}, d:{dens/dens0}. mtot:{mtot/fac}, mmax:{mmax/fac}')

            dens = dens / dens0

            # Re-center all the halos around the LG center, to compute inertia tensors and PCA
            x = np.subtract(this_halos[x_col_ahf].T.values, this_x)
            m = this_halos[m_col].values

            # The eigenvalues here are not ordered!
            inertia_t = t.inertia_tensor(x=x)
            inertia_tw = t.inertia_tensor(x=x, w=m)
            pca_t = t.spatial_pca(data=this_halos, cols=x_col_ahf)

            # Add all the values to a list
            new_row.append(rad)
            new_row.append(dens)
            new_row.append(mtot)
            new_row.append(mmax)

            # Loop over the newly computed eigenvalues
            for tt in inertia_t:
                new_row.append(tt)
            for tw in inertia_tw:
                new_row.append(tw)
            for pt in pca_t:
                new_row.append(pt)

            # Keep track of the lg pair number and simulation
            new_row.append(i)
            new_row.append(run)

            # Make this a row of a pandas dataframe
            this_series = pd.Series(new_row, index = this_df.columns)
            this_df = this_df.append(this_series, ignore_index = True)

        if verbose:
            print(this_df.head(10))

        # Save all the lists
        df_list.append(this_df)
        mf_list.append(masses_r)

    # Now dump all the lists to appropriate pkl files
    print(f'Saving files to {df_list_out} and {mf_list_out}')
    pkl.dump(df_list, open(df_list_out, 'wb'))
    pkl.dump(mf_list, open(mf_list_out, 'wb'))
    print('Done')
    
    return df_list, mf_list


def find_lg_fb(run_min=0, run_max=5):
    """ Find all the local groups in a full box simulation """

    # Choose catalog type
    use_ahf = True; use_rs = False
    #use_ahf = False; use_rs = True

    # Simulation & catalog
    if use_ahf == True:
        file_single='snapshot_054.z0.000.AHF_halos'
        #file_single='snapshot_full_054.z0.000.AHF_halos'
        #base_file_out = 'output/lg_fb_new_'
        base_file_out = 'output/lg_fb_list_'
        base_file_halo_out = 'output/LG_HALOS/halos_around_lg_'
        box_size = 100000.0
        #base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'
        base_path = '/home/edoardo/CLUES/PyRCODIO/output/full.'
        base_suffix = '.snapshot_054.z0.000.AHF_halos.periodic.csv'
        #base_suffix = '.snapshot_054.z0.000.AHF_halos.csv'
        #base_suffix = '.snapshot_054.z0.000.AHF_halos.small.csv'
        #base_path='/media/edoardo/Elements/CLUES/DATA/2048/00_06/'
        sub_runs = cfg.gen_runs(run_min, run_max)
        print(sub_runs)

    elif use_rs == True:
        box_size = 2500000.0
        file_single = '.part'
        base_file_out = 'output/lg_fullbox_rs_'
        base_path = '/home/edoardo/CLUES/DATA/RS/out_79_csv/'
        #base_path = '/z/carlesi/STORE/MultiDark/RockStarCSV/BigMD_3840_Planck1/out_79_csv/'
        sub_runs = []

        n_start, n_parts = 10, 20
        for i in range(n_start, n_parts):
            sub = '%04d' % i
            sub_runs.append(sub)

    lg_models, index = cfg.lg_models()
    this_model = lg_models[index['GENERIC']]
    
    # Here we set the parameters to split the box into sub-volumes
    kpcFac = 1.0e+3
    radius = 13.0 * kpcFac
    side_buffer = 1.0 * kpcFac

    n_sub_x = int(np.ceil(box_size / radius))
    n_sub_y = int(n_sub_x)
    n_sub_z = int(n_sub_y)
    n_tot = np.power(n_sub_x, 3)

    print('Subdivision in ', n_sub_x, ' subcubes per axis, radius: ', radius, ' and side_buffer: ', side_buffer)

    # Loop on all the different realizations
    for run in sub_runs:
        all_lgs = []
        print(run, sub_runs)

        if use_ahf:
            #this_ahf = base_path + run + '/' + file_single
            this_ahf = base_path + run + base_suffix
            print('Reading file: ', this_ahf)
            #halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
            #halo_df = rf.read_ahf_halo(this_ahf, file_mpi=True)
            halo_df = pd.read_csv(this_ahf)
            n_halo = len(halo_df)
            print('Found: ', n_halo, ' objects.')

            x_min, y_min, z_min = 0.0, 0.0, 0.0;         

        elif use_rs:
            this_rs = base_path + run + file_single

            print('Reading file: ', this_rs)
            halo_df = pd.read_csv(this_rs)
            halo_df.columns = t.header_rs2ahf(halo_df.columns)

            # The file is divided in sub files, determine the extension
            x_max = halo_df['Xc(6)'].max() * kpcFac
            y_max = halo_df['Yc(7)'].max() * kpcFac
            z_max = halo_df['Zc(8)'].max() * kpcFac
            x_min = halo_df['Xc(6)'].min() * kpcFac
            y_min = halo_df['Yc(7)'].min() * kpcFac
            z_min = halo_df['Zc(8)'].min() * kpcFac

            n_sub_x = int(np.ceil((x_max - x_min) / radius))
            n_sub_y = int(np.ceil((y_max - y_min) / radius))
            n_sub_z = int(np.ceil((z_max - z_min) / radius))
            n_tot = n_sub_x * n_sub_y * n_sub_z

            print('N subcubes ', n_sub_x, n_sub_y, n_sub_z, ' Ntot: ', n_tot)

        n_count = 0
        old_time = time.time()
        
        for ix in tqdm(range(0, int(n_sub_x))):
            for iy in range(0, n_sub_y):
                new_time = time.time()
                dif_time = '%.3f' % (new_time - old_time)
                percent = '%.3f' % (100.0 * n_count/n_tot)
                #print('Done: ', percent, '% in ', dif_time, ' seconds. Tot LGs: ', len(all_lgs), flush = True)

                for iz in range(0, n_sub_z):
                    n_count += 1
                    this_center = np.array([radius * (0.5 + ix)+x_min, radius * (0.5 + iy)+y_min, radius * (0.5 + iz)+z_min])
                    this_radius = radius * 0.5 + side_buffer
                    these_lgs = hu.find_lg(halo_df, this_center, this_radius, lgmod=this_model, center_cut=True, search='Box', verbose=False)

                    for this_lg in these_lgs:
                        this_lg.code_simu = 'FB'
                        this_lg.code_sub = run

                        all_lgs.append(this_lg)

        this_lg_df = pd.DataFrame(columns = this_lg.header(dump=False))

        for i, lg in enumerate(all_lgs):
            this_row = lg.info(dump=False)
            this_series = pd.Series(this_row, index = this_lg_df.columns)
            this_lg_df = this_lg_df.append(this_series, ignore_index=True)
            this_halos_csv = base_file_halo_out + run + '_' + str(i) + '.csv'
            this_halos_pkl = base_file_halo_out + run + '_' + str(i) + '.pkl'
            lg_center = lg.geo_com()
            this_halos = hu.find_halos(data=halo_df, center=lg_center, radius=radius, search='Box')
            pkl.dump(lg, open(this_halos_pkl, 'wb'))

            # Discard some useless columns to save up space on the disk
            this_halos = this_halos.drop(this_halos.columns[18:], axis=1)
            this_halos = this_halos.drop(this_halos.columns[0:2], axis=1)
            this_halos.to_csv(this_halos_csv)

        this_csv = base_file_out + run + '.csv'
        this_lg_df.drop_duplicates(inplace = True)
        this_lg_df.to_csv(this_csv, float_format='%.3f')

        return this_halos


def find_lg_lgf(res='512', run_min=0, run_max=80):
    """ Given a set of catalogs find the LG-like objects and export the output """

    # Use AHF / csv catalogs
    csvAhf = True; hiResAhf = False
    #csvAhf = False; hiResAhf = True

    # Configure the LG model and subpaths
    if csvAhf == True:
        code_runs = cfg.gen_all_runs(run_min, run_max, 0, 40)
        [model_run, dict_model] = cfg.lg_models()

    elif hiResAhf == True:
        [model_run, dict_model] = cfg.lg_models()
        code_runs = cfg.gen_all_runs(0, 100, 0, 40)

    # Local data path, file names and file format
    data_path = '/home/edoardo/CLUES/PyRCODIO/data/'
    file_ahf = 'snapshot_054.0000.z0.000.AHF_halos'

    # Full dataset base path
    if csvAhf == True:
        #base_path = '/home/edoardo/Elements/CLUES/DATA/CSV/' + res + '/'
        base_path = '/home/edoardo/CLUES/DATA/LGF/' + res + '/CSV/'
        base_path_out = 'output/LG_' + res + '/halos_around_lg_'

    elif hiResAhf == True:
        base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'

    # Select a subsample from the full catalog to look for local groups
    facKpc = 1.0e+3
    radius = 10.0 * facKpc
    center = [50.0 * facKpc] * 3

    if res == '512':
        n_cat_min = 20000
    else:
        n_cat_min = 200

    # Read the | Gyr / z / a | time conversion table
    time = rf.read_time(data_path)

    all_halo_mah = []

    # Output files
    if csvAhf == True:
        out_base_pkl = base_path + 'lg_'
        out_all_lgs_csv = 'output/lg_pairs_' + res + '.csv'

    elif hiResAhf == True:
        out_base_pkl = 'saved/lg_pair_'
        out_all_lgs_csv = 'output/lg_pairs_2048.csv'

    # Use a dummy halo to properly generate the file header
    h = hu.Halo('void')
    cols = hu.LocalGroup(h, h).header(dump=False)
    all_lgs_df = pd.DataFrame(columns=cols)

    # Now loop on all the simulations and gather data
    for code in code_runs:

        # Append all found LGs to this structure
        these_lgs = []

        if csvAhf == True:
            this_ahf = base_path + 'ahf_' + code + '.csv'

        elif hiResAhf == True:
            this_path = base_path + code + '/' + sub + '/'
            this_ahf = this_path + file_ahf

        # Check that file exists
        if os.path.isfile(this_ahf):
            print('Reading AHF file: ', this_ahf)

            if csvAhf == True:
                halo_df = pd.read_csv(this_ahf)
                this_model = model_run[dict_model['GENERIC']]

            elif hiResAhf == True:
                halo_df = rf.read_ahf_halo(this_ahf)
                this_model = model_run[dict_model[code]]

            tot_mass = np.sum(halo_df[m_col].values)
            max_mass = halo_df[m_col].max()
            vol = (100.0e+3) ** 3.0
            #print('Mass Density: ', tot_mass/vol, ' max: ', max_mass/1.0e+14)
            #print(halo_df.head())
            #print(halo_df[m_col])
            #print(halo_df.columns)
            print('Looking for Local Group candidates...')

            if len(halo_df) > n_cat_min:
                these_lgs = hu.find_lg(halo_df, center, radius, lgmod=this_model, center_cut=True, search='Box')
            else:
                #print('Problems with file: ', this_ahf, '. Reading the file results in zero length.')
                print('Skipping catalog: ', this_ahf) 

            # Check if there is at least one LG in the selection
            if len(these_lgs) > 0:

                # Save simu info in LG object
                for i in range(0, len(these_lgs)):
                    these_lgs[i].code_simu = code
                    this_series = pd.Series(these_lgs[i].info(dump=False), index=all_lgs_df.columns)
                    all_lgs_df = all_lgs_df.append(this_series, ignore_index=True)

                    these_halos = hu.find_halos(data=halo_df, center=these_lgs[i].geo_com(), radius=radius, search='Box')
                    #these_halos = these_halos.drop(these_halos.columns[18:], axis=1)
                    #print(f'R={radius}, maxMass={these_halos[m_col].max()/1.0e+14}')
                    out_file_pkl = base_path_out + code + '_' + str(i) + '.pkl'
                    out_file_csv = base_path_out + code + '_' + str(i) + '.csv'

                    print(f'Saving LG output to files {out_file_pkl} and {out_file_csv}')
                    pkl.dump(these_lgs[i], open(out_file_pkl, 'wb'))
                    these_halos.to_csv(out_file_csv)

    # Dump lg list to csv
    print('Dumping to csv file: ', out_all_lgs_csv)
    all_lgs_df.to_csv(out_all_lgs_csv)

    return these_lgs


def plot_halos_around_lg(res='512'):
    """ This function returns the properties of triaxiality/density around the LG for both simulations """

    out_lg_csv = 'output/lg_pairs_' + res + '.csv'
    out_fb_csv = 'output/lg_pairs_FB.csv'
    out_fb = 'output/lg_fb_df.pkl'
    out_lg = 'output/lg_lgf_df' + res + '.pkl'
    print(f'Plotting main properties of the LG environment from files {out_fb}, {out_lg}')
    
    #radii_cols = ['R', 'Dens', 'Mtot', 'Mmax', 'I_a', 'I_b', 'I_c', 'Iw_a',  'Iw_b', 'Iw_c', 'PCA_a', 'PCA_b', 'PCA_c', 'ID', 'simu']
    data_fb = pkl.load(open(out_fb, 'rb'))
    data_lg = pkl.load(open(out_lg, 'rb'))

    df_list_fb = pd.read_csv(out_fb_csv)
    df_list_lg = pd.read_csv(out_lg_csv)

    # Column names, group them
    col_mmax = 'Mmax'
    col_mtot = 'Mtot'
    col_I = ['I_a', 'I_b', 'I_c']
    col_Iw = ['Iw_a', 'Iw_b', 'Iw_c']
    col_pca = ['PCA_a', 'PCA_b', 'PCA_c']

    # Do we want to read all the data and sort it or load from .pkl
    sort_data = False

    # Initialize some variables
    radii = data_lg[0]['R'].values

    file_sorted_fb = 'output/sorted_fb.pkl'
    file_sorted_lg = 'output/sorted_lg' + res + '.pkl'
    center = [50000.0] * 3


    def sort_data(data=None, mvirgo=1.0e+14, select_lg=True):
        """ This procedure will be repeated for LGF and FB halos """

        n_rows = len(data[0])
        good_lgs = 0 

        # Initialize the lists that will contain all the data
        dens = [[] for i in range(0, n_rows)]
        mmax = [[] for i in range(0, n_rows)]
        mtot = [[] for i in range(0, n_rows)]
        i_t = [[] for i in range(0, n_rows)]
        iw_t = [[] for i in range(0, n_rows)]
        pca_t = [[] for i in range(0, n_rows)]
        triax = [[] for i in range(0, n_rows)]
        virgo = np.zeros(n_rows)

        # Loop over the list of dataframes, collect all data
        for ind, df in enumerate(tqdm(data[0:])):

            dist = t.distance(df_list_lg[x_col].iloc[ind], center) 
            if select_lg:
                if dist < 5000.0:
                    check_model = True
                    good_lgs += 1
                else:
                    check_model = False

            else:
                check_model = True

            if check_model:
                for i, row in df.iterrows():

                    # Gather mass and maximum mass at a given radius
                    dens[i].append(row[dens_col])
                    mmax[i].append(row[col_mmax])
                    mtot[i].append(row[col_mtot])

                    if float(mmax[i][len(mmax[i])-1]) > mvirgo:
                        virgo[i] += 1.0

                    # Gather the inertia tensors and triaxialities
                    i_t_ord = np.sort(row[col_I].values)[::-1]
                    iw_t_ord = np.sort(row[col_Iw].values)[::-1]
                    pca_t_ord = np.sort(row[col_pca].values)[::-1]

                    # Compute the triaxialities using three eigenvalues
                    triax_i = t.triaxiality(*i_t_ord)
                    triax_iw = t.triaxiality(*iw_t_ord)
                    triax_pca = t.triaxiality(*pca_t_ord)

                    # Append it all
                    i_t[i].append(list(i_t_ord))
                    iw_t[i].append(iw_t_ord)
                    pca_t[i].append(pca_t_ord)
                    triax[i].append([triax_i, triax_iw, triax_pca])

        # Set some useful variables
        #percentiles = [0, 20, 50, 80, 100]
        percentiles = [25, 50, 75]
        #rho0 = 41.0
        rho0 = 1.0
        #print(good_lgs, 100)
        n_perc = len(percentiles)

        # These contain median, lower percentile and higher percentile
        dens_perc = np.zeros((n_perc, n_rows), dtype=float)
        mmax_perc = np.zeros((n_perc, n_rows), dtype=float)

        # Eigenvalue sets
        i_t_perc = np.zeros((3, n_perc, n_rows), dtype=float)
        iw_t_perc = np.zeros((3, n_perc, n_rows), dtype=float)
        pca_t_perc = np.zeros((3, n_perc, n_rows), dtype=float)

        # Triax is three triax values with three percentile intervals
        triax_perc = np.zeros((n_perc, 3, n_rows), dtype=float)

        # Make this stuff some humanly-manageable matrix
        i_t = np.array(i_t).T
        iw_t = np.array(iw_t).T
        pca_t = np.array(pca_t).T
        triax = np.array(triax).T

        # Loop on each row (i.e. radius)
        for i in range(0, n_rows):

            dens_perc[:, i] = np.percentile(dens[i], q=percentiles) / rho0
            mmax_perc[:, i] = np.percentile(mmax[i], q=percentiles)

            # Loop on a, b, c axes and three triaxialities
            for j in range(0, 3):
                i_t_perc[j, :, i] = np.percentile(i_t[j, :, i], q=percentiles)
                iw_t_perc[j, :, i] = np.percentile(iw_t[j, :, i], q=percentiles)
                pca_t_perc[j, :, i] = np.percentile(pca_t[j, :, i], q=percentiles)
                triax_perc[j, :, i] = np.percentile(triax[j, :, i], q=percentiles)

        return dens_perc, mmax_perc, i_t_perc, iw_t_perc, pca_t_perc, triax_perc, virgo


    def plot_f_r(x=radii, y0=None, err0=False, y1=None, err1=False, y_label=None, fout=None, log_r=False):
        """ Plot a radius-dependent quantity, with or without error """

        # General plot values
        color0 = 'black'
        color1 = 'blue'

        # Line labels
        line0 = 'LGF-L'
        line1 = 'RAND'

        n_xmax = 2

        if log_r:
            #x = np.log10(np.array(x))
            #plt.xlabel(r'$\log_{10} R $')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$R \quad [h^{-1} Mpc]$')
        else:
            x = np.array(x) / 1000.0
            plt.xlabel(r'$R \quad [h^{-1} Mpc]$')

        # First plot densities
        plt.ylabel(y_label)
    
        if err0 and err1:
            plt.plot(x[n_xmax:], y0[1, n_xmax:], color=color0, label=line0)
            plt.plot(x[n_xmax:], y1[1, n_xmax:], color=color1, label=line1)
            plt.fill_between(x[n_xmax:], y0[0, n_xmax:], y0[2, n_xmax:], color=color0, alpha=0.3)
            plt.fill_between(x[n_xmax:], y1[0, n_xmax:], y1[2, n_xmax:], color=color1, alpha=0.3)

        else:
            plt.plot(x, y0, color=color0, label=line0)
            plt.plot(x, y1, color=color1, label=line1)

        plt.legend()
        plt.tight_layout()
        plt.savefig(fout)
        plt.close()
        plt.cla()
        plt.clf()

        return None


    # If we want to load and sort the data from scratch
    if sort_data:
        print('Sorting data before plotting...')

        # Get the values using the above defined function
        d_lg, m_lg, i_lg, iw_lg, pca_lg, t_lg, virgo_lg = sort_data(data=data_lg, select_lg=True)
        d_fb, m_fb, i_fb, iw_fb, pca_fb, t_fb, virgo_fb = sort_data(data=data_fb, select_lg=False)

        sorted_lg = (d_lg, m_lg, i_lg, iw_lg, pca_lg, t_lg, virgo_lg)
        pkl.dump(sorted_lg, open(file_sorted_lg, 'wb'))

        sorted_fb = (d_fb, m_fb, i_fb, iw_fb, pca_fb, t_fb, virgo_fb)
        pkl.dump(sorted_fb, open(file_sorted_fb, 'wb'))

    # Otherwise read the plot-ready files from the output
    else:
        print('Loading sorted data for plotting: {file_sorted_lg} {file_sorted_fb}')

        # Load the pkl sorted values
        sorted_lg = pkl.load(open(file_sorted_lg, 'rb'))
        d_lg, m_lg, i_lg, iw_lg, pca_lg, t_lg, virgo_lg = sorted_lg

        sorted_fb = pkl.load(open(file_sorted_fb, 'rb'))
        d_fb, m_fb, i_fb, iw_fb, pca_fb, t_fb, virgo_fb = sorted_fb

    # Median density as a function of radius
    print('Plotting density...')
    f_out = 'output/lg_fb_dens.png'
    y_label = r'$\Delta$'
    #y_label = r'$\log_{10}\Delta$'
    #d_lg = np.log10(d_lg)
    #d_fb = np.log10(d_fb)
    d_lg = d_lg
    d_fb = d_fb
    plot_f_r(y0=d_lg, y1=d_fb, y_label=y_label, fout=f_out, err0=True, err1=True, log_r=True)

    # Max mass as a function of radius
    print('Plotting maximum mass...')
    f_out = 'output/lg_fb_mmax.png'
    y_label = r'$M_{max}$'
    plot_f_r(y0=np.log10(m_lg), y1=np.log10(m_fb), y_label=y_label, fout=f_out, err0=True, err1=True)

    # Triaxialities
    print('Plotting triaxialities...')
    f_out_base = 'output/lg_fb_triax_'
    y_labels = ['T_I', 'T_Iw', 'T_PCA']

    for i, y_label in enumerate(y_labels):
        f_out = f_out_base + y_label + '.png'
        plot_f_r(y0=t_lg[i], y1=t_fb[i], y_label=y_label, fout=f_out, err0=True, err1=True)

    # PCA evs
    f_out_pca = 'output/lg_fb_pca_'

    # Inertia tensor
    print('Plotting inertia tensor...')
    f_out_base = 'output/lg_fb_inertia_'
    y_labels = ['T_I_a', 'T_I_b', 'T_I_c']

    for i, y_label in enumerate(y_labels):
        f_out = f_out_base + y_label + '.png'
        plot_f_r(y0=i_lg[i], y1=i_fb[i], y_label=y_label, fout=f_out, err0=True, err1=True)


    # Weighter inertia tensor
    print('Plotting weighted inertia tensor...')
    f_out_base = 'output/lg_fb_inertia_'
    y_labels = ['T_Iw_a', 'T_Iw_b', 'T_Iw_c']

    for i, y_label in enumerate(y_labels):
        f_out = f_out_base + y_label + '.png'
        plot_f_r(y0=iw_lg[i], y1=iw_fb[i], y_label=y_label, fout=f_out, err0=True, err1=True)

    # Plot Virgo properties
    #print(virgo_fb)
    #print(virgo_lg)

    print('Done.')

    return None


def plot_lg_densities():
    """ Number of object per Mpc cube """

    data = pd.read_csv('lg_halo_dens.txt')
    
    vol_lg = np.pi * 4.0 / 3.0 * 125.0 * 1000.0
    vol_fb = 100.0 ** 3.0 * 5
    delta = 0.44

    models = [1, 2, 3, 4, 5, 6]
    models_labels = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6']


    def plot_lines_ratio(x=None, y=None, z=None, plottype=None, setrange=False):
        """ """

        ratio1 = np.ones(6)
        ratio2 = x / y
        ratio3 = x / z

        color0 = 'black'
        color1 = 'blue'
        color2 = 'green'
    
        size = 8
        (fig, axs) = plt.subplots(ncols=1, nrows=2, figsize=(size, size), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
        if setrange:
            axs[0].set_ylim([1.e-3, 1.5e-2])

        axs[0].set_yscale('log')
        axs[0].set_ylabel('$N \quad h^{3} Mpc^{-3}$')
        axs[1].set_ylabel('ratio')
        axs[1].set_xlabel('Model')
        axs[0].plot(models, x, label='LGF-L', color=color0)
        axs[0].plot(models, y, label='RAND', color=color1)
        axs[0].plot(models, z, label=r'RAND$_{\Delta}$', color=color2)
        axs[1].plot(models, ratio1, color=color0)
        axs[1].plot(models, ratio2, color=color1)
        axs[1].plot(models, ratio3, color=color2)
        
        fout = 'output/dens_' + plottype + '.png'
        fig.legend(loc = 'center')
        fig.tight_layout()
        fig.savefig(fout)
        plt.close()
        plt.cla()
        plt.clf()
   
        return None


    fac_lgf1 = 1.0
    y_lg = data['lgfI'].values * fac_lgf1 / vol_lg
    y_fb = data['fb'].values / vol_fb
    y_delta = data['fb_delta'].values / (vol_fb * delta)
    plot_lines_ratio(x=y_lg, y=y_fb, z=y_delta, plottype='lgs')

    fac_lgf1 = 0.4
    fac_lgf2 = 4.0
    delta = 0.6
    y_lg = data['halos_lgf'].values * fac_lgf1 / vol_lg
    y_fb = data['halos_fb'].values * fac_lgf2 / vol_fb
    y_delta = data['halos_fb_delta'].values * fac_lgf2 / (vol_fb * delta)
    plot_lines_ratio(x=y_lg, y=y_fb, z=y_delta, plottype='halos', setrange=True)

    return None


# FIXME TODO
def plot_mf_around_lg(res='1024', do_bins=True):
    """ Read and compare the mass function bins """

    data_fb = pkl.load(open('output/lg_fb_mf.pkl', 'rb'))
    data_lg = pkl.load(open('output/lg_lgf_mf' + res + '.pkl', 'rb'))

    print('Loading mass functions')

    # Check wether we re-bin the x from scratch or we have to do it 
    if do_bins == True:

        # How to bin the mass function - settings
        n_bins = 25
        bin_min = 1.5e+9
        bin_max = 5.0e+13

        # These mass bins will be used for all LGs
        x_bins = t.gen_bins(nbins=n_bins, binmax=bin_max, binmin=bin_min, binmode='log')
    else:
        # TODO fix filename
        x_bins = pkl.load(open('output/lg_fb_mf_x.pkl', 'rb'))

    # This one contains the actual x axis points
    x = np.zeros(n_bins-1)
    
    # Take the median point per mass bin
    for i in range(0, n_bins-1):
        x[i] = 0.5 * (x_bins[i+1] + x_bins[i])
    

    def sort_mf(data=None, vol=None):
        """ Get all the mass functions and look for the percentiles at each mass bin """

        if vol == 'Box':
            #vol = 10 **3.0
            vol = (5.5 **3.0) * np.pi * 4.0 / 3.0
        elif vol == 'Sphere':
            vol = 125.0 * np.pi * 4.0 / 3.0

        percentiles = [25, 50, 75]
        n_bins = data[0][0].size
        n_data = len(data)

        y_full = np.zeros((n_bins, n_data))
        y = np.zeros((3, n_bins))

        for i, mf in enumerate(data):
            for j in range(0, n_bins):
                #if mf[0][j] > 1:
                y_full[j, i] = np.log10(mf[0][j] / vol)
                #else:
                #    y_full[j, i] = 1.e-10

        for i in range(0, n_bins):
            y[:, i] = np.percentile(y_full[i, :], q=percentiles)

        return y

    # 
    y0 = sort_mf(data=data_lg, vol='Sphere')
    y1 = sort_mf(data=data_fb, vol='Box')

    # General plot values
    color0 = 'black'
    color1 = 'blue'

    # Line labels
    line0 = 'LGF-H'
    line1 = 'RAND'

    n_min = 1
    n_max = 18

    # Plot mass functions
    plt.xlabel(r'$M \quad [h^{-1} M_{\odot}]$')
    plt.ylabel(r'$ n h^3 Mpc^{-3}$')

    plt.plot(x[n_min:n_max], y0[1, n_min:n_max], color=color0, label=line0)
    plt.plot(x[n_min:n_max], y1[1, n_min:n_max], color=color1, label=line1)

    plt.fill_between(x[n_min:n_max], y0[0, n_min:n_max], y0[2, n_min:n_max], color=color0, alpha=0.3)
    plt.fill_between(x[n_min:n_max], y1[0, n_min:n_max], y1[2, n_min:n_max], color=color1, alpha=0.3)

    fout = 'output/lg_fb_mf_5mpc.png'

    plt.legend()
    plt.tight_layout()
    plt.savefig(fout)
    plt.close()
    plt.cla()
    plt.clf()

    print('Done.')

    return None


# Wrapper for the above functions 
if __name__ == "__main__":

    """
    How to run the script:
    1) Find the lgs
    2) Check the properties of halos around them
    3) Do the plots
    """

    '''
    rads = [5.0e+3, 7.0e+3, 8.0e+3, 9.0e+3, 10.0e+3]

    for r in rads:
        print(f'Rad={r}')
        lg_density(dist=r)
    '''

    #find_lg_fb(run_max=1)

    #def halos_around_lg(verbose=False, simu='fullbox', res='512', run_min=0, run_max=5):
    #find_lg_lgf(res='1024', run_max=80)
    #halos_around_lg(simu='lgf', res='1024', run_max=80)
    #find_lg_fb(run_min=3, run_max=5)
    #halos_around_lg(simu='fullbox', run_min=0, run_max=5)
    plot_halos_around_lg(res='1024')
    #plot_mf_around_lg()
    #plot_lg_densities()

    #halos_around_lg(run_max=1)
    #print(f'FB: {vol_fb}, LGF: {vol_lgf_I}')
    #for i in range(0, 6):
    #print(f'{i} {n_fb[i]} {n_lgf_I[i]/vol_lgf_I}')
    #def frac_from_virgo(m_virgo=0.7e+14):
    #frac_from_virgo()
    #fb_halos_around_lg()
    #halo_density_lgf()
    #halo_density_fb()
