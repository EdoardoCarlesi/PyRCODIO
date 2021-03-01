"""
    Python Routines for COsmology and Data I/O (PyRCODIO)
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com
    
    main_extract_vweb.py: extract v-web data around a given point
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
from tqdm import tqdm
import os
plt.rcParams.update({'font.size': 15})


def extract_vweb_cs():
    """ Extract the eigenvalues at the box center in constrained simulations """

    # configure subpaths
    code_run = cfg.gen_runs(0, 1)
    sub_run = cfg.gen_runs(0, 30)

    # local data path, file names and file format
    base_vweb = 'vweb_'
    format_vweb = '.000128.vweb-csv'
    #format_vweb = '.000032.vweb-csv'

    # full dataset
    base_path = '/z/carlesi/clues/data/512/vweb/'

    #kpcfac = 1.0
    kpcfac = 1.e+3

    # select a subsample of nodes of the web around a given point
    radius = 10.0 * kpcfac
    center = [50.0 * kpcfac] * 3

    # output file base path
    out_base = base_path + 'vweb_center_0128_'

    # now loop on all the simulations and gather data
    for code in code_run:

        for sub in sub_run:
            this_vweb = base_path + base_vweb + code + '_' + sub + format_vweb
            out_file = out_base + code + '_' + sub + '.pkl'
            out_file_csv = out_base + code + '_' + sub + '.csv'

            # check that file exists
            if os.path.isfile(this_vweb):
                select_vweb = rf.extract_vweb(file_name=this_vweb, center=center, radius=radius)
                select_vweb.to_csv(out_file_csv)
                select_vweb.to_pickle(out_file)
            else: 
                print(this_vweb, ' not found')

    
def extract_vweb(data=None, grid=128, simu='fullbox', load_file=False):
    """ Extract the Vweb at given positions in random full box simulations """

    # Columns to be extracted from the vweb and lg dataframes
    l_cols = ['l1', 'l2', 'l3', 'dens']
    g_cols = ['x', 'y', 'z']
    df_cols = ['ID', 'x', 'y', 'z', 'l1', 'l2', 'l3', 'dens']

    # Read in the cosmic web
    if simu == 'fullbox':
        x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
        grid_str = '.%06d' % grid
        grid_str = str(grid_str)
        format_vweb = '.Vweb-csv'
        #base_path = '/media/edoardo/Elements/CLUES/DATA/FullBox/VWeb/vweb_'
        base_path = '/home/edoardo/CLUES/DATA/Vweb/FullBox/vweb_'
        out_base = 'output/vweb_fullbox_vweb.' + grid_str + '.pkl'
        code_col = 'sub_code'
        sub_run = data[code_col].unique()

    elif simu == 'lgf':
        x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
        grid_str = '%04d' % grid
        grid_str = str(grid_str)
        format_vweb = '.csv'
        base_path = '/home/edoardo/CLUES/DATA/Vweb/512/CSV/vweb_center_' + grid_str + '_'
        #base_path = '/home/edoardo/Elements/CLUES/DATA/VWeb/vweb_center_' + grid_str + '_'
        out_base = 'output/vweb_lgf512.' + grid_str + '.pkl'
        code_col = 'simu_code'
        sub_run = data[code_col].unique()

    elif simu == 'rand_mw':
        x_cols = ['Xc(6)', 'Yc(7)', 'Zc(8)']
        grid_str = '.%06d' % grid
        grid_str = str(grid_str)
        format_vweb = '.Vweb-csv'
        #base_path = '/media/edoardo/Elements/CLUES/DATA/FullBox/VWeb/vweb_'
        base_path = '/home/edoardo/CLUES/DATA/Vweb/FullBox/vweb_'
        out_base = 'output/vweb_fullbox_mw.00.pkl'
        code_col = 'sub_code'
        n_data = len(data)
        codes = [0] * n_data
        data[code_col] = np.array(codes)
        sub_run = data[code_col].unique()

    n_lgs = len(data)
    all_evs = [] #np.zeros((3, n_lgs))

    #kpcfac = 1.0
    kpcfac = 1.e+3

    if os.path.exists(out_base) and load_file:
        print('Loading from file: ', out_base)
        all_evs = pkl.load(open(out_base, 'rb'))
        all_evs = np.array(all_evs)

    else:

        # now loop on all the simulations and gather data
        for sub in tqdm(sub_run):
            if simu == 'fullbox' or simu == 'rand_mw':
                sub = '%02d' % sub
                this_vweb = base_path + sub + grid_str + format_vweb

            elif simu == 'lgf':
                this_vweb = base_path + sub + format_vweb
                
            # output file base path
            if os.path.exists(this_vweb):
                vweb = pd.read_csv(this_vweb)
        
                for i, row in data[data[code_col] == int(sub)].iterrows():
                    center = np.reshape(row[x_cols].values, (3,1))
                
                    if simu == 'lgf':
                        vweb['D'] = t.apply_distance(data=vweb, x_col=g_cols, center=center, col='D')
                        d_min = vweb['D'].min()
                        ev_select = vweb[vweb['D'] == d_min][l_cols].values
                        all_evs.append(ev_select)
                    
                    elif simu == 'fullbox' or simu == 'rand_mw':

                        index = t.find_nearest_node_index(x=center, grid=grid, box=100.0e+3)

                        # Some LGs might be in the buffer zone of the periodic boundaries
                        if index < grid  * grid * grid and index > 0:
                            row = vweb.loc[index]
                            #row = vweb[l_cols].iloc[index]
                            ev_select = row[l_cols].values
                            #print(ev_select)
                            all_evs.append(list(ev_select))

        print('Saving to file: ', out_base)
        pkl.dump(all_evs, open(out_base, 'wb'))

    return all_evs


def plot_vweb_fb(grid=128, evs_lg=None, evs_fb=None, add_mw=False, evs_mw=None):
    """ Simple eigenvalue distribution plot """

    l_str = [r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$']
    l_cols = ['l1', 'l2', 'l3']

    data_lg = pd.DataFrame(columns=l_cols)
    data_fb = pd.DataFrame(columns=l_cols)

    if add_mw:
        data_mw = pd.DataFrame(columns=l_cols)  

    percentiles0 = [30, 50, 70]
    percentiles = [25, 50, 75]

    for i, col in enumerate(l_cols):
        data_lg[col] = evs_lg[:, 0, i]
        data_fb[col] = evs_fb[:, i] * 1.15

        if add_mw:
            data_mw[col] = evs_mw[:, i] * 2.25

    color0 = 'black'
    color1 = 'blue'
    color2 = 'green'
    label0 = 'LGF-L'
    label1 = 'RAND'
    label2 = 'MW'
    n_bins = 50
    size = 6

    for i, col in enumerate(l_cols):
        plt.figure(figsize=(size, size))
        x_lg = data_lg[(data_lg[col] != 0.0) & (data_lg[col] != 1.0e-3) & (data_lg[col] < 1.1) & (data_lg[col] > -1.0) & (data_lg[col] != 0.002)][col].values
        x_fb = data_fb[data_fb[col] < 2.0][col].values
        plt.hist(x_lg, bins=n_bins, color=color0, density=True, label=label0, alpha=0.7)
        plt.hist(x_fb, bins=n_bins, color=color1, density=True, label=label1, alpha=0.7)

        if add_mw:
            x_mw = data_mw[data_mw[col] < 4.0][col].values
            plt.hist(x_mw, bins=n_bins, color=color2, density=True, label=label2, alpha=0.7)
            bin_mw = np.percentile(x_mw, q=percentiles)

        plt.xlabel(l_str[i])
    
        bin_lg = np.percentile(x_lg, q=percentiles)
        bin_fb = np.percentile(x_fb, q=percentiles)


        if add_mw:
            lg_str = '$%.3f ^{+%.3f} _{-%.3f}$' % (bin_lg[1], bin_lg[2]-bin_lg[1], bin_lg[1]-bin_lg[0])
            fb_str = '& $%.3f ^{+%.3f} _{-%.3f}$' % (bin_fb[1], bin_fb[2]-bin_fb[1], bin_fb[1]-bin_fb[0])
            mw_str = '& $%.3f ^{+%.3f} _{-%.3f}$ \\\ ' % (bin_mw[1], bin_mw[2]-bin_mw[1], bin_mw[1]-bin_mw[0])
            print(l_str[i], lg_str, fb_str, mw_str)
        else:
            lg_str = '$%.3f ^{+%.3f} _{-%.3f}$' % (bin_lg[1], bin_lg[2]-bin_lg[1], bin_lg[1]-bin_lg[0])
            fb_str = '& $%.3f ^{+%.3f} _{-%.3f}$ \\\ ' % (bin_fb[1], bin_fb[2]-bin_fb[1], bin_fb[1]-bin_fb[0])
            print(l_str[i], lg_str, fb_str)

        file_name = 'output/hist_ev_' + col + '_' + str(grid) + '.png'
        plt.legend()
        plt.tight_layout()
        #print('Saving vweb eigenvalue distribution to: ', file_name)
        plt.savefig(file_name)
        plt.clf()
        plt.clf()
        plt.close()


if __name__ == '__main__':
    """ Wrapper to execute the functions in correct order """

    data_lg = pd.read_csv('output/lg_pairs_1024.csv')
    data_fb = pd.read_csv('output/lg_pairs_FB.csv')
    data_mw = pd.read_csv('output/mw_halos.csv')
    #print(data_lg.head())
    #print(data_fb.head())

    vrad = 0.0
    R = 1300

    data_lg = data_lg[data_lg['Vrad'] < vrad] 
    data_lg = data_lg[data_lg['R'] < R] 
    data_fb = data_fb[data_fb['Vrad'] < vrad] 
    data_fb = data_fb[data_fb['R'] < R] 

    #print(len(data_lg))
    #print(len(data_fb))
    #print(len(select_lg))
    #print(len(select_fb))

    print('Extracting VWeb at LG positions from FullBox.\n'); 
    #print('Plotting VWeb at LG positions from FullBox.\n'); 
    evs_lg = extract_vweb(data=data_lg, simu='lgf', load_file=True)
    evs_fb = extract_vweb(data=data_fb, simu='fullbox', load_file=True)
    evs_mw = extract_vweb(data=data_mw, simu='rand_mw', load_file=True)
    print(evs_fb[:, 0])

    grid = 128
    plot_vweb_fb(grid=grid, evs_lg=evs_lg, evs_fb=evs_fb, add_mw=True, evs_mw=evs_mw)


