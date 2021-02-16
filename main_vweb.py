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
    x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
    l_cols = ['l1', 'l2', 'l3', 'dens']
    g_cols = ['x', 'y', 'z']
    df_cols = ['ID', 'x', 'y', 'z', 'l1', 'l2', 'l3', 'dens']

    # Read in the cosmic web
    if simu == 'fullbox':
        grid_str = '.%06d' % grid
        grid_str = str(grid_str)
        format_vweb = '.Vweb-csv'
        #base_path = '/media/edoardo/Elements/CLUES/DATA/FullBox/VWeb/vweb_'
        base_path = '/home/edoardo/CLUES/DATA/Vweb/FullBox/vweb_'
        out_base = 'output/vweb_fullbox_vweb.' + grid_str + '.pkl'
        code_col = 'sub_code'
        sub_run = [2,3,4]

    elif simu == 'lgf':
        grid_str = '%04d' % grid
        grid_str = str(grid_str)
        format_vweb = '.csv'
        base_path = '/home/edoardo/CLUES/DATA/Vweb/512/CSV/vweb_center_' + grid_str + '_'
        #base_path = '/home/edoardo/Elements/CLUES/DATA/VWeb/vweb_center_' + grid_str + '_'
        out_base = 'output/vweb_lgf512.' + grid_str + '.pkl'
        code_col = 'simu_code'
        sub_run = data[code_col].unique()

    n_lgs = len(data)
    all_evs = [] #np.zeros((3, n_lgs))

    #kpcfac = 1.0
    kpcfac = 1.e+3

    if os.path.exists(out_base) and load_file:
        print('Loading from file: ', out_base)
        all_evs = pkl.load(open(out_base, 'rb'))

    else:

        # now loop on all the simulations and gather data
        for sub in tqdm(sub_run):
            if simu == 'fullbox':
                sub = '%02d' % sub
                this_vweb = base_path + sub + grid_str + format_vweb

            elif simu == 'lgf':
                this_vweb = base_path + sub + format_vweb
                    
            
                    
            # output file base path
            if os.path.exists(this_vweb):
                vweb = pd.read_csv(this_vweb)
                print(vweb.head())
                for i, row in data[data[code_col] == int(sub)][0:100].iterrows():
                    center = np.reshape(row[x_cols].values, (3,1))
                    
                    if simu == 'lgf':
                        vweb['D'] = t.apply_distance(data=vweb, x_col=g_cols, center=center, col='D')
                        d_min = vweb['D'].min()
                        ev_select = vweb[vweb['D'] == d_min][l_cols].values
                        all_evs.append(ev_select)
                    
                    elif simu == 'fullbox':

                        index = t.find_nearest_node_index(x=center, grid=grid, box=100.0e+3)
            
                        # Some LGs might be in the buffer zone of the periodic boundaries
                        if index < grid  * grid * grid:
                            row = vweb.loc[index]
                            #row = vweb[l_cols].iloc[index]
                            ev_select = row[l_cols]
                            all_evs.append(ev_select)

        print('Saving to file: ', out_base)
        pkl.dump(all_evs, open(out_base, 'wb'))

    return np.array(all_evs)


def plot_vweb_fb(grid=128, evs_lg=None, evs_fb=None):
    """ Simple eigenvalue distribution plot """

    l_str = [r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$']
    l_cols = ['l1', 'l2', 'l3']
    #plt.xlim(-1.0, 1.5)

    color0 = 'black'
    color1 = 'blue'
    n_bins = 100

    for i, l in enumerate(l_cols):
        sns.histplot(evs_lg[i, :], bins=n_bins, color=color0)
        sns.histplot(evs_fb[i, :], bins=n_bins, color=color1)
        plt.ylabel(l_str[i])
    
        file_name = 'output/hist_ev_' + l_cols[i] + '_' + str(grid) + '.png'
        plt.tight_layout()
        print('Saving vweb eigenvalue distribution to: ', file_name)
        plt.savefig(file_name)
        plt.clf()
        plt.clf()
        plt.close()


if __name__ == '__main__':
    ''' Wrapper to execute the functions in correct order '''

    #print('Extracting VWeb from CS.\n'); extract_vweb_cs()
    #print('Extracting VWeb at LG positions in CS.\n'); extract_vweb_lg_cs()

    data_lg = pd.read_csv('output/lg_pairs_1024.csv')
    data_fb = pd.read_csv('output/lg_pairs_FB.csv')
    print(data_lg.head())
    print(data_fb.head())

    vrad = 0.0
    R = 1300

    data_lg = data_lg[data_lg['Vrad'] < vrad] 
    data_lg = data_lg[data_lg['R'] < R] 
    data_fb = data_fb[data_fb['Vrad'] < vrad] 
    data_fb = data_fb[data_fb['R'] < R] 

    print(len(data_lg))
    print(len(data_fb))
    #print(len(select_lg))
    #print(len(select_fb))

    print('Extracting VWeb at LG positions from FullBox.\n'); 
    #print('Plotting VWeb at LG positions from FullBox.\n'); 
    evs_lg = extract_vweb(data=data_lg, simu='lgf', load_file=True)
    evs_fb = extract_vweb(data=data_fb, simu='fullbox', load_file=False)
    grid = 128
    plot_vweb_fb(grid=grid, evs_lg=evs_lg, evs_fb=evs_fb)


