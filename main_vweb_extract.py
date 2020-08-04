'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com
    
    main_extract_vweb.py: extract v-web data around a given point
'''

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os

'''
    Extract the eigenvalues at the box center in constrained simulations
'''
def extract_vweb_cs():
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

'''
    Extract the eigenvalues/vectors at LG positions in CS simulations
'''
def extract_vweb_lg_cs():
    out_path = ''
    sub_runs = cfg.gen_runs(0, 100)
    # TODO

    
'''
    Extract the Vweb at given positions in random full box simulations
'''
def extract_vweb_fb():
    # configure subpaths
    sub_run = cfg.gen_runs(4, 5)

    # local data path, file names and file format
    base_vweb = 'vweb_'
    grid_size  = 256
    grid = '%03d' % grid_size
    format_vweb = '.000' + grid + '.Vweb-csv'
    
    # Columns to be extracted from the vweb and lg dataframes
    x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
    l_cols = ['l1', 'l2', 'l3', 'dens']
    g_cols = ['x', 'y', 'z']
    df_cols = ['ID', 'x', 'y', 'z', 'l1', 'l2', 'l3', 'dens']

    # full dataset
    base_path = '/home/edoardo/CLUES/DATA/Vweb/'
    out_base = 'output/lg_fullbox_vweb_' + grid
    lg_base = 'output/lg_fullbox_'

    #kpcfac = 1.0
    kpcfac = 1.e+3

    # now loop on all the simulations and gather data
    for sub in sub_run:

        # output file base path
        this_vweb = base_path + base_vweb + sub + format_vweb
        this_lgs = lg_base + sub + '.csv'
        out_file = out_base + '_' + sub + '.csv'
        out_file_csv = out_base + '_' + sub + '.csv'

        vweb = pd.read_csv(this_vweb)
        lgs = pd.read_csv(this_lgs)
        lg_ev_df = pd.DataFrame(columns=df_cols)

        for ind, row in lgs.iterrows():
            this_x = row[x_cols]/kpcfac
            this_ind = t.find_nearest_node_index(this_x, grid=grid_size, box=100.0)
            this_l = vweb[l_cols].iloc[this_ind]
            this_g = vweb[g_cols].iloc[this_ind]
            this_row = (ind, this_g[0], this_g[1], this_g[2], this_l[0], this_l[1], this_l[2], this_l[3])
            new_row = pd.Series(this_row, index=df_cols)
            lg_ev_df = lg_ev_df.append(new_row, ignore_index=True)

        print('Saving file to csv: ', out_file_csv)
        lg_ev_df.to_csv(out_file_csv)

'''
    Simple eigenvalue distribution plot
'''
def plot_vweb_fb(grids = [32, 64, 128]):
    
    l1 = 'l1_'
    l2 = 'l2_'
    l3 = 'l3_'

    l_cols = ['l1', 'l2', 'l3', 'dens']
    for grid in grids[2:3]:
        this_web = rf.read_lg_vweb(grid_size = grid)
        plt.xlim(-1.0, 1.5)
        sns.distplot(this_web[l_cols[0]], bins=100)
        sns.distplot(this_web[l_cols[1]], bins=100)
        sns.distplot(this_web[l_cols[2]], bins=100)
        print(this_web.info())


    plt.show()


    ############################################################
    #      Main program: Choose what kind of web to extract    #
    ############################################################

#print('Extracting VWeb from CS.\n'); extract_vweb_cs()
#print('Extracting VWeb at LG positions in CS.\n'); extract_vweb_lg_cs()
#print('Extracting VWeb at LG positions from FullBox.\n'); extract_vweb_fb()
print('Plotting VWeb at LG positions from FullBox.\n'); plot_vweb_fb()
