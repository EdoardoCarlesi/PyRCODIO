'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    main_halo_pca_extract.py: extract halo properties and dump them to a separate file
'''

import numpy as np
import read_files as rf
import halo_utils as hu
import config as cfg
import pickle as pkl
import pandas as pd
import tools as t
import os


def plot_density():

    # Local data path
    data_path = '/home/edoardo/CLUES/PyRCODIO/data/'

    # Full dataset
    #base_path = '/media/edoardo/data1/DATA/'
    base_path = '/media/edoardo/Elements/CLUES/DATA/2048/'
    snap_path = 'snapshot_054'

    mpc2kpc = 1.e+3

    # Plot properties
    side_size = 3.0 * mpc2kpc
    thickness = 1.5 * mpc2kpc
    n_files = 1
    frac = 1.00
    units = 'kpc'
    part_types = [4]
    grid_size = 500
    rand_seed = 1
    fig_size = 8
    version = ''
    show_plot = False
    velocity = False
    augment = False
    shift = False
    legend = False
    randomize = True
    vel_components = ['Vx', 'Vy', 'Vz']
    n_min_part = 1000

    # Configure the LG model and subpaths
    code_run = cfg.simu_runs()
    sub_run = cfg.gen_runs(0, 10)

    # Now loop on all the simulations and gather data
    for code in code_run:

        for sub in sub_run:

            this_path = base_path + code + '/' + sub + '/'
            this_snap = this_path + snap_path
            input_all_csv = 'output/clusters_' + code + '_' + sub + '.csv' 

            try:
                data_all = pd.read_csv(input_all_csv)
                com = ['Xc(6)', 'Yc(7)', 'Zc(8)']
                these_com = data_all[com].values

            except:
                print('Error, file: ', input_all_csv, ' could not be found.')

            # Check that file exists
            if os.path.isfile(this_snap) and os.path.isfile(input_all_csv):

                print('Found: ', len(these_com), ' clusters in ', this_snap)

                # Read the full snapshot here
                part_df = rf.read_snap(file_name=this_snap, velocity=velocity, part_types=part_types, n_files=1)

                for i, this_com in enumerate(these_com):

                    for z_axis in range(0, 3):

                        # Select particles for the plot, do a selection first. z_axis is orthogonal to the plane of the projection
                        ax0 = (z_axis + 1) % 3
                        ax1 = (z_axis + 2) % 3

                        if shift == True:
                            'Implement this function'
                            center = t.shift(this_com, side_size * 0.5)
                            this_fout = 'output/cluster_shift_' + code + '_' + sub + '.' + str(i) + '.'
                        else:
                            center = this_com
                            this_fout = 'output/cluster_' + version + code + '_' + sub + '.' + str(i) + '.'
                            this_fout_rand = 'output/cluster_' + version + code + '_' + sub + '.' + str(i) + '.rand.'

                        print('Plot axes, z=', z_axis, ' ax0=', ax0, ' ax1=', ax1, ' center of mass: ', this_com)
                
                        try:
                            # Select a slab around a given axis, this function returns a dataframe
                            slab_part_df = t.find_slab(part_df=part_df, side=side_size, thick=thickness, center=center, reduction_factor=frac, z_axis=z_axis, rand_seed=rand_seed)

                           # Do a plot only if there are enough particles
                            if len(slab_part_df) > n_min_part:
                                # Feed the previously chosen dataframe and plot its 2D density projection
                                pu.plot_density(data=slab_part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, legend=legend,
                                    grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size, velocity=velocity, vel=vel_components)
                        except:
                            print('Could not generate a plot for: ', this_snap, '. Data read error.')

                        
                       # Add some scatter to the plot properties 
                        if randomize == True:
                            print('Randomizing ...')

                            # The randomization of the center needs to be very small - a small % on tens of kpc might shift the LG out of the frame
                            for i, c in enumerate(center):
                                center[i] = c - thickness * 0.25 + 0.5 * np.random.uniform()
        
                            this_fout = this_fout_rand 
                            thickness0 = thickness * (0.9 + 0.4 * np.random.uniform()) 
                            side_size0 = side_size * (0.9 + 0.4 * np.random.uniform()) 
                            grid_size0 = int(grid_size * (0.75 + 0.5 * np.random.uniform()))
                            fig_size0 = fig_size * (0.8 + 0.5 * np.random.uniform())
                            frac0 = frac 

                            try:
                                print('Adding randomized plots around: ', center, ' with size: ', side_size0, ' and frac: ', frac0, ' grid: ', grid_size0, ' fig_size: ', fig_size0)
                                slab_part_df = t.find_slab(part_df=part_df, side=side_size0, thick=thickness0, center=center, 
                                        reduction_factor=frac0, z_axis=z_axis, rand_seed=rand_seed)

                                if len(slab_part_df > n_min_part):
                                    # Feed the previously chosen dataframe and plot its 2D density projection
                                    pu.plot_density(data=slab_part_df, axes_plot=[ax0, ax1], file_name=this_fout, show_plot=show_plot, legend=legend,
                                        grid_size=grid_size, margin=0.1, data_augment=augment, fig_size=fig_size0, velocity=velocity, vel=vel_components)

                            except:
                                print('Could not print slices for file: ', this_fout)


def extract_halo():

    # Resimulation mode
    #mode='2048'
    mode='1024'

    base_path = '/media/edoardo/Elements/CLUES/DATA/LGF/' + mode + '/'

    # Configure the LG model and subpaths
    if mode == '2048':
        code_run = cfg.simu_runs()
        sub_run = cfg.sub_runs()

    elif mode == '1024':
        num_run = cfg.gen_runs(0, 80)
        sub_run = cfg.gen_runs(0, 30)

    # Read csv list of LGs
    lgs_csv = 'output/lg_pairs_' + mode + '.csv'
    df_lgs = pd.read_csv(lgs_csv)
    print('TotLen: ', len(df_lgs))

    x_cols = ['Xc_LG', 'Yc_LG', 'Zc_LG']
    x_cols_ahf = ['Xc(6)', 'Yc(7)', 'Zc(8)']
    box_center = np.array([5e+4] * 3)

    r_max = 10000.0

    '''
    # Do some LG filtering
    v_max = - 0.0
    d_max = 7000.0

    df_lgs = df_lgs[df_lgs['Vrad'] < v_max]
    print('TotLen: ', len(df_lgs))

    # Select LGs depending on their distance from the box center
    '''

    #file_ahf = 'snapshot_full_054.z0.000.AHF_halos'
    #file_ahf = 'snapshot_full_054.z0.001.AHF_halos'
    file_ahf = 'snapshot_054.z0.000.AHF_halos'

    # Loop on all the simulations and gather data
    if mode == '2048':
        for code in code_run:

            for sub in sub_run:
                this_path = base_path + code + '/' + sub + '/'
                this_ahf = this_path + file_ahf

                # Check that file exists
                if os.path.isfile(this_ahf):
                    print('Reading AHF file: ', this_ahf)
                    halo_df = rf.read_ahf_halo(this_ahf, file_mpi=False)
                    halo_select_df = halo_df[halo_df['Mvir(4)'] > m_thresh]
                    out_all_halo_csv = out_base + code + '_' + sub + '.csv'
                    halo_select_df.to_csv(out_all_halo_csv)

    elif mode == '1024':
        
        for ilg, row in df_lgs.iterrows():
            num = str('%02d' % int(row['simu_code']))
            sub = str('%02d' % int(row['sub_code']))

            this_ahf = base_path + num + '_'  + sub + '/' + file_ahf
            this_x = row[x_cols].values
            #print(this_x)

            if os.path.isfile(this_ahf):
                
                halo_df = rf.read_ahf_halo(this_ahf, file_mpi=True)
                halo_df['Dist'] = halo_df[x_cols_ahf].T.apply(lambda x: t.distance(x, this_x))
                halo_df = halo_df[halo_df['Dist'] < r_max]
                #print(halo_df.head())

                f_out = 'output/LG_' + mode + '/lg_center_' + num + '_' + sub + '.' + str(ilg) + '.csv'

                print('Saving to: ', f_out)
                halo_df.to_csv(f_out)

if __name__ == "__main__":
    """ Wrapper """
