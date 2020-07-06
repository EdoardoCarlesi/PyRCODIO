'''
    Python Routines for COsmology and Data I/O
    PyRCODIO Pandas Version
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    plot_utils.py = functions and routines used to plot
'''

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import seaborn as sns
import pandas as pd
import numpy as np
import math
from pygadgetreader import *
from matplotlib import rc


'''
    Find the particles belonging to a slab around a given point in space.
    Slab size, thickness and so on need to be specified.
'''
def find_slab(file_name=None, center=None, side=None, thick=None, velocity=False, rand_seed=69,
        reduction_factor=1.0, z_axis=2, part_type=1, units='kpc', n_files=1):

    # Set some parameters
    kpcThresh = 1.e+4
    kpc2Mpc = 1.e-3

    minima = np.zeros((3))
    maxima = np.zeros((3))

    # TODO: fix for more files (n_files > 1 needs to be implemented)
    if n_files > 1:
        print('Multiple files are not supported yet.')
        return -1       

    # Do we want to read in velocity data as well? Or positions only?
    if velocity == True:
        print('Reading velocity and position particle info...')
        particles = readsnap(file_name, 'pos', part_type)
        velocities = readsnap(file_name, 'vel', part_type)
    
        # Concatenate the columns 
        full_data = np.concatenate((particles, velocities), axis=1)
        cols = ['X', 'Y', 'Z', 'Vx', 'Vy', 'Vz']

    else:
        print('Reading position particle info...')
        particles = readsnap(file_name, 'pos', part_type)
        full_data = particles
        cols = ['X', 'Y', 'Z']

    # Read the snapshot
    n_part = len(full_data)

    print('Found ', n_part, ' particles in total.')
    part_df = pd.DataFrame(data=full_data, columns=cols)

    # Select the two axes for the 2D projection
    ax0 = (z_axis + 1) % 3
    ax1 = (z_axis + 2) % 3
    ax2 = z_axis

    # Column names
    col0 = cols[ax0]
    col1 = cols[ax1]
    col2 = cols[ax2]

    # Sanity check on the units
    half_n = int(n_part * 0.5)
    sum_coord = particles[half_n][0] + particles[half_n][1] + particles[half_n][2]
    
    # Make sure the units are consistent
    if sum_coord < kpcThresh:
        side = side * kpc2Mpc
        center = center * ([kpc2Mpc] *3) 
        thick = thick * kpc2Mpc

    # Set the minima and maxima for the particles to be used in the plot
    minima[ax0] = center[ax0] - side * 0.5
    minima[ax1] = center[ax1] - side * 0.5
    minima[ax2] = center[ax2] - thick * 0.5

    maxima[ax0] = center[ax0] + side * 0.5
    maxima[ax1] = center[ax1] + side * 0.5
    maxima[ax2] = center[ax2] + thick * 0.5

    # Find the particles in the slab
    condition_x = (part_df[col0] > minima[ax0]) & (part_df[col0] < maxima[ax0])
    condition_y = (part_df[col1] > minima[ax1]) & (part_df[col1] < maxima[ax1])
    condition_z = (part_df[col2] > minima[ax2]) & (part_df[col2] < maxima[ax2])
    part_select = part_df[(condition_x) & (condition_y) & (condition_z)]

    print('Found: ', len(part_select), ' particles in the slab')

    # Now select a random subsample of the full particle list
    if reduction_factor < 1.0:
        part_select = part_select.sample(frac=reduction_factor, random_state=rand_seed)

        print('The number of particles to be used has been reduced to: ', len(part_select))

    # Return the selected particles' properties in a dataframe
    return part_select


'''
    This function plots the simple mass density starting from a particle distribution.
    The plot is a 2D projection.
'''
def plot_density(data=None, axes_plot=None, file_name=None, legend=False, show_plot=False, grid_size=100, margin=0.5, data_augment=False, fig_size=10, velocity=False, vel=None):
    print('Plotting density slices...')

    if (velocity == False) and (vel != None):
        print('Velocity is set to false but vel= is different than None! Seting vel to None...')
        vel = None

    # Plot properties
    colorscale = 'inferno'
    #colorscale = 'gray'
    #colorscale = 'hot'
    #colorscale = 'gist_gray'
    #colorscale = 'bwr'
    ax0 = axes_plot[0]
    ax1 = axes_plot[1]
    coord = ['X', 'Y', 'Z']
    axis_label = ['SGX', 'SGY', 'SGZ']

    plt.figure(figsize=(fig_size, fig_size))

    # If we are going to use the images with CNNs then by default legend is set to False
    if legend == True:
        axis_size = fig_size * 2
        plt.rc({'text.usetex': True})
        plt.rc('axes',  labelsize=axis_size)
        plt.rc('xtick', labelsize=axis_size)
        plt.rc('ytick', labelsize=axis_size)
        plt.xlabel(r'$h^{-1}$Mpc')
        plt.ylabel(r'$h^{-1}$Mpc')
        #plt.xlabel(axis_label[ax0]+' '+axis_units)
        #plt.ylabel(axis_label[ax1]+' '+axis_units)

        file_out = file_name + 'density_' + axis_label[ax0] + axis_label[ax1]
    else:
        axis_size = 0
        file_out = file_name + 'rho_no_labels_' + axis_label[ax0] + axis_label[ax1]

    # Find the maxima and minima of the plot, reduce by a small margin to remove border imperfections
    x_min = data[coord[ax0]].min() + margin
    y_min = data[coord[ax1]].min() + margin
    x_max = data[coord[ax0]].max() - margin
    y_max = data[coord[ax1]].max() - margin

    # Set the axis maxima and minima
    plt.axis([x_min, x_max, y_min, y_max])

    # Unset the ticks and set margins to zero if legend is true
    if legend == False:
        plt.margins(0.0) 
        plt.xticks([])
        plt.yticks([])
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    else:
        plt.tight_layout()

    # Do the plot using hexagonal bins
    plt.hexbin(data[coord[ax0]], data[coord[ax1]], gridsize=grid_size, cmap=colorscale, bins='log')
    
    # Actually show the plot, otherwise we will just save it to a .png file
    if show_plot == True:
        plt.show()
    
    # Save the file
    if legend == False:
        plt.tight_layout()

    plt.savefig(file_out + '.png')
    print('File saved to: ', file_out + '.png')

    # Check if we want to plot velocities as well
    if velocity == True:
        if vel == None:
            print('Error: vel keyword is not set and velocity is set to True.')
            return 0
        else:
            # Check if vel is a list, otherwise make it a list
            if isinstance(vel, list) == False:
                vel = [vel]
    
            # Loop on velocity projections
            for vx in vel:
                file_v_out = file_name + 'V_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_' + vx
                plt.hexbin(data[coord[ax0]], data[coord[ax1]], C=data[vx], gridsize=grid_size, cmap=colorscale)
                plt.savefig(file_v_out + '.png')
                print('File saved to: ', file_v_out + '.png')

    # Do some transformations on the data to increase the number of samples
    if data_augment == True:
        print('Data augmentation is set to True, printing additional images...')

        file_out = file_name + 'rho_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_augment0'
        plt.axis([x_max, x_min, y_min, y_max])
        plt.hexbin(data[coord[ax0]], data[coord[ax1]], gridsize=grid_size, cmap=colorscale, bins='log')
        plt.savefig(file_out + '.png')
        #print('File saved to: ', file_out + '.png')

        file_out = file_name + 'rho_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_augment1'
        plt.axis([x_min, x_max, y_max, y_min])
        plt.hexbin(data[coord[ax0]], data[coord[ax1]], gridsize=grid_size, cmap=colorscale, bins='log')
        plt.savefig(file_out + '.png')
        #print('File saved to: ', file_out + '.png')

        file_out = file_name + 'rho_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_augment2'
        plt.axis([x_max, x_min, y_max, y_min])
        plt.hexbin(data[coord[ax0]], data[coord[ax1]], gridsize=grid_size, cmap=colorscale, bins='log')
        plt.savefig(file_out + '.png')
        #print('File saved to: ', file_out + '.png')

        # Print also "augmented" velocity maps in case
        if velocity == True:
            for vx in vel:
                plt.axis([x_max, x_min, y_min, y_max])
                file_v_out = file_name + 'V_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_' + vx + '_augment0'
                plt.hexbin(data[coord[ax0]], data[coord[ax1]], C=data[vx], gridsize=grid_size, cmap=colorscale)
                plt.savefig(file_v_out + '.png')

                plt.axis([x_min, x_max, y_max, y_min])
                file_v_out = file_name + 'V_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_' + vx + '_augment1'
                plt.hexbin(data[coord[ax0]], data[coord[ax1]], C=data[vx], gridsize=grid_size, cmap=colorscale)
                plt.savefig(file_v_out + '.png')

                plt.axis([x_max, x_min, y_max, y_min])
                file_v_out = file_name + 'V_no_labels_' + axis_label[ax0] + axis_label[ax1] + '_' + vx + '_augment2'
                plt.hexbin(data[coord[ax0]], data[coord[ax1]], C=data[vx], gridsize=grid_size, cmap=colorscale)
                plt.savefig(file_v_out + '.png')














    ##################################################################################################################
    # Plot the LG, if plot_pos = "true" then repeat the computation in the PoS moment of inertial eigenvector system #
    ##################################################################################################################

def plot_lg(f_snap, f_out, lg0, lg1, reduce_fac, ptype, plot_pos):
    facMpc = 1000.
    buffPlot = 1.25 * lg0.r # Extra buffer on the edges, largest Rvir of the two halos
    thickn = 1000.

    # Plot properties
    ptsize_lg = 5.0
    col_mw="red"
    col_m31="orange"
    col_lv="black"
    col_lg="lime"

    plot_col = 3

    if plot_pos == "true":
        plot_row = 2
        size_row = 60
        size_col = 40
    else:
        plot_row = 1
        size_row = 60
        size_col = 20

    npt_lg_min = 100
    axis_margins = 1
    axis_size = 50
    axis_units = 'Mpc/h'
    axis_label = []
    axis_label.append('SGX')
    axis_label.append('SGY')
    axis_label.append('SGZ')

    x_ptslg0 = [[] for ix in range(0, 3)]
    x_ptslg1 = [[] for ix in range(0, 3)]
    x_ptslgall = [[] for ix in range(0, 3)]

    # Read the particles of type ptype ONLY (high-res region!!!)
    partDM = readsnap(f_snap, 'pos', ptype)

    # Identify the particles belonging to the different objects
    n_ptlg0 = 0
    n_ptlg1 = 0

    # Identify M31 and MW TODO Just plot all the particles for the moment
    #t1 = time.clock()
    #(x_ptslg0[0], x_ptslg0[1], x_ptslg0[2]) = select_particles_halo(lg0, partDM)
    #(x_ptslg1[0], x_ptslg1[1], x_ptslg1[2]) = select_particles_halo(lg1, partDM)
    #t2 = time.clock()
    #print 'Particles found in %.3f s. ' % (t2-t1)

    lgcom = center_of_mass([lg0.m, lg1.m], [lg0.x, lg1.x])
    d_mwm31 = lg0.distance(lg1.x)

    # Center of plot
    c_plot = lgcom
    side_lg = d_mwm31 + buffPlot
    print('Plot center: ', c_plot)
    print('Plot side  : ', side_lg)

    # This vector contains the minima of the XYZ coordinates for the plot
    min_lg_xyz = [0.0] * 3

    # Plot is centered on HALF of the distance between Virgo and LG, then a buffer is added
    for ix in range(0, 3):
        min_lg_xyz[ix] = (c_plot[ix] - d_mwm31 * 0.5 - buffPlot)

    # These contain the particles in each slab along the XY, XZ, YZ plane for the different resolution levels
    x_plotlg = [[] for ix in range(0, 3)]
    y_plotlg = [[] for ix in range(0, 3)]

    # Find slab of thickness +/- thickn around the axes
    for ix in range(0, 3):
        ixp1 = (ix+1) % 3
        ixp2 = (ix+2) % 3

        t1 = time.clock()

        # High res particles
        (x_plotlg[ixp1], y_plotlg[ixp2]) = find_slab(partDM, ix, lgcom, min_lg_xyz, side_lg, thickn, reduce_fac)
        t2 = time.clock()

        print('Slab (%s, %s) found in %.3f s.' % (axis_label[ixp1], axis_label[ixp2], (t2-t1)))
        plt.ylabel(axis_label[ixp2]+' '+axis_units)


    # General plot settings
    plt.figure(figsize=(size_row,size_col))
    plt.rc('xtick', labelsize=axis_size)
    plt.rc('ytick', labelsize=axis_size)
    plt.rc('axes',  labelsize=axis_size)
    plt.margins(axis_margins)

    for ix in range(0, 3):
        ixp1 = (ix+1) % 3
        ixp2 = (ix+2) % 3

        x_min = (min_lg_xyz[ixp1])
        x_max = (x_min + side_lg)
        y_min = (min_lg_xyz[ixp2])
        y_max = (y_min + side_lg)

        # These plots are in Mpc/h not kpc/h
        x_min /= facMpc
        x_max /= facMpc
        y_min /= facMpc
        y_max /= facMpc

        # Plot settings for each subplot
        plt.subplot(plot_row, plot_col, ix+1)
        plt.axis([x_min, x_max, y_min, y_max])
        plt.xlabel(axis_label[ixp1]+' '+axis_units)
        plt.ylabel(axis_label[ixp2]+' '+axis_units)

        #print len(x_plotlg[ixp1])
        #print len(y_plotlg[ixp2])

        # Actual plot
        plt.scatter(x_plotlg[ixp1][:], y_plotlg[ixp2][:], s=ptsize_lg, c=col_mw)
        #plt.scatter(x_plotlg[ixp1][:], y_plotlg[ixp2][:], s=ptsize_lv, c=col_lg)
        #plt.scatter(x_ptslg0[ixp1][:], x_ptslg0[ixp2][:], s=ptsize_lg, c=col_mw)
        #plt.scatter(x_ptslg1[ixp1][:], x_ptslg1[ixp2][:], s=ptsize_lg, c=col_m31)

    # Save to file
    plt.tight_layout()
    plt.savefig(f_out)


def plot_massfunctions(x_m, y_m, n_mf, f_out, n_bins):
    size_x = 20
    size_y = 20
    lnw = 1.0
    col = 'b'
    axis_margins = 2
#       print 'Plotting massfunctions to file: ', n_mf, f_out, y_max

    #n_bins = 15
    y_bins = [ [] for i in range(0, n_bins-1) ]
    #y_bins = np.zeros((3, n_bins))

    x_min = 1.e+15; x_max = 1.e+7
    y_min = 10000.; y_max = 1.0;

    for im in range(0, n_mf):

        try:
            x_max0 = np.max(x_m[im])
            x_min0 = np.min(x_m[im])
            y_max0 = np.max(y_m[im])
            y_min0 = np.min(y_m[im])

            if x_max0 > x_max:
                x_max = x_max0

            if x_min0 < x_min:
                x_min = x_min0

            if y_max0 > y_max:
                y_max = y_max0

            if y_min0 < y_min:
                y_min = y_min0

        except:
            porco = 0.0

    y_min = 0.0
    #print x_min/1.e+9, ' ', x_max/1.e+9, ' ', y_min, ' ', y_max

    x_bins = np.logspace(np.log10(x_min * 0.99), np.log10(x_max * 1.01), num=n_bins, endpoint=True)
    #x_bins = np.logspace(7, 10.0, num=n_bins, dtype=float, base=10.0, endpoint=True)

    #print x_bins

    #x_min = 5.e+8;         x_max = 5.e+11
    #y_min = 1;     y_max = 50 #max_list(y_n)

    for im in range(0, n_mf):
        n_mm = len(x_m[im])
        #m_mm = len(y_m[im])

        for jm in range(0, n_mm):
            m0 = x_m[im][jm]
            y0 = y_m[im][jm]

            for km in range(0, n_bins-1):
                mbin0 = x_bins[km]
                mbin1 = x_bins[km+1]

                if m0 > mbin0 and m0 < mbin1:
                    y_bins[km].append(y0)

            # At the last step set all the remaining bins above m0 to zero
            if jm == n_mm-1:

                for km in range(0, n_bins-1):
                    mbin0 = x_bins[km]
                    mbin1 = x_bins[km+1]

                    if m0 > mbin0 and m0 < mbin1:
                        y_bins[km].append(0.0)

        #if x_m[im][n_mm-1] < 0.7e+11
        #               y_bins[].append(0)

    mf_poisson = [ [] for i in range(0, 3)]
    mf_median = []
    mf_max = []
    mf_min = []
    mass_bins = []
    nmedStep = 0
    nminStep = 0
    nmaxStep = 0

    for km in range(0, n_bins-1):
        mbin0 = x_bins[km]
        mbin1 = x_bins[km+1]
        mmed0 = 0.5 * (mbin0 + mbin1)

        if km == n_bins-2:
            mmed0 = mbin0

        if mbin0 > 1.e+10:
            n_thresh = 0
        else:
            n_thresh = 3

        if len(y_bins[km]) > n_thresh:
            nmax0 = np.percentile(y_bins[km], 100)
            nmin0 = np.percentile(y_bins[km], 0)
            #nmax0 = np.percentile(y_bins[km], 80)
            #nmin0 = np.percentile(y_bins[km], 20)
            nmed0 = np.mean(y_bins[km])

            if km == 0:
                nmedStep = nmed0
                nmaxStep = nmax0
                nminStep = nmin0
            else:
                if nmed0 > nmedStep:
                    nmed0 = nmedStep
                else:
                    nmedStep = nmed0

                if nmin0 > nminStep:
                    nmin0 = nminStep
                else:
                    nminStep = nmin0

                if nmax0 > nmaxStep:
                    nmax0 = nmaxStep
                else:
                    nmaxStep = nmax0

            nmax = nmax0
            nmin = nmin0
            nmed = nmed0
            mmed = np.log10(mmed0)

            mass_bins.append(mmed)
            mf_median.append(nmed)
            mf_min.append(nmin)
            mf_max.append(nmax)
            mf_poisson[0].append(np.sqrt(nmed))
            mf_poisson[1].append(nmed + np.sqrt(nmed))
            mf_poisson[2].append(nmed - np.sqrt(nmed))

    y_max = max(mf_poisson[1])
    x_max = np.log10(x_max)
    x_min = np.log10(x_min)
    (fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(6, 6))

    plt.rc({'text.usetex': True})
    plt.xlabel(r'$log_{10}M_{\odot} h^{-1}$')

    #plt.ylabel('$N(>M)$')
    axs.set_yscale('log'); 
    #plt.ylabel(r'$log_{10}N(>6.25\times 10^7 M_{\odot})$'); y_min=1.0
    plt.ylabel(r'$log_{10}N(>5\times 10^8 M_{\odot})$'); y_min=1.0

    y_max = 1.8 * y_max
    axs.axis([x_min, x_max, y_min, y_max])

    #usecolor='blue'
    #usecolor='red'
    usecolor='grey'
    pois_col='red'

    axs.plot(mass_bins, mf_median, linewidth=3, color='black')
    axs.plot(mass_bins, mf_poisson[1], linewidth=2, dashes=[2, 5], color=pois_col)
    axs.plot(mass_bins, mf_poisson[2], linewidth=2, dashes=[2, 5], color=pois_col)
    axs.fill_between(mass_bins, mf_min, mf_max, facecolor=usecolor)

    #print mf_poisson[0]
    #if n_mf > 1:
    #       for im in range(0, n_mf):
    #               axs.plot(x_m[im], y_n[im], linewidth=lnw, color=col)

    plt.tight_layout()
    plt.savefig(f_out)
    plt.clf()
    plt.cla()
    plt.close()





def plot_mass_accretion(time, mah, f_out):
    size_x = 20
    size_y = 20
    lnw = 1.0
    col = 'b'

    x_min = np.min(time); x_max = np.max(time)
    y_min = np.min(mah); y_max = np.max(mah)
    plt.yscale('log')
    plt.axis([x_min, x_max, y_min, y_max])

    plt.plot(time, mah, linewidth=lnw, color=col)

    plt.savefig(f_out)
    plt.clf()
    plt.cla()
    plt.close()



def plot_mass_accretions(time, mahs, f_out):
    size_x = 20
    size_y = 20
    lnw = 1.0
    col = 'b'

    n_plots = len(mahs[:,0])
    x_min = np.min(time); x_max = np.max(time)
    #y_min = np.log10(np.min(mahs)); y_max = np.log10(np.max(mahs) * 1.5)
    #y_min = np.min(mahs); y_max = np.max(mahs) * 1.5
    y_min = np.min(1.e+9); y_max = np.max(mahs) * 1.5

    if y_min < 1.e+5:
        y_min = y_max / 200.

    n_steps = len(time)

    all_mahs = [[] for i in range(0, n_steps)]

    (fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

    plt.yscale('log')
    axs.axis([x_min, x_max, y_min, y_max])

    for istep in range(0, n_steps):
        for iplot in range(0, n_plots):
            this_mahs = mahs[iplot, istep]
            #all_mahs[istep].append(np.log10(this_mahs))
            all_mahs[istep].append(this_mahs)

    med_mah = [[] for i in range(0, n_steps)];
    min_mah = [[] for i in range(0, n_steps)];
    max_mah = [[] for i in range(0, n_steps)];


    for istep in range(0, n_steps):
        med_mah[istep] = np.percentile(all_mahs[istep], 50)
        min_mah[istep] = np.percentile(all_mahs[istep], 80)
        max_mah[istep] = np.percentile(all_mahs[istep], 20)

    #plt.margins(axis_margins)
    #axs.set_xscale('log')
    plt.rc({'text.usetex': True})
    plt.xlabel('GYr')
    plt.ylabel('M')
    #axs.set_yscale('log')
    axs.axis([x_min, x_max, y_min, y_max])

    #print mf_poisson[0]
    #print mf_poisson[1]

    #print med_mah

    axs.plot(time, med_mah, color='black')
    #axs.plot(mass_bins, min_mah[1], linewidth=4, dashes=[2, 5], color='black')
    #axs.plot(mass_bins, mf_poisson[2], linewidth=4, dashes=[2, 5], color='black')
    #axs.fill_between(time, min_mah, max_mah, facecolor='azure')
    axs.fill_between(time, min_mah, max_mah, facecolor='grey')

    #for iplot in range(0, n_plots):
        #mah = mahs[iplot]
        #axs.plot(time, mah, linewidth=lnw, color=col)

    plt.tight_layout()
    plt.savefig(f_out)
    plt.clf()
    plt.cla()
    #plt.close(fig)
    #plt.gcf().show()

