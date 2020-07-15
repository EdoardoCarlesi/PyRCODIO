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
from matplotlib import rc


'''
    This function plots the simple mass density starting from a particle distribution.
    The plot is a 2D projection.
'''
def plot_density(data=None, axes_plot=None, file_name=None, legend=False, show_plot=False, grid_size=100, margin=0.5, data_augment=False, 
            hex_plot = True, fig_size=10, velocity=False, vel=None, colorscale=None):
    print('Plotting density slices...')

    if (velocity == False) and (vel != None):
        print('Velocity is set to false but vel= is different than None! Seting vel to None...')
        vel = None

    # If the colorscale is undefined by the user then set a default one among those below
    if colorscale == None:        
        #colorscale = 'inferno'
        colorscale = 'gray'
        #colorscale = 'hot'
        #colorscale = 'gist_gray'
        #colorscale = 'bwr'


    # Other plot properties
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

    if hex_plot == True:
        # Do the plot using hexagonal bins
        plt.hexbin(data[coord[ax0]], data[coord[ax1]], gridsize=grid_size, cmap=colorscale, bins='log')
    else:
        plt.hist2d(data[coord[ax0]], data[coord[ax1]], bins=grid_size, cmap=colorscale)

    # Actually show the plot, otherwise we will just save it to a .png file
    if show_plot == True:
        plt.show()
    
    # Save the file
    if legend == True:
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

    # Close all figures and clear stuff
    plt.clf()
    plt.cla()
    plt.close()



'''
    Plot the MAH of a single halo
'''
def plot_mass_accretion(time, mah, f_out, size=10, scale='lin'):

    # Plot properties
    lnw = 1.0
    col = 'b'

    # Max & min x/y scale
    x_min = np.min(time); x_max = np.max(time)
    y_min = np.min(mah); y_max = np.max(mah)

    plt.figure(figsize=(size, size))
    plt.yscale('log')

    # If scale is set to log then also to the x axis in log scale
    if scale == 'log':
       plt.xscale('log')

    plt.axis([x_min, x_max, y_min, y_max])

    # Do the actual plot
    plt.plot(time, mah, linewidth=lnw, color=col)

    # Save the figure and clean plt
    plt.savefig(f_out)
    plt.clf()
    plt.cla()
    plt.close()


'''
    Plot the MAH for an ensamble of halos
'''
def plot_mass_accretions(time, mahs, f_out, percentiles=[25, 50, 75], size=10, scale='lin'):

    # Plot properties
    line_width = 1.0
    plt.figure(figsize=(size, size))
    plt.yscale('log')
    plt.rc({'text.usetex': True})
    plt.xlabel('GYr')
    plt.ylabel('M')

    # Set plot max and min 
    x_min = np.min(time); x_max = np.max(time)
    y_min = np.min(mahs); y_max = np.max(mahs)

    # Set axes
    axs.axis([x_min, x_max, y_min, y_max])

    # Change x axis scale
    if scale == 'log':
        axs.set_xscale('log')

    # TODO: what is this??
    if y_min < 1.e+5:
        y_min = y_max / 200.

    # Find number of steps and plots
    n_steps = len(time)
    n_plots = len(mahs[:,0])

    # Initialize empty list
    all_mahs = [[] for i in range(0, n_steps)]

    # Initialize all_mah list
    for istep in range(0, n_steps):
        for iplot in range(0, n_plots):
            this_mahs = mahs[iplot, istep]
            all_mahs[istep].append(this_mahs)

    # Initialize some empty lists
    med_mah = [[] for i in range(0, n_steps)];
    min_mah = [[] for i in range(0, n_steps)];
    max_mah = [[] for i in range(0, n_steps)];

    # For each step of the MAH find the percentiles (min, med, max) corresponding to that mass bin
    for istep in range(0, n_steps):
        min_mah[istep] = np.percentile(all_mahs[istep], percentile[0])
        med_mah[istep] = np.percentile(all_mahs[istep], percentile[1])
        max_mah[istep] = np.percentile(all_mahs[istep], percentile[2])

    # Plot the median with a solid line + the gray shaded area
    axs.plot(time, med_mah, color='black')
    axs.fill_between(time, min_mah, max_mah, facecolor='grey')

    # Save the figure and clean plt
    plt.tight_layout()
    plt.savefig(f_out)
    plt.clf()
    plt.cla()


'''
    Plot a given set of halo mass functions, with variance and Poissonian error bars
'''
def plot_massfunctions(x_m, y_m, n_mf, f_out, n_bins=10, percentiles=[25, 50, 75], Poisson=True, size=10):

    # Initialize some plot properties
    line_width = 1.0
    color='black'
    fill_color='grey'
    poisson_color='red'

    plt.figure(figsize=(size, size))

    # This has to be initialized as a list as these bins will contain different numbers of entries (some MFs might not have haloes at large masses)
    y_bins = [ [] for i in range(0, n_bins-1) ]

    # Initialize to very high/low values so they can be reset
    x_min = 1.e+15; x_max = 1.e+7; y_max = 1.0;

    # Find the max / min values among all the mass functions
    for im in range(0, n_mf):
    
        # Try - because some bins might be empty and getting a zero might cause issues
        try:
            x_max0 = np.max(x_m[im])
            x_min0 = np.min(x_m[im])
            y_max0 = np.max(y_m[im])
            y_min0 = np.min(y_m[im])

            # Set the new maxima and minima
            if x_max0 > x_max:
                x_max = x_max0

            if x_min0 < x_min:
                x_min = x_min0

            if y_max0 > y_max:
                y_max = y_max0

        except:
            'Do nothing just avoid crashing'

    # The mass function minimum is always set to zero
    y_min = 0.0

    # Shift the maximum and minimum about an eps factor
    eps = 0.1
    fac_min = 1.0 - eps
    fac_max = 1.0 + eps

    # Now generate a new set of mass bins in log space
    x_bins = np.logspace(np.log10(x_min * fac_min), np.log10(x_max * fac_max), num=n_bins, endpoint=True)
    
    for im in range(0, n_mf):
        n_mm = len(x_m[im])

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

    if Poisson == True:
        mf_poisson = np.zeros((2, n_bins)) 

    mf_numeric = np.zeros((3, n_bins))
    mass_bins = []

    nmedStep = 0;     nminStep = 0;     nmaxStep = 0

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
            nmax0 = np.percentile(y_bins[km], percentiles[2])
            nmin0 = np.percentile(y_bins[km], percentiles[0])
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

            if Poisson == True:
               mf_poisson[1].append(nmed + np.sqrt(nmed))
               mf_poisson[0].append(nmed - np.sqrt(nmed))

    y_max = max(mf_poisson[1])
    x_max = np.log10(x_max)
    x_min = np.log10(x_min)
    (fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(size, size))

    plt.rc({'text.usetex': True})
    plt.xlabel(r'$log_{10}M_{\odot} h^{-1}$')
    #plt.ylabel('$N(>M)$')
    axs.set_yscale('log'); 
    #plt.ylabel(r'$log_{10}N(>6.25\times 10^7 M_{\odot})$'); y_min=1.0
    plt.ylabel(r'$log_{10}N(>5\times 10^8 M_{\odot})$'); y_min=1.0

    #y_max = 1.8 * y_max
    axs.axis([x_min, x_max, y_min, y_max])

    # Do the actual plot
    axs.plot(mass_bins, mf_median, linewidth=3, color=color)
    axs.fill_between(mass_bins, mf_min, mf_max, facecolor=fill_color)

    # Add poissonian error bars
    if Poisson == True:
        axs.plot(mass_bins, mf_poisson[1], linewidth=2, dashes=[2, 5], color=poisson_col)
        axs.plot(mass_bins, mf_poisson[0], linewidth=2, dashes=[2, 5], color=poisson_col)

    # Save file and clean plt
    plt.tight_layout()
    plt.savefig(f_out)
    plt.clf()
    plt.cla()
    plt.close()



