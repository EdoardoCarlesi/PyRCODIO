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
    axis_data (False by default) prints also units on the axes and the legend.
'''
def find_slab(file_name=None, center=None, side=None, thick=None, velocity=False, rand_seed=69,
        reduction_factor=1.0, z_axis=2, part_type=1, units='kpc', n_files=1, axis_data=False):

    # Set some parameters
    kpcThresh = 1.e+4
    kpc2Mpc = 1.e-3

    minima = np.zeros((3))
    maxima = np.zeros((3))

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


#def plot_density(center, side_size, f_out, nbins, f_rescale, thickn, units, slab, bw_smooth, slab_type):

def plot_density(data=None, axes_plot=None, file_name=None, legend=False):
    print('Plotting density slices...')

    # Plot properties
    axis_margins = 1
    fig_size = 10
    ax0 = axes_plot[0]
    ax1 = axes_plot[1]
    coord = ['X', 'Y', 'Z']

    # If we are going to use the images with CNNs then by default legend is set to False
    if legend == True:
        axis_size = 12
        axis_units = 'Mpc/h'
        axis_label = ['SGX', 'SGY', 'SGZ']
        plt.title('Matter density')
        plt.xlabel(axis_label[ax0]+' '+axis_units)
        plt.ylabel(axis_label[ax1]+' '+axis_units)
    else:
        axis_size = 0
        
    # General plot settings
    plt.figure(figsize=(fig_size,fig_size))
    plt.rc('xtick', labelsize=axis_size)
    plt.rc('ytick', labelsize=axis_size)
    plt.rc('axes',  labelsize=axis_size)
    plt.margins(axis_margins)

    # Find the maxima and minima of the plot. Use a small reduction factor to make sure the picture is well centered
    eps = 0.001
    fac_min = 1.00 + eps
    fac_max = 1.00 - eps

    x_min = data[coord[ax0]].min() * fac_min
    y_min = data[coord[ax1]].min() * fac_min
    x_max = data[coord[ax0]].max() * fac_max
    y_max = data[coord[ax1]].max() * fac_max

    print('XMin: ', x_min, 'YMin: ', y_min)
    print('XMax: ', x_max, 'YMax: ', y_max)


    #sns.kdeplot(data = data[coord[ax0]], data2 = data[coord[ax0]])
    #sns.jointplot(coord[ax0], coord[ax1], data=data, kind='kde')
    sns.jointplot(coord[ax0], coord[ax1], data=data, kind='hex', bins='log')
    plt.show()
    
    '''

    plt.axis([x_min, x_max, y_min, y_max])
    colorscale = 'rainbow'

    levels = [np.percentile(zi, 75), np.percentile(zi, 90), np.percentile(zi, 99)]

    print('Printing contour and color mesh on ', len(levels), ' levels ...')
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=colorscale, shading='gouraud') 
    plt.contour(xi, yi, zi.reshape(xi.shape), levels=levels, linewidth=[10.0], colors='black')
    f_out = f_out + '_dm.png'
    print('Done.')

    # Save to file
    plt.tight_layout()
    plt.savefig(f_out)
    '''
    



def simple_plot_rho(center, side_size, f_out, nbins, f_rescale, thickn, units, slab, bw_smooth, ptype):
    print('Plotting density slices for snapshot: ', slab, ' particle type= ', ptype)

    # Plot properties
    axis_margins = 1
    axis_size = 12
    
    # Select factor margin: select particles slightly outside of the side_size restriction, for nicer binning
    sf = 1.1

    if units == 'Mpc':
        axis_units = 'Mpc/h'; facMpc = 1.
    elif units == 'kpc':
        axis_units = 'Mpc/h'; facMpc = 1000.

    axis_label = ['SGX', 'SGY']

    f_slab = open(slab[0], 'rb')
    (x_plot, y_plot) = pickle.load(f_slab)

    # Overplot stars to gas
    if ptype == 0:
        f_slab4 = open(slab[1], 'rb')
        x_s, y_s = pickle.load(f_slab4)
        print(f_slab4, ' has ', len(x_s), ' star particles.')

    n_x = len(x_plot)
    print('N Part in slab: ', n_x)
    print('Slab (%s, %s) with %d particles found.' % (axis_label[0], axis_label[1], n_x))

    plt.xlabel(axis_label[0]+' '+axis_units)
    plt.ylabel(axis_label[1]+' '+axis_units)

    # General plot settings
    plt.figure(figsize=(8,8))
    plt.rc('xtick', labelsize=axis_size)
    plt.rc('ytick', labelsize=axis_size)
    plt.rc('axes',  labelsize=axis_size)
    plt.margins(axis_margins)

    x_min = -side_size; x_max = side_size
    y_min = -side_size; y_max = side_size

    # These plots are in Mpc/h not kpc/h
    x_min /= facMpc
    x_max /= facMpc
    y_min /= facMpc
    y_max /= facMpc

    print('XMin: ', x_min, ' XMax: ', x_max)
    print('New particles number: ', n_x, ' n bins: ', nbins)
    data_x = [];     data_y = []
    datas_x = [];     datas_y = []
    
    if ptype == 0:
        for iss in range(0, len(x_s)):
            xs = (x_s[iss] - center[0])/facMpc
            ys = (y_s[iss] - center[1])/facMpc
            datas_x.append(xs)
            datas_y.append(ys)


    for ip in range(0, n_x):
        x = (x_plot[ip] - center[0])/facMpc
        y = (y_plot[ip] - center[1])/facMpc

        # Select only particles within the area
        if (x > x_min*sf and x < x_max*sf and y > y_min*sf and y < y_max*sf):
            data_x.append(x) 
            data_y.append(y) 

    data = np.zeros((len(data_x), 2), dtype='float')
    for ip in range (0, len(data_x)):
            data[ip, 0] = data_x[ip]
            data[ip, 1] = data_y[ip]

    # Smooth gas density
    if ptype == 0:
        xi, yi = np.mgrid[min(data_x):max(data_x):nbins*1j, min(data_y):max(data_y):nbins*1j]
        k = kde.gaussian_kde(data.T, bw_method=bw_smooth)
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.axis([x_min, x_max, y_min, y_max])
    plt.xlabel(axis_label[0]+' '+axis_units)
    plt.ylabel(axis_label[1]+' '+axis_units)

    if ptype == 0:
#        plt.title('Gas + Stars')
        plt.title('Gas')
        colorscale = 'rainbow'
        #colorscale = 'viridis'
        #plt.contour(xi, yi, zi.reshape(xi.shape)) #, levels=[1.0])
        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=colorscale, shading='gouraud') 
        plt.contour(xi, yi, zi.reshape(xi.shape), levels=[1.0, 2.0, 5.0], linewidth=10.0, colors='black')
        f_out = f_out + '_gas.png'

    elif ptype == 1:
        plt.title('DM')
        colorscale = 'inferno'
        plt.hexbin(data_x, data_y, gridsize=nbins, cmap=colorscale, bins='log')
        #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=colorscale) #, shading='gouraud')
        f_out = f_out + '_dm.png'

    '''
    elif ptype == 4:
        plt.title('')
        colorscale = 'hot'
        plt.scatter(data_x, data_y, s=0.15, color='black') 
        f_out = f_out + '_stars.png'
    '''

    # Save to file
    plt.tight_layout()
    plt.savefig(f_out)



def plot_rho(f_snap, center, side_size, f_out, nbins, f_rescale, thickn, units, n_files):
    print('Plotting density slices for snapshot: ', f_snap)

    # Plot properties
    ptsize_lv = 2.0
    col_lv = 'red'
    plot_col = 3
    plot_row = 1

    npt_lg_min = 100
    axis_margins = 1
    axis_size = 12

    if units == 'Mpc':
        axis_units = 'Mpc/h'; facMpc = 1.
    elif units == 'kpc':
        axis_units = 'Mpc/h'; facMpc = 1000.
        #axis_units = 'kpc/h'; facMpc = 1000.

    axis_label = []
    axis_label.append('SGX')
    axis_label.append('SGY')
    axis_label.append('SGZ');
    parts1 = readgadget(f_snap, 'pos', 1, n_files)
    parts2 = []
    parts3 = []
    parts4 = []

    #parts = [parts0, parts1, parts2]
    #print parts.max()
    #print parts.min()

    # Identify the particles belonging to the different objects
    i_type = 0

    x_plotlv = [[] for ix in range(0, 3)]
    y_plotlv = [[] for ix in range(0, 3)]

    minx = center[0] - side_size;   miny = center[1] - side_size;   minz = center[2] - side_size
    minima = [minx, miny, minz];    #print minima

    # Find slab of thickness +/- thickn around the axes
    for ix in range(0, 3):
        ixp1 = (ix+1) % 3
        ixp2 = (ix+2) % 3

        t1 = time.clock()
        (x_plot_tmp1, y_plot_tmp1) = find_slab(parts1, ix, center, minima, side_size, thickn, f_rescale * 512.0, units)
        (x_plot_tmp2, y_plot_tmp2) = find_slab(parts2, ix, center, minima, side_size, thickn, f_rescale * 64.0, units)
        n_tmp1 = len(x_plot_tmp1);              n_tmp2 = len(x_plot_tmp2)

        print('N Part1 in slab: ', n_tmp1)
        print('N Part2 in slab: ', n_tmp2)

        for ijk in range(0, n_tmp1):
            x_plotlv[ixp1].append(x_plot_tmp1[ijk])
            y_plotlv[ixp2].append(y_plot_tmp1[ijk])

        for ijk in range(0, n_tmp2):
            x_plotlv[ixp1].append(x_plot_tmp2[ijk])
            y_plotlv[ixp2].append(y_plot_tmp2[ijk])

        if (units == 'Mpc' and side_size > 4.0) or (units == 'kpc' and side_size > 4000.0):
            print('Selecting additional slabs')
            (x_plot_tmp3, y_plot_tmp3) = find_slab(parts3, ix, center, minima, side_size, thickn, f_rescale * 8.0, units)
            (x_plot_tmp4, y_plot_tmp4) = find_slab(parts4, ix, center, minima, side_size, thickn, f_rescale * 1.0, units)
            n_tmp3 = len(x_plot_tmp3);              n_tmp4 = len(x_plot_tmp4)

            for ijk in range(0, n_tmp3):
                x_plotlv[ixp1].append(x_plot_tmp3[ijk])
                y_plotlv[ixp2].append(y_plot_tmp3[ijk])

            for ijk in range(0, n_tmp4):
                x_plotlv[ixp1].append(x_plot_tmp4[ijk])
                y_plotlv[ixp2].append(y_plot_tmp4[ijk])
        else:
            print('Zoom mode')
            x_plot_tmp3 = []; x_plot_tmp4 = [];
            n_tmp3 = 0; n_tmp4 = 0;

        t2 = time.clock()

        print('Slab (%s, %s) with found in %.3f s.' % (axis_label[ixp1], axis_label[ixp2], (t2-t1)))
        print('Selected a total of ', n_tmp1 + n_tmp2 + n_tmp3 + n_tmp4, ' particles.')
        plt.ylabel(axis_label[ixp2]+' '+axis_units)

    # General plot settings
    plt.figure(figsize=(24,8))
    plt.rc('xtick', labelsize=axis_size)
    plt.rc('ytick', labelsize=axis_size)
    plt.rc('axes',  labelsize=axis_size)
    plt.margins(axis_margins)

    #fig, axes = plt.subplots(ncols=6, nrows=1, figsize=(21, 5))

    for ix in range(0, 3):
        ixp1 = (ix+1) % 3
        ixp2 = (ix+2) % 3

        x_min = -side_size; x_max = side_size
        y_min = -side_size; y_max = side_size

        # These plots are in Mpc/h not kpc/h
        x_min /= facMpc
        x_max /= facMpc
        y_min /= facMpc
        y_max /= facMpc

        print('XMin: ', x_min, ' XMax: ', x_max)

        # Plot settings for each subplot
        plt.subplot(plot_row, plot_col, ix+1)
        plt.axis([x_min, x_max, y_min, y_max])
        plt.xlabel(axis_label[ixp1]+' '+axis_units)
        plt.ylabel(axis_label[ixp2]+' '+axis_units)

        this_x = x_plotlv[ixp1][:]
        this_y = y_plotlv[ixp2][:]

        n_x = len(x_plotlv[ixp1])
        data_xy = np.zeros((2, n_x), dtype='float')

        # Convert units to Mpc
        for ip in range(0, len(this_x)):
            data_xy[0, ip] = (this_x[ip] - center[ixp1])/facMpc
            data_xy[1, ip] = (this_y[ip] - center[ixp2])/facMpc

        #colorscale = 'inferno'
        colorscale = 'rainbow'
        #(counts, xbins, ybins) = np.histogram2d(data_xy[0, :], data_xy[1, :], bins=nbins)
        #(counts, xbins, ybins, image) = plt.hist2d(data_xy[0, :], data_xy[1, :], bins=nbins) #, cmap=plt.cm.BuGn_r)

        #print counts
        #print this_x

        #smoothed = gaussian_filter(counts, sigma=2)
        #print smoothed
        #plt.pcolormesh(xbins, ybins, smoothed, cmap=plt.cm.BuGn_r)
        #plt.pcolormesh(xbins, ybins, smoothed, norm=colors.LogNorm(vmin=smoothed.min(), vmax=smoothed.max()), cmap=plt.cm.viridis)
        #plt.pcolormesh(xbins, ybins, smoothed, norm=colors.LogNorm(vmin=smoothed.min(), vmax=smoothed.max()), cmap=plt.cm.rainbow)
        plt.hexbin(data_xy[0, :], data_xy[1, :], gridsize=nbins, cmap=colorscale, bins='log')
        #, bins=nbins) #, cmap=plt.cm.BuGn_r)

        '''
        print 'Estimating gaussian kernel... '
        k = kde.gaussian_kde(data_xy)
        xi, yi = np.mgrid[data_xy[0].min():data_xy[0].max():nbins*1j, data_xy[1].min():data_xy[1].max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        print 'Done.'

        # plot a density
        #axes[3].set_title('Calculate Gaussian KDE')
        #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.BuGn_r)
        #axes[3].pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.BuGn_r)

        plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)
        plt.contour(xi, yi, zi.reshape(xi.shape) )

        #axes[2].set_title('2D Histogram')
        #(counts, ybins, xbins, image) = hist2d(this_x, this_y, nbins)
        #(counts, xbins, ybins, image) = plt.hist2d(this_x, this_y) #, bins=nbins, cmap=plt.cm.BuGn_r)
        #(counts, ybins, xbins, image) = plt.hist2d(this_x, this_y, gridsize=nbins, bins='log', cmap=plt.cm.BuGn_r)

        xv = [];        yv = [];        zv = []

        # Set up a regular grid of interpolation points
        #xi, yi = np.linspace(x_min, x_max, nbins), np.linspace(y_min, y_max, n_bins)
        #xi, yi = np.meshgrid(xi, yi)

        # Interpolate; there's also method='cubic' for 2-D data such as here
        #zi = scipy.interpolate.griddata((x, y), rho, (xi, yi), method='linear')

        for ib in range(0, len(counts)):
                new_x = 0.5 * (xbins[ib] + xbins[ib+1])
                new_y = 0.5 * (ybins[ib] + ybins[ib+1])
                new_z = counts[ib]

                if (new_x < x_max) and (new_x > x_min) and (new_y < y_max) and (new_y > y_min):
                        xv.append(new_x)
                        yv.append(new_y)
                        zv.append(new_z)

        #print len(counts)
        #print xv, yv
        print len(xv)
        print len(yv)
        print len(zv)

        plt.imshow(zv, extent=[x_min, x_max, y_min, y_max])

        '''

        #locations = [0.1, 0.5, 1.0, 5.0, 10.0, 15.0, 20.0]
        #plt.contour(xv, yv, np.transpose(zv), 5, fontsize=5, colors='black')#, manual=locations)
        #plt.contour(xv, yv, zv, 5, fontsize=5, colors='black')

        # plot a density
        #plt.pcolormesh(xv, yv, zv, cmap=plt.cm.BuGn_r)
        #print 'Plot edges: %.3f, %.3f, %.3f, %.3f\n' % (x_min, x_max, y_min, y_max)

    # Save to file
    plt.tight_layout()
    plt.savefig(f_out)


def plot_lglv(f_snap, h_ahf, f_out, lg0, lg1, x_virgo, reduce_fac, n_types):
    facMpc = 1000.
    buffPlot = 4000.        # Extra buffer on the edges
    thickn = 2500.
    r_virgo = 1500.

    # This samples 8 times less particles in the high - res region
    reduce_factors = [0] * n_types

    for i_red in range(0, n_types):
        reduce_factors[i_red] = pow(i_red+1,3) * reduce_fac

    print('Plotting LG & LV slices for snapshot: ', f_snap)
#       print reduce_factors

    # Plot properties
    ptsize_lg = 150.0
    ptsize_lv = 10.0
    col_mw="red"
    col_m31="orange"
    col_lv="black"
    col_lg="lime"

    plot_col = 3
    plot_row = 2

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
    x_ptsvirgo = [[] for ix in range(0, 3)]

    partDM = []

    # Read the particles
    for i_type in range (0, n_types):
        parts = readsnap(f_snap, 'pos', i_type+1)
        partDM.append(parts)

    # Identify the particles belonging to the different objects
    n_ptlg0 = 0
    n_ptlg1 = 0
    i_type = 0

    while (n_ptlg0 < npt_lg_min):

        t1 = time.clock()
        if n_ptlg0 < npt_lg_min:
            (x_ptslg0[0], x_ptslg0[1], x_ptslg0[2]) = select_particles_halo(lg0, partDM[i_type])

        if n_ptlg1 < npt_lg_min:
            (x_ptslg1[0], x_ptslg1[1], x_ptslg1[2]) = select_particles_halo(lg1, partDM[i_type])
        t2 = time.clock()

        n_ptlg0 = len(x_ptslg0[0])
        n_ptlg1 = len(x_ptslg1[0])

        if n_ptlg0 > npt_lg_min:
            print('M31 type:%d particles found in %.3f s. ' % (i_type, t2-t1))
        else:
            print('M31 type:%d particles not found after %.3f s, looking for particle types: %d ' % (t2-t1, i_type, i_type))

        if n_ptlg1 > npt_lg_min:
            print('MW  type:%d particles found in %.3f s. ' % (i_type, t2-t1))
        else:
            print('MW  type:%d particles not found after %.3f s, looking for particle types: %d ' % (t2-t1, i_type, i_type))

        if (n_ptlg0 < npt_lg_min) or (n_ptlg0 < npt_lg_min):
            i_type += 1

#       t1 = time.clock()
#       (x_ptsvirgo[0], x_ptsvirgo[1], x_ptsvirgo[2]) = select_particles(x_virgo, partDM2, r_virgo)
#       t2 = time.clock()
#       print 'Virgo particles found in %.3f s. ' % (t2-t1)

    lgcom = center_of_mass([lg0.m, lg1.m], [lg0.x, lg1.x])
    d_lgvirgo = distance(lgcom, x_virgo)
    d_mwm31 = distance(lg0.x, lg1.x)
    print('Distance LG to Virgo: %.3f ' % d_lgvirgo)

    # Center of plot
    c_plot = center_of_mass([1.0, 1.0], [lgcom, x_virgo])
    side_lv = d_lgvirgo + buffPlot * 2.0
    side_lg = d_mwm31 * 2
    print('Plot center: ', c_plot)
    print('Plot side  : ', side_lg)

    # This vector contains the minima of the XYZ coordinates for the plot
    min_lv_xyz = [0.0] * 3
    min_lg_xyz = [0.0] * 3

    # Plot is centered on HALF of the distance between Virgo and LG, then a buffer is added
    for ix in range(0, 3):
        min_lv_xyz[ix] = (c_plot[ix] - d_lgvirgo * 0.5 - buffPlot)
        min_lg_xyz[ix] = (lgcom[ix] - d_mwm31)
    #       print ix, min_lg_xyz[ix]

    # These contain the particles in each slab along the XY, XZ, YZ plane for the different resolution levels
    x_plotlg = [[] for ix in range(0, 3)]
    y_plotlg = [[] for ix in range(0, 3)]
    x_plotlv = [[] for ix in range(0, 3)]
    y_plotlv = [[] for ix in range(0, 3)]

    # Find slab of thickness +/- thickn around the axes
    for ix in range(0, 3):
        ixp1 = (ix+1) % 3
        ixp2 = (ix+2) % 3

        t1 = time.clock()
        # Lowest res particles - very numerous. Ignore the intermediate layers!
        (x_plotlv[ixp1], y_plotlv[ixp2]) = find_slab(partDM[n_types-1], ix, lgcom, min_lg_xyz, side_lv, thickn, reduce_factors[0])

        # High res particles
        (x_plotlg[ixp1], y_plotlg[ixp2]) = find_slab(partDM[0], ix, lgcom, min_lg_xyz, side_lv, thickn, reduce_factors[n_types-1])
        t2 = time.clock()

        print('Slab (%s, %s) found in %.3f s.' % (axis_label[ixp1], axis_label[ixp2], (t2-t1)))
        plt.ylabel(axis_label[ixp2]+' '+axis_units)


    # General plot settings
    plt.figure(figsize=(60,40))
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

        # Background high-res particles
        plt.scatter(x_plotlg[ixp1][:], y_plotlg[ixp2][:], s=ptsize_lg, c=col_lg)

        # Actual plot
        plt.scatter(x_ptslg0[ixp1][:], x_ptslg0[ixp2][:], s=ptsize_lg, c=col_mw)
        plt.scatter(x_ptslg1[ixp1][:], x_ptslg1[ixp2][:], s=ptsize_lg, c=col_m31)

    for ix in range(0, 3):
        ixp1 = (ix+1) % 3
        ixp2 = (ix+2) % 3

        x_min = (min_lv_xyz[ixp1])
        x_max = (x_min + side_lv)
        y_min = (min_lv_xyz[ixp2])
        y_max = (y_min + side_lv)

        x_min /= facMpc
        x_max /= facMpc
        y_min /= facMpc
        y_max /= facMpc

        print('Plot edges: %.3f, %.3f, %.3f, %.3f\n' % (x_min, x_max, y_min, y_max))

        plt.subplot(plot_row, plot_col, 3+ix+1)
        plt.axis([x_min, x_max, y_min, y_max])
        plt.scatter(x_plotlg[ixp1][:], y_plotlg[ixp2][:], s=ptsize_lv, c=col_lg)
        plt.scatter(x_plotlv[ixp1][:], y_plotlv[ixp2][:], s=ptsize_lv, c=col_lv)

        # Re-plot the local groups on the larger picture
        plt.scatter(x_ptslg0[ixp1][:], x_ptslg0[ixp2][:], s=ptsize_lg, c=col_mw)
        plt.scatter(x_ptslg1[ixp1][:], x_ptslg1[ixp2][:], s=ptsize_lg, c=col_m31)

        plt.xlabel(axis_label[ixp1]+' '+axis_units)
        plt.ylabel(axis_label[ixp2]+' '+axis_units)

    # Save to file
    plt.tight_layout()
    plt.savefig(f_out)



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

