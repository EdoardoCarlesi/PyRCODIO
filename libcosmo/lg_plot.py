import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from scipy.ndimage.filters import gaussian_filter

from matplotlib import rc
import time
import numpy as np
import math
from libcosmo.utils import *
from pygadgetreader import *
from particles import *
from scipy.stats import kde



def plot_rho(f_snap, center, side_size, f_out, nbins, f_rescale, thickn, units):
	print 'Plotting density slices for snapshot: ', f_snap

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
	axis_label.append('SGZ')

	# Read the particles
	parts1 = readsnap(f_snap, 'pos', 1)
	#parts2 = readsnap(f_snap, 'pos', 2)
	#parts3 = readsnap(f_snap, 'pos', 3)
	#parts4 = readsnap(f_snap, 'pos', 4)

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

	minx = center[0] - side_size;	miny = center[1] - side_size;	minz = center[2] - side_size
	minima = [minx, miny, minz];	print minima

	# Find slab of thickness +/- thickn around the axes
	for ix in range(0, 3):
		ixp1 = (ix+1) % 3
		ixp2 = (ix+2) % 3

		t1 = time.clock()
		(x_plot_tmp1, y_plot_tmp1) = find_slab(parts1, ix, center, minima, side_size, thickn, f_rescale * 512.0, units) 
		(x_plot_tmp2, y_plot_tmp2) = find_slab(parts2, ix, center, minima, side_size, thickn, f_rescale * 64.0, units) 
		n_tmp1 = len(x_plot_tmp1); 		n_tmp2 = len(x_plot_tmp2)

		print 'N Part1 in slab: ', n_tmp1
		print 'N Part2 in slab: ', n_tmp2

		for ijk in range(0, n_tmp1):
			x_plotlv[ixp1].append(x_plot_tmp1[ijk])
			y_plotlv[ixp2].append(y_plot_tmp1[ijk])

		for ijk in range(0, n_tmp2):
			x_plotlv[ixp1].append(x_plot_tmp2[ijk])
			y_plotlv[ixp2].append(y_plot_tmp2[ijk])

		if (units == 'Mpc' and side_size > 4.0) or (units == 'kpc' and side_size > 4000.0):
			print 'Selecting additional slabs'
			(x_plot_tmp3, y_plot_tmp3) = find_slab(parts3, ix, center, minima, side_size, thickn, f_rescale * 8.0, units) 
			(x_plot_tmp4, y_plot_tmp4) = find_slab(parts4, ix, center, minima, side_size, thickn, f_rescale * 1.0, units) 
			n_tmp3 = len(x_plot_tmp3); 		n_tmp4 = len(x_plot_tmp4)
	
			for ijk in range(0, n_tmp3):
				x_plotlv[ixp1].append(x_plot_tmp3[ijk])
				y_plotlv[ixp2].append(y_plot_tmp3[ijk])

			for ijk in range(0, n_tmp4):
				x_plotlv[ixp1].append(x_plot_tmp4[ijk])
				y_plotlv[ixp2].append(y_plot_tmp4[ijk])
		else:
			print 'Zoom mode'
			x_plot_tmp3 = []; x_plot_tmp4 = []; 
			n_tmp3 = 0; n_tmp4 = 0;

		t2 = time.clock()
		
		print 'Slab (%s, %s) with found in %.3f s.' % (axis_label[ixp1], axis_label[ixp2], (t2-t1))
		print 'Selected a total of ', n_tmp1 + n_tmp2 + n_tmp3 + n_tmp4, ' particles.'
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

		print 'XMin: ', x_min, ' XMax: ', x_max

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
	
		colorscale = 'inferno'
		#colorscale = 'rainbow'
		#(counts, xbins, ybins) = np.histogram2d(data_xy[0, :], data_xy[1, :], bins=nbins)
		#(counts, xbins, ybins, image) = plt.hist2d(data_xy[0, :], data_xy[1, :], bins=nbins) #, cmap=plt.cm.BuGn_r)
		
		#print counts
		#print this_x

		#smoothed = gaussian_filter(counts, sigma=2)
		#print smoothed
		#plt.pcolormesh(xbins, ybins, smoothed, cmap=plt.cm.BuGn_r)
		#plt.pcolormesh(xbins, ybins, smoothed, norm=colors.LogNorm(vmin=smoothed.min(), vmax=smoothed.max()), cmap=plt.cm.viridis)
		#plt.pcolormesh(xbins, ybins, smoothed, norm=colors.LogNorm(vmin=smoothed.min(), vmax=smoothed.max()), cmap=plt.cm.rainbow)
		plt.hexbin(data_xy[0, :], data_xy[1, :], gridsize=nbins, cmap=colorscale, bins='log') #, bins=nbins) #, cmap=plt.cm.BuGn_r)

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

		xv = []; 	yv = [];	zv = []

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
	buffPlot = 4000.	# Extra buffer on the edges
	thickn = 2500.
	r_virgo = 1500.

	# This samples 8 times less particles in the high - res region
	reduce_factors = [0] * n_types

	for i_red in range(0, n_types):
		reduce_factors[i_red] = pow(i_red+1,3) * reduce_fac

	print 'Plotting LG & LV slices for snapshot: ', f_snap
#	print reduce_factors

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
			print 'M31 type:%d particles found in %.3f s. ' % (i_type, t2-t1)
		else:
			print 'M31 type:%d particles not found after %.3f s, looking for particle types: %d ' % (t2-t1, i_type, i_type)
		
		if n_ptlg1 > npt_lg_min:
			print 'MW  type:%d particles found in %.3f s. ' % (i_type, t2-t1)
		else:
			print 'MW  type:%d particles not found after %.3f s, looking for particle types: %d ' % (t2-t1, i_type, i_type)

		if (n_ptlg0 < npt_lg_min) or (n_ptlg0 < npt_lg_min):
			i_type += 1

#	t1 = time.clock()
#	(x_ptsvirgo[0], x_ptsvirgo[1], x_ptsvirgo[2]) = select_particles(x_virgo, partDM2, r_virgo)
#	t2 = time.clock()
#	print 'Virgo particles found in %.3f s. ' % (t2-t1)

	lgcom = center_of_mass([lg0.m, lg1.m], [lg0.x, lg1.x])
	d_lgvirgo = distance(lgcom, x_virgo)
	d_mwm31 = distance(lg0.x, lg1.x)
	print 'Distance LG to Virgo: %.3f ' % d_lgvirgo

	# Center of plot
	c_plot = center_of_mass([1.0, 1.0], [lgcom, x_virgo])
	side_lv = d_lgvirgo + buffPlot * 2.0
	side_lg = d_mwm31 * 2
	print 'Plot center: ', c_plot
	print 'Plot side  : ', side_lg

	# This vector contains the minima of the XYZ coordinates for the plot
	min_lv_xyz = [0.0] * 3 
	min_lg_xyz = [0.0] * 3 

	# Plot is centered on HALF of the distance between Virgo and LG, then a buffer is added 
	for ix in range(0, 3):
		min_lv_xyz[ix] = (c_plot[ix] - d_lgvirgo * 0.5 - buffPlot) 
		min_lg_xyz[ix] = (lgcom[ix] - d_mwm31) 
	#	print ix, min_lg_xyz[ix]

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

		print 'Slab (%s, %s) found in %.3f s.' % (axis_label[ixp1], axis_label[ixp2], (t2-t1))
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

		print 'Plot edges: %.3f, %.3f, %.3f, %.3f\n' % (x_min, x_max, y_min, y_max)

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
	buffPlot = 1.25 * lg0.r	# Extra buffer on the edges, largest Rvir of the two halos
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
	print 'Plot center: ', c_plot
	print 'Plot side  : ', side_lg

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

		print 'Slab (%s, %s) found in %.3f s.' % (axis_label[ixp1], axis_label[ixp2], (t2-t1))
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


	##################################################################################
	# One LG type only either MW or M31, per ONE realisation and N sub-realisations	 #				
	##################################################################################

def bin_lg_sub(lgs, n_lgs):
	n_tot = len(lgs)

	m_lgs = np.zeros((n_lgs, n_tot))
	n_sub = np.zeros((n_lgs, n_tot))

	bin_m = np.zeros((n_lgs, 5))
	bin_n = np.zeros((n_lgs, 5))
	
	perc_min0 = 5.
	perc_max0 = 100. - perc_min0
	perc_min1 = 32.
	perc_max1 = 100. - perc_min1

	for ilg in range(0, n_tot):
		m_lgs[0][ilg] = lgs[ilg].LG1.m
		m_lgs[1][ilg] = lgs[ilg].LG2.m		
		n_sub[0][ilg] = lgs[ilg].LG1.nsub
		n_sub[1][ilg] = lgs[ilg].LG2.nsub

	for ilg in range(0, n_lgs):
		bin_m[ilg][0] = np.median(m_lgs[ilg][:])
		bin_m[ilg][1] = np.percentile(m_lgs[ilg][:], perc_min0)
		bin_m[ilg][2] = np.percentile(m_lgs[ilg][:], perc_max0)
		bin_m[ilg][3] = np.percentile(m_lgs[ilg][:], perc_min1)
		bin_m[ilg][4] = np.percentile(m_lgs[ilg][:], perc_max1)

		bin_n[ilg][0] = np.median(n_sub[ilg][:])
		bin_n[ilg][1] = np.percentile(n_sub[ilg][:], perc_min0)
		bin_n[ilg][2] = np.percentile(n_sub[ilg][:], perc_max0)
		bin_n[ilg][3] = np.percentile(n_sub[ilg][:], perc_min1)
		bin_n[ilg][4] = np.percentile(n_sub[ilg][:], perc_max1)

	return [bin_m, bin_n]



def plot_lg_bins(x_bins, y_bins, f_out):
	size_p = 5
	size_x = 10
	size_y = 5

	x_label = 'Nsub'
	y_label = 'Msun/h'
	title_size = 10
	axis_size = 10
	axis_margins = 0.5
	tickSize = 0.4

	plt.rc({'text.usetex': True})
	#plt.rcParams['text.usetex'] = True
	#plt.rcParams['text.latex.unicode'] = True

	#normx = 1.e+12
	normx = 1.0000
	normy = 1.0000

	x_min0 = 0.6e+12 / normx; x_max0 = 4.2e+12 / normx
	y_min0 = 1.0 / normy; 	y_max0 = 140. / normy
	x_min1 = 0.5e+12 / normx; x_max1 = 2.1e+12 / normx
	y_min1 = 1.0 / normy; 	y_max1 = 60. / normy

	#print 'X Axis min=%.3f  max=%.3f\n' % (x_min, x_max)
	#print 'Y Axis min=%.3f  max=%.3f\n' % (y_min, y_max)

	n_bins = len(x_bins)
	print 'Nbins = ', n_bins

	xs = np.zeros(n_bins)
	ys = np.zeros(n_bins)
	x_err = np.zeros((2, n_bins))
	y_err = np.zeros((2, n_bins))

	#(fig, axs) = plt.subplots(ncols=3, nrows=1, figsize=(12, 4))
	#plt.rc(labelsize=axis_size)    
	plt.rc('xtick', labelsize=axis_size)    
	plt.rc('ytick', labelsize=axis_size)    

	(fig, axs) = plt.subplots(ncols=2, nrows=1, figsize=(size_x, size_y)) #, sharey = True)

	#axs[0].yaxis.set_ticks_position('left')
	axs[0].set_xlabel('M / $10^{12}M_{\odot} $', fontsize=axis_size)
	axs[1].set_xlabel('M / $10^{12}M_{\odot} $', fontsize=axis_size)
	axs[0].set_ylabel('$N_{sub} (M > 3\\times 10^{8} M_{\odot})$', fontsize=axis_size) 
	axs[1].set_ylabel('$N_{sub} (M > 3\\times 10^{8} M_{\odot})$', fontsize=axis_size) 
	axs[0].axis([x_min0, x_max0, y_min0, y_max0])
	axs[1].axis([x_min1, x_max1, y_min1, y_max1])
	axs[0].set_title('M31'); 	axs[1].set_title('MW')
	#axs[0].set_xscale('log'); 	axs[1].set_xscale('log')	
	#axs[0].set_yscale('log'); 	axs[1].set_yscale('log')	

	#axs[0].xaxis.set_major_locator(plt.MaxNLocator(4))
	#axs[1].xaxis.set_major_locator(plt.MaxNLocator(4))
	#plt.margins(axis_margins)		
	#plt.setp(axs[0].get_xticklabels(), visible = True)
	#plt.setp(axs[0].get_yticklabels(), visible = True)
	#plt.setp(axs[1].get_xticklabels(), visible = True)

	# Do two loops - M31 and MW
	for ilg in range(0, 2):
		#plt.subplot(1, 2, ilg+1)		

		# Now loop for each halo on a 
		for ibin in range(0, n_bins):
			xs[ibin] = float(x_bins[ibin, ilg, 0])
			x_err[0][ibin] = float(xs[ibin] - x_bins[ibin][ilg][3]) / normx
			x_err[1][ibin] = float(x_bins[ibin][ilg][4] - xs[ibin]) / normx
			#x_err[0][ibin] = float(xs[ibin] - x_bins[ibin][ilg][1]) / normx
			#x_err[1][ibin] = float(x_bins[ibin][ilg][2] - xs[ibin]) / normx
			xs[ibin] /= normx
			ys[ibin] = float(y_bins[ibin][ilg][0])
			y_err[0][ibin] = float(ys[ibin] - y_bins[ibin][ilg][3]) / normy
			y_err[1][ibin] = float(y_bins[ibin][ilg][4] - ys[ibin]) / normy
			#y_err[0][ibin] = float(ys[ibin] - y_bins[ibin][ilg][1]) / normy
			#y_err[1][ibin] = float(y_bins[ibin][ilg][2] - ys[ibin]) / normy
			ys[ibin] /= normy

		axs[ilg].errorbar(xs, ys, yerr=[y_err[0], y_err[1]], xerr=[x_err[0], x_err[1]], markersize=size_p, fmt='o')
	
	#plt.xticks
	#axs[0].xaxis.set_major_locator(ticker.MultipleLocator(tickSize))
	#axs[1].xaxis.set_major_locator(ticker.MultipleLocator(tickSize))
	#axs[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%7.2f"))
	#axs[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%7.2f"))
	#axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

	#axs[0].set_title(r'MW')#, size=title_size)
	#axs[1].set_title(r'M31')#, size=title_size)
	
	print 'Saving png to file: ', f_out
	plt.savefig(f_out)
	plt.close()
	plt.clf()
	plt.cla()


def plot_anisotropies(anisotropies, i_main, n_sub, n_snap, f_out):
	size_x = 20
	size_y = 20
	lnw = 1.0
	col = 'b'
	
	x_min = 0; 	x_max = 54 * 0.25 # GYrs
	y_min = 0; 	y_max = 1.0

	print 'Plotting anisotropies to file: ', f_out

	(fig, axs) = plt.subplots(ncols=3, nrows=1, figsize=(12, 4))
	
	for iax in range(0, 3):
		axs[iax].axis([x_min, x_max, y_min, y_max])

	x_m = np.zeros((n_snap))
	y_n = np.zeros((n_snap))
	z_n = np.zeros((n_snap))
	w_n = np.zeros((n_snap))

	e1 = np.zeros((n_snap))
	e2 = np.zeros((n_snap))
	e3 = np.zeros((n_snap))

	for i_sub in range(0, n_sub):
		
		e1 = anisotropies[i_main, i_sub, :, 0]
		e2 = anisotropies[i_main, i_sub, :, 1]
		e3 = anisotropies[i_main, i_sub, :, 2]

		for i_snap in range(0, n_snap):
			x_m[i_snap] = (n_snap - i_snap) * 0.25	# FIXME this assumes a fixed 250 myr timestep
			y_n[i_snap] = e1[i_snap] / e3[i_snap]
			z_n[i_snap] = e2[i_snap] / e3[i_snap]
			w_n[i_snap] = (e2[i_snap]-e1[i_snap]) / e3[i_snap]

		for i_y in range(0, len(e1)):
			if math.isnan(y_n[i_y]):
				y_n[i_y] = 0.0

			if math.isnan(z_n[i_y]):
				z_n[i_y] = 0.0

			if math.isnan(w_n[i_y]):
				w_n[i_y] = 0.0

		# TODO plot median in shaded region!!!!!!
		axs[0].set_title("a/c")
		axs[1].set_title("b/c")
		axs[2].set_title("(b - c) / a")

		axs[0].plot(x_m, y_n, linewidth=lnw, color=col)
		axs[1].plot(x_m, z_n, linewidth=lnw, color=col)
		axs[2].plot(x_m, w_n, linewidth=lnw, color=col)

	plt.savefig(f_out)
	plt.clf()
	plt.cla()
	plt.close()


def plot_massfunctions(x_m, y_m, n_mf, f_out, n_bins):
	size_x = 20
	size_y = 20
	lnw = 1.0
	col = 'b'
	axis_margins = 2	
#	print 'Plotting massfunctions to file: ', n_mf, f_out, y_max

	#n_bins = 15
	y_bins = [ [] for i in range(0, n_bins-1) ]
	#y_bins = np.zeros((3, n_bins))

	x_min = 1.e+15;	x_max = 1.e+7
	y_min = 10000.;	y_max = 1.0;

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

	#x_min = 5.e+8; 	x_max = 5.e+11
	#y_min = 1; 	y_max = 50 #max_list(y_n)

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
		#		y_bins[].append(0)				

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
			nmax0 = np.percentile(y_bins[km], 80)
			nmin0 = np.percentile(y_bins[km], 20)
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
	plt.xlabel('$log_{10}M_{\odot} h^{-1}$')

	#plt.ylabel('$N(>M)$')
	axs.set_yscale('log'); plt.ylabel('$log_{10}N(>M)$'); y_min=1.0

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
	#	for im in range(0, n_mf):
	#		axs.plot(x_m[im], y_n[im], linewidth=lnw, color=col)

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


def plot_trajectory(all_x, all_y, label_x, label_y, f_out):
	size = 20
	n_plots = len(all_x)

	axis_margins = 1
	axis_size = 50
	axis_units = 'ckpc/h'
	axis_label = []
	axis_label.append(label_x)
	axis_label.append(label_y)

	plt.figure(figsize=(size,size))
	plt.rc('xtick', labelsize=axis_size)    
	plt.rc('ytick', labelsize=axis_size)    
	plt.rc('axes',  labelsize=axis_size)    
	plt.margins(axis_margins)		

	x_min = min_list(all_x)
	x_max = max_list(all_x)
	y_min = min_list(all_y)
	y_max = max_list(all_y)

	plt.axis([x_min, x_max, y_min, y_max])
	plt.xlabel(axis_label[0]+' '+axis_units)
	plt.ylabel(axis_label[1]+' '+axis_units)

	for iplot in range(0, n_plots):
		#plt.lines.Line2D(all_x[iplot], all_y[iplot])
		plt.plot(all_x[iplot], all_y[iplot])

	plt.savefig(f_out)
	plt.clf()
	plt.gcf().show()

