import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib import rc
import time
import numpy as np
import math
from libcosmo.utils import *
from pygadgetreader import *
from particles import *



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
def bin_lg_sub(lgs):
	n_tot = len(lgs)

	m_lgs = np.zeros((2, n_tot))
	n_sub = np.zeros((2, n_tot))

	bin_m = np.zeros((2, 5))
	bin_n = np.zeros((2, 5))
	
	perc_min0 = 5.
	perc_max0 = 100. - perc_min0
	perc_min1 = 25.
	perc_max1 = 100. - perc_min1

	for ilg in range(0, n_tot):
		m_lgs[0][ilg] = lgs[ilg].LG1.m
		m_lgs[1][ilg] = lgs[ilg].LG2.m		
		n_sub[0][ilg] = lgs[ilg].LG1.nsub
		n_sub[1][ilg] = lgs[ilg].LG2.nsub

	for ilg in range(0, 2):
		#bin_m[ilg][0] = np.amin(m_lgs[ilg][:])
		bin_m[ilg][0] = np.percentile(m_lgs[ilg][:], perc_min0)
		bin_m[ilg][1] = np.median(m_lgs[ilg][:])
		bin_m[ilg][2] = np.percentile(m_lgs[ilg][:], perc_max0)
		bin_m[ilg][3] = np.percentile(m_lgs[ilg][:], perc_min1)
		bin_m[ilg][4] = np.percentile(m_lgs[ilg][:], perc_max1)
		bin_n[ilg][0] = np.percentile(n_sub[ilg][:], perc_min0)
		bin_n[ilg][1] = np.median(n_sub[ilg][:])
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
			xs[ibin] = float(x_bins[ibin, ilg, 1])
			x_err[0][ibin] = float(xs[ibin] - x_bins[ibin][ilg][0]) / normx
			x_err[1][ibin] = float(x_bins[ibin][ilg][2] - xs[ibin]) / normx
			xs[ibin] /= normx
			ys[ibin] = float(y_bins[ibin][ilg][1])
			y_err[0][ibin] = float(ys[ibin] - y_bins[ibin][ilg][0]) / normy
			y_err[1][ibin] = float(y_bins[ibin][ilg][2] - ys[ibin]) / normy
			ys[ibin] /= normy

		#	print xs[ibin], ys[ibin]

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


def plot_massfunctions(x_m, y_n, n_mf, f_out):
	size_x = 20
	size_y = 20
	lnw = 1.0
	col = 'b'
	axis_margins = 2	

	x_min = 5.e+8; 	x_max = 5.e+11
	y_min = 1; 	y_max = 50 #max_list(y_n)

	print 'Plotting massfunctions to file: ', n_mf, f_out, y_max

	(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

	plt.rc({'text.usetex': True})
	plt.margins(axis_margins)		
	plt.xlabel('$M_{\odot}$')
	plt.ylabel('N')
	axs.set_xscale('log')
	axs.set_yscale('log')
	axs.axis([x_min, x_max, y_min, y_max])
	
	if n_mf > 1:
		for im in range(0, n_mf):
			axs.plot(x_m[im], y_n[im], linewidth=lnw, color=col)

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
	y_min = np.min(mahs); y_max = np.max(mahs) * 1.85

	if y_min < 1.e+5:
		y_min = y_max / 200.

	(fig, axs) = plt.subplots(ncols=1, nrows=1, figsize=(4, 4))

	plt.yscale('log')
	axs.axis([x_min, x_max, y_min, y_max])

	for iplot in range(0, n_plots):
		mah = mahs[iplot]
		axs.plot(time, mah, linewidth=lnw, color=col)

	plt.savefig(f_out)
	plt.clf()
	plt.cla()
	plt.close(fig)
	plt.gcf().show()


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

