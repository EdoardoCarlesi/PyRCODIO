#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
#import plotly.plotly as py
import numpy as np
import os

from libio.read_ascii import *
from config import *
from libcosmo.utils import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.lg_plot import *
import pickle

resolution='2048'
simuruns = simu_runs()

# Which realisation
this_simu = 0
simu_init = 0
simu_end = 10

# Number of subrun
sub_init = 0
sub_end = 10

# Snapshots
snap_init = 0
snap_end = 54

base_path = '/home/eduardo/CLUES/PyRCODIO/'
outp_path = 'output/'
save_path = 'saved/'

file_web = 'zoom_2048_054.000064.Vweb-ascii'; grid_web = 64.0; box_web = 100.0

file_m_bins = save_path + 'm_bins.pkl'
file_n_bins = save_path + 'n_bins.pkl'
file_masses = save_path + 'masses.pkl'
file_n_subs = save_path + 'n_subs.pkl'

m_sub_min = 3.e+8
min_part = 30
stepMyr = 0.25

'''
	SELECT WHAT KIND OF ANALYSIS NEEDS TO BE DONE
'''
# Plot mass accretion histories and evolution of satellite anisotropy
#do_evolution = True
do_evolution = False

# General combined statistics of all LG realisations - needed to gather the data
do_all_lgs = False
#do_all_lgs = True

# Only do some post-post processing plots
#do_plots_only = False
do_plots_only = True

# Subhalo trajectories
#do_trajectories = True
do_trajectories = False

# Plot mass accretion functions
#do_plot_mfs = True
do_plot_mfs = False

do_subs = True
#do_subs = False

#do_vweb = False
do_vweb = True

#simurun = simuruns[this_simu]

lg_names = ['MW', 'M31']
n_lg = len(lg_names)

mass_histories = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init))
anisotropies = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init, 3))

all_lgs = []

# Total mass functions of ALL subhaloes
mfmw_z0_x = []; 		mfmw_z0_y = []
mfmw_max_x = []; 		mfmw_max_y = []
mfm31_z0_x = []; 		mfm31_z0_y = []
mfm31_max_x = []; 		mfm31_max_y = []
n_simu = 0
sub_skip = 0
	
n_sub_stat = 8		# Properties being tracked: mass, mass at infall, number of Rvir crossings etc.
subs_stats_m31 = []; subs_stats_mw = []

time = np.zeros((snap_end-snap_init))

for i_time in range(0, snap_end-snap_init):
	time[i_time] = (snap_end - i_time) * stepMyr

# This skips the following loop on the halo evolution
if do_evolution == False:
	next_sub_init = sub_init
	next_sub_end = sub_end
	sub_init = 0
	sub_end = 0

all_init =  simu_init * sub_init
all_end  = (simu_end - simu_init) * (sub_end - sub_init)

for i_all in range(all_init, all_end):

	this_simu = simu_init + int(i_all/sub_end)
	simurun = simuruns[this_simu]
	i_sub = sub_init + (i_all % sub_end)
	sub_skip = 0
	subrun = '%02d' % i_sub

	s_fname='saved/'+simurun+'_'+subrun+'_sats.pkl'
	m_fname='saved/'+simurun+'_'+subrun+'_mains.pkl'
	
	# Try to load the pre-saved pickle format binary output
	try:
		hand_main = open(m_fname, 'r')
		hand_sats = open(s_fname, 'r')
		main = pickle.load(hand_main)
		sats = pickle.load(hand_sats)

		n_main = len(main)
		r_sub = 1.5
		n_lg = 2
		n_simu += 1
	except:	
		n_lg = 0

	# The first two are the main LG members
	for i_lg in range(0, n_lg):
		this_lg = main[i_lg]
		this_center = this_lg.x_t()
		this_run = '%02d' % i_lg
		lg_r = this_lg.halo[0].r
		this_mt = this_lg.m_t()

		# Plot mass accretion history of the main halo
		out_mah = outp_path + lg_names[i_lg] + '_main_' + simurun + '_' + subrun + '_mah.png'

		try:
			mass_histories[i_lg, i_sub] = this_mt
		except:		
			sub_skip += 0.5
			print 'MT not found, skipping ', sub_skip

		# For each LG member we are appending subhalo positions and histories
		subs = []; 		subs_stats = [];	n_subs = 0
		subs_xt = []; 		subs_yt = []; 		subs_zt = []
		masses_max = []; 	masses_z0 = []

		# This is a loop on all the other haloes stored at z=0
		for i_main in range(n_lg, n_main):
			this_main = main[i_main]
			this_x = this_main.halo[0].x
			this_d = this_lg.halo[0].distance(this_x) 
			dist_rvir = this_d / lg_r

			if this_d < lg_r * r_sub:
				this_xt = main[i_main].x_t_center(this_center)
				this_mt = main[i_main].m_t()

				this_sub_z = SubHaloThroughZ(snap_end-snap_init)

				this_sub_z.host = this_lg
				this_sub_z.assign_halo_z(this_main)
				acc_time = this_sub_z.accretion_time()
				m_z0 = this_sub_z.halo[0].m
				n_cross = len(acc_time)
	
				if n_cross > 0:
					n_subs += 1
					first_cross = acc_time[n_cross-1] / 1000.
					(r_min, m_rmin, istep) = this_sub_z.r_min(this_center)
					m_max = this_sub_z.m_max()

				#	print 'Rvir crossing times: ', this_sub_z.accretion_time()
				#	print 'Subh %d in %s has m(z=0): %e m_max: %e, crossRvir %d, first %3.3f. D(z=0) = %3.3f, minD=%3.3f, m=%e' % \
				#		(i_main, lg_names[i_lg], m_z0, m_max, n_cross, first_cross, dist_rvir, r_min/lg_r, m_rmin)
					masses_max.append(m_max)
					masses_z0.append(m_z0)

					subs.append(this_sub_z)
					subs_xt.append(this_xt[0, :])
					subs_yt.append(this_xt[1, :])
					subs_zt.append(this_xt[2, :])

					if m_max > m_sub_min:
						sub_stat = np.zeros((n_sub_stat))
						sub_stat[0] = m_z0
						sub_stat[1] = m_max
						sub_stat[2] = float(n_cross)
						sub_stat[3] = first_cross
						sub_stat[4] = dist_rvir
						sub_stat[5] = r_min/lg_r
						sub_stat[6] = m_rmin
						sub_stat[7] = lg_r
			
						subs_stats.append(sub_stat)
	
		(this_mfz0_x, this_mfz0_y) = mass_function(masses_z0)
		(this_mfmax_x, this_mfmax_y) = mass_function(masses_max)
	
		if i_lg == 0:
			mfmw_z0_x.append(this_mfz0_x); 		mfmw_z0_y.append(this_mfz0_y) 
			mfmw_max_x.append(this_mfmax_x); 	mfmw_max_y.append(this_mfmax_y) 

		elif i_lg == 1:
			mfm31_z0_x.append(this_mfz0_x); 	mfm31_z0_y.append(this_mfz0_y) 
			mfm31_max_x.append(this_mfmax_x); 	mfm31_max_y.append(this_mfmax_y) 

		if do_trajectories == True:
			file_trajectory = outp_path + 'trajectories_lg' + this_run + '_' + simurun + '_' + subrun + '.png'
			plot_trajectory(subs_xt, subs_yt, 'x', 'y', file_trajectory)

		#Planes of satellites in substructure
		for i_snap in range(0, snap_end-snap_init):

			# Only computes the inertia tensor considering satellites above N=min_part particles
			if sats[i_lg][i_snap].n_sub > 5:
				try:
					(evals, red_evals, evecs, red_evecs) = sats[i_lg][i_snap].anisotropy('part', min_part)
					anisotropies[i_lg, i_sub, i_snap] = evals
				except:
					evals = (0, 0, 0)
		'''
			Save some informations on most massive subhalo statistics
		'''
		sub_fname='saved/sub_stats_'+lg_names[i_lg]+'_'+simurun+'_'+subrun+'_mains.pkl'

		print 'Saving subhalo statistics to file: ', sub_fname
		f_subs = open(sub_fname, 'w')
		pickle.dump(subs_stats, f_subs)

	# Plot all the subhalo mass function for a given LG run
	if do_plot_mfs == True and i_sub == sub_end-1:

		out_mwz0 = 'output/mf_' + simurun + '_' + lg_names[0] + '_mz0_subs.png'
		out_m31z0 = 'output/mf_' + simurun + '_' + lg_names[1] + '_mz0_subs.png'
		out_mwmax = 'output/mf_' + simurun + '_' + lg_names[0] + '_mmax_subs.png'
		out_m31max = 'output/mf_' + simurun + '_' + lg_names[1] + '_mmax_subs.png'
	
		plot_massfunctions(mfmw_z0_x, mfmw_z0_y, n_simu, out_mwz0)
		plot_massfunctions(mfm31_z0_x, mfm31_z0_y, n_simu, out_m31z0)
		plot_massfunctions(mfmw_max_x, mfmw_max_y, n_simu, out_mwmax)
		plot_massfunctions(mfm31_max_x, mfm31_max_y, n_simu, out_m31max)
	
		fout_mwz0 = 'saved/stats_' + simurun + '_' + lg_names[0] + '_mz0_subs.pkl'
		fout_m31z0 = 'saved/stats_' + simurun + '_' + lg_names[1] + '_mz0_subs.pkl'
		fout_mwmax = 'saved/stats_' + simurun + '_' + lg_names[0] + '_mmax_subs.pkl'
		fout_m31max = 'saved/stats_' + simurun + '_' + lg_names[1] + '_mmax_subs.pkl'

		sub_skip = int(sub_skip)

	# Plot mass accretion histories and evolution of satellite anisotropy
	if do_evolution == True and i_sub == sub_end-1:

		for i_lg in range(0, 2):
			out_fname = 'output/anisotropy_' + simurun + '_' + lg_names[i_lg] + '_' + subrun + '_' + this_run + '.png'
			plot_anisotropies(anisotropies, i_lg, sub_end-sub_skip, snap_end, out_fname)

		print 'Plotting mass accretions...'
		out_fname = 'output/mah_' + simurun + '_' + lg_names[0] + '.png'
		plot_mass_accretions(time, mass_histories[0, :, :], out_fname)
		out_fname = 'output/mah_' + simurun + '_' + lg_names[1] + '.png'
		plot_mass_accretions(time, mass_histories[1, :, :], out_fname)

	# We reset the lists
	if i_sub == sub_end - 1:
		mfmw_z0_x = []; 		mfmw_z0_y = []
		mfmw_max_x = []; 		mfmw_max_y = []
		mfm31_z0_x = []; 		mfm31_z0_y = []
		mfm31_max_x = []; 		mfm31_max_y = []
		subs = []; 		subs_stats = [];	n_subs = 0
		subs_xt = []; 		subs_yt = []; 		subs_zt = []
		masses_max = []; 	masses_z0 = []
		n_simu = 0
		print 'Cleaning loop...'

'''
	This section computes some GLOBAL statistics, i.e. taking into account all of the LG simulations
'''
if do_all_lgs == True:

	print 'Do all lgs only.'
	n_bins = np.zeros((simu_end-simu_init, 2, 3))
	m_bins = np.zeros((simu_end-simu_init, 2, 3))
	all_masses = [[] for ijk in range(0, 2)] 
	all_n_subs = [[] for ijk in range(0, 2)] 

	#for simurun in simuruns: 
	for i_simu in range(simu_init, simu_end): 
		simurun = simuruns[i_simu]
		these_lgs = []
		these_masses = [[] for ijk in range(0, 2)] 
		these_n_subs = [[] for ijk in range(0, 2)] 

		print 'Loading .pkl files for ', simurun

		for i_sub in range(next_sub_init, next_sub_end):
			this_lg = LocalGroup(simurun)
			subrun = '%02d' % i_sub

			s_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_sats.pkl'
			m_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_mains.pkl'
	
			if do_subs == True:
				sub_fname_mw='saved/sub_stats_'+lg_names[0]+'_'+simurun+'_'+subrun+'_mains.pkl'
				sub_fname_m31='saved/sub_stats_'+lg_names[1]+'_'+simurun+'_'+subrun+'_mains.pkl'

				try:
					f_subs_mw = open(sub_fname_mw, 'r')
					subs_stats_mw = pickle.load(f_subs_mw)
					f_subs_m31 = open(sub_fname_m31, 'r')
					subs_stats_m31 = pickle.load(f_subs_m31)
				except:
					print f_subs_mw
					print f_subs_m31

			# Try to load the pre-saved pickle format binary output
			try:
			#	print 'Loading ', m_fname
				hand_main = open(m_fname, 'r')
				main = pickle.load(hand_main)
				this_lg.init_halos(main[0].halo[0], main[1].halo[0])
				these_lgs.append(this_lg)

				for ijk in range(0, 2):
					these_masses[ijk].append(main[ijk].halo[0].m)
					these_n_subs[ijk].append(main[ijk].halo[0].nsub)

				n_main = len(main)
				r_sub = 1.35
				n_lg = 2
			except:	
				n_lg = 0

		for ijk in range(0, 2):
			all_masses[ijk].append(these_masses[ijk])
			all_n_subs[ijk].append(these_n_subs[ijk])

		(m_bin, n_bin) = bin_lg_sub(these_lgs)
		
		n_bins[i_simu, :, :] = n_bin
		m_bins[i_simu, :, :] = m_bin

	file_m = open(file_m_bins, 'w')
	file_n = open(file_n_bins, 'w')
	file_ms = open(file_masses, 'w')
	file_ns = open(file_n_subs, 'w')

	pickle.dump(m_bins, file_m)
	pickle.dump(n_bins, file_n)
	pickle.dump(all_masses, file_ms)
	pickle.dump(all_n_subs, file_ns)

'''
	PLOT THE THE MASS vs. SATELLITE FLUCTUATIONS 
'''
if do_plots_only == True:
	file_m = open(file_m_bins, 'r')
	file_n = open(file_n_bins, 'r')
	m_bins = pickle.load(file_m)
	n_bins = pickle.load(file_n)

	tot_bin = 9
	
	var_m = np.zeros((tot_bin, 2, 2))
	var_n = np.zeros((tot_bin, 2, 2))
	
	for ibin in range(0, tot_bin):
		for ilg in range(0, 2):
			med_m = m_bins[ibin, ilg, 1]
			min_m = m_bins[ibin, ilg, 0]
			max_m = m_bins[ibin, ilg, 2]

			var_m[ibin, ilg, 0] = abs(med_m - min_m)/med_m
			var_m[ibin, ilg, 1] = abs(med_m - max_m)/med_m

			med_n = n_bins[ibin, ilg, 1]
			min_n = n_bins[ibin, ilg, 0]
			max_n = n_bins[ibin, ilg, 2]

			var_n[ibin, ilg, 0] = abs(med_n - min_n)/med_n
			var_n[ibin, ilg, 1] = abs(med_n - max_n)/med_n

	f_out = 'output/sat_n_bins.png'
	plot_lg_bins(m_bins, n_bins, f_out)

	'''
		NOW PLOT THE HISTOGRAMS OF THE MASS & SATELLITE FLUCTUATIONS WRT THE MEDIAN
	'''
	file_ms = open(file_masses, 'r')
	file_ns = open(file_n_subs, 'r')
	masses = pickle.load(file_ms)
	n_subs = pickle.load(file_ns)

	delta_m = [[] for ijk in range(0, 2)] 
	delta_n = [[] for ijk in range(0, 2)] 

	n_runs = len(n_subs[0])


	for ijk in range(0, 2):
		for i_run in range(0, n_runs):
			#this_m_med = np.mean(masses[ijk][i_run])
			#this_n_med = np.mean((n_subs[ijk][i_run]))
			this_m_med = np.median(masses[ijk][i_run])
			this_n_med = np.median((n_subs[ijk][i_run]))
			n_lgs = len(masses[ijk][i_run])

			for ilg in range(0, n_lgs):
				#print i_run, n_lgs, this_m_med
				this_ms = abs(masses[ijk][i_run][ilg] - this_m_med) / this_m_med
				this_ns = abs(float(n_subs[ijk][i_run][ilg]) - float(this_n_med)) / this_n_med

				delta_m[ijk].append(this_ms)
				delta_n[ijk].append(this_ns)

	for ijk in range(0, 2):
		print np.median(delta_m[ijk]), np.mean(delta_m[ijk])
		print np.median(delta_n[ijk]), np.mean(delta_n[ijk])

		n_bins = 30

		# Plot Msub difference histograms
		(mny, mbins, mpat) = plt.hist(delta_m[ijk], n_bins)
		fout_m = 'output/hist_' + lg_names[ijk] + '_delta_m.png'
	        plt.rc({'text.usetex': True})
		plt.xlabel('$\Delta M$')
		plt.ylabel('N')
		plt.title(lg_names[ijk])
		plt.savefig(fout_m)
	        plt.clf()
        	plt.cla()
		plt.close()
	
		# Plot Nsub differences histograms
		(nny, nbins, npat) = plt.hist(delta_n[ijk], n_bins)
		fout_n = 'output/hist_' + lg_names[ijk] + '_delta_n.png'	
	        plt.rc({'text.usetex': True})
		plt.xlabel('$\Delta N_{sub}$')
		plt.ylabel('N')
		plt.title(lg_names[ijk])
		plt.savefig(fout_n)
	        plt.clf()
        	plt.cla()
	        
'''
	Analysis of the cosmic web
'''
if do_vweb == True:
		
	#TODO look for the XX_XX modes

	for i_sub in range(0, 10):
		this_sub = '%02d' % i_sub
		lg_fname = 'saved/lg_' + main_run + this_sub + '.pkl'

		try:
			this_lg = pickle.load(lg_fname)
			search_web = True
		else:
			print 'LG data not available at this step.'
			search_web = False
		
		if search_web == True:
			this_com = this_lg 




