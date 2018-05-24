#!/usr/bin/python

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
this_simu = 2
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

all_m_bins = save_path + 'm_bins.pkl'
all_n_bins = save_path + 'n_bins.pkl'

min_part = 30
stepMyr = 0.25

'''
	SELECT WHAT KIND OF ANALYSIS NEEDS TO BE DONE
'''

do_evolution = True
#do_evolution = False

# General combined statistics of all LG realisations
do_all_lgs = False
#do_all_lgs = True

do_plots_only = False
#do_plots_only = True

#do_trajectories = True
do_trajectories = False

#do_subs = True
do_subs = False

simurun = simuruns[this_simu]

lg_names = ['MW', 'M31']
n_lg = len(lg_names)

mass_histories = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init))
anisotropies = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init, 3))

all_lgs = []

time = np.zeros((snap_end-snap_init))

for i_time in range(0, snap_end-snap_init):
	time[i_time] = (snap_end - i_time) * stepMyr

# This skips the following loop on the halo evolution
if do_evolution == False:
	next_sub_init = sub_init
	next_sub_end = sub_end
	sub_init = 0
	sub_end = 0

for i_sub in range(sub_init, sub_end):

	subrun = '%02d' % i_sub

	s_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_sats.pkl'
	m_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_mains.pkl'
	
	# Try to load the pre-saved pickle format binary output
	try:
		hand_main = open(m_fname, 'r')
		hand_sats = open(s_fname, 'r')
		main = pickle.load(hand_main)
		sats = pickle.load(hand_sats)

		n_main = len(main)
		r_sub = 1.5
		n_lg = 2
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
		#plot_mass_accretion(time, this_mt, out_mah)
		mass_histories[i_lg, i_sub] = this_mt

		# For each LG member we are appending subhalo positions and histories
		subs = []
		subs_xt = []
		subs_yt = []
		subs_zt = []

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
				n_cross = len(acc_time)
	
				if n_cross > 0:
					first_cross = acc_time[n_cross-1] / 1000.
					#print 'Rvir crossing times: ', this_sub_z.accretion_time()
					print 'Subhalo %d in host %s has m(z=0): %e and m_max: %e, crossed Rvir %d times first at %3.3f. D(z=0) = %3.3f' % \
						(i_main, lg_names[i_lg], this_sub_z.halo[0].m, this_sub_z.m_max(), n_cross, first_cross, dist_rvir)

					subs.append(this_sub_z)
					subs_xt.append(this_xt[0, :])
					subs_yt.append(this_xt[1, :])
					subs_zt.append(this_xt[2, :])

		
		if do_trajectories == True:
			file_trajectory = outp_path + 'trajectories_lg' + this_run + '_' + simurun + '_' + subrun + '.png'
			plot_trajectory(subs_xt, subs_yt, 'x', 'y', file_trajectory)

		if do_subs == True:
			n_sat = len(subs)
			for i_sat in range(0, n_sat):	
				m0 = subs[i_sat].halo[0].m
				m1 = subs[i_sat].m_max()
				print '%d) Subhalo %d has m(z=0): %e and m_max: %e ' % (i_sub, i_sat, m0, m1)
				#print m1 #m1

		#Planes of satellites in substructure
		for i_snap in range(0, snap_end-snap_init):

			# Only computes the inertia tensor considering satellites above N=min_part particles
			if sats[i_lg][i_snap].n_sub > 5:
				(evals, red_evals, evecs, red_evecs) = sats[i_lg][i_snap].anisotropy('part', min_part)			
				anisotropies[i_lg, i_sub, i_snap] = evals
				#anisotropies[i_lg, i_sub, i_snap] = red_evals
	
if do_evolution == True:
	for i_lg in range(0, 2):
		out_fname = 'output/anisotropy_' + simurun + '_' + lg_names[i_lg] + '_' + subrun + '_' + this_run + '.png'
		plot_anisotropies(anisotropies, i_lg, sub_end, snap_end, out_fname)

	print 'Plotting mass accretions...'
	out_fname = 'output/mah_' + simurun + '_' + lg_names[0] + '.png'
	plot_mass_accretions(time, mass_histories[0, :, :], out_fname)
	out_fname = 'output/mah_' + simurun + '_' + lg_names[1] + '.png'
	plot_mass_accretions(time, mass_histories[1, :, :], out_fname)

'''
	This section computes some GLOBAL statistics, i.e. taking into account all of the LG simulations
'''

if do_all_lgs == True:

#	n_bins = []
#	m_bins = []

	n_bins = np.zeros((simu_end-simu_init, 2, 3))
	m_bins = np.zeros((simu_end-simu_init, 2, 3))

	for i_simu in range(simu_init, simu_end):
		simurun = simuruns[i_simu]
		these_lgs = []

		for i_sub in range(next_sub_init, next_sub_end):
	
			subrun = '%02d' % i_sub
			this_lg = LocalGroup(simurun)

		#	s_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_sats.pkl'
			m_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_mains.pkl'
	
			# Try to load the pre-saved pickle format binary output
			try:
			#	print 'Loading ', m_fname
				hand_main = open(m_fname, 'r')
				main = pickle.load(hand_main)
				this_lg.init_halos(main[0].halo[0], main[1].halo[0])
				these_lgs.append(this_lg)
			#	print this_lg.info()
			#	hand_sats = open(s_fname, 'r')
			#	sats = pickle.load(hand_sats)
				n_main = len(main)
				r_sub = 1.35
				n_lg = 2
			except:	
				n_lg = 0

		(m_bin, n_bin) = bin_lg_sub(these_lgs)
		
		n_bins[i_simu, :, :] = n_bin
		m_bins[i_simu, :, :] = m_bin

	file_m = open(all_m_bins, 'w')
	file_n = open(all_n_bins, 'w')

	pickle.dump(m_bins, file_m)
	pickle.dump(n_bins, file_n)

	print n_bins
	print m_bins

if do_plots_only == True:

	file_m = open(all_m_bins, 'r')
	file_n = open(all_n_bins, 'r')
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

	
	#for ibin in range(0, tot_bin):
	#	print m_bins[ibin, :, :]
	#	print n_bins[ibin, :, :]
	#for ibin in range(0, tot_bin):
	#	print '-'
	#	print var_m[ibin, 0, :]
	#	print var_n[ibin, 0, :]
	#print var_n[:, :, :]
	
	#print np.median(var_n)	
	#print np.median(var_n)	
	#print np.median(var_m)	

	f_out = 'output/sat_n_bins.png'
	plot_lg_bins(m_bins, n_bins, f_out)



# Add some V-Web stuff
#if 




