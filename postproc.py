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

file_web = 'zoom_2048_054.000064.Vweb-ascii'; grid_web = 64; box_web = 100000.0

file_m_bins = save_path + 'm_bins.pkl'
file_n_bins = save_path + 'n_bins.pkl'
file_masses = save_path + 'masses.pkl'
file_n_subs = save_path + 'n_subs.pkl'

m_sub_min = 1.e+8
min_part = 30
stepMyr = 0.25

'''
	SELECT WHAT KIND OF ANALYSIS NEEDS TO BE DONE
'''
# Plot mass accretion histories and evolution of satellite anisotropy
do_evolution = True
#do_evolution = False

# General combined statistics of all LG realisations (gathers data on mass/nsub bins etc.)
do_all_lgs = True
#do_all_lgs = False

# Do plots of the mass - subhalo scatter for MW and M31
do_plots_mass_sub = True
#do_plots_mass_sub = False

# Subhalo trajectories
#do_trajectories = True
do_trajectories = False

# Plot sub halo mass functions and anisotropies
do_sub_plots = True
#do_sub_plots = False

# Plot mass accretion functions
do_plot_mfs = True
#do_plot_mfs = False

# Save informations about subhalos at different steps
#do_subs = True
do_subs = False

# Compute the cosmic web at each LG position
#do_vweb = True
do_vweb = False

# Use the pre-computed cosmic web to 
#do_subs_web = True
do_subs_web = False

lg_names = ['MW', 'M31']
n_lg = len(lg_names)

mass_histories = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init))
anisotropies = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init, 3))

all_lgs = []

n_sub_stat = 8		# Properties being tracked: mass, mass at infall, number of Rvir crossings etc.

n_dyn_stat = 6		# Velocity and position (histories!) in the MW / M31 reference frame

time = np.zeros((snap_end-snap_init))

for i_time in range(0, snap_end-snap_init):
	time[i_time] = (snap_end - i_time) * stepMyr

# This skips the following loop on the halo evolution
true_sub_init = sub_init
true_sub_end = sub_end

if do_evolution == False:
	sub_init = 0
	sub_end = 0

for this_simu in range(simu_init, simu_end):
	simurun = simuruns[this_simu]

	# Total mass functions of ALL subhaloes
	mfmw_z0_x = []; 		mfmw_z0_y = []
	mfmw_max_x = []; 		mfmw_max_y = []
	mfm31_z0_x = []; 		mfm31_z0_y = []
	mfm31_max_x = []; 		mfm31_max_y = []
	n_simu = 0
	sub_skip = 0
	
	for i_sub in range(sub_init, sub_end):
		subs_stats_m31 = []; true_subs_stats_mw = []

		sub_skip = 0
		subrun = '%02d' % i_sub

		s_fname='saved/'+simurun+'_'+subrun+'_sats.pkl'
		#m_fname='saved/'+simurun+'_'+subrun+'_mains_all.pkl'
		m_fname='saved/'+simurun+'_'+subrun+'_mains_all.pkl'
		print 'Loading file... ', s_fname
	
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
			masses_max = []; 	masses_z0 = [];		dyns_stats = []

			# This is a loop on all the other haloes stored at z=0
			for i_main in range(n_lg, n_main):
				this_main = main[i_main]
				this_x = this_main.halo[0].x
				this_d = this_lg.halo[0].distance(this_x) 
				dist_rvir = this_d / lg_r

				# This halo is within the range from the virial radius
				if this_d < lg_r * r_sub:
					this_xt = main[i_main].x_t_center(this_center)
					this_mt = main[i_main].m_t()
					n_steps = main[i_main].n_steps
			
					this_sub_z = SubHaloThroughZ(snap_end-snap_init)

					this_sub_z.host = this_lg
					this_sub_z.assign_halo_z(this_main)
					acc_time = this_sub_z.accretion_time()
					m_z0 = this_sub_z.halo[0].m
					n_cross = len(acc_time)
					#print i_main, m_z0/m_sub_min
	
					# This halo actually crossed the main halo's virial radius at some point
					if n_cross > 0:
						n_subs += 1
						first_cross = acc_time[n_cross-1] / 1000.
						(r_min, m_rmin, istep) = this_sub_z.r_min(this_center)
						m_max = this_sub_z.m_max()

						masses_max.append(m_max)
						masses_z0.append(m_z0)

						subs.append(this_sub_z)
						subs_xt.append(this_xt[0, :])
						subs_yt.append(this_xt[1, :])
						subs_zt.append(this_xt[2, :])

						# Check for the number of particle
						if m_max > m_sub_min:
			
						# Stat: (0 = mz0), (1=mmax), (2=n cross), (3=first cross), (4=dist rvir), (5=r_min/rvir), (6=m_rmin), (7=rvir)
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
						
							dyns_stats.append(this_xt)

			# End of the loop on i_main, still looping on i_lg
			(this_mfz0_x, this_mfz0_y) = mass_function(masses_z0)
			(this_mfmax_x, this_mfmax_y) = mass_function(masses_max)
	
			'''	
				Save informations on subhalo positions through time
			'''
			if i_lg == 0:
				mfmw_z0_x.append(this_mfz0_x); 		mfmw_z0_y.append(this_mfz0_y) 
				mfmw_max_x.append(this_mfmax_x); 	mfmw_max_y.append(this_mfmax_y) 
				dynmw_fname='saved/dyn_stats_'+lg_names[i_lg]+'_'+simurun+'_'+subrun+'.pkl'
				f_dyns_mw = open(dynmw_fname, 'w')
				print 'Storing informations on subhalo trajectories for MW to: ', dynmw_fname
				pickle.dump(dyns_stats, f_dyns_mw)

			elif i_lg == 1:
				mfm31_z0_x.append(this_mfz0_x); 	mfm31_z0_y.append(this_mfz0_y) 
				mfm31_max_x.append(this_mfmax_x); 	mfm31_max_y.append(this_mfmax_y) 
				dynm31_fname='saved/dyn_stats_'+lg_names[i_lg]+'_'+simurun+'_'+subrun+'.pkl'
				f_dyns_m31 = open(dynm31_fname, 'w')
				print 'Storing informations on subhalo trajectories for M31 to: ', dynm31_fname
				pickle.dump(dyns_stats, f_dyns_m31)

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

	# The loop on the subrun is over
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
	if i_sub == sub_end-1 and do_sub_plots == True:

		for i_lg in range(0, 2):
			out_fname = 'output/anisotropy_' + simurun + '_' + lg_names[i_lg] + '_' + subrun + '_' + this_run + '.png'
			plot_anisotropies(anisotropies, i_lg, sub_end-sub_skip, snap_end, out_fname)

		print 'Plotting mass accretions...'
		out_fname = 'output/mah_' + simurun + '_' + lg_names[0] + '.png'
		plot_mass_accretions(time, mass_histories[0, :, :], out_fname)
		out_fname = 'output/mah_' + simurun + '_' + lg_names[1] + '.png'
		plot_mass_accretions(time, mass_histories[1, :, :], out_fname)

'''
	This section computes some GLOBAL statistics, i.e. taking into account all of the LG simulations
'''
if do_all_lgs == True:
	print 'Do all lgs only.'
	n_bins = np.zeros((simu_end-simu_init, 2, 5))
	m_bins = np.zeros((simu_end-simu_init, 2, 5))
	all_masses = [[] for ijk in range(0, 2)] 
	all_n_subs = [[] for ijk in range(0, 2)] 

	#for simurun in simuruns: 
	for i_simu in range(simu_init, simu_end): 
		simurun = simuruns[i_simu]
		these_lgs = []
		these_masses = [[] for ijk in range(0, 2)] 
		these_n_subs = [[] for ijk in range(0, 2)] 

		print 'Loading .pkl files for ', simurun

		for i_sub in range(true_sub_init, true_sub_end):
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
if do_plots_mass_sub == True:
	file_m = open(file_m_bins, 'r')
	file_n = open(file_n_bins, 'r')
	m_bins = pickle.load(file_m)
	n_bins = pickle.load(file_n)

	tot_bin = 9
	
	var_m = np.zeros((tot_bin, 2, 2))
	var_n = np.zeros((tot_bin, 2, 2))
	
	for ibin in range(0, tot_bin):
		for ilg in range(0, 2):
			med_m = m_bins[ibin, ilg, 3]
			min_m = m_bins[ibin, ilg, 0]
			max_m = m_bins[ibin, ilg, 4]

			var_m[ibin, ilg, 0] = abs(med_m - min_m)/med_m
			var_m[ibin, ilg, 1] = abs(med_m - max_m)/med_m

			med_n = n_bins[ibin, ilg, 3]
			min_n = n_bins[ibin, ilg, 0]
			max_n = n_bins[ibin, ilg, 4]

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
		print 'MedianDeltaM= %.3f, MeanDeltaM= %.3f' % (np.median(delta_m[ijk]), np.mean(delta_m[ijk]))
		print 'MedianDeltaN= %.3f, MeanDeltaN= %.3f' % (np.median(delta_n[ijk]), np.mean(delta_n[ijk]))

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
		ANALYSIS OF THE COSMIC WEB
'''
# Output / input filenames for alignment angles et al.
angles_fname = 'saved/all_angles.pkl'

if do_vweb == True:
#	read_web_ascii = True
	read_web_ascii = False

	all_e3_angles = np.zeros((2, simu_end))
	all_evals = np.zeros((3, 3, simu_end))
	evals_evecs = np.zeros((4, 3))

	for i_simu in range(simu_init, simu_end):	
		main_run = simuruns[i_simu]
		this_e3_angle = []
		this_e3_deg = []
		this_evals = []

		#TODO look for the XX_XX modes
		for i_sub in range(0, 10):
			
			this_sub = '%02d' % i_sub
			lg_fname = 'saved/lg_' + main_run + '_' + this_sub + '.pkl'
			web_fname = 'saved/web_' + main_run + '_' + this_sub + '.pkl'

			try:
				lg_file = open(lg_fname, 'r')
				this_lg = pickle.load(lg_file)
				this_com = this_lg.get_com()
				search_web = True
			except:
			#	print 'LG data not available at this step.'
				search_web = False
		
			# Read the v-web from the original file
			if read_web_ascii == True and search_web == True:

				# This is triggered only if the v-web file is found
				if search_web == True:
					file_web='/home/eduardo/CLUES/DATA/2048/' + main_run + '/' + this_sub + '/zoom_2048_054.000064.Vweb-ascii'; 
					this_grid = read_vweb(file_web, grid_web, box_web)
					vel = this_grid.velPos(this_com)
					evals = this_grid.evalPos(this_com)
					evecs = this_grid.evecPos(this_com)
					evals_evecs[0, :] = evals
					evals_evecs[1:4, :] = evecs

					print 'Saving evals/evecs at LG position to: ', web_fname
					web_file = open(web_fname, 'w')
					pickle.dump(evals_evecs, web_file)

			# Not reading from an ascii file, loading pkl format 
			else:
				
				try:
					web_file = open(web_fname, 'r')
					evals_evecs = pickle.load(web_file)
					evals = evals_evecs[0, :]
					evecs = evals_evecs[1:4, :]
					search_web = True
				except:	
				#	print 'Webfile not found.'	
					search_web = False
			
			if search_web == True:
				virgo_x = [48000.0, 61000.0, 49000.0]
				to_virgo = vec_subt(virgo_x, this_com)
				virgo_n = vec_norm(to_virgo)
				
				this_evals.append(evals)
				e3 = evecs[2, :]; ang3 = abs(angle(virgo_n, e3)); 
				this_e3_angle.append(ang3)

		ang_med = np.median(this_e3_angle)
		ang_std = np.std(this_e3_angle)
		all_e3_angles[0, i_simu] = ang_med
		all_e3_angles[1, i_simu] = ang_std

		print 'Median: %.3f, Stddev: %.3f, Degrees= %.3f' % (ang_med, ang_std, np.arccos(ang_med) * rad2deg())

	f_angles = open(angles_fname, 'w')
	pickle.dump(all_e3_angles, f_angles)

	#print np.median(np.arccos(all_e3_angles[0, 0:6-8])) * rad2deg()

'''
	Analyze subhalo properties in the LSS environment
'''		
if do_subs_web == True:

	m_min = 1.e+8
	t_max = 10.0
	n_rand_runs = 500

	# Anisotropies above / below expectation
	n_e3_rej = 0
	n_e3_acc = 0

	# Which snapshot to select 
	step_num = 0
	this_step = '%02d' % step_num

	fout_anis_m31 = 'output/hist_percentiles_satellite_anisotropy_M31'+this_step+'.png'	
	fout_anis_mw = 'output/hist_percentiles_satellite_anisotropy_MW'+this_step+'.png'	
	anis_pkl = 'saved/percentiles_satellite_anisotropy_'+this_step+'_LG.pkl'	
	angles_pkl = 'saved/angles_align_satellite_anisotropy_'+this_step+'_LG.pkl'	

#	select_mass = True
	select_mass = False
	select_accr = True
#	select_accr = False

	lg_perc_anis = [[] for i in range(0, 2)]
	lg_angles_anis = [[] for i in range(0, 3)]

	for i_simu in range(simu_init, simu_end):	
		main_run = simuruns[i_simu]

		for i_sub in range(true_sub_init, true_sub_end):

			for i_lg in range(0, 2):		

				subrun = '%02d' % i_sub
				lg_fname = 'saved/lg_' + main_run + '_' + subrun + '.pkl'
				web_fname = 'saved/web_' + main_run + '_' + subrun + '.pkl'
	
				sub_fname='saved/sub_stats_'+lg_names[i_lg]+'_'+simurun+'_'+subrun+'_mains.pkl'
				dyn_fname='saved/dyn_stats_'+lg_names[i_lg]+'_'+simurun+'_'+subrun+'.pkl'

				try:
					print 'Loading pkl files at ', i_lg, i_sub, main_run
					f_lg = open(lg_fname, 'r')
					f_web = open(web_fname, 'r')
					f_sub = open(sub_fname, 'r')
					f_dyn = open(dyn_fname, 'r')

					# Stat: (0 = mz0), (1=mmax), (2=n cross), (3=first cross), (4=dist rvir), (5=r_min/rvir), (6=m_rmin), (7=rvir)
					sats = pickle.load(f_sub)
					dyns = pickle.load(f_dyn)
					vweb = pickle.load(f_web)

					nsats = len(sats)
					nselect = 0	

					positions = np.zeros((1,3))			
					masses = np.zeros((1))
					all_mass = np.zeros((nsats))
					all_pos	= np.zeros((nsats, 3))

					for isat in range(0, nsats):
						selected = False
						this_sat =  sats[isat]
						this_pos = dyns[isat]
						this_mz0 = this_sat[0]
						this_mmax = this_sat[1]
						this_tacc = this_sat[3]
						this_rvir = this_sat[7]
					
						all_mass[isat] = 1.0
						all_pos[isat, :] = this_pos[:, step_num]

						if select_mass and this_mmax > m_min:
							nselect += 1
							selected = True

						elif select_accr and this_tacc < t_max:
							# print 'selected ', this_tacc
							nselect += 1
							selected = True

						if selected == True:
							if nselect == 1:	
								positions[nselect-1, :] = this_pos[:, step_num]
								masses[nselect-1] = 1.0
							else:
								positions.resize((nselect, 3))
								positions[nselect-1, :] = this_pos[:, step_num]
								masses.resize((nselect))
								masses[nselect-1] = 1.0

					# Compute the moment of inertia of all the satellites or only of those accreted before tmax
					(sel_evals, sel_evecs) = moment_inertia(positions, masses)
					(all_evals, all_evecs) = moment_inertia(all_pos, all_mass)

					ev1 = vweb[1, :]
					ev3 = vweb[3, :]
					#(sel_rand_eval, sel_rand_disp, sel_rand_ang) = random_triaxialities_and_angles(nselect, n_rand_runs, ev3)
					(sel_rand_eval, sel_rand_disp, perc_eval) = random_triaxialities(nselect, n_rand_runs, sel_evals[0]/sel_evals[2])

					#print perc_eval
					lg_perc_anis[i_lg].append(perc_eval)

					new_positions = np.zeros((nselect, 3))

					#print sel_evecs

					for isat in range(0, nselect):
						old_pos = positions[isat, :]
						new_pos = change_basis(old_pos, all_evecs)
						new_positions[isat, :] = new_pos

					#print 'Dispersions, satellites vs. Rand:'
					#print 'I_web  I_IT  ev * evec  vweb  IM evals'
					for ix in range(1, 4):
						ev = vweb[ix, :]
						for iy in range(0, 1):
							this_angle = abs(angle(ev, sel_evecs[iy, :]))
							#print ix, iy, this_angle, vweb[0, ix-1], sel_evals[ix-1]/sel_evals[2]
							lg_angles_anis[ix-1].append(this_angle)		
					#print lg_angles_anis

						#print ix, np.std(abs(new_positions[:, ix])) / this_rvir, sel_rand_disp[ix, :]
					#print 'Main eigenvectors, satellites vs. Rand:'
					#for ix in range(0, 3):
					#print sel_evals[ix]/sel_evals[2], sel_rand_eval[ix, :]
					#this_ang = abs(angle(ev3, sel_evecs[2, :]))
					#rand_ang = abs(angle(ev3, sel_rand_evecs[2, :]))

					#if sel_evals[0]/sel_evals[2] < sel_rand_eval[0, 0]:
					#	n_e3_acc += 1
					#else:
					#	n_e3_rej += 1

					#print 'Vweb ev3 -smaller IM eigenvalues align: '
					#print this_ang, sel_rand_ang
					#	print np.arccos(this_ang) * 180.0 / math.pi, sel_rand_ang
				
					#for imh in range(0, 1):
					#	mi_sel = sel_evecs[imh, :]
					#	mi_all = all_evecs[imh, :]
						#print imh, iweb, angle(ev, mi_all)
						#for iweb in range(3, 4):	
							#ev = vweb[iweb, :]
							#print imh, iweb, angle(ev, mi_all), angle(ev, mi_sel)


				except:
					print 'Skipping step ', i_sub, ' / ', main_run	

# Plot Nsub differences histograms
#print 'Saving anisotropies to file ', anis_pkl
#f_anis_pkl = open(anis_pkl, 'w')
#pickle.dump(lg_perc_anis, f_anis_pkl)

	f_angles_pkl = open(angles_pkl, 'w')
	pickle.dump(lg_angles_anis, f_angles_pkl)

	n_bins = 20
	plt.rc({'text.usetex': True})
	plt.xlabel('$Percentile$')
	plt.ylabel('N')
	plt.title('MW')
	plt.hist(lg_perc_anis[0], n_bins)
	plt.savefig(fout_anis_mw)
	plt.clf()
	plt.cla()

	plt.rc({'text.usetex': True})
	plt.xlabel('$Percentile$')
	plt.ylabel('N')
	plt.title('M31')
	plt.hist(lg_perc_anis[1], n_bins)
	plt.savefig(fout_anis_m31)
	plt.clf()
	plt.cla()


#print lg_perc_anis[1]
