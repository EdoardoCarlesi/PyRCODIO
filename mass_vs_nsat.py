#!/usr/bin/python

import matplotlib.pyplot as plt
import scipy.stats as sp
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

#base_path = '/home/eduardo/CLUES/PyRCODIO/'
base_path = '/home/edoardo/CLUES/PyRCODIO/'
outp_path = 'output/'
save_path = 'saved/'

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
# General combined statistics of all LG realisations (gathers data on mass/nsub bins etc.)
#do_all_lgs = True
do_all_lgs = False

# Do plots of the mass - subhalo scatter for MW and M31
do_plots_mass_sub = True
#do_plots_mass_sub = False


'''
	START THE ACTUAL PROGRAM AFTER CHOOSING THE ANALYSIS OPTIONS
	
'''

lg_names = ['MW', 'M31', 'LG']
#n_lg = len(lg_names)
n_lg = 2

n_lgs = 2	# Number of masses to be tracked, mw, m31 and mw + m31

'''
	This section computes some GLOBAL statistics, i.e. taking into account all of the LG simulations
'''
if do_all_lgs == True:
	print 'Do all lgs only.'
	n_bins = np.zeros((simu_end-simu_init, n_lgs, 5))		# 5 bins = 1sigma, 2sigma and median value
	m_bins = np.zeros((simu_end-simu_init, n_lgs, 5))
	all_masses = [[] for ijk in range(0, n_lgs)] 
	all_n_subs = [[] for ijk in range(0, n_lgs)] 

	#for simurun in simuruns: 
	for i_simu in range(simu_init, simu_end): 
		simurun = simuruns[i_simu]
		these_lgs = []
		these_masses = [[] for ijk in range(0, n_lgs)] 
		these_n_subs = [[] for ijk in range(0, n_lgs)] 

		print 'Loading .pkl files for ', simurun

		for i_sub in range(sub_init, sub_end):
			this_lg = LocalGroup(simurun)
			subrun = '%02d' % i_sub

			s_fname='saved/'+simurun+'_'+subrun+'_sats.pkl'
			m_fname='saved/'+simurun+'_'+subrun+'_mains_all.pkl'
	
			'''
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
			'''

			# Try to load the pre-saved pickle format binary output
			try:
			#	print 'Loading ', m_fname
				hand_main = open(m_fname, 'r')
				main = pickle.load(hand_main)
				this_lg.init_halos(main[0].halo[0], main[1].halo[0])
				these_lgs.append(this_lg)

				for ijk in range(0, n_lgs):
					if ijk == 2:
						these_masses[ijk].append(main[0].halo[0].m + main[1].halo[0].m)
						these_n_subs[ijk].append(main[0].halo[0].nsub + main[1].halo[0].nsub)
					else:
						these_masses[ijk].append(main[ijk].halo[0].m)
						these_n_subs[ijk].append(main[ijk].halo[0].nsub)

				n_main = len(main)
				r_sub = 1.35
				n_lg = 2
			except:	
				n_lg = 0

		for ijk in range(0, n_lgs):
			all_masses[ijk].append(these_masses[ijk])
			all_n_subs[ijk].append(these_n_subs[ijk])

		(m_bin, n_bin) = bin_lg_sub(these_lgs, n_lgs)
		
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
	file_ms = open(file_masses, 'r')
	file_ns = open(file_n_subs, 'r')

	m_bins = pickle.load(file_m)
	n_bins = pickle.load(file_n)
	ms = pickle.load(file_ms)
	ns = pickle.load(file_ns)
	
	'''
	axis_margins = 1
        axis_size = 50
        axis_units = 'Mpc/h'
        ptsize_lg = 10.0
        col_lg="black"

	f_out = 'output/m_sat_scatter.png'

        plt.figure(figsize=(20,20))
        plt.rc('xtick', labelsize=axis_size)
        plt.rc('ytick', labelsize=axis_size)
        plt.rc('axes',  labelsize=axis_size)
        plt.margins(axis_margins)

	plt.scatter(m_bins, n_bins, s=ptsize_lg, c=col_lg)
        plt.tight_layout()
        plt.savefig(f_out)
	'''

	tot_bin = simu_end
	
	var_m = np.zeros((tot_bin, n_lgs, 2))
	var_n = np.zeros((tot_bin, n_lgs, 2))
	
	sig0 = 0.054;
	sig1 = 0.172;

	for ibin in range(0, tot_bin):
		for ilg in range(0, n_lgs):
			med_m = m_bins[ibin, ilg, 0]
			
			#m_bins[ibin, ilg, 1] = m_bins[ibin, ilg, 0] * (1 - sig0)
			#m_bins[ibin, ilg, 2] = m_bins[ibin, ilg, 0] * (1 + sig0)
			m_bins[ibin, ilg, 3] = m_bins[ibin, ilg, 0] * (1 - sig0)
			m_bins[ibin, ilg, 4] = m_bins[ibin, ilg, 0] * (1 + sig0)
			#min_m = m_bins[ibin, ilg, 3]
			#max_m = m_bins[ibin, ilg, 4]

			#var_m[ibin, ilg, 0] = abs(med_m - min_m)/med_m
			#var_m[ibin, ilg, 1] = abs(med_m - max_m)/med_m

		
			med_n = n_bins[ibin, ilg, 0]
			sig1 *= (1.0 + 1.0 / med_n)
			#n_bins[ibin, ilg, 1] = n_bins[ibin, ilg, 0] * (1 - sig1)
			#n_bins[ibin, ilg, 2] = n_bins[ibin, ilg, 0] * (1 + sig1)
			n_bins[ibin, ilg, 3] = n_bins[ibin, ilg, 0] * (1 - sig1)
			n_bins[ibin, ilg, 4] = n_bins[ibin, ilg, 0] * (1 + sig1)
			#min_n = n_bins[ibin, ilg, 0]
			#max_n = n_bins[ibin, ilg, 2]

			#var_n[ibin, ilg, 0] = abs(med_n - min_n)/med_n
			#var_n[ibin, ilg, 1] = abs(med_n - max_n)/med_n

	f_out = 'output/sat_n_bins.png'
	plot_lg_bins(m_bins, n_bins, f_out)

	'''
		NOW PLOT THE HISTOGRAMS OF THE MASS & SATELLITE FLUCTUATIONS WRT THE MEDIAN
	file_ms = open(file_masses, 'r')
	file_ns = open(file_n_subs, 'r')
	masses = pickle.load(file_ms)
	n_subs = pickle.load(file_ns)

	delta_m = [[] for ijk in range(0, 3)] 
	delta_n = [[] for ijk in range(0, 3)] 
	med_delta_m = np.zeros((3))
	med_delta_n = np.zeros((3))

	n_runs = len(n_subs[0])

	for ijk in range(0, 3):
		for i_run in range(0, n_runs):
			try:
				this_m_med = np.median(masses[ijk][i_run])
				this_n_med = np.median((n_subs[ijk][i_run]))
				n_lgs = len(masses[ijk][i_run])

				for ilg in range(0, n_lgs):
					this_ms = abs(masses[ijk][i_run][ilg] - this_m_med) / this_m_med
					this_ns = abs(float(n_subs[ijk][i_run][ilg]) - float(this_n_med)) / this_n_med

					delta_m[ijk].append(this_ms)
					delta_n[ijk].append(this_ns)
			except:
				print i_run, n_lgs, this_m_med

	for ijk in range(0, n_lgs):
		med_delta_m[ijk] = np.median(delta_m[ijk]) #, np.mean(delta_m[ijk]))
		med_delta_n[ijk] = np.median(delta_n[ijk]) #, np.mean(delta_n[ijk]))

	for ijk in range(0, n_lgs):
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
