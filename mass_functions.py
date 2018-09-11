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


simuruns = simu_runs()
n_subrun = 10
snap_end = 54
snap_init = 0


for i_simu in range(1, 10):
	simu = simuruns[i_simu]
	n_simu = 0

	mfmw_z0_x = []; 		mfmw_z0_y = []
	mfmw_max_x = []; 		mfmw_max_y = []
	mfm31_z0_x = []; 		mfm31_z0_y = []
	mfm31_max_x = []; 		mfm31_max_y = []

	for i_sub in range(0, n_subrun):
		subrun = '%02d' % i_sub
		m_fname='saved/'+simu+'_'+subrun+'_mains_all.pkl'
		s_fname='saved/'+simu+'_'+subrun+'_sats.pkl'
		n_lg = 0

		try:
			hand_main = open(m_fname, 'r')
			hand_sats = open(s_fname, 'r')
			main = pickle.load(hand_main)
			sats = pickle.load(hand_sats)

			print 'Found: ', s_fname
			print 'Found: ', m_fname
			#print main[0].halo
			n_main = len(main)
			#print n_main
			r_sub = 1.5
			n_lg = 2
			n_simu += 1
		except:	
			print '' #Nothing found'

		masses_max = []; 	masses_z0 = [];		
		
		for i_lg in range(0, n_lg):
			this_lg = main[i_lg]
			this_center = this_lg.x_t()
			this_run = '%02d' % i_lg
			lg_r = this_lg.halo[0].r
			this_mt = this_lg.m_t()

			#print i_lg

			# This is a loop on all the other haloes stored at z=0
			for i_main in range(n_lg, n_main):
				this_main = main[i_main]
				this_x = this_main.halo[0].x
				this_d = this_lg.halo[0].distance(this_x) 
				dist_rvir = this_d / lg_r

				# This halo is within the range from the virial radius
				if this_d < lg_r * r_sub:
					#acc_time = this_sub_z.accretion_time()
					this_sub_z = SubHaloThroughZ(snap_end-snap_init)
					this_sub_z.host = this_lg
					this_sub_z.assign_halo_z(this_main)
					acc_time = this_sub_z.accretion_time()
					m_z0 = this_sub_z.halo[0].m
					m_max = this_sub_z.m_max()
					masses_max.append(m_max)
					masses_z0.append(m_z0)

			# End of the loop on i_main, still looping on i_lg
			(this_mfz0_x, this_mfz0_y) = mass_function(masses_z0)
			(this_mfmax_x, this_mfmax_y) = mass_function(masses_max)

			'''	
				Save informations on subhalo positions through time
			'''
			if i_lg == 0:
				mfmw_z0_x.append(this_mfz0_x); 		mfmw_z0_y.append(this_mfz0_y) 
				mfmw_max_x.append(this_mfmax_x); 	mfmw_max_y.append(this_mfmax_y) 

			elif i_lg == 1:
				mfm31_z0_x.append(this_mfz0_x); 	mfm31_z0_y.append(this_mfz0_y) 
				mfm31_max_x.append(this_mfmax_x); 	mfm31_max_y.append(this_mfmax_y) 

	print 'Loop on ', simu, ' finished. Computing mass functions...'

	f_mwz0 = simu + '_mf_z0_MW.png'
	plot_massfunctions(mfmw_z0_x, mfmw_z0_y, n_simu, f_mwz0)

	f_m31z0 = simu + '_mf_z0_M31.png'
	plot_massfunctions(mfm31_z0_x, mfm31_z0_y, n_simu, f_m31z0)

	f_mwmax = simu + '_mf_max_MW.png'
	plot_massfunctions(mfmw_max_x, mfmw_max_y, n_simu, f_mwmax)

	f_m31max = simu + '_mf_max_M31.png'
	plot_massfunctions(mfm31_max_x, mfm31_max_y, n_simu, f_m31max)




