#!/usr/bin/python

import numpy as np
import os

from libio.read_ascii import *

from libcosmo.utils import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.lg_plot import *
import pickle

resolution='2048'

sub_init = 0
sub_end = 1

snap_init = 0
snap_end = 54

base_path = '/home/eduardo/CLUES/PyRCODIO/'
outp_path = 'output/'
save_path = 'saved/'

n_part = 50

stepMyr = 0.25

simurun = '00_06'
#simurun = '01_12'

#do_plots = "true"
#do_plots = "false"

n_main = 2

mass_histories = np.zeros((n_main, sub_end - sub_init, snap_end-snap_init))
anisotropies = np.zeros((n_main, sub_end - sub_init, snap_end-snap_init, 3))

time = np.zeros((snap_end-snap_init))

for i_time in range(0, snap_end-snap_init):
	time[i_time] = (snap_end - i_time) * stepMyr

print time

for i_sub in range(sub_init, sub_end):

	subrun = '%02d' % i_sub

	s_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_sats.pkl'
	m_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_mains.pkl'

	hand_main = open(m_fname, 'r')
	hand_sats = open(s_fname, 'r')

	main = pickle.load(hand_main)
	sats = pickle.load(hand_sats)

	'''
		Analyze the trajectories and properties of the main LG halos
	'''

	n_main = len(main)
	m31_main = 1
	mw_main = 0
	r_sub = 1.5

	mw_subs = []
	mw_subs_xt = []
	mw_subs_yt = []
	mw_subs_zt = []
	m31_subs = []
	m31_subs_xt = []

	x_mw0 = main[mw_main].halo[0].x
	r_mw0 = main[mw_main].halo[0].r
	xt_mw = main[mw_main].x_t()

	x_m310 = main[m31_main].halo[0].x
	r_m310 = main[mw_main].halo[0].r
	xt_m31 = main[m31_main].x_t()

	center = 0.5 * (xt_m31 + xt_mw)
	#print center

	for i_main in range(0, n_main):
		this_sub = main[i_main].halo[0]
		this_run = '%02d' % i_main

		if this_sub.distance(x_mw0) < r_mw0 * r_sub:
			#this_xt = main[i_main].x_t()
			this_xt = main[i_main].x_t_center(center)

			this_mt = main[i_main].m_t()

			out_mah = 'output/' + simurun + '_' + subrun + '_mah_' + this_run + '.png'
			plot_mass_accretion(time, this_mt, out_mah)

			this_sub_z = SubHaloThroughZ(snap_end-snap_init)

			this_sub_z.host = main[mw_main]
			this_sub_z.assign_halo_z(main[i_main])

			print this_sub_z.accretion_time()

	#		print i_main, this_xt
			mw_subs_xt.append(this_xt[0, :])
			mw_subs_yt.append(this_xt[1, :])

	out_fname = 'output/' + simurun + '_' + subrun + '_all_trajectories_' + this_run + '.png'
	plot_trajectory(mw_subs_xt, mw_subs_yt, 'x', 'y', out_fname)

	'''
		Planes of satellites in substructure
	for i_snap in range(0, snap_end-snap_init):

		if sats[this_main][i_snap].n_sub > 5:
			(evals, red_evals) = sats[this_main][i_snap].anisotropy('part', n_part)			
			#print i_snap, evals
			#print i_sub, i_snap, (evals[1] - evals[0])/evals[2]
			anisotropies[this_main, i_sub, i_snap] = evals 

	#print 'LMM: ', main[this_main].last_major_merger()
	#print 'FT : ', main[this_main].formation_time()
	#print this_mah
f_out = 'output/anis_test_e1.png'
plot_anisotropies(anisotropies, this_main, sub_end, snap_end, f_out)

	'''
