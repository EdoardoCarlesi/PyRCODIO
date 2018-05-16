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

# Number of subrun
sub_init = 0
sub_end = 10

# Snapshots
snap_init = 0
snap_end = 54

base_path = '/home/eduardo/CLUES/PyRCODIO/'
outp_path = 'output/'
save_path = 'saved/'

min_part = 30
stepMyr = 0.25

#simurun = '00_06'
#simurun = '01_12'
#simurun = '45_17'
#simurun = '09_18'
#simurun = '17_10'
#simurun = '17_13'
#simurun = '55_02'
#simurun = '34_13'
simurun = '64_14'

#do_plots = "true"
#do_plots = "false"

lg_names = ['MW', 'M31']
n_lg = len(lg_names)

mass_histories = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init))
anisotropies = np.zeros((n_lg, sub_end - sub_init, snap_end-snap_init, 3))
satellite_number = 

time = np.zeros((snap_end-snap_init))

for i_time in range(0, snap_end-snap_init):
	time[i_time] = (snap_end - i_time) * stepMyr

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
		r_sub = 1.35
		n_lg = 2
	except:	
		n_lg = 0

	# The first two are the main LG members
	for i_lg in range(0, n_lg):
		this_lg = main[i_lg]
		this_center = this_lg.x_t()
		lg_r = this_lg.halo[0].r
		this_mt = this_lg.m_t()

		# Plot mass accretion history of the main halo
		out_mah = 'output/' + lg_names[i_lg] + '_main_' + simurun + '_' + subrun + '_mah.png'
		#plot_mass_accretion(time, this_mt, out_mah)
		mass_histories[i_lg, i_sub] = this_mt

		# For each LG member we are appending subhalo positions and histories
		subs = []
		subs_xt = []
		subs_yt = []
		subs_zt = []

		# This is a loop on all the other haloes stored at z=0
		for i_main in range(n_lg, n_main):
			this_main = main[i_main].halo[0]
			this_run = '%02d' % i_main
			this_x = this_main.x

			if this_lg.halo[0].distance(this_x) < lg_r * r_sub:
				this_xt = main[i_main].x_t_center(this_center)
				this_mt = main[i_main].m_t()

				this_sub_z = SubHaloThroughZ(snap_end-snap_init)

				this_sub_z.host = this_lg
				this_sub_z.assign_halo_z(this_lg)

				subs_xt.append(this_xt[0, :])
				subs_yt.append(this_xt[1, :])

		#Planes of satellites in substructure
		for i_snap in range(0, snap_end-snap_init):

			# Only computes the inertia tensor considering satellites above N=min_part particles
			if sats[i_lg][i_snap].n_sub > 5:
				(evals, red_evals, evecs, red_evecs) = sats[i_lg][i_snap].anisotropy('part', min_part)			
				anisotropies[i_lg, i_sub, i_snap] = evals
				#anisotropies[i_lg, i_sub, i_snap] = red_evals

for i_lg in range(0, 2):
	out_fname = 'output/anisotropy_' + simurun + '_' + lg_names[i_lg] + '_' + subrun + '_' + this_run + '.png'
	plot_anisotropies(anisotropies, i_lg, sub_end, snap_end, out_fname)

print 'Plotting mass accretions...'
out_fname = 'output/mah_' + simurun + '_' + lg_names[0] + '_' + this_run + '.png'
plot_mass_accretions(time, mass_histories[0, :, :], out_fname)
out_fname = 'output/mah_' + simurun + '_' + lg_names[1] + '_' + this_run + '.png'
plot_mass_accretions(time, mass_histories[1, :, :], out_fname)

