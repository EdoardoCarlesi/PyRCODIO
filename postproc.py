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
sub_end = 7

snap_init = 0
snap_end = 54

base_path = '/home/eduardo/CLUES/PyRCODIO/'
outp_path = 'output/'
save_path = 'saved/'

n_part = 75

#simurun = '00_06'
simurun = '01_12'

#do_plots = "true"
#do_plots = "false"

n_main = 2
this_main = 0

mass_histories = np.zeros((n_main, sub_end - sub_init, snap_end-snap_init))
anisotropies = np.zeros((n_main, sub_end - sub_init, snap_end-snap_init, 3))

for i_sub in range(sub_init, sub_end):

	subrun = '%02d' % i_sub

	s_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_sats.pkl'
	m_fname='/home/eduardo/CLUES/PyRCODIO/saved/'+simurun+'_'+subrun+'_mains.pkl'

	hand_main = open(m_fname, 'r')
	hand_sats = open(s_fname, 'r')

	main = pickle.load(hand_main)
	sats = pickle.load(hand_sats)

	#this_mah = main[this_main].m_t() 
	#this_mah[:] /= this_mah[this_main]
	#mass_histories[this_main, i_sub] = this_mah


	for i_snap in range(snap_init, snap_end):

		if sats[this_main][i_snap].n_sub > 5:
			(evals, red_evals) = sats[this_main][i_snap].anisotropy('part', n_part)			
			#print i_snap, evals
			print i_sub, i_snap, (evals[1] - evals[0])/evals[2]
			anisotropies[this_main, i_sub, i_snap] = evals 

	#print 'LMM: ', main[this_main].last_major_merger()
	#print 'FT : ', main[this_main].formation_time()
	#print this_mah
f_out = 'output/anis_test_e1.png'
plot_anisotropies(anisotropies, this_main, sub_end, snap_end, f_out)

