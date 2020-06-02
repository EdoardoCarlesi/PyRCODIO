#!/usr/bin/python

'''
	This script is a wrapper for computing the following:
		- Find an LG candidate at each step including nearby satellites
		- For each halo (main LGs and satellites) trace its merging history
		- At each step, re-identify the subhalos for each candidate (simply taking all the halos within a factor * Rvir)
		- Use the subhalo distribution at each step to re-compute the anisotropy & inertia tensor eigenvalues
		- Do not do any plot, simply dump all the stuff to .txt - the analysis will be done separately
'''

import numpy as np
import os
import pickle

from config import *
from libio.read_ascii import *
from libcosmo.utils import *
from libcosmo.halos import *
from libcosmo.track_halos import *
from libcosmo.find_halos import *
from libcosmo.lg_plot import *

resolution='2048'

run_init = 0
run_end = 10

subrun_init = 0
subrun_end = 10

ini_snap = 0
end_snap = 54

base_path = '/home/eduardo/CLUES/'
outp_path = 'output/'
env_type = 'zoom'
snap_base = 'snapshot_'

save_path = 'saved/'
save_ext = '.pkl'

hubble = 0.67		

# Local group selection parameters
center = [50000., 50000., 50000.]

# Allocate LG Models
all_runs = simu_runs()
(all_lg_models, hash_run) = lg_models()

lg_dummy = LocalGroup(all_runs[0])
file_lg_header = lg_dummy.header()

# General halo settings
settings = Settings(base_path, outp_path, env_type, resolution, snap_base)
settings.box_center = center

# Subhalo identification criterion
fac_r = 1.2	# This is used for the global R around the LG as well as for Rvir around each halo
min_common = 15
part_min = 1000.
mmin = 1.e+8	# Track haloes above this threshold at z=0

# Loop on the different base - realisations, i.e. different LGs
for run_j in range(run_init, run_end):
	base_run = all_runs[run_j]
	this_run = hash_run[base_run]
	lg_model = all_lg_models[this_run]

	print base_run

	# Loop on the different small scale realisations of the LGs
	for subrun_i in range(subrun_init, subrun_end):
	
		# First find a Local Group candidate
		run_num = '%02d' % subrun_i
		settings.init_files(base_run, run_num)
		settings.re_init()
	
		# Look for LGs in the last snapshot first, use only one MPI task, initialize other useful variables
		snap_last = 54
		this_task = 0
		good_lgs = 0
		ids_sub = []
		main_ids =[]
		main_parts = []

		(this_file_part, this_file_halo) = settings.get_ahf_files(snap_last, this_task)

		print 'Reading files: ', this_file_halo, this_file_part
		halos = read_ahf(this_file_halo)
		(ids, parts) = read_particles(this_file_part) 

		this_lg = find_lg(halos, lg_model)
		n_lgs = int(len(this_lg))

		print 'Found a total of %d LG pairs.' % (n_lgs)

		if n_lgs > 0:
			# If there are more candidates we need to find the right one
			rating = 1000

			for ind in range(0, n_lgs):
				lg = this_lg[ind]
				lg.c_box = settings.box_center

				if lg.rating() < rating and (lg.LG1.npart > part_min) and (lg.LG2.npart > part_min):
					print 'Old rating: %f new rating %f this index %d' % (rating, lg.rating, ind)
					good_lgs += 1
					rating = lg.rating
					best_lg = lg

			if good_lgs > 0:
				print 'Best LG: ', best_lg.info()

				# Identify the main halos to be tracked at all zs
				com = best_lg.get_com()
				rad = best_lg.r_halos() * fac_r
				main_halos = find_halos_mass_radius(com, halos, rad, mmin)

				(main_ids, main_parts) = find_ids(main_halos, ids, parts)

				# Select only the two main LG candidates to keep track of all the satellites
				ids_sub.append(best_lg.LG1.ID)
				ids_sub.append(best_lg.LG2.ID)
				track_subs = True
	                        '''
				# NOW DO THE FULL BLOWN IDENTIFICATION OF THE LGs and their satellites at all steps!
				(mains, sats) = merger_tree(end_snap, ini_snap, min_common, main_halos, main_parts, \
							main_ids, halos, settings, track_subs, ids_sub, fac_r)
                                
				file_mains = save_path + base_run + '_' + run_num + '_mains_all' + save_ext
				file_sats = save_path + base_run + '_' + run_num + '_sats' + save_ext

				filehand_mains = open(file_mains, 'w')
				filehand_sats = open(file_sats, 'w')

				pickle.dump(mains, filehand_mains)
				pickle.dump(sats, filehand_sats)
                                '''
