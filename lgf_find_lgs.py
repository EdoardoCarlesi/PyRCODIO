#!/usr/bin/python

from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from config import *
import time
import pickle
import os.path
from libcosmo.grid import *

# Simulation & catalog 
box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/LGF/1024/'
file_single0 = 'snapshot_054.z0.000.AHF_halos'
#file_single0 = 'snapshot_054.0000.z0.000.AHF_halos'
file_single1 = 'snapshot_054.0000.z0.001.AHF_halos'

# Search distance from the box center
r_select = 8500.
m_select = 5.e+11
box_center = np.zeros((3))

for ix in range(0, 3):
	box_center[ix] = box_size * 0.5

# Local group selection parameters
iso_radius = 2000.
r_max = 1300.
r_min = 350. 
m_min = 5.e+11  
m_max = 5.0e+12 
ratio_max = 4.0
vrad_max = 10.0

# Initialzie a local group model
lg_model = LocalGroupModel(iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)

i_ini = 0
i_end = 80
g_ini = 0
g_end = 20
n_tot_lgs = 0

for i_dir in range(i_ini, i_end):
	i_dir_str = '%02d' % i_dir
	
	for g_dir in range(g_ini, g_end):
		g_dir_str = '%02d' % g_dir
		sub_path = i_dir_str + '_' + g_dir_str
		ahf_file0 = base_path + sub_path + '/' + file_single0
		ahf_file1 = base_path + sub_path + '/' + file_single1
		this_ahf_file = None

		if os.path.isfile(ahf_file0):
			this_ahf_file = ahf_file0
		elif os.path.isfile(ahf_file1):
			this_ahf_file = ahf_file1

		if this_ahf_file != None:
			#print('Reading file : %s' % this_ahf_file)
			all_halos = read_ahf(this_ahf_file)

			n_halos = len(all_halos)
			#print('%d halos have been read in, finding Local Groups...' % n_halos)

			# Filter by halo mass
			m_halos = []

			#print('Removing halos below %e and outside of %f Mpc/h of the box center.' % (m_min, r_select))

			for halo in all_halos:
				this_r = halo.distance(box_center)

				if halo.m > m_select and this_r < r_select:
					m_halos.append(halo)

			print('Found %d halos above mass threshold within the search radius ' % len(m_halos))
			all_lgs = find_lg(m_halos, lg_model)

			for one_lg in all_lgs:
				print(one_lg.LG1.distance(box_center), one_lg.info())

			n_lgs = len(all_lgs)
			print('Found %s LG candidates.' % n_lgs)

			if n_lgs > 0:
				n_tot_lgs += 1
				out_lgs = 'saved/lgs_' + sub_path + '.pkl'
				f_out_lgs = open(out_lgs, 'wb')
				print('Saving LGs of %s to file %s.' % (sub_path, out_lgs))
				pickle.dump(all_lgs, f_out_lgs)
