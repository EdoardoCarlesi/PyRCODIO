#!/usr/bin/python

from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.particles import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from config import *
import time
import pickle
from libcosmo.grid import *

# Simulation & catalog
file_single='snapshot_054.z0.000.AHF_halos'
#file_single='snapshot_054.0000.z0.000.AHF_halos'
box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/'
sub_path = '00'

# Local group selection parameters
iso_radius = 2000.
radius = 1000. 
r_max = 1500.
r_min = 250. 
m_min = 4.e+11  
m_max = 5.0e+12 
ratio_max = 10.0
vrad_max = 100.0

i_ini = 0
i_end = 5

for i_dir in range(i_ini, i_end):
	sub_path = '%02d' % i_dir

	lg_model = LocalGroupModel(iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
	this_ahf_file = base_path + sub_path + '/' + file_single

	print('Reading file : %s' % this_ahf_file)
	all_halos = read_ahf(this_ahf_file)
	n_halos = len(all_halos)

	print('%d halos have been read in, finding Local Groups...' % n_halos)

	# Filter by halo mass
	print 'Removing halos below ', m_min
	m_halos = []

	for halo in all_halos:
		if halo.m > m_min:
			m_halos.append(halo)

	print 'Found ', len(m_halos), ' above mass threshold'

	all_lgs = find_lg(m_halos, lg_model)
	n_lgs = len(all_lgs)
	print 'Found ', n_lgs

	out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
	f_out_lgs = open(out_lgs, 'w')
	pickle.dump(all_lgs, f_out_lgs)
