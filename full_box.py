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

file_single='snapshot_054.z0.000.AHF_halos'

box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/'
sub_path = '00'
root_file = 'snapshot_054.'
suff_file = '.z0.000.AHF_halos'

tot_files = 120

m_max = 2.e+15
m_min = 1.e+11

# Local group selection parameters
iso_radius = 2000.
radius = 1000. 
r_max = 1300.
r_min = 350. 
m_min = 2.e+11  
m_max = 4.0e+12 
ratio_max = 4.
vrad_max = 10.0

lg_model = LocalGroupModel(radius, iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
this_file = base_path + sub_path + '/' + file_single

all_halos = read_ahf_mass_range(this_file, m_max, m_min)
all_lgs = find_lg(all_halos, lg_model)

out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
f_out_lgs = open(out_lgs, 'w')

pickle.dump(f_out_lgs, all_lgs)


'''
all_lgs = []
all_halos = []
for i_cpu in range(0, tot_files):
	this_cpu_num = '%04d' % i_cpu
	this_file = base_path + sub_path + root_file + this_cpu_num + suff_file
	this_ahf = read_ahf_mass_range(this_file, m_max, m_min)
	all_halos.append(this_ahf)

all_lgs = find_lg(all_halos, lg_model)
print 'Found a total LG candidates = ', len(all_lgs)

this_lgs = find_lg(this_ahf, lg_model)
n_lgs = len(this_lgs)
	
	print 'Found ', n_lgs
	if n_lgs > 0:
		all_lgs.append(this_lgs)
		for ilg in range(0, n_lgs):
			print this_lgs[ilg].info()
'''
