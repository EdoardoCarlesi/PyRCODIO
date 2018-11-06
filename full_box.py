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
#file_single='snapshot_054.0000.z0.000.AHF_halos'

box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/'
sub_path = '02'
root_file = 'snapshot_'
suff_halo = '.z0.000.AHF_halos'
suff_part = '.z0.000.AHF_particles'
this_ahf_file = base_path + sub_path + '/' + file_single

tot_files = 1
use_files = 1

m_max = 2.e+15
m_min = 1.e+11

# Local group selection parameters
iso_radius = 2500.
radius = 1000. 
r_max = 1200.
r_min = 450. 
m_min = 5.e+11  
m_max = 5.0e+12 
ratio_max = 3.
vrad_max = 0.0

resolution = '1024'
env_type = 'std'
out_dir = 'output/'

min_common = 15

end_snap = 55
ini_snap = 54

settings = Settings(base_path, out_dir, env_type, resolution, root_file)
settings.base_file_chunk = base_path + sub_path + '/' + root_file

#settings.ahf_path = base_path + sub_path + '/snapshot_' + '*.0000.*' 
settings.ahf_path = base_path + sub_path #+ '/snapshot_'

lg_model = LocalGroupModel(radius, iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
this_root = base_path + sub_path + '/' + root_file + '054.'

#all_halos = read_ahf_chunks_mass_range(this_root, suff_halo, tot_files, m_max, m_min)
all_halos = read_ahf(this_ahf_file)
#all_halos = read_ahf_chunks(this_root, suff_halo, use_files)
n_halos = len(all_halos)

print n_halos, 'read in, finding Local Groups...'

# Filter by halo mass
m_halos = []

print 'Removing halos below ', m_min

for ih in range(0, n_halos):
	if all_halos[ih].m > m_min:
		m_halos.append(all_halos[ih])

print 'Found ', len(m_halos), ' above mass threshold'

#(all_ids, all_parts) = read_particles_chunks(this_root, suff_part, use_files, n_halos)
all_lgs = find_lg(m_halos, lg_model)

#print all_ids

n_lgs = len(all_lgs)
	
print 'Found ', n_lgs

main_halos = []

'''
if n_lgs > 0:
	for ilg in range(0, n_lgs):
		print all_lgs[ilg].info()
		main_halos.append(all_lgs[ilg].LG1)
		main_halos.append(all_lgs[ilg].LG2)


(main_ids, main_parts) = find_ids(main_halos, all_ids, all_parts)

lg_halos_z = merger_tree_chunks(end_snap, ini_snap, tot_files, use_files, min_common, main_halos, main_parts, main_ids, settings)


#n_halos = 15171
#(parts, ids) = read_particles_chunks(this_root, suff_part, tot_files, n_halos)


out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
f_out_lgs = open(out_lgs, 'w')
#pickle.dump(all_lgs, f_out_lgs)

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
