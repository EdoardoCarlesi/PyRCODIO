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

#file_single='snapshot_054.z0.000.AHF_halos'
#file_single='snapshot_054.0000.z0.000.AHF_halos'
file_single='snapshot_176.0000.z0.000.AHF_halos'

box_size = 100000.0
#base_path = '/home/eduardo/CLUES/DATA/FullBox/'
base_path = '/home/edoardo/CLUES/DATA/SIMULATIONS/LGF/1024/'
#sub_path = '00'
sub_path = '00_06'
root_file = 'snapshot_'
suff_halo = '.z0.000.AHF_halos'
suff_part = '.z0.000.AHF_particles'
this_ahf_file = base_path + sub_path + '/' + file_single

tot_files = 1
use_files = 1

m_max = 2.e+15
m_min = 1.e+11

# Local group selection parameters
iso_radius = 2200.
radius = 1000. 
r_max = 1500.
r_min = 350. 
m_min = 5.e+11  
m_max = 5.0e+12 
ratio_max = 5.0
vrad_max = 0.0

resolution = '1024'
env_type = 'std'
out_dir = 'output/'

end_snap = 176
ini_snap = 177

settings = Settings(base_path, out_dir, env_type, resolution, root_file)
settings.base_file_chunk = base_path + sub_path + '/' + root_file
settings.ahf_path = base_path + sub_path 

lg_model = LocalGroupModel(radius, iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)
this_root = base_path + sub_path + '/' + root_file + '176.'

all_halos = read_ahf(this_ahf_file)
n_halos = len(all_halos)

print('read in %d halos, finding Local Groups...', n_halos)

# Filter by halo mass
m_halos = []

print('Removing halos below ', m_min)

for ih in range(0, n_halos):
	if all_halos[ih].m > m_min:
		m_halos.append(all_halos[ih])

print('Found ', len(m_halos), ' above mass threshold')

all_lgs = find_lg(m_halos, lg_model)

n_lgs = len(all_lgs)
	
print('Found ', n_lgs)

#out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
out_lgs = 'saved/test_lgs_' + sub_path + '.pkl'
f_out_lgs = open(out_lgs, 'w')
pickle.dump(all_lgs, f_out_lgs)
