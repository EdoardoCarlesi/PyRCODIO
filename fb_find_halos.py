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
box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'

# Local group selection parameters: very generous. Select better afterwards.
iso_radius = 2000.0
r_max = 1500.0
r_min = 250.0
m_min = 9.0e+11  
m_max = 5.0e+13 
ratio_max = 20.0
vrad_max = 1000.0

i_ini = 0
i_end = 4

for i_dir in range(i_ini, i_end):
        sub_path = '%02d' % i_dir
        this_ahf_file = base_path + sub_path + '/' + file_single

        print('Reading file : %s' % this_ahf_file)
        all_halos = read_ahf(this_ahf_file) 
        n_halos = len(all_halos)
        print('%d halos have been read in, finding suitable halos...' % n_halos)

	# Filter by halo mass
        print('Removing halos below ', m_min, ' and above ', m_max)
        m_halos = []

        for halo in all_halos:
            if halo.m > m_min and halo.m < m_max:
                    m_halos.append(halo)

        print('Found ', len(m_halos), ' within mass limits')

        out_halos = 'saved/rand_halos_' + sub_path + '.pkl'
        f_out_halos = open(out_halos, 'wb')
        pickle.dump(m_halos, f_out_halos)
