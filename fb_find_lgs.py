from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halos import *
from libcosmo.find_halo import *
from libcosmo.particles import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from config import *
import numpy as np
import time
import pickle
from libcosmo.grid import *

# Simulation & catalog
file_single='snapshot_054.z0.000.AHF_halos'
box_size = 100000.0
base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'

# Local group selection parameters: very generous.
# Model 0 - root

#refineSelection = False
refineSelection = True
n_models = 5
select_params = np.zeros((n_models, 7), dtype='float')

iso_radius = 2000.0

ms_min = [4.0e+11, 5.0e+11, 6.0e+11, 7.0e+11, 7.5e+11, 8.0e+11]
ms_max = [5.0e+12, 4.0e+12, 3.0e+12, 2.5e+12, 2.0e+12, 1.5e+12]
n_dens = [0, 0, 0, 0, 0, 0]
n_loc_dens = [0, 0, 0, 0, 0, 0]


# TODO implement the six LG models for the paper

if refineSelection == True:
    i_ini = 0
    i_end = 0
else:
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
    print('Removing halos below ', m_min)
    m_halos = []

    for im in range(0, 6):
        n_loc_dens[im] = 0

    for halo in all_halos:
        if halo.m > m_min:
            m_halos.append(halo)
'''
# THIS COMMENTED BLOCK IS USED TO COMPUTE THE DENSITY OF HALOS IN DIFFERENT MASS RANGES
        for im in range(0, 6):
                        mmin = ms_min[im]
                        mmax = ms_max[im]

                        if halo.m > mmin and halo.m < mmax:
                            n_dens[im] = n_dens[im] + 1
                            n_loc_dens[im] = n_loc_dens[im] + 1

        for im in range(0, 6):
            print(im, '] Found ', n_loc_dens[im], ' above mass threshold, in total: ', n_dens[im])

for im in range(0, 6):
    print(im, ' ', float(n_dens[im]) / 5.e+6)
'''

if refineSelection == False:
    # Find suitable pairs in the selected halos
    all_lgs = find_lg(m_halos, lg_model)
    n_lgs = len(all_lgs)
    print('Found ', n_lgs)

    # Save all the identified LG pairs to a binary file
#        out_lgs = 'saved/rand_lgs_' + sub_path + '.pkl'
#        out_lgs = 'saved/rand_web_lgs_' + sub_path + '.pkl'
    out_lgs = 'saved/rand_all_lgs_' + sub_path + '.pkl'
#        out_lgs = 'saved/rand_select_lgs_' + sub_path + '.pkl'
    f_out_lgs = open(out_lgs, 'wb')
    pickle.dump(all_lgs, f_out_lgs)


i_ini = 0
i_end = 5
vol = 100.0 * 100.0 * 100.0

if refineSelection == True:
    for im in range(0, n_models):
        i_lgs = 0
        lg_tot = 0

        for i_dir in range(i_ini, i_end):
            sub_path = '%02d' % i_dir
            out_lgs = 'saved/rand_select_lgs_' + sub_path + '.pkl'
            f_lgs = open(out_lgs, 'rb')
            these_lgs = pickle.load(f_lgs)

            lg_tot = lg_tot + len(these_lgs)
#            print(out_lgs, len(these_lgs))

            for lg in these_lgs:
                condition1 = (lg.LG1.m > select_params[im, 4] and lg.LG1.m < select_params[im, 3])
                condition2 = (lg.LG2.m > select_params[im, 4] and lg.LG2.m < select_params[im, 3])
                condition3 = (lg.r_halos() > select_params[im, 2] and lg.r_halos() < select_params[im, 1])
                condition4 = (lg.v_radial() < select_params[im, 6])
                condition5 = (lg.LG1.m / lg.LG2.m < select_params[im, 5])

                if (condition1 and condition2 and condition3 and condition4 and condition5):
                    i_lgs = i_lgs + 1

#        print('N LGs: ', i_lgs, ' condition: ', im, ' tot LGS: ', lg_tot)

#        print('LG density: ', i_lgs / vol, ' per cubic Mpc/h')
        print(i_lgs / vol / 5.0)

print(lg_tot / vol / 5.0)

print('LG density: ', lg_tot / vol / 5.0, ' per cubic Mpc/h')
