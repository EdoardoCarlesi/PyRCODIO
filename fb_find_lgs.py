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
iso_radius = 2000.0
r_max = 1500.0
r_min = 250.0
m_min = 4.0e+11 
m_max = 5.0e+12 
ratio_max = 4.0
vrad_max = 25.0

refineSelection = False
#refineSelection = True
n_models = 5
select_params = np.zeros((n_models, 7), dtype='float')

iso_radius = 2000.0

ms_min = [4.0e+11, 5.0e+11, 6.0e+11, 7.0e+11, 7.5e+11, 8.0e+11]
ms_max = [5.0e+12, 4.0e+12, 3.0e+12, 2.5e+12, 2.0e+12, 1.5e+12]
n_dens = [0, 0, 0, 0, 0, 0]
n_loc_dens = [0, 0, 0, 0, 0, 0]

if refineSelection == True:
    # Model 1 
    r_max = 1300.0
    r_min = 300.0
    m_min = 4.5e+11 
    m_max = 4.0e+12 
    ratio_max = 3.0
    vrad_max = 0.0

    select_params[0, 0] = iso_radius
    select_params[0, 1] = r_max
    select_params[0, 2] = r_min
    select_params[0, 3] = m_max
    select_params[0, 4] = m_min 
    select_params[0, 5] = ratio_max
    select_params[0, 6] = vrad_max 

    # Model 2 
    r_max = 1000.0
    r_min = 350.0
    m_min = 5.0e+11 
    m_max = 3.0e+12 
    ratio_max = 3.0
    vrad_max = -25.0

    select_params[1, 0] = iso_radius
    select_params[1, 1] = r_max
    select_params[1, 2] = r_min
    select_params[1, 3] = m_max
    select_params[1, 4] = m_min 
    select_params[1, 5] = ratio_max
    select_params[1, 6] = vrad_max 

    # Model 3 
    r_max = 900.0
    r_min = 400.0
    m_min = 5.5e+11 
    m_max = 2.5e+12 
    ratio_max = 2.5
    vrad_max = -50.0

    select_params[2, 0] = iso_radius
    select_params[2, 1] = r_max
    select_params[2, 2] = r_min
    select_params[2, 3] = m_max
    select_params[2, 4] = m_min 
    select_params[2, 5] = ratio_max
    select_params[2, 6] = vrad_max 

    # Model 4 
    r_max = 800.0
    r_min = 450.0
    m_min = 6.0e+11 
    m_max = 2.0e+12 
    ratio_max = 2.5
    vrad_max = -75.0

    select_params[3, 0] = iso_radius
    select_params[3, 1] = r_max
    select_params[3, 2] = r_min
    select_params[3, 3] = m_max
    select_params[3, 4] = m_min 
    select_params[3, 5] = ratio_max
    select_params[3, 6] = vrad_max 

    # Model 5
    r_max = 700.0
    r_min = 500.0
    m_min = 6.5e+11 
    m_max = 1.5e+12 
    ratio_max = 2.0
    vrad_max = -100.0

    select_params[4, 0] = iso_radius
    select_params[4, 1] = r_max
    select_params[4, 2] = r_min
    select_params[4, 3] = m_max
    select_params[4, 4] = m_min 
    select_params[4, 5] = ratio_max
    select_params[4, 6] = vrad_max 


'''
TEST PARAMETERS
iso_radius = 2000.0
r_max = 1500.0
r_min = 250.0
m_min = 1.0e+11 
m_max = 5.0e+13 
ratio_max = 20.0
vrad_max = 1000.0
'''

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
    
    select_params[2, 0] = iso_radius
    select_params[2, 1] = r_max
    select_params[2, 2] = r_min
    select_params[2, 3] = m_max
    select_params[2, 4] = m_min 
    select_params[2, 5] = ratio_max
    select_params[2, 6] = vrad_max 
'''

i_ini = 0
i_end = 5


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
            print(out_lgs, len(these_lgs))

            for lg in these_lgs:
                condition1 = (lg.LG1.m > select_params[im, 4] and lg.LG1.m < select_params[im, 3])                
                condition2 = (lg.LG2.m > select_params[im, 4] and lg.LG2.m < select_params[im, 3])                
                condition3 = (lg.r_halos() > select_params[im, 2] and lg.r_halos() < select_params[im, 1])                
                condition4 = (lg.v_radial() < select_params[im, 6])
                condition5 = (lg.LG1.m / lg.LG2.m < select_params[im, 5])

                if (condition1 and condition2 and condition3 and condition4 and condition5):
                    i_lgs = i_lgs + 1

        print('N LGs: ', i_lgs, ' condition: ', im, ' tot LGS: ', lg_tot)

        vol = 100.0 * 100.0 * 100.0
        print('LG density: ', i_lgs / vol, ' per cubic Mpc/h')

