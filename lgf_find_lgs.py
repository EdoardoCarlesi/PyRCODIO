from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halos import *
from libcosmo.find_halo import *
#from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from config import *
import time
import pickle
import os.path
import numpy as np
from libcosmo.grid import *

# Simulation & catalog
box_size = 100000.0

#base_path = '/home/eduardo/CLUES/DATA/LGF/1024/'
#file_single0 = 'snapshot_054.z0.000.AHF_halos'
file_single0 = 'snapshot_054.0000.z0.000.AHF_halos'
#file_single0 = 'snapshot_054.0000.z0.000.AHF_halos'
file_single1 = 'snapshot_054.0000.z0.001.AHF_halos'
'snapshot_054.0000.z0.001.AHF_halos'

#all_dirs=['37_11']
all_dirs=simu_runs()

#simu='37_11'
#res='4096'
#res='2048'
res='512'

#print(all_dirs)
#dir_path = '/home/eduardo/CLUES/DATA/' + res + '/'
#dir_path = '/z/carlesi/STORE/LGF/SNAPS/' + res + '_Hestia_LSS/'
dir_path = '/z/carlesi/STORE/LGF/SNAPS/' + res + '/'
#dir_path = '/home/eduardo/CLUES/DATA/' + res + '/'

base_path = dir_path

# Search distance from the box center
r_select = 5000.
box_center = np.zeros((3))

for ix in range(0, 3):
    box_center[ix] = box_size * 0.5

doMoreSelect = False
doMoreSelect = True

# Local group selection parameters
r_max = 1500.
r_min = 250.
m_min = 4.0e+11
m_max = 5.0e+12
ratio_max = 4.0
vrad_max = 100.0

iso_radius = 2000.
m_select = m_min

mMinStr = str(m_min)
rStr = str(r_select)

# Initialzie a local group model
lg_model = LocalGroupModel(iso_radius, r_max, r_min, m_max, m_min, ratio_max, vrad_max)

i_ini = 0
i_end = 10
g_ini = 0
g_end = 30

# initialize some counters
n_files = 0
n_tot_lgs = 0
n_tot = 0

file_lg = 'info_lg.txt'
file_mw = 'info_mw.txt'
file_m31 = 'info_m31.txt'
file_vweb = 'info_vweb.txt'

f_lg = open(file_lg, 'w')
f_mw = open(file_mw, 'w')
f_m31 = open(file_m31, 'w')

mw_head = '#SimuCode(0)  ' + Halo().header(1) + '\n'
lg_head = '' + LocalGroup('').header() + '\n'
f_lg.write(lg_head)
f_mw.write(mw_head)
f_m31.write(mw_head)

#for dir_str in all_dirs:       ### 2048 & 4096 mode
for i_dir in range(i_ini, i_end):
    #simu = dir_str
    #base_path = dir_path + dir_str + '/'
    i_dir_str = '%02d' % i_dir

    for g_dir in range(g_ini, g_end):
        g_dir_str = '%02d' % g_dir
        sub_path = i_dir_str + '_' + g_dir_str  ### USE THIS FOR 512 / 1024

        ahf_file0 = base_path + sub_path + '/' + file_single0
        ahf_file1 = base_path + sub_path + '/' + file_single1
        ahf_file2 = base_path + sub_path + '/' + 'HESTIA_100Mpc_512_' + sub_path + '.127.z0.000.AHF_halos'
        this_ahf_file = None

        if os.path.isfile(ahf_file0):
            this_ahf_file = ahf_file0
        elif os.path.isfile(ahf_file1):
            this_ahf_file = ahf_file1
        elif os.path.isfile(ahf_file2):
            this_ahf_file = ahf_file2

        isPkl = False
        this_pkl = 'saved/lgs_r_' + rStr + '_mMin' + mMinStr + '_' + res + '_' + sub_path + '.pkl'

        if os.path.isfile(this_pkl):
            isPkl = True

        if this_ahf_file != None:
            n_files = n_files + 1

            all_halos = read_ahf(this_ahf_file)

            for halo in all_halos:
                this_r = halo.distance(box_center)

                if halo.m > m_min and this_r < r_select and halo.m < m_max:
                    m_halos.append(halo)

            n_tot = n_tot + len(m_halos)

vol = np.power(r_select/1000., 3) * 4.0 / 3.0 * 3.14
dens = n_tot / (vol * n_files)
print('nTot = ', n_tot, ' dens= ', dens, ' nFiles: ', n_files)

'''
        if this_ahf_file != None and isPkl == False:
            print('Reading file : %s' % this_ahf_file)
            all_halos = read_ahf(this_ahf_file)

            n_halos = len(all_halos)

        # Filter by halo mass
            m_halos = []

            for halo in all_halos:
                    this_r = halo.distance(box_center)

                    if halo.m > m_select and this_r < r_select:
                            m_halos.append(halo)

            all_lgs = find_lg(m_halos, lg_model)

            for one_lg in all_lgs:
                    one_lg.code = sub_path
                    f_lg.write(one_lg.info() + '\n')
                    f_mw.write(sub_path + '\t' + one_lg.LG2.line)
                    f_m31.write(sub_path + '\t' + one_lg.LG1.line)

            n_lgs = len(all_lgs)
            #print('Found %s LG candidates.' % n_lgs)

            if n_lgs == 0:
                    out_lgs = this_pkl
                    f_out_lgs = open(out_lgs, 'wb')
                    void_lg = [LocalGroup('EMPTY'), LocalGroup('EMPTY')]
                    pickle.dump(void_lg, f_out_lgs)

            elif n_lgs > 0:
                    n_tot_lgs = n_tot_lgs + n_lgs
                    sub_halos = find_halos_mass_radius(all_lgs[0].geo_com(), all_halos, iso_radius, 0.0)
                    #out_lgs = 'saved/lgs_' + res + '_' + simu + '_' + sub_path + '.pkl'
                    #out_subs = 'saved/subs_' + res + '_' + simu + '_' +  sub_path + '.pkl'
                    #out_lgs = 'saved/lgs_' + res + '_' + sub_path + '.pkl'
                    #out_subs = 'saved/subs_' + res + '_' + sub_path + '.pkl'
                    out_lgs = this_pkl #'saved/lgs_cv_' + res + '_' + sub_path + '.pkl'
                    out_subs = 'saved/subs_cv_' + res + '_' + sub_path + '.pkl'
                    f_out_lgs = open(out_lgs, 'wb')
                    f_out_subs = open(out_subs, 'wb')
                    print('Saving LGs of %s to file %s.' % (sub_path, out_lgs))
                    #print('Saving SUBs of %s to file %s.' % (sub_path, out_subs))

                    pickle.dump(sub_halos, f_out_subs)
                    f_out_subs.close()

            lgs_pkl = all_lgs

            pickle.dump(lgs_pkl, f_out_lgs)
            f_out_lgs.close()

    if isPkl == True:
            f_pkl = open(this_pkl, 'rb')
            lg_pkl = pickle.load(f_pkl)

            #print('Found: ', this_pkl, len(lg_pkl), ' - ', lg_pkl[0].code)
            #print(lg_pkl)
            if lg_pkl[0].code == 'EMPTY':
                    'do nothing'
            else:
                    n_lgs = len(lg_pkl)
                    n_tot_lgs = n_tot_lgs + n_lgs

f_lg.close()
f_mw.close()

vol = 4.3 / 3.14 * np.power(r_select/1000.0, 3.0)
density = n_tot_lgs / (n_files * vol)

print('LG CS density: ', density, ' totLGs found: ', n_tot_lgs, ' out of ', n_files)

if doMoreSelect == True:

    #select = 'zero'
    #select = 'one'
    select = 'two'
    #select = 'three'
    #select = 'four'
    #select = 'five'

    if select == 'zero':
            iso_radius = 2000.
            r_max = 1500.
            r_min = 250.
            m_min = 4.0e+11
            m_max = 5.0e+12
            ratio_max = 4.0
            vrad_max = 10.0


    if select == 'one':
            iso_radius = 2000.
            r_max = 1200.
            r_min = 350.
            m_min = 5.0e+11
            m_max = 4.0e+12
            ratio_max = 3.0
            vrad_max = 0.0

    if select == 'two':
            r_max = 1200.
            r_min = 350.
            m_min = 6.0e+11
            m_max = 4.0e+12
            ratio_max = 3.0
            vrad_max = -25.0

    if select == 'three':
            r_max = 1000.0
            r_min = 400.0
            m_min = 7.0e+11
            m_max = 3.0e+12
            ratio_max = 2.5
            vrad_max = -50.0

    if select == 'four':
            r_max = 900.0
            r_min = 450.0
            m_min = 7.5e+11
            m_max = 2.5e+12
            ratio_max = 2.5
            vrad_max = -75.0

    if select == 'five':
            r_max = 750.0
            r_min = 500.0
            m_min = 8.0e+11
            m_max = 2.0e+12
            ratio_max = 2.0
            vrad_max = -100.0

    n_lgs = 0

    for i_dir in range(i_ini, i_end):
            i_dir_str = '%02d' % i_dir

            for g_dir in range(g_ini, g_end):
                    g_dir_str = '%02d' % g_dir
                    sub_path = i_dir_str + '_' + g_dir_str

                    this_pkl = 'saved/lgs_r_' + rStr + '_mMin' + mMinStr + '_' + res + '_' + sub_path + '.pkl'

                    if os.path.isfile(this_pkl):
                            isPkl = True

                            f_pkl = open(this_pkl, 'rb')
                            lg_pkl = pickle.load(f_pkl)

                            if lg_pkl[0].code == 'EMPTY':
                                    'do nothing'
                            else:
                                    #print(this_pkl, 'FOUND')
                                    for lg in lg_pkl:
                                            condition1 = (lg.LG1.m > m_min and lg.LG1.m < m_max)
                                            condition2 = (lg.LG2.m > m_min and lg.LG2.m < m_max)
                                            condition3 = (lg.r_halos() > r_min and lg.r_halos() < r_max)
                                            condition4 = (lg.v_radial() > vrad_max and lg.r_halos() < r_max)
                                            condition5 = True #(lg.LG1.m / lg.LG2.m < ratio_max)
                                    if (condition1 and condition2 and condition3 and condition4 and condition5):
                                            n_lgs = n_lgs + 1

density = n_lgs / (n_files * vol)
print('LG CS refined density: ', density, ' totLGs found: ', n_lgs, ' out of ', n_files)
'''
