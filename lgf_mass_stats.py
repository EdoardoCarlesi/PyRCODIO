#!/usr/bin/python

from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halos import *
from libcosmo.find_halo import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *
from config import *
import time
import pickle
import os.path
from libcosmo.grid import *

import matplotlib.pyplot as plt

# Simulation & catalog
box_size = 100000.0
'/home/eduardo/CLUES/DATA1/LGF/1024/'
base_path = '/home/eduardo/CLUES/DATA1/LGF/1024/'
#base_path = '/home/eduardo/CLUES/DATA/1024/'
#base_path = '/home/eduardo/CLUES/DATA/4096/37_11/'
#file_single0 = 'snapshot_054.z0.000.AHF_halos'
file_single1 = 'snapshot_055.0000.z0.000.AHF_halos'
#file_single0 = 'snapshot_054.0000.z0.000.AHF_halos'
file_single0 = 'snapshot_054.z0.000.AHF_halos'
#file_single1 = 'snapshot_054.0000.z0.001.AHF_halos'

# Search distance from the box center
radii = [10000.0, 9000.0, 8000.0, 7000.0, 6000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0]
radii_eps = [10100.0, 9100.0, 8100.0, 7100.0, 6100.0, 5100.0, 4100.0, 3100.0, 2100.0, 1100.0]

i_ini = 0
i_end = -1
g_ini = 0
g_end = 30

i_lg = 0
n_lgs_max = 500
rhos = np.zeros((len(radii), n_lgs_max), dtype='float')

rho0 = 45.0

for i_dir in range(i_ini, i_end):
    i_dir_str = '%02d' % i_dir

    for g_dir in range(g_ini, g_end):
        g_dir_str = '%02d' % g_dir
        sub_path = i_dir_str + '_' + g_dir_str
        ahf_file0 = base_path + sub_path + '/' + file_single0
        ahf_file1 = base_path + sub_path + '/' + file_single1
        this_ahf_file = None

        #print(ahf_file0)
        if os.path.isfile(ahf_file0):
            this_ahf_file = ahf_file0
        elif os.path.isfile(ahf_file1):
            this_ahf_file = ahf_file1

        if this_ahf_file != None:
            #print('Reading file : %s' % this_ahf_file)
            all_halos = read_ahf(this_ahf_file)

            n_halos = len(all_halos)

            out_lgs = 'saved/lgs_' + sub_path + '.pkl'

            if os.path.isfile(out_lgs):
                f_out_lgs = open(out_lgs, 'rb')

                #print('Reading LGs from to file %s.' % (out_lgs))
                all_lgs = pickle.load(f_out_lgs)
    
                lg = all_lgs[0]
                i_r = 0

                for rad in radii:
                    if i_r == 0:
                        halos = find_halos_point(lg.geo_com(), all_halos, rad)
                    else:
                        halos = find_halos_point(lg.geo_com(), halos, rad)

                    m_tot_loc = 0.0
                    masses = []

                    for hl in halos:
                        m_tot_loc += hl.m
                        masses.append(hl.m)

                    thisVol = 3.14 * rad * rad * rad * 4. / 3.
                    thisRho = m_tot_loc / thisVol

                    rhos[i_r, i_lg] = thisRho
                    i_r = i_r + 1

                print(i_dir, i_lg, ') Rad: ', rad, ' n_halos: ', len(halos), ' localRho: ', thisRho, thisRho / rho0)
                i_lg = i_lg + 1

if i_end > 0:
    out_rhos = 'saved/lgf_rhos_all.pkl'
    f_out_rhos = open(out_rhos, 'wb')
    pickle.dump(rhos, f_out_rhos)
    f_out_rhos.close()
else:
    out_rhos = 'saved/lgf_rhos_all.pkl'
    f_out_rhos = open(out_rhos, 'rb')
    n_rhos = pickle.load(f_out_rhos)
    f_out_rhos.close()

    i_lg = 0
    while n_rhos[0, i_lg] > 0.0:
        i_lg = i_lg +1

    print(i_lg)
    rhos = np.zeros((len(radii), i_lg), dtype='float')

    for i in range(0, i_lg):
        for j in range(0, len(radii)):
            rhos[j, i] = n_rhos[j, i]

    
    med_vals = np.zeros((len(radii), 3))
    err_vals = np.zeros((2,len(radii)))

    for i_r in range(0, len(radii)):
        this_r = rhos[i_r, :]
        med_vals[i_r, 0] = np.percentile(this_r, 25)
        med_vals[i_r, 1] = np.percentile(this_r, 50)
        med_vals[i_r, 2] = np.percentile(this_r, 65)
        err_vals[0, i_r] = (med_vals[i_r, 1] - med_vals[i_r, 0])/rho0
        err_vals[1, i_r] = (med_vals[i_r, 2] - med_vals[i_r, 1])/rho0
        med_vals[i_r, 1] = np.percentile(this_r, 50)/rho0

        print(i_r, med_vals[i_r, 0]/rho0, med_vals[i_r, 1]/rho0, med_vals[i_r, 2]/rho0)

print(err_vals)
out_rhos = 'saved/rand_rhos_lgs_all_medians.pkl'
f_out_rhos = open(out_rhos, 'rb')

all_med = pickle.load(f_out_rhos)
all_err = np.zeros((2, len(radii)))

rs = np.zeros((len(radii)))
rs_eps = np.zeros((len(radii)))

for i in range(0, len(radii)):
    all_err[0, i] = all_med[i, 1] - all_med[i, 0]
    all_err[1, i] = all_med[i, 2] - all_med[i, 1]
    rs[i] = radii[i] / 1000.0
    rs_eps[i] = radii_eps[i] / 1000.0

print(all_err)

plt.rc({'text.usetex': True})
plt.xlabel(r'$R_{smooth} \quad h^{-1}$ Mpc')
plt.ylabel(r'$<{\rho / \bar{\rho}}>$')
plt.axis([1.5, 10.5, 0.0, 5.0])
plt.errorbar(rs_eps, med_vals[:, 1], yerr=err_vals, fmt='o', color='black')
plt.errorbar(rs, all_med[:, 1], yerr=all_err, fmt='o', color='red')

plt.tight_layout()
plt.savefig('ProjectsPlots/delta_rho_cs_rand.png')






