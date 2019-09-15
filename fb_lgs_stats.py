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
import time
import pickle
from libcosmo.grid import *

# Simulation & catalog
base_path = '/home/eduardo/CLUES/DATA/FullBox/catalogs/'
file_single='snapshot_054.z0.000.AHF_halos'
box_size = 100000.0

radii = [10000.0, 9000.0, 8000.0, 7000.0, 6000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0]

i_ini = 0
i_end = -1

rho0 = 40.0 # Units of Msun over kpc ^ 3
print('Mean density: ', rho0)

volume = box_size * box_size * box_size

rand_mf_x = []
rand_mf_y = []

for i_dir in range(i_ini, i_end):
        m_tot = 0.0
        sub_path = '%02d' % i_dir
        out_lgs = 'saved/rand_web_lgs_' + sub_path + '.pkl'
        print('Loading file: ', out_lgs)
        f_out_lgs = open(out_lgs, 'rb')
        all_lgs = pickle.load(f_out_lgs)
        n_lgs = len(all_lgs)
        print('Found ', n_lgs, ' LGs.') 

        this_ahf_file = base_path + sub_path + '/' + file_single

        rhos = np.zeros((len(radii), n_lgs), dtype='float')

        print('Reading file : %s' % this_ahf_file)
        all_halos = read_ahf(this_ahf_file) 
        n_halos = len(all_halos)
        print(n_halos, ' halos have been read in, computing densities...') 

        i_lg = 0
        for lg in all_lgs:
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
            i_lg = i_lg + 1; 

        out_rhos = 'saved/rand_rhos_lgs_' + sub_path + '.pkl'
        f_out_rhos = open(out_rhos, 'wb')
        pickle.dump(rhos, f_out_rhos)
        f_out_rhos.close()

i_ini = 0
i_end = 5

med_vals = np.zeros((i_end, len(radii), 3), dtype='float')

for i_dir in range(i_ini, i_end):
    sub_path = '%02d' % i_dir
    out_rhos = 'saved/rand_rhos_lgs_' + sub_path + '.pkl'
    f_out_rhos = open(out_rhos, 'rb')
    rhos = pickle.load(f_out_rhos)

    for i_r in range(0, len(radii)):
        this_r = rhos[i_r, :]
        med_vals[i_dir, i_r, 0] = np.percentile(this_r, 25)
        med_vals[i_dir, i_r, 1] = np.percentile(this_r, 50)
        med_vals[i_dir, i_r, 2] = np.percentile(this_r, 75)
        #print(i_dir, i_r, med_vals[i_dir, i_r, 1]/rho0)

sub_path = '%02d' % i_dir
out_rhos = 'saved/rand_rhos_lgs_all_medians.pkl'
f_out_rhos = open(out_rhos, 'wb')

all_med = np.zeros((len(radii), 3), dtype='float')

for i_r in range(0, len(radii)):
    rad = radii[i_r]
    all_med[i_r, 0] = np.average(med_vals[:, i_r, 0]/rho0)
    all_med[i_r, 1] = np.average(med_vals[:, i_r, 1]/rho0)
    all_med[i_r, 2] = np.average(med_vals[:, i_r, 2]/rho0)
    print(rad, np.average(med_vals[:, i_r, 0]/rho0), np.average(med_vals[:, i_r, 1]/rho0), np.average(med_vals[:, i_r, 2]/rho0))
    #print(rad, (med_vals[:, i_r, 1]/rho0))

pickle.dump(all_med, f_out_rhos)


