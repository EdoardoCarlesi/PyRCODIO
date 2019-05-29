'''
        This script is a wrapper for computing the following:
                - Find an LG candidate at each step including nearby satellites
                - For each halo (main LGs and satellites) trace its merging history
                - At each step, re-identify the subhalos for each candidate (simply taking all the halos within a factor * Rvir)
                - Use the subhalo distribution at each step to re-compute the anisotropy & inertia tensor eigenvalues
                - Do not do any plot, simply dump all the stuff to .txt - the analysis will be done separately
'''

import numpy as np
import os
import pickle

from config import *
from libio.read_ascii import *
from libcosmo.utils import *
from libcosmo.track_halos import *
from libcosmo.find_halo import *
#from libcosmo.lg_plot import *

hubble = 0.67
# Local group selection parameters
center = [50000., 50000., 50000.]

# Allocate LG Models
#resolution = '2048'
resolution = '4096'
#this_code = '00_06'
#this_code = '37_11'
#this_code = '09_18'
this_code = '45_17'
lg_dummy = LocalGroup(this_code)

base_path = '/z/carlesi/HestiaNoam/RE_SIMS/'+resolution+'/GAL_FOR/'+this_code+'/AHF_output/'

#snap_file = 'HESTIA_100Mpc'+this_code+'_127.z0.000.AHF_halos'
snap_file = 'HESTIA_100Mpc_'+resolution+'_'+this_code+'.127.z0.000.AHF_halos'

# Subhalo identification criterion
fac_r = 1.5     # This is used for the global R around the LG as well as for Rvir around each halo
min_common = 15
part_min = 1000.
r_select = 7500.
mmin = 1.e+6    # Track haloes above this threshold at z=0
mcut = 0.5e+10

n_lg_good = 0
main_ids =[]; this_file_halo = base_path + '/' + snap_file

print('Reading files: ', this_file_halo)
(all_lg_models, hash_run) = lg_models()
halos = read_ahf(this_file_halo)
select_halos = find_halos_mass_radius(center, halos, r_select, mmin)

lg_model = all_lg_models[hash_run[this_code]]
this_lg = find_lg(select_halos, lg_model)
n_lgs = int(len(this_lg))

print('Found a total of %d LG pairs.' % (n_lgs))
good_lgs = 0
if n_lgs > 0:
# If there are more candidates we need to find the right one
    rating = 1000.0

    for ind in range(0, n_lgs):
        lg = this_lg[ind]
        lg.c_box = center #settings.box_center

        if lg.rating() < rating and (lg.LG1.npart > part_min) and (lg.LG2.npart > part_min):
            good_lgs += 1
            rating = lg.rating
            best_lg = lg

        if good_lgs > 0:
            print('Best LG: ', best_lg.info());
            n_lg_good += 1;
	    run_num = this_code

            # Identify the main halos to be tracked at all zs
            com = best_lg.get_com()
            rad_lg  = best_lg.r_halos() * fac_r
            rad_mw  = best_lg.LG1.r * fac_r
            rad_m31 = best_lg.LG2.r * fac_r
            lg_halos  = find_halos_mass_radius(com, halos, rad_lg, mmin)
            mw_halos  = find_halos_mass_radius(best_lg.LG1.x, halos, rad_mw, mmin)
            m31_halos = find_halos_mass_radius(best_lg.LG2.x, halos, rad_m31, mmin)
                
            # Print the subhalo list
	    fname_mw = '_MW_'
	    fname_m31 = '_M31_'
	    fname_lg = '_LG_'
            print_subhalos(com, 10. * mcut, lg_halos,  run_num, fname_lg)
            print_subhalos(best_lg.LG1.x, mcut, mw_halos,  run_num, fname_mw)
            print_subhalos(best_lg.LG2.x, mcut, m31_halos, run_num, fname_m31)

print('Found ', n_lg_good, 'viable local group pairs.')
