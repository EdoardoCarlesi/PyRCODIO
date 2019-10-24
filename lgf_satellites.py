import glob
import numpy as np
import os

from config import *
from libio.read_ascii import *
from libcosmo.utils import *
from libcosmo.halos import *
from libcosmo.find_halo import *
from libcosmo.lg_plot import *

#resolution='1024'
#resolution='2048'
resolution='4096'

subrun_init = 0
subrun_end = 10

base_path = '/home/eduardo/CLUES/DATA/trees/'
z_file = '/home/eduardo/CLUES/DATA/4096/zs.txt'

nSteps = 54
#do_plots = "true"
#do_plots = "false"

# Allocate LG Models
if resolution == '2048':
    all_runs = simu_runs()
elif resolution == '4096':
    all_runs = ['37_11']

run_init = 0
run_end = len(all_runs)

sub_path = base_path + resolution + '/'

for run_j in range(run_init, run_end):
    x_mf1 = [];     y_mf1 = []
    x_mf2 = [];     y_mf2 = []
    n_mf = 0

    this_code = all_runs[run_j]

    for subrun_i in range(subrun_init, subrun_end):
        run_num = '%02d' % subrun_i
        this_path = sub_path + this_code + '/' + run_num + '/'

        this_m31 = sub_path + this_code + '/lgs_' + resolution + '_' + this_code + '_' + run_num + '_subs_M31.txt'
        this_m31_pkl = 'saved/lgs_' + resolution + '_' + this_code + '_' + run_num + '_subs_M31.pkl'
        this_mw = sub_path + this_code + '/lgs_' + resolution + '_' + this_code + '_' + run_num + '_subs_M31.txt'
        this_mw_pkl = 'saved/lgs_' + resolution + '_' + this_code + '_' + run_num + '_subs_M31.pkl'
        this_lg = sub_path + this_code + '/lgs_' + resolution + '_' + this_code + '_' + run_num + '_subs_M31.txt'
        this_lg_pkl = 'saved/lgs_' + resolution + '_' + this_code + '_' + run_num + '_subs_M31.pkl'
        
        if os.path.exists(this_m31):
            print('Reading in AHF file: ', this_m31)
            head = read_header(this_m31)
            headID = str(head[3].rstrip('\n'))
            #            print(ID)
            #halos_m31 = read_ahf(this_m31)
           #print('Found ', len(halos_m31), ' subhalos.')
            thisTree = this_path + 'halo_' + headID + '.allinfo'
            hostM31 = HaloThroughZ(nSteps)
            hostM31.load_files(thisTree, z_file)
        
            for halo in halos_m31[0:2]:
                thisTree = this_path + 'halo_' + str(halo.ID) + '.allinfo'
                thisHalo = HaloThroughZ(nSteps)

                if os.path.exists(thisTree):
                    print(thisTree, ' found.')
                    thisHalo.load_files(thisTree, z_file)
#                    print(thisHalo.m_t())
                else:
                    print(thisTree, ' not found.')
    
            '''
            print(halos_m31[0].info())

            #these_sub2 = find_halos(best_lg.LG2, ahf_all, fac_r * best_lg.LG2.r)
            #subs2 = SubHalos(best_lg.LG2, these_sub2)

            (x_m, y_n) = subs1.mass_function()
            x_mf1.append(x_m);      y_mf1.append(y_n)

            (x_m, y_n) = subs2.mass_function()
            x_mf2.append(x_m);      y_mf2.append(y_n)
            n_mf += 1

    # Do some plots at the end
    if do_plots == "true":
            plot_lg(this_file_gad, file_png_name, best_lg.LG1, best_lg.LG2, reduce_fac, 1, plot_pos)

            if (subrun_i == subrun_end-1):
                    plot_massfunctions(x_mf1, y_mf1, n_mf, file_png_mfs1)
                    plot_massfunctions(x_mf2, y_mf2, n_mf, file_png_mfs2)
                    del x_m ;       del y_n
                    del x_mf1;      del x_mf2
                    del y_mf1;      del y_mf2
                    n_mf = 0
            '''
