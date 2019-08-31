#!/usr/bin/python

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
from libcosmo.grid import *
from config import *
import time
import pickle
import os.path

ini_num=0
end_num=0

# Should we read the gadget file and export the particles, or just read the already exported particles?
#doReadSlab = True
doReadSlab = False

# When reading and storing particle data, reduce by this factor
#f_rescale = 1.0 
f_rescale = 4.0 
#f_rescale = 16.0 

#all_lg_base = simu_runs()
#all_lg_base=['00_06']
#all_lg_base=['17_10']
#all_lg_base=['34_13']
all_lg_base=['09_18']

#subruns = ['00']
subruns = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']

max_range=0

for isub in range(0, max_range):
    sub_str = '%02d' % isub
    subruns.append(sub_str)

#box = '500'; resolution = '512'; n_files = 8; withGas = False
#box = '100'; resolution = '4096'; n_files = 1; withGas = False
box = '100'; resolution = '2048'; n_files = 1; withGas = False
#box = '100'; resolution = '1024'; n_files = 1; withGas = True

#box_size = 100.0e+3; plot_side = 10.0e+3; thickn = 5000.0; units = 'kpc'
box_size = 100.0; plot_side = 1.00; thickn = 2.5; units = 'Mpc'

snap_name='snapshot_054'

#fn_lg = 'saved/lgs_00.pkl'
fn_lg = 'saved/lgs_'+resolution+'_'+all_lg_base[0]+'_00.pkl'
f_lg = open(fn_lg, 'rb')
this_lg = pickle.load(f_lg)

#print(this_lg[0].geo_com())
#box_center = [46.6, 50.7, 47.8]
#box_center= this_lg[0].geo_com()
box_center= this_lg.geo_com()

for ip in range(0, 3): 
    box_center[ip] = box_center[ip] / 1000.0

print('Rescaling the plot around the LG position: ', box_center)


for subrun in subruns:
    base_path = '/home/eduardo/CLUES/DATA/'+resolution+'/'+all_lg_base[0]+'/' + subrun + '/'
    #base_path = '/home/oem/CLUES/DATA/'+resolution+'/'+all_lg_base[0]+'/' + subrun + '/'

    # File Names for the output slabs
    fn_0 = base_path + 'slab_xy0_' + str(f_rescale) + '.pkl'
    fn_1 = base_path + 'slab_xy1_' + str(f_rescale) + '.pkl'
    fn_4 = base_path + 'slab_xy4.pkl'

    if doReadSlab == True:
        # Double the side of the slab just in case
        plot_side = plot_side * 2
        print(base_path + snap_name)

        if withGas:
            [x0, y0] = return_slab(base_path + snap_name, 2, box_center, plot_side, thickn, n_files, f_rescale, units, 0)
            # Never rescale star particles
            [x4, y4] = return_slab(base_path + snap_name, 2, box_center, plot_side, thickn, n_files, 1.0, units, 4)
            f_0 = open(fn_0, 'wb')
            pickle.dump([x0, y0], f_0)
            f_4 = open(fn_4, 'wb')
            pickle.dump([x4, y4], f_4)
            print('Dumping to files: ', fn_0, fn_4)

        [x1, y1] = return_slab(base_path + snap_name, 2, box_center, plot_side, thickn, n_files, f_rescale, units, 1)
        f_1 = open(fn_1, 'wb')
        pickle.dump([x1, y1], f_1)

        print('Dumping to files: ', fn_0, fn_1, fn_4)

        f_1.close()

        if withGas:
            f_0.close()
            f_4.close()

    else:
        # DM
        slab = [fn_1]; ptype = 1
        bw_smooth = 0.25; nbins = 750
        f_out = 'ProjectsPlots/lg_rhos_' + resolution + '_' + all_lg_base[0] + '_' + subrun + '_grid' + str(nbins)
        simple_plot_rho(box_center, plot_side, f_out, nbins, f_rescale, thickn, units, slab, bw_smooth, ptype)

        # gas and stars
        if withGas:
            slab = [fn_0, fn_4]; ptype = 0
            #bw_smooth = 0.025; nbins = 256
            bw_smooth = 0.1; nbins = 64
            simple_plot_rho(box_center, plot_side, f_out, nbins, f_rescale, thickn, units, slab, bw_smooth, ptype)
