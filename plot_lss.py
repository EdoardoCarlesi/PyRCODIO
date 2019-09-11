import numpy as np
import pylab
import matplotlib
import os
from matplotlib import ticker, cm
from libcosmo.grid import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
import pickle

# First read the original GADGET file then 
#loadFile=True
loadFile=False

'''
simutype='CF3'

#codes = ['70_00','71_00','72_00']
codes = ['70_00']
simnames = ['CF3_YH_h78', 'CF3_RG_1500', 'CF3_RG_1000', 'CF3_YH_v4']
snapname='snapshot_127'; kpcU=1.0; c=250.0*kpcU; n_files=8; box = '500'


simutype='CF2P5'
simnames = ['CF2P5_YH']
codes = ['59050', '58151', '57252', '56353', '55454']
snapname='snapshot_019'; kpcU=1.e+3; c=335.0*kpcU; n_files=8; box = '670'
'''

#simutype='CF3'
#simnames = ['CF3_YH_BGc']
#codes = ['59050']
#snapname='snapshot_019'; kpcU=1.e+3; c=500.0*kpcU; n_files=8; box = '1000'

simutype='CF3'
simnames = ['CF3_BGc_256']
codes = ['59050', '58151', '57252', '56353'] #, '55454']
snapname='snapshot_019'; kpcU=1.e+3; c=250.0*kpcU; n_files=8; box = '500'

# Compress information: take one out of reduce_fac particles
#reduce_fac = 32    
#reduce_fac = 16    
reduce_fac = 2

# Plot Properties
#gridSize = 512
#gridSize = 256
gridSize = 128
#gridSize = 64
#gridSize = 32

# Smoothing length
#smooth = 0.05
smooth = 0.05

# Thickness
thickn = 5.0 * kpcU

# Only plot a sub-square of the full slide
#side = 200.0
#side = 150.0
side = 250.0 * kpcU

center=[c, c, c]; 
#center=[47, 50, 47]; 
units='kpc'
#units='Mpc'


if loadFile:
    for code in codes:
        for simname in simnames:

            # Reset variables
            i_tmp = 0
            data_x = []; data_y = []
            print('Looping on ', code, simname)
            
            base_dir = '/home/eduardo/CLUES/DATA/' + simutype + '/' + box + '/' + simname + '/' + code + '/'
            slab_x_fname = base_dir + 'slab_x_fac'+str(reduce_fac)+'.pkl'
            slab_y_fname = base_dir + 'slab_y_fac'+str(reduce_fac)+'.pkl'
            snap = base_dir + snapname

            if os.path.exists(slab_x_fname) and os.path.exists(slab_y_fname):
                loop_n_files = -1   # Skip the loop on reading the files
                print('File: ', slab_x_fname, ' found, skipping to read and compress the next one...')
            else:
                # Do the loop normally
                loop_n_files = n_files

            for i in range (0, loop_n_files):
                snap_loc = snap + '.' + str(i)
                print(snap_loc)

                # Here we put the side = c ---> no zoom, we save all the particles in the slab
                tmp_x, tmp_y = return_slab(snap_loc, 2, center, c, thickn, 1, reduce_fac, units, 1)

                for i in range(0, len(tmp_x)):
                    x = tmp_x[i]; y = tmp_y[i]
                    x = x - c
                    y = y - c

                    if abs(x) < c and abs(y) < c:
                        x = x/kpcU
                        y = y/kpcU
                        data_x.append(x)
                        data_y.append(y)
                        i_tmp += 1

                print('Found a total of ', i_tmp, ' particles. ')

            print('Saving files: ', slab_x_fname)
            x_out = open(slab_x_fname, 'wb')
            y_out = open(slab_y_fname, 'wb')
            pickle.dump(data_x, x_out)
            pickle.dump(data_y, y_out)
            print('Done.')

# Use the other routines for plotting
else:
    for code in codes:
        for simname in simnames:

            print('Looping on ', code, simname)
            base_dir = '/home/eduardo/CLUES/DATA/' + simutype + '/' + box + '/' + simname + '/' + code + '/'
            
            slab_x_fname = base_dir + 'slab_x_fac' + str(reduce_fac) + '.pkl'
            slab_y_fname = base_dir + 'slab_y_fac' + str(reduce_fac) + '.pkl'
            snap = base_dir + snapname

            print('Loading files: ', slab_x_fname)
            slabs = [slab_x_fname, slab_y_fname]

            figname = 'ProjectsPlots/' + snapname + '_' + simname + '_' + code + '_' \
                    + str(gridSize) + '_sm' + str(smooth) + '_fac' + str(reduce_fac) 

            print('Dumping to file: ', figname + '_dm.png')

            try:
                simple_plot_dm(center, side, figname, gridSize, reduce_fac, thickn, units, slabs, smooth, 2)
            except:
                print('File: ', slabs, ' not found.')
