#import yt
#import yt.units as units
#from lenstools.simulations import Gadget2Snapshot
import numpy as np
import pylab
import matplotlib
from matplotlib import ticker, cm
from libcosmo.grid import *
from libcosmo.lg_plot import *
from libio.read_ascii import *
import pickle
import seaborn as sns

# Should we read the GADGET original file or not?
#loadFile=True
loadFile=False
oldLoad=False


#simname='_LG_'
#base_dir='/home/eduardo/CLUES/DATA/01/'; snapname='snapshot_035'; kpcU=1.0; c = 50.0 * kpcU; n_files=1
#base_dir='/home/oem/CLUES/DATA/1024/37_11/00/'; snapname='snapshot_054'; kpcU=1.0; c = 50.0 * kpcU; n_files=4

#simname='_YH_cf2p5_'
#code='59050';base_dir='/home/eduardo/CLUES/DATA/CF2P5/670/512/'+code+'/';snapname='snapshot_019';kpcU=1.e+3;c=335.0*kpcU;n_files=8
#code='58151';base_dir='/home/eduardo/CLUES/DATA/CF2P5/670/512/'+code+'/';snapname='snapshot_019';kpcU=1.e+3;c=335.0*kpcU;n_files=8
#code='57252';base_dir='/home/eduardo/CLUES/DATA/CF2P5/670/512/'+code+'/';snapname='snapshot_019';kpcU=1.e+3;c=335.0*kpcU;n_files=8
#code='56353';base_dir='/home/eduardo/CLUES/DATA/CF2P5/670/512/'+code+'/';snapname='snapshot_019';kpcU=1.e+3;c=335.0*kpcU;n_files=8
#code='55454';base_dir='/home/eduardo/CLUES/DATA/CF2P5/670/512/'+code+'/';snapname='snapshot_019';kpcU=1.e+3;c=335.0*kpcU;n_files=8

#simname='_YH_h74_'
#code='70_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_v4/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='71_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_v4/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='72_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_v4/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='73_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_v4/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='74_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_v4/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8

simname='_YH_h78_'
code='70_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_h78/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='71_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_h78/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='72_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_YH_h78/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8

#simname='_RG1500_'
#code='70_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_RG_1500/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='71_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_RG_1500/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='72_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_RG_1500/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8

#simname='_RG1000_'
#code='70_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_RG_1000/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='71_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_RG_1000/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8
#code='72_00';base_dir='/home/eduardo/CLUES/DATA/CF3/500/CF3_RG_1000/'+code+'/';snapname='snapshot_127';kpcU=1.0;c=250.0*kpcU;n_files=8

#simname = '_LGss_'
#code='00';base_dir='/home/eduardo/CLUES/DATA/1024/'+code+'/';snapname='snapshot_054';kpcU=1.0;c=50.0*kpcU;n_files=4


snap=base_dir + snapname

# take one out of reduce_fac particles
#reduce_fac = 128
#reduce_fac = 64   
reduce_fac = 32    
#reduce_fac = 16    
#reduce_fac = 8

slab_x_fname = base_dir + 'slab_x_fac'+str(reduce_fac)+'.pkl'
slab_y_fname = base_dir + 'slab_y_fac'+str(reduce_fac)+'.pkl'
#slab_x_fname = base_dir + 'slab_x.pkl'
#slab_y_fname = base_dir + 'slab_y.pkl'

# Plot Properties
#gridSize = 512
#gridSize = 256
gridSize = 128
#gridSize = 64
#gridSize = 32
#smooth = 0.05
smooth = 0.05

thickn = 5.0 * kpcU

# Only plot a sub-square of the full slide
#side = 200.0
#side = 150.0
side = 250.0 * kpcU
#side = c / kpcU

doBar=False
#doBar=True
n_levels = 10

# Small box settings
#contour = np.logspace(-6.0, -1.0, num=20)

# Full box
contour = np.logspace(-6.5, -4.0, num=8)
#[1.e-8, 1.e-7, 1.e-6, 1.e-5, 1.e-4, 1.e-3, 1.e-2]
#print(contour)


center=[c, c, c]; 
#center=[47, 50, 47]; 
units='kpc'
#units='Mpc'

data_x = []; data_y = []
i_tmp = 0

if loadFile:
    for i in range (0, n_files):
        snap_loc = snap + '.' + str(i)
        print(snap_loc)

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
    
    side = c

elif oldLoad:
    print(snap)
    print('Loading files: ', slab_x_fname)
    x_out = open(slab_x_fname, 'rb')
    y_out = open(slab_y_fname, 'rb')
    all_data_x = pickle.load(x_out)
    all_data_y = pickle.load(y_out)

    data_x = [];     data_y = []

    for ix in range(0, len(all_data_x)):
        dx = all_data_x[ix]
        dy = all_data_y[ix]

        if abs(dx) < side and abs(dy) < side:
            data_x.append(dx)
            data_y.append(dy)

    print('Done. Selected ', len(data_x), ' particles out of ', len(all_data_x))

    figname=snapname+simname+'_'+code+'_' + str(gridSize) + '_sm' + str(smooth) + '_fac' + str(reduce_fac) + '_side' + str(side) + '.png'

    plt.figure(figsize=(20, 20))
    matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel("SGY (Mpc/h)")
    plt.xlabel("SGX (Mpc/h)")

    print('Smoothing out plot on a ', gridSize, ' grid with a ', smooth, ' Gaussian kernel.')
    sns.kdeplot(data_x, data_y, cmap="coolwarm", shade=True, shade_lowest=True, gridsize=gridSize,
        bw=smooth, levels=n_levels, cbar=doBar)
    # SCATTERPLOT just in case
    #plt.plot(data_x, data_y, linestyle='', marker='o', markersize=0.1)
    #plt.savefig('test_scatter.png', dpi=100)
    plt.tight_layout()
    print('Saving output to file: ', figname)
    plt.savefig(figname, dpi=300)

# Use the other routines for plotting
else:
    print(snap)
    print('Loading files: ', slab_x_fname)
    slabs = [slab_x_fname, slab_y_fname]
    figname='ProjectsPlots/'+snapname+simname+'_'+code+'_' + str(gridSize) + '_sm' + str(smooth) + '_fac' + str(reduce_fac) 
    
    simple_plot_dm(center, side, figname, gridSize, reduce_fac, thickn, units, slabs, smooth, 2)
