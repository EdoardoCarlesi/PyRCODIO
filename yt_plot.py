import yt
import numpy as np
import yt.units as units
import pylab

#fname = '/home/eduardo/CLUES/DATA/1024/00/snapshot_054.0'
fname = '/home/eduardo/CLUES/DATA/HESTIA/2048/37_11/snapshot_127.0.hdf5'
#fname = 'GadgetDiskGalaxy/snapshot_200.hdf5'

unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}

#bbox_lim = 100.0 #kpc
bbox_lim = 100.0 #Mpc/h

'''
bbox = [[-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim]]
'''

bbox = [[0,bbox_lim],
        [0,bbox_lim],
        [0,bbox_lim]]

ds = yt.load(fname,unit_base=unit_base,bounding_box=bbox, index_ptype="PartType0")
#ds = yt.load(fname,unit_base=unit_base,bounding_box=bbox, index_ptype="Gas")
px = yt.ProjectionPlot(ds, 'x', ('gas', 'density'))
px.show()
#ds = yt.load(fname)
#ds.index
#ad= ds.all_data()

