from libio.read_ascii import *
from libcosmo.find_halo import *
from libcosmo.halos import *
from libcosmo.grid import *
import pickle
import numpy as np
import os

runNum = '00'
file_ahf = '/home/eduardo/CLUES/DATA/FullBox/catalogs/' + runNum + '/snapshot_054.z0.000.AHF_halos'

box=100.0
nodes=20
radius=5.0

mpc2kpc=1.e+3
mMin = 1.e+8
mMax = 1.e+16

thrRho = 0.9

print('Reading: ', file_ahf)
#allHalos = read_ahf(file_ahf)
ahfGrid = Grid(nodes, box * mpc2kpc)
pklGrid = 'saved/grid_' + str(nodes) + '_' + str(mMin) + '.pkl'

doLgs = True

if os.path.isfile(pklGrid):
    allHalos = []
    fGrid = open(pklGrid, 'rb')
    ahfGrid = pickle.load(fGrid)
    fGrid.close()
    mTot = 4.027e+16   
else:
    allHalos = read_ahf_mass_range(file_ahf, mMax, mMin)
    print('Done. Found', len(allHalos), ' halos.')

    mTot = 0
    for h in allHalos:
        mTot = mTot + h.m

        thisIndex = ahfGrid.phys2grid(h.x)
        ahfGrid.rho[thisIndex[0], thisIndex[1], thisIndex[2]] = ahfGrid.rho[thisIndex[0], thisIndex[1], thisIndex[2]] + h.m

rho0 = mTot / np.power(box, 3)
print('Halo mass density: ', mTot, ' rho0: ', rho0)

fGrid = open(pklGrid, 'wb')
pickle.dump(ahfGrid, fGrid)

nodeVol = np.power(box/nodes, 3)
nUnder = 0; nOver = 0;
for ix in range(0, nodes):
    for iy in range(0, nodes):
        for iz in range(0, nodes):
            ahfGrid.rho[ix, iy, iz] = ahfGrid.rho[ix, iy, iz] / nodeVol / rho0
            thisRho = ahfGrid.rho[ix, iy, iz] 

            if thisRho < thrRho:
                nUnder = nUnder + 1
            else:
                nOver = nOver + 1

print('nUnder = ', nUnder, ' nOver = ', nOver, ' nTot = ', np.power(nodes, 3), ' %: ', float(nOver)/np.power(nodes, 3))

nUnder = 0; nOver = 0;
if doLgs == True:
    #lgPkl = 'saved/rand_lgs_' + runNum + '.pkl'
    lgPkl = 'saved/rand_select_lgs_' + runNum + '.pkl'
    fLG = open(lgPkl, 'rb')
    allLGs = pickle.load(fLG)
    nLGs = len(allLGs)

    for lg in allLGs:
        xyz = lg.geo_com()
        #print(xyz)
        ixyz = ahfGrid.phys2grid(xyz)
        
        thisRho = ahfGrid.rho[ixyz[0], ixyz[1], ixyz[2]]
        if thisRho < thrRho:
            nUnder = nUnder + 1
        else:
            nOver = nOver + 1

print('Overd: ', nOver, ', Underd: ', nUnder, ' nTot: ', nLGs)

#        print(xyz, ixyz)
'''
        WARNING: this should be the exact way of computing it, except it's too computationally expensive.
        The one above is a good approximation

nodeVol = np.power(radius, 3) * 3.14 * 4.0 / 3.0
thisNode = np.zeros((3))

for ix in range(0, nodes):
    thisNode[0] = ix * ahfGrid.cell * mpc2kpc

    for iy in range(0, nodes):
        thisNode[1] = iy * ahfGrid.cell * mpc2kpc

        for iz in range(0, nodes):
            thisNode[2] = iz * ahfGrid.cell * mpc2kpc

            nodeHalos = find_halos_mass_radius(thisNode, allHalos, radius * mpc2kpc, 0.0)

            mNode = 0
            for nH in nodeHalos:
                mNode = nH.m + mNode

            dNode = (mNode / nodeVol)
            dNode0 = dNode / rho0

            ahfGrid.rho[ix, iy, iz] = dNode0

            print('Node(', ix, iy, iz, ') dens: ', dNode0, \
                    ' nHalos: ', len(nodeHalos), ' mass: ', mNode, ' node: ', thisNode, ' nVol: ', nodeVol)
'''
