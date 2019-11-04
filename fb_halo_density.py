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
nodes=10

mpc2kpc=1.e+3
mMin = 1.e+8
mMax = 1.e+16

#thrRho = 0.78
thrRho = 0.5

hMmin = 0.4e+12
hMmax = 5.0e+12

print('Reading: ', file_ahf)
pklGrid = 'saved/grid_' + str(nodes) + '_' + str(mMin) + '.pkl'
pklHalo = 'saved/halo_' + str(nodes) + '_' + str(mMin) + '.pkl'

doLgs = True
#doLgs = False
doHalos = True

if os.path.isfile(pklGrid):
#if False:
    allHalos = []
    fGrid = open(pklGrid, 'rb')
    ahfGrid = pickle.load(fGrid)
    fGrid.close()
    mTot = 4.027e+16   
else:
    allHalos = read_ahf_mass_range(file_ahf, mMax, mMin)
    ahfGrid = Grid(nodes, box * mpc2kpc)
    #allHalos = read_ahf(file_ahf)

    print('Done. Found', len(allHalos), ' halos.')

    mTot = 0
    nTot = 0
    mHalos = []
    for h in allHalos:
        mTot = mTot + h.m

        if h.m > hMmin and h.m < hMmax:
            nTot = nTot + 1
            mHalos.append(h)

        thisIndex = ahfGrid.phys2grid(h.x)
        ahfGrid.rho[thisIndex[0], thisIndex[1], thisIndex[2]] = ahfGrid.rho[thisIndex[0], thisIndex[1], thisIndex[2]] + h.m
    
    print('Saving to file ', pklHalo, ' mHalos: ', len(mHalos))
    fHalo = open(pklHalo, 'wb')
    pickle.dump(mHalos, fHalo)
    fHalo.close()

    fGrid = open(pklGrid, 'wb')
    pickle.dump(ahfGrid, fGrid)
    fGrid.close()

rho0 = mTot / np.power(box, 3)
print('Halo mass density: ', mTot, ' rho0: ', rho0)
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

fracVol = float(nOver)/np.power(nodes, 3)
print('Tot nUnder = ', nUnder, ' nOver = ', nOver, ' nTot = ', np.power(nodes, 3), ' fracVol: ', fracVol) 

mMin = 0.4e+12; mMax = 5.0e+12; v_max = 25; r_min = 250; r_max = 1500
#mMin = 0.45e+12; mMax = 4.0e+12; v_max = 0; r_min = 300; r_max = 1300
#mMin = 0.5e+12; mMax = 3.0e+12; v_max = -25; r_min = 350; r_max = 1000
#mMin = 0.55e+12; mMax = 2.5e+12; v_max = -50; r_min = 400; r_max = 900
#mMin = 0.6e+12; mMax = 2.0e+12; v_max = -75.; r_min = 450; r_max = 800
#mMin = 0.65e+12; mMax = 1.5e+12; v_max = -100; r_min = 500; r_max = 700

nUnder = 0; nOver = 0;
if doLgs == True:
    #lgPkl = 'saved/rand_lgs_' + runNum + '.pkl'
    lgPkl = 'saved/rand_select_lgs_' + runNum + '.pkl'
    fLG = open(lgPkl, 'rb')
    allLGs = pickle.load(fLG)
    nLGs = len(allLGs)

    for lg in allLGs:
        xyz = lg.geo_com()
        ixyz = ahfGrid.phys2grid(xyz)

        thisRho = ahfGrid.rho[ixyz[0], ixyz[1], ixyz[2]]
        condition = (lg.LG1.m > mMin and lg.LG1.m < mMax and lg.LG2.m > mMin 
                    and lg.LG2.m < mMax and lg.v_radial() < v_max and lg.r_halos() > r_min and lg.r_halos() < r_max) 
        #print(xyz)
        #condition = True 
        #        print(xyz, ixyz)
        #if thisRho < thrRho :
            #print(thisRho)

        if (thisRho > thrRho) and condition:
            nOver = nOver + 1
        else:
            nUnder = nUnder + 1

#    fracVol = float(nOver)/np.power(nodes, 3)
    subVol = fracVol * np.power(box, 3)
    print('LG Overd: ', nOver, ', Underd: ', nUnder, ' nTot: ', nLGs, ' newDens: ', nOver / subVol )





mMin = 0.4e+12; mMax = 5.0e+12; v_max = 25; r_min = 250; r_max = 1500
mMin = 0.45e+12; mMax = 4.0e+12; v_max = 0; r_min = 300; r_max = 1300
#mMin = 0.5e+12; mMax = 3.0e+12; v_max = -25; r_min = 350; r_max = 1000
#mMin = 0.55e+12; mMax = 2.5e+12; v_max = -50; r_min = 400; r_max = 900
#mMin = 0.6e+12; mMax = 2.0e+12; v_max = -75.; r_min = 450; r_max = 800
#mMin = 0.65e+12; mMax = 1.5e+12; v_max = -100; r_min = 500; r_max = 700

nUnder = 0; nOver = 0;
if doHalos == True:
    print('Opening file ', pklHalo)
    fHalos = open(pklHalo, 'rb')
    allHalos = pickle.load(fHalos)
    nHalos = len(allHalos)
    print('halos: ', nHalos)

    for h in allHalos:
        xyz = h.x
        ixyz = ahfGrid.phys2grid(xyz)

        thisRho = ahfGrid.rho[ixyz[0], ixyz[1], ixyz[2]]
        condition = (h.m > mMin and h.m < mMax)
    
        if (thisRho > thrRho) and condition:
            nOver = nOver + 1
        else:
            nUnder = nUnder + 1

#    fracVol = float(nOver)/np.power(nodes, 3)
    subVol = fracVol * np.power(box, 3)
    print('AllHalos Overd: ', nOver, ', Underd: ', nUnder, ' nTot: ', nHalos, ' newDens: ', nOver / subVol )


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
