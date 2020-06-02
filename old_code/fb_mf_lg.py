from libio.read_ascii import *
from libcosmo.find_halo import *
from libcosmo.halos import *
from libcosmo.grid import *
import pickle
import numpy as np
import os

radMF = 5.e+3

print('File read in, looking for MFs')
for i in range(1, 5):
    runNum = '%02d' % i
    file_ahf = '/home/eduardo/CLUES/DATA/FullBox/catalogs/' + runNum + '/snapshot_054.z0.000.AHF_halos'
    pklMFs = 'saved/rand_mfs_' + runNum + '.pkl'
    lgPkl = 'saved/rand_select_lgs_' + runNum + '.pkl'

    print('Reading: ', file_ahf)
    allHalos = read_ahf(file_ahf)

    allMFs = []
    fLG = open(lgPkl, 'rb')
    allLGs = pickle.load(fLG)
    nLGs = len(allLGs)

    print('Will write to:', pklMFs)

    for lg in allLGs:
        xyz = lg.geo_com()
        radSubs = find_halos_mass_radius(xyz, allHalos, radMF, 0.0)
        
        masses = []
        for h in radSubs:
            masses.append(h.m)

        thismf = mass_function(masses)
        allMFs.append(thismf)

    print('Saving mass functions to file: ', pklMFs)
    fPkl = open(pklMFs, 'wb')
    pickle.dump(allMFs, fPkl)
#    print(allMFs)


