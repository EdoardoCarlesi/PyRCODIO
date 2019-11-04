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

mpc2kpc=1.e+3
mMin = 1.e+8
mMax = 1.e+16
radMF = 5.e+3

hMmin = 0.4e+12
hMmax = 5.0e+12

print('Reading: ', file_ahf)
allHalos = read_ahf(file_ahf)

pklMFs = 'saved/rand_mfs_' + runNum + '.pkl'
doLGs = True
allMFs = []

print('File read in, looking for MFs')
if doLGs == True:
    lgPkl = 'saved/rand_select_lgs_' + runNum + '.pkl'
    fLG = open(lgPkl, 'rb')
    allLGs = pickle.load(fLG)
    nLGs = len(allLGs)

    for lg in allLGs:
        xyz = lg.geo_com()
        radSubs = find_halos_mass_radius(xyz, allHalos, radMF, 0.0)
        
        masses = []
        for h in radSubs:
            masses.append(h.m)

        thismf = mass_function(masses)
#        print(len(thismf[0]))

    allMFs.append(thismf)

print('Saving mass functions to file: ', pklMFs)
fPkl = open(pklMFs, 'wb')
pickle.dump(allMFs, fPkl)
