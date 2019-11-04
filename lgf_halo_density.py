from libcosmo.halo import *
import numpy as np
import pickle
import os

iIni=0
iEnd=80
gIni=0
gEnd=30

rhoZero=4.027e+16 / np.power(100.0, 3)
volZero=125.0 * 4.0 * 3.14 / 3.0

allRhos = []

for i in range(iIni, iEnd):
    iStr = "%02d" % i
    for g in range(gIni, gEnd):
        gStr = "%02d" % g
        runStr = iStr + '_' + gStr

        #pklFile = 'saved/lgs_r_5000.0_mMin4e+11_512_' + runStr + '.pkl'
        pklFile = 'saved/lgs_center_5000.0_' + runStr + '.pkl'
        #'saved/lgs_center_5000.0_65_01.pkl'

        if os.path.isfile(pklFile):
            #fPkl = open(pklFile, 'r') #, encoding='latin1')
            #halos = pickle.load(fPkl, encoding='latin1')
            with open(pklFile, 'rb') as fPkl:
                halos = pickle.load(fPkl, encoding='latin1')

            mTot = 0
            for h in halos:
                mTot = mTot + h.m

            thisRho = mTot / volZero / rhoZero
            allRhos.append(thisRho)
            #print('Found: ', pklFile, ', n: ', len(halos))
            
rhoMed = np.median(allRhos)
rhoMea = np.mean(allRhos)
rhoSig = np.std(allRhos)

print('Med: ', rhoMed, ' rhoMea: ', rhoMea, ' rhoSig: ', rhoSig)

