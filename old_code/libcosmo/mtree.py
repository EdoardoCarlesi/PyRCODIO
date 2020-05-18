import numpy as np
import pickle
import random
import math 
import sys
from operator import *
#from .find_halos import *
#from .utils import *
#from .units import *

class MergerTree:
    # Major merger
    MM = 0

    # Formation time
    FT = 0

    # Number of steps to track the tree
    nSteps = 0

    # Number of particles per step
    nPart = np.zeros((0), dtype=int)

    # Normalized to z=0
    nPartNorm = np.zeros((0))

    # Smoothed mah normalized to z=0 
    nPartSmooth = np.zeros((0))

    # List all halo IDs of the main progenitor
    IDs = []

    # Vmax 
    Vmax = []

    # Keep track of its original subdir
    subDir = ''

    def __init__(self, nSteps, readPart, readIDs, subDir):
        self.MM = 0
        self.FT = 0
        self.IDs = []
        self.Vmax = []
        self.nSteps = nSteps
        self.nPart = np.zeros((nSteps), dtype=int)
        self.nPartNorm = np.zeros((nSteps))
        self.nPartSmooth = np.zeros((nSteps))
        self.subDir = subDir

        # Initialize and compute additional tree properties
        self.init_parts(readPart)
        self.init_ids(readIDs)
        self.norm_parts()
        self.smooth_tree()

    def init_parts(self, parts):
        iN = 0
        for part in parts:
            self.nPart[iN] = int(part)
            iN += 1


    # We are assuming that N_Particles(z = 0) is the first element, nPart[0]
    def last_major_merger(self, smooth):
        MajorMerg = 0.1
    
        if smooth == True:
            self.smooth()
        else:
            self.unsmooth()

        #thesePart = self.nPartNorm
        thesePart = self.nPart

        #if smooth == True:
        #    thesePart = self.nPartSmooth
        #else:

        for iMM in range(0, self.nSteps-1):
            nNP0 = float(thesePart[iMM])
            nNP1 = float(thesePart[iMM+1])

            dNP = abs(nNP1 - nNP0) / nNP0
            #dNP = abs(nNP1 - nNP0) / nNP1

            if dNP > MajorMerg:
                #print(iMM, nNP1, nNP0)
                return (iMM)

    def formation_time(self, smooth):
        mHalf = 0.5

        '''
        if smooth == True:
            self.smooth()
        else:
            self.unsmooth()
        '''

        thesePart = self.nPartNorm

        indexForm = np.where(thesePart < mHalf)
        thisIndex = indexForm[0]

#        print(thisIndex, thesePart)

        return (thisIndex[0])

    def init_ids(self, allIDs):
        for thisID in allIDs:
           self.IDs.append(thisID)

    def norm_parts(self):
        nNorm = float(self.nPart[0])
        iN = 0
        for thisNPart in self.nPart:
            try:
                self.nPartNorm[iN] = (1.0 * thisNPart) / nNorm
                iN += 1
            except:
                dummy = 0.0

    def smooth(self):
        iN = 1
        for thisN in self.nPart[1:(self.nSteps-2)]:
            self.nPart[iN] = int( (0.5 * self.nPart[iN-1] + thisN + 0.5 *self.nPart[iN+1])/2.0 )
            #self.nPart[iN] = int( (self.nPart[iN-1] + thisN/2 ))
            iN = iN + 1


        '''
        for thisN in self.nPart[0:(self.nSteps-1)]:
            if self.nPart[iN+1] < thisN:
                fact = 0.5 / (1.0 + np.sqrt(iN))
                self.nPart[iN+1] += int(fact * (thisN - self.nPart[iN+1])) 
                #self.nPart[iN+1] = self.nPart[iN+1]
        '''
        self.norm_parts()


    def unsmooth(self):
        np.random.seed()
        thisN = 10; nextN = 0;        delta = -10

        iN = 0
        for thisN in self.nPart[0:(self.nSteps-1)]:
            nextN = self.nPart[iN+1]
            #extra = random.uniform(0, nextN)
            #extra = np.random.normal(nextN, 10)
            extra = np.random.randn() * np.sqrt(nextN)/1.25
            #delta = int(abs(int(thisN) - int(nextN)) * thr)
            #delta = int(abs(int(thisN) - int(nextN)) * 0.1)
            #delta = abs(int(thisN) - int(nextN)) 
            delta = abs(int(thisN) - int(nextN)) 
            self.nPart[iN] -= int(delta/3) + int(extra) #int(extra/(1 + iN)))
            #self.nPart[iN] = int((1.0 + thr) * thisN)
            #if thisN > 10000000:
            #print(iN, self.nPart[iN], self.nSteps) 
            #print(iN, delta, thisN, nextN, self.nSteps)
            #print(iN, thisN, nextN, delta, self.nSteps)
            iN += 1
        '''
        '''
        self.norm_parts()

    def smooth_tree(self):
        nPtsAvg = 1

        for iPts in range(0, nPtsAvg):
            self.nPartSmooth[iPts] = self.nPartNorm[iPts]

        for iPts in range(0, nPtsAvg):
            self.nPartSmooth[self.nSteps-iPts-1] = self.nPartNorm[self.nSteps-iPts-1]

        for iPts in range(nPtsAvg, self.nSteps-nPtsAvg):
            for iSmooth in range(iPts-nPtsAvg, iPts+nPtsAvg+1):
                if iSmooth < iPts-1:
                    thisNorm = self.nPartSmooth[iSmooth] 
                else:
                    thisNorm = self.nPartNorm[iSmooth] 

                if thisNorm > 1.0:
                    self.nPartSmooth[iPts] += thisNorm * abs(1-0.75*(thisNorm - 1.0))
                else:
                    self.nPartSmooth[iPts] += thisNorm

                self.nPartSmooth[iPts] = self.nPartSmooth[iPts] / float(2*nPtsAvg)

        return self.nPartSmooth

    def info(self):
        print('MTree with: %d steps' % self.nSteps)
        print('DirCode   : %s' % self.subDir)
