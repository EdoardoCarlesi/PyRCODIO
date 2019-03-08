import numpy as np
import pickle
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
    nPart = np.zeros((0))

    # Normalized to z=0
    nPartNorm = np.zeros((0))

    # Normalized to z=0 
    nPartSmooth = np.zeros((0))

    # Keep track of its original subdir
    subDir = ''

    def __init__(self, nSteps, readPart, subDir):
        self.MM = 0
        self.FT = 0
        self.nSteps = nSteps
        self.nPart = readPart #np.zeros((nSteps))
        self.nPartNorm = np.zeros((nSteps))
        self.nPartSmooth = np.zeros((nSteps))
        self.subDir = subDir

        nNorm = float(self.nPart[0])

        iN = 0
        for thisNPart in self.nPart:
            try:
                self.nPartNorm[iN] = (1.0 * thisNPart) / nNorm
                iN += 1
            except:
                dummy = 0.0

        # Automatically smooth the tree when reading in
        self.smooth_tree()

    # We are assuming that N_Particles(z = 0) is the first element, nPart[0]
    def last_major_merger(self, smooth):
        MajorMerg = 0.1 

        if smooth == True:
            thesePart = self.nPartSmooth
        else:
            thesePart = self.nPartNorm


        for iMM in range(0, self.nSteps-1):
            nNP0 = thesePart[iMM]
            nNP1 = thesePart[iMM+1]

            dNP = abs(nNP1 - nNP0) / nNP1

            if dNP > MajorMerg:
                return iMM

    def formation_time(self, smooth):
        mHalf = 0.5
        mZero = 2.0

        if smooth == True:
            thesePart = self.nPartSmooth
        else:
            thesePart = self.nPartNorm

        indexForm = np.where(thesePart < 0.5)
        thisIndex = indexForm[0]

        return thisIndex[0]

    def smooth_tree(self):
        nPtsAvg = 2

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

                self.nPartSmooth[iPts] /= float(2*nPtsAvg+1)

        return self.nPartSmooth

    def info(self):
        print('MTree with: %d steps' % self.nSteps)
