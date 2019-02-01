import math 
import pickle
import numpy as np
import find_halos as fh
import sys
from scipy import interpolate
from operator import *
from utils import *
from units import *

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

	def __init__(self, nSteps, readPart):
		self.MM = 0
		self.FT = 0
		self.nSteps = nSteps
		self.nPart = readPart #np.zeros((nSteps))
		self.nPartNorm = np.zeros((nSteps))

		nNorm = float(self.nPart[0])

		iN = 0
		for thisNPart in self.nPart:
			self.nPartNorm[iN] = (1.0 * thisNPart) / nNorm
			iN += 1


	# We are assuming that N_Particles(z = 0) is the first element, nPart[0]
	def last_major_merger(self):
		MajorMerg = 0.1 

		for iMM in range(0, self.nSteps-1):
			nNP0 = self.nPartNorm[iMM]
			nNP1 = self.nPartNorm[iMM+1]
			
			dNP = abs(nNP1 - nNP0) / nNP1

			if dNP > MajorMerg:
				print(self.nPartNorm[iMM-2:iMM+3])
				return iMM
	
	def formation_time(self):
		mHalf = 0.5
		mZero = 2.0

		indexForm = np.where(self.nPartNorm < 0.5)
		thisIndex = indexForm[0]

		#print(self.nPartNorm[thisIndex[0]-2:thisIndex[0]+3])
		return thisIndex[0]

		'''
		iFT = 0


		while mZero > mHalf and iFT < self.nSteps-1:
			if 
			
			iFT += 1
		'''

		return self.f_time

	def info(self):
		print('MTree with: %d steps' % self.nSteps)



'''	
	def smooth_history(self):
		return smooth
'''
