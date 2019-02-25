import numpy as np
import sys, os
from .mtree import *


class ReadSettings:
	baseFile = ''
	suffFile = ''
	nChunk = 0
	nSnaps = 0 
	nSteps = 0
	treeFiles = []
	allTrees = []
	treeIDz0 = []

	def __init__(self, baseFile, suffFile, nChunk, nSnaps, nSteps):
		print('ReadSettings, initialize with baseFile, suffFile, nChunk, nSteps')
		self.nChunk = nChunk
		self.nSnaps = nSnaps
		self.nSteps = nSteps
		self.baseFile = baseFile
		self.suffFile = suffFile
		self.treeFiles = []
		self.allTrees = []
		self.treeIDz0 = []

		self.gen_file_list()	

	def gen_file_list(self):		
		for iStep in range(0, self.nSteps):			
			strStep = "%03d." % (self.nSnaps - iStep)
			theseChunks = []

			for iChunk in range(0, self.nChunk):
				strChunk = "%d." % iChunk 			
				thisFile = self.baseFile + strStep + strChunk + self.suffFile
				theseChunks.append(thisFile)

			self.treeFiles.append(theseChunks)

	# All the trees files will be read and the trees properly built
	def read_trees(self):
	
		self.allTrees = []
 
		# Here we store the progenitor index for the next step
		allOldIndex = dict()
		allTmpIndex = dict()

		for iStep in range(0, self.nSteps):
				
			# We have to merge all the chunks-per-snapshot into a single dictionary		
			allChunkDescProg = dict()

			for iChunk in range(0, self.nChunk):
				chunkDescProg = read_metrocpp_tree(self.treeFiles[iStep][iChunk])
		
				for idDesc, descProg in chunkDescProg.items():
					allChunkDescProg[idDesc] = descProg

					# At the first step, initialize all the merger trees
					if iStep == 0:
						thisTree = MTree(self.nSteps, idDesc)
						thisTree.update_desc_prog(iStep, descProg)
						iTree = len(self.allTrees)
						allOldIndex[descProg.progID] = iTree
						self.allTrees.append(thisTree) 
						
						del thisTree

			# Print for consistency check 
			# print("%d) New dictionary has %d elements, old dictionary has %d elements" % (iStep, len(allChunkDescProg), len(allOldIndex)))

			# Do this once all the chunks have been gathered
			if iStep > 0:

				# Match the progenitor at iStep-1 with the descendant at iStep
				# allOldIndex is initialized at iStep-1, allChunkDP is initialized at this step
				for idProg, indTree in allOldIndex.items():	
					if idProg in allChunkDescProg.keys():	
						thisDescProg = allChunkDescProg[idProg]
						self.allTrees[indTree].update_desc_prog(iStep, thisDescProg)
						allTmpIndex[thisDescProg.progID] = indTree

				allOldIndex.clear()

				# Manually copy everything, just to be on the safe side...
				for idProg, indTree in allTmpIndex.items():	
					allOldIndex[idProg] = indTree

				allTmpIndex.clear()


		return self.allTrees		
			


def read_metrocpp_tree(fileMetroCpp):
	
	print('Reading file: %s in MetroCPP format.' % fileMetroCpp)

	# IDs are already initialized as strings as they will be stored into dictionaries
	descID = '123123123123123123123123'
	progID = '123123123123123123123123'

	descHalo = np.full((3), 0) 
	progHalo = np.full((2), 0)

	# Create a dictionary that connects the descendant IDs to their main progenitor 
	descProg = dict()

	with open(fileMetroCpp) as fileIn:
		allLines = fileIn.read().splitlines()
		
	iLineProg = 0

	for thisLine in allLines:
		
		splitLine = thisLine.split()

		# First line contains: main halo, number of particles, number of progenitors and 0/1 value for token halos
		if len(splitLine) == 4:
			[descID, descHalo[0], descHalo[1], descHalo[2]] = thisLine.split()

			# Reset the line number of progenitor halos
			iLineProg = 0
			
			'''
			# Consistency check
			if descHalo[1] > 3:
				printLine = "ID:{0}, Np:{1}, Nd:{2}".format(descID, int(descHalo[0]), int(descHalo[1]))
				print(printLine)
			'''

		# Following Nprog lines contain: common particles, progenitor ID and number of particles per progenitor
		if len(splitLine) == 3:
			iLineProg += 1
			[progHalo[0], progID, progHalo[1]] = thisLine.split()
			
			if iLineProg == 1:
				# Here we save the number of progenitor/descendand id, the number of particles of both and the number of progenitors per descendant
				thisDescProg = DescendantProgenitor(progID, progHalo[1], descID, descHalo[0], descHalo[1])
				descProg[descID] = thisDescProg

	return descProg


def read_ahf_tree(fileAhf):

	return 0



def read_z(fileZ):

	nLines = 0 # Read the number of file lines
	z = np.zeros((nLines))
	
	return z


def read_a(fileA):

	nLines = 0 # Read the number of file lines
	a = np.zeros((nLines))
	
	return a


def read_a_t(fileAT):

	nLines = 0 # Read the number of file lines
	a_t = np.zeros((2, nLines))
	
	return a_t


def read_z_t(fileZT):

	nLines = 0 # Read the number of file lines
	z_t = np.zeros((2, nLines))
	
	return z_t

