import numpy as np
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


class DescendantProgenitor:
	'A wrapper/class (structure-like) to save a few properties of the main branch halos. Useful for dictionaries.'

	progID = '123123123123123123'
	descID = '123123123123123123'
	progNPart = 0
	descNPart = 0
	descNProg = 0
	
	def __init__(self, progID, progNPart, descID, descNPart, descNProg):
		self.progID = progID
		self.descID = descID
		self.progNPart = progNPart
		self.descNPart = descNPart
		self.descNProg = descNProg



class MTree:
	'The Merger tree class'
	
	# How many snapshots is our merger tree composed of (upper limit)
	nSteps = 0
	
	# Main halo ID at z = 0. We use strings as we will be using dictionaries all along to connect these objects
	mainID = '123123123123123123123'

	# Time, expansion factor and redshift
	t = np.zeros((0))
	a = np.zeros((0))
	z = np.zeros((0))
	
	# Main branch halo: number of particles, main halo ID and number of descendants
	mainBranchID = [] 
	mainBranchNPart = np.full((0), 0)
	mainBranchNProg = np.full((0), 0)
	 
	# Keep track of the steps at which the halo was lost and replaced by a token instead
	isToken = np.full((0), False)

	# Initialize the Merger Tree by its (expected) number of steps
	def __init__(self, n, ID):
		self.mainID = ID
		self.nSteps = n
		self.z = np.zeros((n))
		self.a = np.zeros((n))
		self.t = np.zeros((n))
		self.mainBranchID = [''] * n
		self.mainBranchNPart = np.full((n), 0)
		self.mainBranchNProg = np.full((n), 0)
		self.isToken = np.full((n), False)

	# 
	def fill_mass_id(self, nps, ids):
		for iM in range(0 , self.nSteps):
			self.mainBranchNPart[iM] = nps[iM]
			self.mainBranchID[iM] = ids[iM]

	# 
	def get_mass_id(self):
		return [self.mainBranchNPart, self.mainBranchID]

	# Print the number of particles and mass ID corresponding 
	def print_mass_id(self):
		for iM in range(0, self.nSteps):	
			print("%s %d" % (self.mainBranchID[iM], self.mainBranchNPart[iM]))
	
	# Return the normalized (z=0) mass across all history
	def norm_mass(self):
		norm_mass = []
		m_zero = self.mainBranchNPart[0]
		
		return self.mainBranchNPart[:]/m_zero
	
	# Smooth the mass accretion history of the object	TODO
	def smooth_mass(self):

		smooth = []

		return smooth

	# Add the descendant progenitor pairs to the next step
	def update_desc_prog(self, iStep, dP):
		self.mainBranchID[iStep] = dP.descID
		self.mainBranchNPart[iStep] = dP.descNPart
		self.mainBranchNProg[iStep] = dP.descNProg

	# Read & implement the steps from an expansion factor file
	def init_a(self, file_name):
		self.a = read_a(file_name)

	# Read & implement the steps from a redshift file
	def init_z(self, file_name):
		self.z = read_z(file_name)

	# When did the halo experience its last major merger
	def last_major_merger(self):
		lmm = 0.0
		return lmm

	# When did the halo form
	def formation_time(self):
		ft = 0.0
		return ft

