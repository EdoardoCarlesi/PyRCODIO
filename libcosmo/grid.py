import math
import numpy as np


class Grid:
	'Generic class to store grids'	
	
	size = None
	box = None
	cell = None
	
	# Base grid only contains info about density and velocity
	rho = None
	vel = None

	def __init__(self, size, box):
		self.size = size
		self.box = box
		self.cell = float(box) / float(size)
		self.rho = np.zeros((size, size, size))
		self.vel = np.zeros((3, size, size, size))

	# This converts a generic position given in physical units to grid units
	def phys2grid(self, pos):
		ijk = np.zeros((3), dtype=int)

		for ix in range(0, 3):
			ijk[ix] = int(pos[ix] / self.cell)

		return ijk
	
	def rhoPos(self, pos):
		ijk = self.phys2grid(pos)

		return self.rho[ijk[0], ijk[1], ijk[2]]
	
	def velPos(self, pos):
		ijk = self.phys2grid(pos)

		return self.vel[:, ijk[0], ijk[1], ijk[2]]

	def reverse_index(self, ind):
		i = ind % self.size
		j = (ind / self.size) % self.size
		k = (ind / (self.size * self.size) ) % self.size
		return [i, j, k]

	def index(self, i, j, k):
		ind = i + self.size * j + self.size * self.size * k
		return ind

	def grid2phys(self, ijk):
		halfcell = 0.5 * self.cell
		pos = np.zeros((3))

		for ix in range(0, 3):
			pos[ix] = ijk[ix] * self.cell + halfcell

		return pos

class VWeb(Grid):

	evals = None
	evecs = None

	def __init__(self, size, box):
		self.size = size
		self.box = box
		self.cell = float(box) / float(size)
		self.rho = np.zeros((size, size, size))
		self.vel = np.zeros((3, size, size, size))

		# At each size node save three eigenvalues, plus three coordinates for each eigenvector
		self.evals = np.zeros((3, size, size, size))
		self.evecs = np.zeros((3, 3, size, size, size))
	
	def evalPos(self, pos):
		ijk = self.phys2grid(pos)

		return self.evals[:, ijk[0], ijk[1], ijk[2]]
	
	def evecPos(self, pos):
		ijk = self.phys2grid(pos)

		return self.evecs[:, :, ijk[0], ijk[1], ijk[2]]




