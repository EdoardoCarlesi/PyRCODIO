import math 

class Halo:
	ID = 1234567890123456789
	m  = 0.0
	x  = [0.0, 0.0, 0.0]
	v  = [0.0, 0.0, 0.0]
	r = 0.0

	def __init__(self, ind, mass, pos, vel, rad, n_sub, n_part):
		self.ID = ind
		self.m  = mass
		self.x  = pos 
		self.v  = vel
		self.r  = rad
		self.nsub = n_sub	
		self.npart = n_part

	def assign(self, ind, mass, pos, vel, rad, n_sub, n_part):
		self.ID = ind
		self.m  = mass 
		self.x  = pos 
		self.v  = vel
		self.r  = rad 
		self.nsub = n_sub
		self.npart = n_part

	def distance(self, rad):
		dist = math.sqrt( pow((self.x[0] - rad[0]), 2) + pow((self.x[1] - rad[1]), 2) + pow((self.x[2] - rad[2]), 2) )
		return dist

	def info(self):
		return "Halo %ld - Mass: %.3e - Pos: %s - Vel: %s\n" % (self.ID, self.m, self.x[:], self.v[:])

