import math
import numpy as np
from numpy import linalg as la

def distance(vec1, vec2):
	dist = math.sqrt(pow((vec1[0] - vec2[0]), 2) + pow((vec1[1] - vec2[1]), 2) + pow((vec1[2] - vec2[2]), 2))
	return dist

def module(vec):
	mod = math.sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2))
	return mod

def center_of_mass(m,x):
	n = len(m)
	com = [0.0] * 3 
	
	for j in range(0, 3):
		mtot = 0.0

		for i in range(0, n):
			mtot += m[i] 
			com[j] += x[i][j] * m[i]
		
		com[j] /= mtot

	return com

def vec_subt(x1, x2):
	vs = [x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]]
	return vs

def dot_prod(x1, x2):
	dp = x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] 
	return dp

def vec_divide(x, fac):
	n = len(x)
	nx = [0] * n

	for i in range(0, n):
		nx[i] = x[i] / fac

	return nx
	
def vec_module(v):
	vmod = dot_prod(v, v)
	return math.sqrt(vmod)

def vec_distance(x, y):
	vsub = vec_subt(x, y)
	vdis = vec_module(vsub)
	return vdis

def vec_norm(x):
	n = x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
	vn = [x[0]/n, x[1]/n, x[2]/n]
	return vn

def vel_radial(x1, x2, v1, v2):
	x12 = vec_subt(x2, x1)	
	n12 = math.sqrt(dot_prod(x12, x12))
	r12 = dot_prod(vec_subt(v2, v1), x12) 
	nr12 = r12 / n12
	return nr12

def moment_inertia(coord, masses):
	n_parts = len(masses)
	m_in = [[0 for ix in range(0, 3)] for iy in range(0, 3)]
	
	for ip in range(0, n_parts):
		# Diagonal elements
		m_in[0][0] += masses[ip] * (pow(coord[ip][1], 2) + pow(coord[ip][2], 2))
		m_in[1][1] += masses[ip] * (pow(coord[ip][0], 2) + pow(coord[ip][2], 2))
		m_in[2][2] += masses[ip] * (pow(coord[ip][1], 2) + pow(coord[ip][0], 2))

		# Off-diagonal
		#m_in[1][0] += -masses[ip] * (coord[ip][1] * coord[ip][0])
		#m_in[1][2] += -masses[ip] * (coord[ip][1] * coord[ip][2])
		#m_in[0][2] += -masses[ip] * (coord[ip][0] * coord[ip][2])

	m_in[0][1] = m_in[1][0]
	m_in[2][1] = m_in[1][2]
	m_in[2][0] = m_in[0][2]

	#print m_in[0][0]	
	#print m_in[1][1]	
	#print m_in[2][2]	

	(e_val, e_vec) = la.eig(m_in)
	e_max = np.amax(e_val)
	print e_val/e_max
	#print e_vec



