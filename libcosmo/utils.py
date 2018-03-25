import random
import math
import numpy as np
from numpy import linalg as la

def distance(vec1, vec2):
	dist = math.sqrt(pow((vec1[0] - vec2[0]), 2) + pow((vec1[1] - vec2[1]), 2) + pow((vec1[2] - vec2[2]), 2))
	return dist

def module(vec):
	mod = math.sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2))
	return mod

def max_list(vec):
	n_y = len(vec)
	v0 = vec[0][0]

	for iy in range(0, n_y):
		n_x = len(vec[iy])
		for ix in range(0, n_x):
			if vec[iy][ix] > v0:
				v0 = vec[iy][ix]

	return v0

def angle(v1, v2):
	mv1 = [0.0] * 3	
	mv2 = [0.0] * 3
	mod1 = module(v1)
	mod2 = module(v2)

	for i in range(0, 3):
		mv1[i] = v1[i] / mod1
		mv2[i] = v2[i] / mod2

	v12 = dot_prod(mv1, mv2)
	return v12

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

def abs_val(x):
	return math.sqrt(x * x)

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

def moment_inertia_reduced(coord, masses):
	n_parts = len(masses)
	m_in = np.zeros((3,3))
	center = [0.0] * 3		


	for ip in range(0, n_parts):
		# Diagonal elements
		R = distance(coord[ip], center)
		
		if R != 0.0:
			m_in[0][0] += masses[ip] * (pow(coord[ip][1], 2) + pow(coord[ip][2], 2)) / R
			m_in[1][1] += masses[ip] * (pow(coord[ip][0], 2) + pow(coord[ip][2], 2)) / R
			m_in[2][2] += masses[ip] * (pow(coord[ip][1], 2) + pow(coord[ip][0], 2)) / R

			# Off-diagonal
			m_in[1][0] += -masses[ip] * (coord[ip][1] * coord[ip][0]) / R
			m_in[1][2] += -masses[ip] * (coord[ip][1] * coord[ip][2]) / R
			m_in[0][2] += -masses[ip] * (coord[ip][0] * coord[ip][2]) / R

	for ix in range(0, 3):
		for iy in range(0, 3):
			m_in[ix][iy] /= n_parts

	m_in[0][1] = m_in[1][0]
	m_in[2][1] = m_in[1][2]
	m_in[2][0] = m_in[0][2]

	(e_val, e_vec) = la.eig(m_in)
	e_max = np.amax(e_val)
	print e_val/e_max

	return (e_val, e_vec)

def moment_inertia(coord, masses):
	n_parts = masses.size
	m_in = np.zeros((3,3))

	for ip in range(0, n_parts):
		# Diagonal elements
		m_in[0][0] += masses[ip] * (pow(coord[ip][1], 2) + pow(coord[ip][2], 2))
		m_in[1][1] += masses[ip] * (pow(coord[ip][0], 2) + pow(coord[ip][2], 2))
		m_in[2][2] += masses[ip] * (pow(coord[ip][1], 2) + pow(coord[ip][0], 2))

		# Off-diagonal
		m_in[1][0] += -masses[ip] * (coord[ip][1] * coord[ip][0])
		m_in[1][2] += -masses[ip] * (coord[ip][1] * coord[ip][2])
		m_in[0][2] += -masses[ip] * (coord[ip][0] * coord[ip][2])

	for ix in range(0, 3):
		for iy in range(0, 3):
			m_in[ix][iy] /= n_parts

	m_in[0][1] = m_in[1][0]
	m_in[2][1] = m_in[1][2]
	m_in[2][0] = m_in[0][2]

	(e_val, e_vec) = la.eig(m_in)
	e_max = np.amax(e_val)
	print e_val/e_max
	#print e_vec

	return (e_val, e_vec)

def mass_function(masses):
	n_m = len(masses)
	y_n = [0 for im in range(0, n_m)]
	x_m = sorted(masses)

	for im in range(0, n_m):
		y_n[im] = n_m - im

	return (x_m, y_n)

def rand_points_sphere(n_pts, center, radius):
	xyz_pts = np.zeros((n_pts, 3))
	x_rand = [0.0] * 3	
	x_min = [0.0] * 3	
	x_max = [0.0] * 3	
	i_pts = 0
	i_rej = 0
	i_tot = 0

#	for ix in range(0, 3):
#		x_min[ix] = center[ix] - radius
#		x_max[ix] = center[ix] + radius

	while (i_pts < n_pts):
		i_tot += 1
		nx = random.uniform(0.0, 1.0)
		ny = random.uniform(0.0, 1.0)
		nz = random.uniform(0.0, 1.0)
		#print nx, ny, nz
		dx = radius * (2 * nx - 1.0)
		dy = radius * (2 * ny - 1.0)
		dz = radius * (2 * nz - 1.0)
		#print dx, dy, dz

		x_rand = [dx, dy, dz]
		d_xyz = distance(x_rand, center)	

		if d_xyz < radius:
			xyz_pts[i_pts][:] = x_rand
			#print x_rand
			i_pts += 1
		else:
			i_rej +=1

	#print 'Rejected:%d, Accepted:%d, Total:%d ' % (i_rej, i_pts, i_tot)
	return xyz_pts

