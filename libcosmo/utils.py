import random
import math
import numpy as np
from numpy import linalg as la

def distance(vec1, vec2):
	dist = math.sqrt(pow((vec1[0] - vec2[0]), 2) + pow((vec1[1] - vec2[1]), 2) + pow((vec1[2] - vec2[2]), 2))
	return dist

def module(vec):
	n_v = len(vec)	
	elem = 0.0

	for i in range(0, n_v):
		elem += pow(vec[i], 2)

	return math.sqrt(elem)

def is_there(val, vec):
	is_there = False
	n_vec = len(vec)
	entries = []

	for i_vec in range(0, n_vec):
		if str(vec[i_vec]) == str(val):
			is_there = True
			entries.append(i_vec)

	return (is_there, entries)
	
def max_list(vec):
	n_y = len(vec)
	v0 = vec[0][0]
	#v0 = None

	for iy in range(0, n_y):
		n_x = len(vec[iy])
		for ix in range(0, n_x):
			if vec[iy][ix] > v0:
				v0 = vec[iy][ix]

	return v0

def min_list(vec):
	n_y = len(vec)
	v0 = vec[0][0]
	#v0 = vec[0]

	for iy in range(0, n_y):
		n_x = len(vec[iy])
		for ix in range(0, n_x):
			if vec[iy][ix] < v0:
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

def change_basis(vec, new_base):
	new_vec = np.zeros((3))
	base_change = np.zeros((3, 3))

	I = np.zeros((3, 3))
	
	for iv in range(0, 3):
		I[iv, iv] = 1.0

	# Find out the matrix of base change:
	for iv in range(0, 3):
		base_change[iv, :] = np.linalg.solve(new_base, I[iv, :])

	for iv in range(0, 3):
		new_vec[iv] = dot_prod(vec, base_change[iv, :])

	return new_vec

def abs_val(x):
	absval = math.sqrt(x * x)
	return absval

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
	n = len(x)
	norm = 0.0
	vn = []

	for ix in range(0, n):
		norm += x[ix] * x[ix]
	norm = math.sqrt(norm)	

	for ix in range(0, n):
		vn.append(x[ix] / norm)

	return vn

def vel_radial(x1, x2, v1, v2):
	x12 = vec_subt(x2, x1)	
	n12 = math.sqrt(dot_prod(x12, x12))
	r12 = dot_prod(vec_subt(v2, v1), x12) 
	nr12 = r12 / n12
	return nr12

def inertiaTensor(x, y, z):
	I=[]
	for index in range(9):
        	I.append(0)
	
	#n_x = len(x)
	#for ix in range(0, n_x):
	#	I[0] += y[ix] * y[ix] + z[ix] * z[ix]

	I[0] = np.sum(y*y+z*z) 
	I[1] = np.sum(-y*x)    
	I[2] = np.sum(-x*z)    
	I[3] = np.sum(-y*x)    
	I[4] = np.sum(x*x+z*z) 
	I[5] = np.sum(-y*z)    
	I[6] = np.sum(-z*x)    
	I[7] = np.sum(-z*y)    
	I[8] = np.sum(x*x+y*y) 

	tensor = np.array([(I[0:3]), (I[3:6]), (I[6:9])])
	vals, vects = np.linalg.eig(tensor)  # they come out unsorted, so the command below is needed
	eig_ord = np.argsort(vals)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.
	ord_vals = vals[eig_ord]
	ord_vects = vects[:, eig_ord].T

	TriaxParam = (ord_vals[2]*ord_vals[2]-ord_vals[1]*ord_vals[1])/(ord_vals[2]*ord_vals[2]-ord_vals[0]*ord_vals[0])
	AxisRatio = ord_vals[0]/ord_vals[2]

	return ord_vals, ord_vects, TriaxParam, AxisRatio


def moment_inertia(coord, masses):
	n_parts = masses.size
	m_in = np.zeros((3,3))

	x = coord[:, 0]
	y = coord[:, 1]
	z = coord[:, 2]

	m_in[0][0] = np.sum(y*y + z*z) 
	m_in[1][1] = np.sum(x*x + z*z)
	m_in[2][2] = np.sum(x*x + y*y)

	m_in[1][0] = np.sum(-y * x)
	m_in[1][2] = np.sum(-y * z)
	m_in[0][2] = np.sum(-x * z)

	m_in[0][1] = m_in[1][0]
	m_in[2][1] = m_in[1][2]
	m_in[2][0] = m_in[0][2]

	(e_val, e_vec) = la.eig(m_in)
	e_max = np.amax(e_val)

	eig_ord = np.argsort(e_val)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.
	ord_vals = e_val[eig_ord]
	ord_vecs = e_vec[:, eig_ord].T

	return (ord_vals, ord_vecs)

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

	m_in[0][1] = m_in[1][0]
	m_in[2][1] = m_in[1][2]
	m_in[2][0] = m_in[0][2]

	(e_val, e_vec) = la.eig(m_in)

	eig_ord = np.argsort(e_val)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.
	ord_vals = e_val[eig_ord]
	ord_vecs = e_vec[:, eig_ord].T

	return (ord_vals, ord_vecs)


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

