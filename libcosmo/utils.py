import random
import math
import numpy as np
import scipy.stats as sp
import pickle
#from .halos import *
#from halos import *
from numpy import linalg as la
from pygadgetreader import *


def readgadget(f_snap, read_type, p_type, n_files):
    if n_files == 1:
        print('Read snap: ', f_snap, ' read type: ', read_type, ' ptype: ', p_type)
        parts = readsnap(f_snap, read_type, p_type)
    else:
        parts = np.zeros((1, 3), dtype=float)
        size_part = 0
        for i_file in range(0, n_files):
            f_tmp = f_snap + '.' + str(i_file)
            parts_tmp = readsnap(f_tmp, read_type, p_type)

            old_size = size_part
            size_part += len(parts_tmp)
            parts.resize((size_part, 3))
            parts[old_size:size_part][:] = [xyz for xyz in parts_tmp[:][:]]

            '''
            for i_part in range(0, len(parts_tmp)):
                parts[old_size + i_part, :] = parts_tmp[i_part][:]
#                i_x = 2
#                print(parts_tmp[i_part][i_x])
#                print(parts[old_size + i_part][i_x])
#                for i_x in range(0, 3):
#                    parts[old_size + i_part, i_x] = parts_tmp[i_part][i_x]
                    #print(parts_tmp[i_part][i_x])
                    #print(parts[old_size + i_part][i_x])
            '''
    print('Found a total of ', np.shape(parts), ' particles.')
    return parts


def KS_randBootstrap(sample1, n_boot, n_iter):
    n1 = len(sample1) # Smaller sample
    cdf1 = cdf(sample1)
    n_stat = 0
    
    #print('SAMPLE1: ', sample1)

    # Generate n_iter subsamples (bootstrapping)
    for it in range(0, n_iter):
        new_sample = np.zeros((n_boot))
        
        # Generate
        for ip in range(0, n_boot):
            new_index = random.randint(0, n_boot-1)
            new_sample[ip] = sample1[new_index]

        cdf2 = cdf(new_sample)               
        kst = sp.ks_2samp(cdf1[0], cdf2[0])

        #bmark = 2.4e-6
        bmark = 0.05 #2.4e-6
        thisks = kst[1] * n1

        if thisks < bmark:
            n_stat += 1

        #print(it, kst, kst[1] * n1)

    print('Bootstrap Stat: ', float(n_stat)/float(n_iter))


def KS_randSubSample(sample1, sample2, n_iter):
    n1 = len(sample1) # Smaller sample
#    n1 = 300
    n2 = len(sample2) # Larger sample
    cdf1 = cdf(sample1)
    n_stat = 0
    
    #print('SAMPLE1: ', sample1)

    # Generate n_iter subsamples (bootstrapping)
    for it in range(0, n_iter):
        new_sample = np.zeros((n1))
        
        # Generate
        for ip in range(0, n1):
            new_index = random.randint(0, n2-1)
            new_sample[ip] = sample2[new_index]

        cdf2 = cdf(new_sample)               
        kst = sp.ks_2samp(cdf1[0], cdf2[0])

        bmark = 0.95
        thisks = kst[1] * n1

        if thisks > bmark:
            n_stat += 1

        #print(it, kst, kst[1] * n1)

    #print('Stat: ', n_stat)
    print('Bootstrap Stat: ', float(n_stat)/float(n_iter))

def mahFit(z, a, b):
        zb = np.power(1+z, b)
        za = np.exp(-a * (np.sqrt(1 + z) - 1))
        return (zb * za)

def makeHaloDictionary(halos):
    haloDict = dict()

    iHalo = 0
    for halo in halos:
        haloDict[int(halo.ID)] = iHalo
        iHalo += 1

    return haloDict 

def distance(vec1, vec2):
	dist = math.sqrt(pow((vec1[0] - vec2[0]), 2) + pow((vec1[1] - vec2[1]), 2) + pow((vec1[2] - vec2[2]), 2))
	return dist

def rad2deg():
	coeff = 180.0 / math.pi
	return coeff

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
	try:
		n_parts = masses.size
	except:
		n_parts = len(masses)

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


def rand_points_sphere(n_pts):
	xyz_pts = np.zeros((n_pts, 3))
	x_rand = [0.0] * 3	
	x_min = [0.0] * 3	
	x_max = [0.0] * 3	
	i_pts = 0
	i_rej = 0
	i_tot = 0

	# Set the default center 
	radius = 1.0
	center = [1.0] * 3

	diameter = 2 * radius

	while (i_pts < n_pts):
		i_tot += 1
		nx = random.uniform(0.0, diameter) - center[0]
		ny = random.uniform(0.0, diameter) - center[1]
		nz = random.uniform(0.0, diameter) - center[2]
		x_rand = [nx, ny, nz]

		d_xyz = distance(x_rand, center)	

		#if d_xyz < radius:
		xyz_pts[i_pts, :] = x_rand
		i_pts += 1
		#else:
		#	i_rej +=1

	#print 'Rejected:%d, Accepted:%d, Total:%d ' % (i_rej, i_pts, i_tot)
	return xyz_pts

# The additional vector defines the e3 eigenvector of the VWeb
def random_triaxialities_and_angles(n_pts, n_trials, vector):
	
	med_e = np.zeros((3, 3))
	med_d = np.zeros((3, 3))
	med_c = np.zeros((3, 3))
	masses = [1.0] * n_pts

	rand_evals = np.zeros((3, n_trials))
	rand_disp = np.zeros((3, n_trials))
	rand_cos = np.zeros((3, n_trials))

	percMin = 5.0
	percMax = 100.0 - percMin
	new_positions = np.zeros((3, n_pts))

	for i_trial in range(0, n_trials):
		these_pts = rand_points_sphere(n_pts)
		(evals, evecs) = moment_inertia(these_pts, masses)
		rand_evals[:, i_trial] = evals / evals[2]

		#print evecs
		#print evecs
		#print i_trial, rand_cos[:, i_trial]
		
		rand_cos[0, i_trial] = abs(angle(vector, evecs[0, :]))
		rand_cos[1, i_trial] = abs(angle(vector, evecs[1, :]))
		rand_cos[2, i_trial] = abs(angle(vector, evecs[2, :]))

		#rad2deg = 180.0 / math.pi
		#rand_cos[0, i_trial] = np.arccos(abs(angle(vector, evecs[0, :]))) * rad2deg
		#rand_cos[1, i_trial] = np.arccos(abs(angle(vector, evecs[1, :]))) * rad2deg
		#rand_cos[2, i_trial] = np.arccos(abs(angle(vector, evecs[2, :]))) * rad2deg

		for i_sat in range(0, n_pts):
			new_pos = change_basis(these_pts[i_sat, :], evecs)
			new_positions[:, i_sat] = new_pos			

		rand_disp[0, i_trial] = np.std(abs(new_positions[0, :]))
		rand_disp[1, i_trial] = np.std(abs(new_positions[1, :]))
		rand_disp[2, i_trial] = np.std(abs(new_positions[2, :]))
	
	for ie in range(0, 3):
		med_e[ie, 1] = np.median(rand_evals[ie, :])
		med_e[ie, 0] = np.percentile(rand_evals[ie, :], percMin)
		med_e[ie, 2] = np.percentile(rand_evals[ie, :], percMax)

		med_d[ie, 1] = np.median(rand_disp[ie, :])
		med_d[ie, 0] = np.percentile(rand_disp[ie, :], percMin)
		med_d[ie, 2] = np.percentile(rand_disp[ie, :], percMax)

		med_c[ie, 1] = np.median(rand_cos[ie, :])
		med_c[ie, 0] = np.percentile(rand_cos[ie, :], percMin)
		med_c[ie, 2] = np.percentile(rand_cos[ie, :], percMax)

		#print ie, med_e[ie, 0], med_e[ie, 1], med_e[ie, 2] 
		#print ie, med_d[ie, 0], med_d[ie, 1], med_d[ie, 2] 
		#print ie, med_c[ie, :]

	return (med_e, med_d, med_c)




def random_triaxialities(n_pts, n_trials, compare_value):
	
	med_e = np.zeros((3, 3))
	med_d = np.zeros((3, 3))
	masses = [1.0] * n_pts

	rand_evals = np.zeros((3, n_trials))
	rand_disp = np.zeros((3, n_trials))

	percMin = 5.0
	percMax = 100.0 - percMin
	new_positions = np.zeros((3, n_pts))

	for i_trial in range(0, n_trials):
		these_pts = rand_points_sphere(n_pts)
		(evals, evecs) = moment_inertia(these_pts, masses)
		rand_evals[:, i_trial] = evals / evals[2]

		for i_sat in range(0, n_pts):
			new_pos = change_basis(these_pts[i_sat, :], evecs)
			new_positions[:, i_sat] = new_pos			

		rand_disp[0, i_trial] = np.std(abs(new_positions[0, :]))
		rand_disp[1, i_trial] = np.std(abs(new_positions[1, :]))
		rand_disp[2, i_trial] = np.std(abs(new_positions[2, :]))
	
	for ie in range(0, 3):
		med_e[ie, 1] = np.median(rand_evals[ie, :])
		med_e[ie, 0] = np.percentile(rand_evals[ie, :], percMin)
		med_e[ie, 2] = np.percentile(rand_evals[ie, :], percMax)

		med_d[ie, 1] = np.median(rand_disp[ie, :])
		med_d[ie, 0] = np.percentile(rand_disp[ie, :], percMin)
		med_d[ie, 2] = np.percentile(rand_disp[ie, :], percMax)

		#print ie, med_e[ie, 0], med_e[ie, 1], med_e[ie, 2] 
		#print ie, med_d[ie, 0], med_d[ie, 1], med_d[ie, 2] 
	
	#print rand_evals[0, :]
	score_value = sp.percentileofscore(rand_evals[0, :], compare_value)

	return (med_e, med_d, score_value)



def random_table_triaxialities(n_pts, n_trials, read_table):
		
	points = '%05d' % n_pts
	trials = '%05d' % n_trials
	triax_name = 'saved/rand_triax_' + points + '_' + trials + '.pkl'

	if read_table == True:

		try:
			print('Reading pre-stored table values from: ', triax_name)
			f_evals = open(triax_name, 'r')
			rand_evals = pickle.load(f_evals)
			return rand_evals
		except:
			print('File not found: ', triax_name )
	
	# Generate tables
	else:
		masses = [1.0] * n_pts
		rand_evals = np.zeros((3, n_trials))

		for i_trial in range(0, n_trials):
			these_pts = rand_points_sphere(n_pts)
			(evals, evecs) = moment_inertia(these_pts, masses)
			rand_evals[:, i_trial] = evals / evals[2]

		f_evals = open(triax_name, 'w')
		print('Saving table to: ', triax_name)
		pickle.dump(rand_evals, f_evals)

	#print rand_evals[0, :]
	#score_value = sp.percentileofscore(rand_evals[0, :], compare_value)

	#return (med_e, med_d, score_value)

def cdf(values):
	n_values = len(values)
	cdf = np.zeros((2, n_values))
	
	sv = np.sort(values)
	#norm = np.max(values)
	#norm = values[n_values-1]	

	for i in range(0, n_values):
		cdf[0][i] = sv[i] 
		cdf[1][i] = float(i) / float(n_values)
		#cdf[1][i] = float(n_values - i) / float(n_values)

	return cdf






