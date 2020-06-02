'''
    Python Routines for COsmology and Data I/O (PyRCODIO)
    Pandas Upgrade
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    tools.py: various utilities and simple computational routines used throughout the code
'''

import pandas as pd
import numpy as np
import sys

sys.path.append('pygadgetreader/')
from pygadgetreader import *


'''
    Compute the Euclidean distance between two points in space
'''
def distance(x, center):
    dist = 0;

    for i in range(0, len(x)):
        dist += (x[i] - center[i])**2.0

    dist = np.sqrt(dist)

    return dist

'''
    Below is a series of simple functions copied from the old utils.py
'''

def rad2deg():
	coeff = 180.0 / math.pi
	return coeff

def module(vec):
	n_v = len(vec)	
	elem = 0.0

	for i in range(0, n_v):
		elem += pow(vec[i], 2)

	return math.sqrt(elem)

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

def mass_function(masses):
	n_m = len(masses)
	y_n = [0 for im in range(0, n_m)]
	x_m = sorted(masses)

	for im in range(0, n_m):
		y_n[im] = n_m - im

	return (x_m, y_n)

'''
    Read multiple gadget files (one snap split into several)
'''
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






