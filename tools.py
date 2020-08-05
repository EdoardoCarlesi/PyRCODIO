'''
    Python Routines for COsmology and Data I/O (PyRCODIO)
    Pandas Upgrade
    Edoardo Carlesi 2020
    ecarlesi83@gmail.com

    tools.py: various utilities and simple computational routines used throughout the code
'''

import read_files as rf
import pandas as pd
import numpy as np
import random

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
    Given a set of coordinates, randomly shift them by a maximum of 'r' amount
'''
def shift(center, r):
    new_c = []

    for c in center:
        eps = random.randrange(-r, r)
        c = c + eps
        new_c.append(c)

    return np.array(new_c)


'''
    Below is a series of simple functions copied from the old utils.py
'''


def rad2deg():
	coeff = 180.0 / np.pi
	return coeff

def module(vec):
	n_v = len(vec)	
	elem = 0.0

	for i in range(0, n_v):
		elem += pow(vec[i], 2)

	return np.sqrt(elem)


def find_nearest_node_index(x, grid=None, box=None):
    
    cell = box / grid
    
    ix = np.floor(x[0] / cell)
    iy = np.floor(x[1] / cell)
    iz = np.floor(x[2] / cell)

    index = int(ix + grid * iy + grid * grid * iz)

    return index
    


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
	return np.sqrt(vmod)

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
	norm = np.sqrt(norm)	

	for ix in range(0, n):
		vn.append(x[ix] / norm)

	return vn

def vel_radial(x1, x2, v1, v2):
	x12 = vec_subt(x2, x1)	
	n12 = np.sqrt(dot_prod(x12, x12))
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
    Find the center of mass of a particle distribution
'''
def particles_com(part_df, cols=['X', 'Y', 'Z'], mass_types=1):
    com = [0.0] * 3

    if mass_types > 1:
        print('Error. Only available for ONE particle mass for all the distribution')
        return 0

    for i, col in enumerate(cols):
        com[i] = part_df[col].mean()

        # Convert to kpc in case
        if com[i] < 1.e+4:
            com[i] = com[i] * 1.e+3


    return np.array(com)

'''
    Find the particles belonging to a slab around a given point in space.
    Slab size, thickness and so on need to be specified.
'''
def find_slab(part_df=None, center=None, side=None, thick=None, rand_seed=69, reduction_factor=1.0, z_axis=2, units='kpc', cols=['X', 'Y', 'Z']):

    # Set some parameters
    kpcThresh = 1.e+4
    kpc2Mpc = 1.e-3

    minima = np.zeros((3))
    maxima = np.zeros((3))

    # Select the two axes for the 2D projection
    ax0 = (z_axis + 1) % 3
    ax1 = (z_axis + 2) % 3
    ax2 = z_axis

    # Column names
    col0 = cols[ax0]
    col1 = cols[ax1]
    col2 = cols[ax2]

    n_part = len(part_df)

    # Sanity check on the units
    half_n = int(n_part * 0.5)
    sum_coord = part_df[col0].iloc[half_n] + part_df[col1].iloc[half_n] + part_df[col2].iloc[half_n] 

    # Make sure the units are consistent
    if sum_coord < kpcThresh:
        side = side * kpc2Mpc
        center = center * ([kpc2Mpc] *3) 
        thick = thick * kpc2Mpc
        
        #print(part_df[part_df['Type'] == 4.0].head())
        #print(sum_coord, center, thick)

    # Set the minima and maxima for the particles to be used in the plot
    minima[ax0] = center[ax0] - side * 0.5
    minima[ax1] = center[ax1] - side * 0.5
    minima[ax2] = center[ax2] - thick * 0.5

    maxima[ax0] = center[ax0] + side * 0.5
    maxima[ax1] = center[ax1] + side * 0.5
    maxima[ax2] = center[ax2] + thick * 0.5

    # Find the particles in the slab
    condition_x = (part_df[col0] > minima[ax0]) & (part_df[col0] < maxima[ax0])
    condition_y = (part_df[col1] > minima[ax1]) & (part_df[col1] < maxima[ax1])
    condition_z = (part_df[col2] > minima[ax2]) & (part_df[col2] < maxima[ax2])
    part_select = part_df[(condition_x) & (condition_y) & (condition_z)]

    print('Found: ', len(part_select), ' particles in the slab')

    # Now select a random subsample of the full particle list
    if reduction_factor < 1.0:
        part_select = part_select.sample(frac=reduction_factor, random_state=rand_seed)

        print('The number of particles to be used has been reduced to: ', len(part_select))

    # Return the selected particles' properties in a dataframe
    return part_select


'''
    vweb is a DataFrame containing all the web information
'''
def smooth_web(vweb, x_point=None, smooth_length=1.5, smooth_type='avg'):
    x_col = ['x', 'y', 'z']
    new_col = 'Distance'

    
    # Take the simple average of all points within a smoothing_length distance
    if smooth_type == 'avg':
        '''
        vweb[x_col].apply(lambda x: distance(x_point))
        smooth_length
        '''

    return smooth


def check_units(data=None, cols=None):
    n_pts = int(len(data) * 0.5)

    vals = data[cols].iloc[n_pts]

    # If this is true, then the units are 
    if np.sum(vals) < 1.e+4:
        factor = 1.0e+3
    else:
        factor = 1.0

    data[cols] = data[cols].apply(lambda x: x * factor)
    #print(data.head())

    return factor

