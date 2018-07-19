#!/usr/bin/python

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde
 

def index(i, j, k, n):
	return (i + j * n + k * n * n);

def find_indexes_in_slab(box, ngrid, data, min_x, max_x):
	cell = box / float(grid)

	dim = max_x[0] - min_x[0]
	slab_x = np.zeros(dim); 	slab_y = np.zeros(dim);		slab_z = np.zeros(dim)
	slab_d = np.zeros((dim, dim))

	for i in range(min_x[0], max_x[0]):
		for j in range(min_x[1], max_x[1]):
			for k in range(min_x[2], max_x[2]):
				ind = index(i, j, k, ngrid)

				ii = i - min_x[0]
				jj = j - min_x[1]
				slab_x[ii] = i * cell
				slab_y[jj] = j * cell
				slab_d[ii, jj] += data[ind]

				'''
				value = data[ind]
				slab_x.append(i * cell) 
				slab_y.append(j * cell)
				slab_z.append(k * cell)
				slab_d.append(value)
				'''

	return (slab_x, slab_y, slab_z, slab_d)

def read_cic(file_base, format_type):
	data = []
	f_in = open(file_name, 'r')

	if format_type == 'ice':
		for line in f_in.readlines():
			values = line.split(' ')

			for value in values:
				if value != '' and value != '\n':
					value = float(value)
					data.append(value)	

	if format_type == 'std':
		for line in f_in.readlines():
			values = line.split('\t')
			value = float(values[3])
			data.append(value)

	#		print value	


	return data





		#################################
		#				#
		#	START THE PROGRAM	#
		#				#
		#################################

grid = 256
box = 400

# Create a figure with 3 plot areas
fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(5, 5))

base_path = '/home/eduardo/CLUES/DATA/CF2YH/286393/'
#base_path = '/home/eduardo/CLUES/DATA/'
#file_base = 'CF2_grouped_YH_201611.co_128_400.000_WF_Dx.txt'
#file_base = 'CF2_grouped_YH_201611.co_256_400.000_Dx_rza.txt.10'
#file_base = 'CF2_grouped_YH_201611.co_256_400.000_WF_Dx.txt'
#file_base = 'CF2_grouped_YH_201611.co_256_400.000_Dx.txt.1'
#file_base = 'CF2_grouped_YH_201611.co_256_400.000_Dx_rza.txt'
#file_base='cf2gvpecc1pt5_256_400.000_Dx_rza.txt'
#file_base = 'cic_128.txt'
file_base = 'cic_256.txt'

file_name = base_path + file_base 

data = read_cic(file_base, 'std')

min_x = np.zeros((3), dtype=int)
max_x = np.zeros((3), dtype=int)

step = 0
thick = 16

min_x[0] = step; max_x[0] = grid - step
min_x[1] = step; max_x[1] = grid - step
min_x[2] = 0.5 * (grid - thick); max_x[2] = 0.5 * (grid + thick)

(sl_x, sl_y, sl_z, sl_d) = find_indexes_in_slab(box, grid, data, min_x, max_x)

plt.contour(sl_x, sl_y, sl_d)

outname = 'dens_' + file_base + '_' + str(step) + '.png'

plt.savefig(outname)

'''
# Thus we can cut the plotting window in several hexbins
nbins = 20
axes[1].set_title('Hexbin')
axes[1].hexbin(x, y, gridsize=nbins, cmap=plt.cm.BuGn_r)
 
# 2D Histogram
axes[2].set_title('2D Histogram')
axes[2].hist2d(x, y, bins=nbins, cmap=plt.cm.BuGn_r)
 
# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
k = kde.gaussian_kde(data.T)
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# plot a density
axes[3].set_title('Calculate Gaussian KDE')
axes[3].pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.BuGn_r)
 
# add shading
axes[4].set_title('2D Density with shading')
axes[4].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)
 
# contour
axes[5].set_title('Contour')
axes[5].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)
axes[5].contour(xi, yi, zi.reshape(xi.shape) )


'''


