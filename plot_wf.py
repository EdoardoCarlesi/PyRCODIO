import numpy as np
import matplotlib.pyplot as plt


# Wiener Filter file in IceCore format
path_wf = '/z/carlesi/CLUES/IceCore/output/'
name_wf = path_wf + 'mock_CF3_z_brent_EDO-format_25_0_256_500.000_WF_Dx.txt'
grid = 256
box = 500.0
cell = box / float(grid)

half_grid = grid / 2

print 'Opening file: ', name_wf, ', box: ', box, ', grid: ', grid, ', cell: ', cell

file_wf = open(name_wf, 'r')
data_wf = file_wf.read().splitlines()
rhos_wf = np.zeros((grid * grid * grid))

this_index = 0

for line_wf in data_wf:
	this_line = line_wf.split()

	for this_col in this_line:
		rhos_wf[this_index] = this_col
		this_index +=1 

xy_ind = []

for ix in range (half_grid - 1, half_grid +1):
	for iy in range(0, grid):
		for iz in range(0, grid):
			this_ind = ix + iy * grid + iz * grid * grid
			xy_ind.append(this_ind)

#print rhos_wf[xy_ind]
plt.contour(rhos_wf[xy_ind], colors='black')
plt.tight_layout()
plt.savefig('test_ft.png')
plt.clf()
plt.cla()


#print this_index
#print grid * grid * grid
