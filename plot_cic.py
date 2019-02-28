from scipy.ndimage import gaussian_filter
import scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

grid = 128; sigma = 0.75
#grid = 256; sigma = 0.5


path_cic='/home/eduardo/CLUES/DATA/CF3/500/'
#model='CF3_YH_h78'
#model='CF3_YH_h76nocmb'
#model='CF3_RG_1000'
#model='CF3_RG_1200'
model='CF3_RG_1500'
#model='CF3_Brent_v0'
#model='CF3_YH_v4'
subrun='71_00'
name_cic='cic_snapshot_127_'+ str(grid)+'.txt'

name_cic = path_cic + model + '/' + subrun + '/' + name_cic
file_out = 'cic_' + model + '_' + subrun + '.png'

box = 500.0
cell = box / float(grid)

half_grid = grid / 2
print 'Opening file: ', name_cic, ', box: ', box, ', grid: ', grid, ', cell: ', cell

file_cic = open(name_cic, 'r')
data_cic = file_cic.read().splitlines()
rhos_cic = np.zeros((grid * grid * grid))

this_index = 0

x_coord = np.zeros((grid))
y_coord = np.zeros((grid))
z_value = np.zeros((grid, grid))

full_data = np.zeros((grid, grid, grid))

for i in range(0, grid):
	x_coord[i] = i * cell - box * 0.5
	y_coord[i] = i * cell - box * 0.5

index = 0; ix = 0; iy = 0; iz = 0
for line_cic in data_cic:

	if iz < grid +2:
		this_line = line_cic.split()
		line = []

		for col in this_line:
			line.append(float(col))
		
		ix = int (line[0] / cell)
		iy = int (line[1] / cell)
		iz = int (line[2] / cell)

		full_data[ix][iy][iz] = np.log10(line[3])
		#full_data[ix][iy][iz] = line[3]

		'''
		print line[0], line[1], line[2], ix, iy, iz, line[3]

		index += 1
	
		iz = index / (grid * grid)
		iy = index / grid - iz * grid
		ix = index - iy * grid - iz * grid * grid

		if iz == half_grid:
			x_coord[ix] = line[0] - box * 0.5
			y_coord[iy] = line[1] - box * 0.5
			z_value[ix][iy] = line[3]
		'''
	else:
		print 'Breaking...'
		break

#gauss_data = scipy.ndimage.filters.gaussian_filter(full_data[:][half_grid][:], sigma)
gauss_data = gaussian_filter(full_data[:][half_grid][:], sigma)
#gauss_data = full_data[:][half_grid][:]

#levels = np.linspace(0.0, 500.0, num=10)
#levels = np.logspace(-1.0, 3.0, num=10)

#levels = [-10.0, -5.0, 0.0, 1.0, 5.0, 10.0, 50.0, 200.0]
#levels = [0.0, 0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 250.0, 500.0]
levels = [0.0, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 1000.0]
level = [1.0] #, 1.0, 200.0]
#level = [1.0, 10.0, 100.0]
limit = 150

ax = plt.gca()
ax.set_ylim([-limit, limit])
ax.set_xlim([-limit, limit])

#plt.contourf(x_coord, y_coord, z_value, levels=levels, cmap="jet")
#plt.contourf(x_coord, y_coord, full_data[:][half_grid-1][:], locator=ticker.LogLocator(), cmap="jet")
plt.contourf(x_coord, y_coord, gauss_data, cmap="jet")
#plt.contourf(x_coord, y_coord, gauss_data, locator=ticker.LogLocator(), cmap="jet")
#plt.contourf(x_coord, y_coord, full_data[:][half_grid-1][:], levels=levels, cmap="jet")
plt.title(model)
#plt.colorbar(ticks=levels)
#plt.colorbar()
#plt.contour(x_coord, y_coord, gauss_data, colors='black', levels=level, linewidths=[0.5, 1.5, 1.0])
plt.contour(x_coord, y_coord, gauss_data, colors='black', levels=level, linewidths=[1.5])
plt.tight_layout()
plt.savefig(file_out)
plt.clf()
plt.cla()


