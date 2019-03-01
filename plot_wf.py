import numpy as np
import matplotlib.pyplot as plt


# Wiener Filter file in IceCore format
#path_wf = '/z/carlesi/CLUES/IceCore/output/'
path_wf = '/home/eduardo/CLUES/DATA/WF/'
#root_wf = 'mock_CF3_z_jenny_EDO-format_50_0_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_25_0_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_0_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_1_256_500.000'; do_mock = True
root_wf = 'mock_CF3_z_brent_EDO-format_50_2_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_3_256_500.000'; do_mock = True

#root_wf='CF3_RG_1300_256_500.000'; do_mock = False
#root_wf='CF3_RG_1200_256_500.000'; do_mock = False
#root_wf='CF3_BT_cmb_256_500.000'; do_mock = False
#root_wf='CF3_BT_nocmb_256_500.000'; do_mock = False
#root_wf='CF3_Brent_v0_256_500.000'; do_mock = False
#root_wf='CF3_YH_h74_256_500.000'; do_mock = False
#root_wf='CF3_YH_h75_256_500.000'; do_mock = False
#root_wf='CF3_YH_h76nocmb_256_500.000'; do_mock = False
#root_wf='CF3_YH_h78_256_500.000'; do_mock = False

name_wf = path_wf + root_wf + '_WF_Dx.txt'
file_out = 'wf_' + root_wf + '.png'

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

x_coord = np.zeros((grid))
y_coord = np.zeros((grid))
z_value = np.zeros((grid, grid))

for ix in range (0, grid):
	x_coord[ix] = ix * cell - half_grid * cell
	for iy in range(0, grid):
		y_coord[iy] = iy * cell - half_grid * cell
	
		# Z axis is fixed 
		iz = half_grid

		this_ind = ix + iy * grid + iz * grid * grid
		z_value[ix][iy] = rhos_wf[this_ind]

if do_mock == False:
	levels = [-5.0, -1.0, 0.0, 0.5, 1.0, 2.5, 5.0]
	level = [-5.0, 1.0, 5.0]
	limit = 150
else:
#	levels = [-10.0, -5.0, -1.0, 0.0, 0.5, 1.0, 2.5, 5.0, 10]
	levels = [-5.0, -1.0, 0.0, 0.5, 1.0, 2.5, 5.0]
	level = [-5.0, 1.0, 5.0]
#	levels = [-50.0, -10.0, -5.0, -1.0, 0.0, 1.0, 5.0, 10.0, 50.0]
#	level = [-10.0, 1.0, 10.0]
	limit = 150

ax = plt.gca()
ax.set_ylim([-limit, limit])
ax.set_xlim([-limit, limit])
plt.contourf(x_coord, y_coord, z_value, levels=levels, cmap="jet")
plt.title(root_wf)
plt.colorbar(ticks=levels)
plt.contour(x_coord, y_coord, z_value, colors='black', levels=level, linewidths=[0.5, 1.5, 1.0])
plt.tight_layout()
plt.savefig(file_out)
plt.clf()
plt.cla()


