import numpy as np
import matplotlib.pyplot as plt

def struct_wf_data(data_wf, grid):

    this_index = 0
    rhos_wf = np.zeros((grid * grid * grid))

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

            this_ind = int(ix + iy * grid + iz * grid * grid)
            z_value[ix][iy] = rhos_wf[this_ind]
#            print(ix, iy, z_value[ix][iy])

    return [x_coord, y_coord, z_value]  





# Wiener Filter file in IceCore format
#path_wf = '/z/carlesi/CLUES/IceCore/output/'
path_wf = '/home/eduardo/CLUES/DATA/WF/'
#root_wf = 'mock_CF3_z_jenny_EDO-format_50_0_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_25_0_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_0_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_1_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_2_256_500.000'; do_mock = True
#root_wf = 'mock_CF3_z_brent_EDO-format_50_3_256_500.000'; do_mock = True

root_wf1='CF3_RG_1200_256_500.000'; do_mock = False
root_wf2='CF3_RG_1200_256_500.000'; do_mock = False
#root_wf2='CF3_BT_cmb_256_500.000'; do_mock = False
root_wf3='CF3_YH_h75_256_500.000'; do_mock = False

#root_wf='CF3_RG_1300_256_500.000'; do_mock = False
#root_wf='CF3_RG_1200_256_500.000'; do_mock = False
#root_wf='CF3_BT_cmb_256_500.000'; do_mock = False
#root_wf='CF3_BT_nocmb_256_500.000'; do_mock = False
#root_wf='CF3_Brent_v0_256_500.000'; do_mock = False
#root_wf='CF3_YH_h74_256_500.000'; do_mock = False
#root_wf='CF3_YH_h75_256_500.000'; do_mock = False
#root_wf='CF3_YH_h76nocmb_256_500.000'; do_mock = False
#root_wf='CF3_YH_h78_256_500.000'; do_mock = False

name_wf1 = path_wf + root_wf1 + '_WF_Dx.txt'
name_wf2 = path_wf + root_wf2 + '_WF_Dx.txt'
name_wf3 = path_wf + root_wf3 + '_WF_Dx.txt'
file_out1 = 'diff_wf_' + root_wf1 + '.png'
file_out2 = 'diff_wf_' + root_wf2 + '.png'
file_out3 = 'diff_wf_' + root_wf3 + '.png'

grid = 256
box = 500.0
cell = box / float(grid)

half_grid = grid / 2

print('Opening file: ', name_wf1, ', box: ', box, ', grid: ', grid, ', cell: ', cell)
print('Opening file: ', name_wf2, ', box: ', box, ', grid: ', grid, ', cell: ', cell)
print('Opening file: ', name_wf3, ', box: ', box, ', grid: ', grid, ', cell: ', cell)

#file_wf1 = open(name_wf1, 'r')
#data_wf1 = file_wf1.read().splitlines()

file_wf2 = open(name_wf2, 'r')
data_wf2 = file_wf2.read().splitlines()

file_wf3 = open(name_wf3, 'r')
data_wf3 = file_wf3.read().splitlines()

#[x_coord, y_coord, z_value1] = struct_wf_data(data_wf1, grid)
[x_coord, y_coord, z_value2] = struct_wf_data(data_wf2, grid)
[x_coord, y_coord, z_value3] = struct_wf_data(data_wf3, grid)

z_value = np.zeros((grid, grid))

for ix in range(0, grid):
    for iy in range(0, grid):
        z_value[ix][iy] = (z_value3[ix][iy] - z_value2[ix][iy]) #* 10.0   
#        print(z_value[ix][iy], z_value3[ix][iy], z_value2[ix][iy])

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
plt.title(root_wf3+' - '+root_wf2)
plt.colorbar(ticks=levels)
plt.contour(x_coord, y_coord, z_value, colors='black', levels=level, linewidths=[0.5, 1.5, 1.0])
plt.tight_layout()
plt.savefig(file_out2)
plt.clf()
plt.cla()


