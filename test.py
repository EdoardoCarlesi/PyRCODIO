#!/usr/bin/python

from libcosmo.track_halos import *
from libcosmo.utils import *
from libcosmo.units import *
from libcosmo.halo import *
from libcosmo.find_halos import *
from libcosmo.particles import *
from libcosmo.std_lg_plot import *
from libio.read_ascii import *
from libio.find_files import *
from pygadgetreader import *

#import matplotlib.pyplot as plt
import time

#base_path = '/home/eduardo/CLUES/DATA/1024/00_06/LV/'
#snap_file = base_path + 'snapshot_020'
#ahf_file = base_path + 'snapshot_020.0000.z0.000.AHF_halos.pttype1'

base_path = '/home/edoardo/CLUES/DATA/SIMULATIONS/LGF/2048/00_06/00/'

zs = ahf_redshifts(base_path)
ss = ahf_snapshots(base_path)

print ss

#base_path = '/home/eduardo/CLUES/DATA/'
#ahf_file = '/home/eduardo/CLUES/DATA/1024/00_06/02/snapshot_054.0000.z0.000.AHF_halos'
#ahf_file = '/home/eduardo/CLUES/DATA/1024/00_06/02/snapshot_054.0000.z0.000.AHF_halos'
#ahf_file = base_path + '/512/17_01/snapshot_054.0000.z0.000.AHF_halos'
#snap_file = base_path + '/512/17_01/snapshot_054'
#snap_file = base_path + '512/00_06_LV/snapshot_019'
#ahf_file = base_path + '512/00_06_LV/snapshot_019.0000.z0.000.AHF_halos'
ahf_file = base_path + '/2048/00_06/01/snapshot_054.0000.z0.000.AHF_halos'
snap_file = base_path + '/2048/00_06/01/snapshot_054'
#ascii_file = '/home/eduardo/CLUES/DATA/output/lg_candidates_2048_00_06.txt'
#ahf_file = base_path + '/512/17_01/snapshot_054.AHF_halos'
#ahf_file = base_path + '/512/17_00/snapshot_054.z0.000.AHF_halos'
#ahf_file = base_path + 'snapshot_020.0000.z0.000.AHF_halos'

x = [50000., 50000., 50000.]
v = [400., -200., 100]

dMyrs = 1000.0

nx = new_x(x, v, dMyrs)
d = distance(nx, x)

print(nx)
print(d)

'''
ahf = read_ahf(ahf_file)

center = [50000., 50000., 50000.]
d_max = 7000.
r_iso = 2000.
r_max = 1500.
r_min = 250.
ratio_max = 2.5
vrad_max = 0.0
m_min = 3.e+11
m_max = 5.e+12


reduce_fac = 8
lg_model = LocalGroupModel(d_max, r_iso, r_max, r_min, m_max, m_min, ratio_max, vrad_max)

#hs = []
#hs.append(ahf[0])
#hs.append(ahf[2])

#hs = find_halos_point(center, ahf, radius)
hs = find_lg(ahf, lg_model)
subs = find_halos(hs[0], ahf, radius)


(x0, m0, mt) = locate_virgo(ahf)

print x0, m0, mt

for ih in range(0, len(hs)):
	print hs[ih].info()

plot_lglv(snap_file, ahf, 'output_00_06_LV.png', hs[2], hs[3], x0, 8, 1)

center = [50000.0] * 3
rads = 7000.

halos = find_halos_point(center, ahf, rads)

nh = len(halos)

for ih in range(0, nh):
	if halos[ih].m > 4.e+11 :
		print halos[ih].info()

n_pts = 10

points = rand_points_sphere(n_pts, [0.0, 0.0, 0.0], 5.0)
masses = np.zeros((n_pts))

for im in range(0, n_pts):
	masses[im] = 1.0

(evals, evecs) = moment_inertia(points, masses)

print evals
nh = len(hs)

nlg = nh / 2

print 'Nhalos in center: ', nh
m512 = simu_pmass(100, 2048)
print 'hires pmass: %e' % m512

for ih in range(0, nlg):
	lg1 = hs[2 * ih]
	lg2 = hs[2 * ih + 1]
	lgcom = center_of_mass([lg1.m, lg2.m], [lg1.x, lg2.x])
	d12 = distance(lg1.x, lg2.x)
	dcom = distance(center, lgcom)
	rate = rate_lg_pair(lg1, lg2, center)
	print lg1.info()
	print lg2.info()

	print 'Com     : ', lgcom
	print 'Centerbo: ', dcom
	print 'Distance: ', d12 / 0.67
	print 'V_radial: ',  vel_radial(lg1.x, lg2.x, lg1.v, lg2.v) + d12 * 0.067
	print '%.3e, %d, %.5e, %.3f' % (lg1.m, lg1.npart, (lg1.m/lg1.npart), lg1.distance(center))
	print '%.3e, %d, %.5e, %.3f' % (lg2.m, lg2.npart, (lg2.m/lg2.npart), lg2.distance(center))
	print 'Rated: %f\n' % rate

plot_pos = "false"
ptype = 1

lgs = read_lgs(ascii_file)

(bin_m, bin_n) = bin_lg_sub(lgs)

print bin_m
print bin_n

#plot_lg(snap_file, 'output_lg.png', hs[0], hs[1], reduce_fac, ptype, plot_pos)
#lgcom = center_of_mass([hs[0].m, hs[1].m], [hs[0].x, hs[1].x])

#print hs[0].info()

#ptypes = 4

#(xv, mv, mvtot) = locate_virgo(ahf)

#plot_lglv(snap_file, ahf, "output.png", hs[0], hs[1], xv, reduce_fac, ptypes)

#red = readhead(snap_file, 'redshift')
#time = readhead(snap_file, 'time')
#dmPart1 = readsnap(snap_file, 'pos', 1)
#dmPart2 = readsnap(snap_file, 'pos', 2)

array = np.zeros((3, 100), dtype='float')

rType1 = 5000.
rType2 = 15000.
shellZ = 2000.

print len(dmPart1)

npymp = 12

ps1 = select_particles(lgcom, dmPart1, rType1, facMpc)

#ps2 = select_particles(lgcom, dmPart2, rType2, facMpc)

(xa, xb) = find_slice(ps1, 2, hs[0].x[2], shellZ)
t3 = time.clock()
print 'Time past: ', (t3-t2)

xamax = np.amax(xa)
xamin = np.amin(xa)
xbmax = np.amax(xb)
xbmin = np.amin(xb)

shellXY = 1500.

xamax = lgcom[0] + shellXY
xamin = lgcom[0] - shellXY
xbmax = lgcom[1] + shellXY
xbmin = lgcom[1] - shellXY

facX = 10000.
x_min = facX * math.floor(xamin / facX)
x_max = facX * math.ceil(xamax / facX)
y_min = facX * math.floor(xbmin / facX)
y_max = facX * math.ceil(xbmax / facX)

x_min = xamin 
x_max = xamax 
y_min = xbmin
y_max = xbmax

print xamax, xamin
print xbmax, xbmin

ptsize = 0.1

plt.axis([x_min, x_max, y_min, y_max])
plt.scatter(xa, xb, s=ptsize, c=2)
plt.show()


#print ps1
#h0 = ahf[0]
#radius = h0.r * 2.
#hs = find_halos(h0, ahf, radius)
#hs = find_lg(ahf, center, radius, m_min, r_min, r_max)
n_m = 6
masses = [1.0] * n_m

coords = [[0] * 3] * n_m

#print coords
#print masses

for ip in range(0, n_m):
	coords[ip][0] = ip
	coords[ip][1] = ip * 2
	coords[ip][2] = ip * 3
	

print coords

moment_inertia(coords, masses)
#print masses

#nh = len(hs)
#for i in range(0, nh):
#	print hs[i].info()

#print 'found subhalos: %d' % len(hs)
#nu.particles.ellipticities(0,0)

#r = [50000, 50000, 50000]
#print ahf[0].distance(r)
#print "Halo mass: %.3f \n" % h.m
'''
