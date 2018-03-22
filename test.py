#!/usr/bin/python

from libcosmo.utils import *
from libcosmo.halo import Halo
from libcosmo.read_ahf import read_ahf
from libcosmo.find_halos import *
from libcosmo.particles import *

from pygadgetreader import *
from libgsr.gsr import *
import matplotlib.pyplot as plt

import time

#from stdlib.lib.python.gsn_numpy_util import * 
#import gsn_numpy_util as nu
#import * 
#as nu
#import stdlib.lib.python.gsn_numpy_util as nu

base_path = '/home/eduardo/CLUES/DATA/'
#ahf_file = '/home/eduardo/CLUES/DATA/1024/00_06/02/snapshot_054.0000.z0.000.AHF_halos'
#ahf_file = '/home/eduardo/CLUES/DATA/1024/00_06/02/snapshot_054.0000.z0.000.AHF_halos'
ahf_file = base_path + '/512/17_00/snapshot_054.0000.z0.000.AHF_halos'
snap_file = base_path + '/512/17_00/snapshot_054'

ahf = read_ahf(ahf_file)

center = [50000., 50000., 50000.]
radius = 8000.
r_iso = 2000.

r_max = 1500.
r_min = 250.
m_min = 5.e+11
facMpc = 1000.
shellSize = 2000.	# plus or minus
hs = find_lg(ahf, center, radius, r_iso, m_min, r_min, r_max)
lgcom = center_of_mass([hs[0].m, hs[1].m], [hs[0].x, hs[1].x])

print 'Reading gadget snapshot: ', snap_file
#snap = Snapshot(snap_file)

#head = readhead(snap_file, 'header')
#red = readhead(snap_file, 'redshift')
#time = readhead(snap_file, 'time')

dmPart1 = readsnap(snap_file, 'pos', 1)
dmPart2 = readsnap(snap_file, 'pos', 2)

array = np.zeros((3, 100), dtype='float')

rType1 = 5000.
rType2 = 15000.
shellZ = 2000.

print len(dmPart1)

npymp = 12

t1 = time.clock()
ps1 = select_particles(lgcom, dmPart1, rType1, facMpc)
t2 = time.clock()

print 'Time past: ', (t2-t1)
#ps2 = select_particles(lgcom, dmPart2, rType2, facMpc)

(xa, xb) = find_slice(ps1, 2, hs[0].x[2], shellZ)
t3 = time.clock()
print 'Time past: ', (t3-t2)

'''
xamax = np.amax(xa)
xamin = np.amin(xa)
xbmax = np.amax(xb)
xbmin = np.amin(xb)
'''

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
'''
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
'''
#print masses

#nh = len(hs)
#for i in range(0, nh):
#	print hs[i].info()

#print 'found subhalos: %d' % len(hs)
#nu.particles.ellipticities(0,0)

#r = [50000, 50000, 50000]
#print ahf[0].distance(r)
#print "Halo mass: %.3f \n" % h.m
