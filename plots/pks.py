#!/usr/bin/python

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
import matplotlib.pyplot as plt
import numpy as np

file1='../data/planck.dat'


f1 = open(file1, 'r')
line1 = f1.readline()

npart1 = 256.0000000
npart2 = 2. * npart1
npart3 = 2. * npart2

box = 100.000

# Draw vertical lines at three different ks

k1 = 3.14 * (npart1) / box
k2 = 3.14 * (npart2) / box
k3 = 3.14 * (npart3) / box

kMin = 3.14 * (2.0) / box
kMax = 3.14 * (2.0 * npart3) / box

#print k1, k2, k3

k = []; pk = []

k_min = 0.01
k_max = 100.

y_min = 0.00035
y_max = 30000.0

plt.axis([k_min, k_max, y_min, y_max])
plt.xscale('log')#, basex=10)	
plt.yscale('log')#, basex=10)	

plt.rc({'text.usetex': True})
plt.xlabel('$k \\quad [$Mpc$^{-1} h]$')
plt.ylabel('$P(k) \\quad [h^{-3}$ Mpc$^3$]')

while line1:
	try:
		line1 = f1.readline()
		char1 = line1[0]
		
		if char1 != '#':
			column1 = line1.strip()
			column1 = line1.split()
			#print float(column1[1])
			#print (column1[1])
			k.append(float(column1[0]))
			pk.append(float(column1[1]))
		
	except:
		print '\n'

vx1 = [k1, k1, k1]; vy1 = [0.00001, 1.0, 100000.] 
vx2 = [k2, k2, k2]; vy2 = [0.00001, 1.0, 100000.] 
vx3 = [k3, k3, k3]; vy3 = [0.00001, 1.0, 100000.] 
vxMin = [kMin, kMin, kMin]; vyMin = [0.00001, 1.0, 100000.] 
vxMax = [kMax, kMax, kMax]; vyMax = [0.00001, 1.0, 100000.] 

lnw = 2
lnw2 = 5
lnw3 = 7
col1='black'

plt.plot(k, pk, linewidth=lnw, color=col1)
plt.plot(vxMin, vyMin, linewidth=lnw3) 
#plt.plot(vx1, vy1, linewidth=lnw2); f_out='../output/pks_01.png'
#plt.plot(vx2, vy2, linewidth=lnw2); f_out='../output/pks_02.png'
plt.plot(vx3, vy3, linewidth=lnw2); f_out='../output/pks_03.png'
plt.plot(vxMax, vyMax, linewidth=lnw3) 

#f_out='../output/pks_allcuts.png'

plt.savefig(f_out)
plt.clf()
plt.cla()
plt.close()


