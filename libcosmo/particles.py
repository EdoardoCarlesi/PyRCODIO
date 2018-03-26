from utils import *
from units import *
from halo import *
import math
import numpy as np


def next_ip_slab(ip, r, h):
	skip1 = 150
	skip2 = 20
	skip3 = 1

	d1 = r
	d2 = 2. * r
	d3 = 4. * r

	# We skip reading clusters of particles faaaaaar away from the desired point
	if h > d3:
		ip += skip1

	elif h > d2 and h < d3: 
		ip += skip2

	else:
		ip += skip3

	return ip


def next_ip_sphere(ip, r, d):
	skip1 = 250
	skip2 = 20
	skip3 = 1

	d1 = r
	d2 = 10. * r
	d3 = 100. * r

	# We skip reading clusters of particles faaaaaar away from the desired point
	if d > d3:
		ip += skip1

	elif d > d2 and d < d3: 
		ip += skip2

	else:
		ip += skip3

	return ip


def select_particles_halo(halo, x_p):
	n_p = len(x_p)
	n_h = halo.npart
	x_com = halo.x
	r = halo.r
	x_a = []
	x_b = []
	x_c = []

	print 'Selecting %d particles within %f kpc/h out of %d' % (n_h, r, n_p)

	units = x_p[n_p-1][0] + x_p[n_p-1][1] + x_p[n_p-1][2]

	if units > 5000. :
		facMpc = 1.0
	else:
		facMpc = 1000.0
	
	ip = 0
	ih = 0
	loc_x = [0.0] * 3

	while (ih < (n_h) and ip < (n_p-1)):
		loc_x[0] = x_p[ip][0] * facMpc
		loc_x[1] = x_p[ip][1] * facMpc
		loc_x[2] = x_p[ip][2] * facMpc

		d = distance(x_com, loc_x)

		if d < r:
#			print ip, d, r, x_com, x_p[ip]
			x_a.append(x_p[ip][0])
			x_b.append(x_p[ip][1])
			x_c.append(x_p[ip][2])
			ih += 1

		ip = next_ip_sphere(ip, r, d)

	print 'Found %d particles within the sphere, read %d particles in total.' % (len(x_a), ip)
	return (x_a, x_b, x_c)


def select_particles(x_com, x_p, r):
	n_p = len(x_p)
	x_a = []
	x_b = []
	x_c = []

	print 'Selecting particles within %f kpc/h out of %d' % (r, n_p)

	units = x_p[n_p-1][0] + x_p[n_p-1][1] + x_p[n_p-1][2]

	print units, x_com

	if units > 5000. :
		facMpc = 1.0
	else:
		facMpc = 1000.0
	
	loc_x = [0.0] * 3

	while ip < (n_p-1):

		loc_x[0] = x_p[ip][0] * facMpc
		loc_x[1] = x_p[ip][1] * facMpc
		loc_x[2] = x_p[ip][2] * facMpc
	
		d = distance(x_com, loc_x)

		if d < r:
			x_a.append(x_p[ip][0])
			x_b.append(x_p[ip][1])
			x_c.append(x_p[ip][2])

		ip = next_ip_sphere(ip, r, d)

	print 'Found %d particles within the sphere.' % len(x_a)
	return (x_a, x_b, x_c)


def find_slab(x_p, axis, center, min_ab, side, thick, reduce_fac):
	halfside = 0.5 * side
	n_p = len(x_p)
	x_a = []
	x_b = []

	units = x_p[n_p-1][0] + x_p[n_p-1][1] + x_p[n_p-1][2]
	
	if units > 5000. :
		facMpc = 1.0
	else:
		facMpc = 1000.0
	
	loc_x = [0.0] * 3
	
	# These are the 2-D coordinates
	min_a = min_ab[0]
	min_b = min_ab[1] 
	max_a = min_ab[0] + side
	max_b = min_ab[1] + side

	print 'Selecting particles within %.3f side and %.3f height around %s' % (side, thick, center)
	
	a = (axis+1) % 3
	b = (axis+2) % 3
	
	ireduce = 0	
	ip = 0

	while ip < (n_p-1):
		loc_x[0] = x_p[ip][0] * facMpc
		loc_x[1] = x_p[ip][1] * facMpc
		loc_x[2] = x_p[ip][2] * facMpc

		delta_c = math.fabs(loc_x[axis] - center[axis])

		if delta_c < thick:
			delta_a = math.fabs(loc_x[axis] - center[axis])
			delta_b = math.fabs(loc_x[axis] - center[axis])

			if delta_b < halfside and delta_a < halfside:
					x_a.append(x_p[ip][a])	
					x_b.append(x_p[ip][b])

		ip = next_ip_slab(ip, thick, delta_c)
		ip += reduce_fac

	print 'Found %d particles within the slab, reduced by a factor of %d.' % (len(x_a), reduce_fac)
	return (x_a, x_b)


def simu_pmass(box, npart):
	pmass0 = 5.26083e+09
	box0 = 100.0
	npart0 = 256

	pmass = pmass0 * (box/box0) * pow((float(npart0)/float(npart)), 3)

	return pmass


def find_slice(x_p, axis, center, shell):
	n_p = len(x_p)
	x_a = []
	x_b = []

	print 'Selecting particles within %f distance out of %d' % (shell, n_p)
	
	a = (axis+1) % 3
	b = (axis+2) % 3
	
	for ip in range(0, n_p):
		delta_shell = math.fabs(x_p[ip][axis] - center)
		if delta_shell < shell:
			x_a.append(x_p[ip][a])	
			x_b.append(x_p[ip][b])	

	print 'Found %d particles within the shell.' % len(x_a)
	return (x_a, x_b)


'''
import pymp
def find_slab_omp(x_p, axis, center, min_ab, side, thick, reduce_fac, npymp):
	halfside = 0.5 * side
	n_p = len(x_p)
	x_a = []
	x_b = []

	units = x_p[n_p-1][0] + x_p[n_p-1][1] + x_p[n_p-1][2]
	
	if units > 5000. :
		facMpc = 1.0
	else:
		facMpc = 1000.0
	
	loc_x = [0.0] * 3
	
	# These are the 2-D coordinates
	min_a = min_ab[0]
	min_b = min_ab[1] 
	max_a = min_ab[0] + side
	max_b = min_ab[1] + side

	print 'Selecting particles within %.3f side and %.3f height around %s' % (side, thick, center)
	
	a = (axis+1) % 3
	b = (axis+2) % 3
	
	ireduce = 0	

	with pymp.Parallel(npymp):

		for ip in p.range(0, n_p):
			loc_x[0] = x_p[ip][0] * facMpc
			loc_x[1] = x_p[ip][1] * facMpc
			loc_x[2] = x_p[ip][2] * facMpc

			delta_c = math.fabs(loc_x[axis] - center[axis])
		
			if delta_c < thick:

				delta_a = math.fabs(loc_x[axis] - center[axis])
				delta_b = math.fabs(loc_x[axis] - center[axis])

				if delta_b < halfside and delta_a < halfside:
					if ireduce == reduce_fac:
						x_a.append(x_p[ip][a])	
						x_b.append(x_p[ip][b])
						ireduce = 0
					else:
						ireduce += 1
			

	print 'Found %d particles within the slab, reduced by a factor of %d.' % (len(x_a), reduce_fac)
	return (x_a, x_b)



def select_particles_omp(x_c, x_p, r, facMpc, npymp):
	n_p = len(x_p)
	x_s = pymp.shared.list()

	print 'Selecting particles within %f kpc/h out of %d' % (r, n_p)

	with pymp.Parallel(npymp) as p:
		for ip in p.range(0, n_p):
			x_p[ip][0] *= facMpc
			x_p[ip][1] *= facMpc
			x_p[ip][2] *= facMpc
			d = vec_distance(x_c, x_p[ip])

		#		with p.lock:
			if d < r:
				x_s.append(x_p[ip])

	print 'Found %d particles within the sphere.' % len(x_s)
	print x_s[0], x_s[100]
	return x_s


def find_slice_omp(x_p, axis, center, shell, facMpc, npymp):
	n_p = len(x_p)
	x_a = pymp.shared.list()
	x_b = pymp.shared.list()
	x_s = pymp.shared.list()

	print 'Selecting particles within %f distance out of %d' % (shell, n_p)
	
	a = (axis+1) % 3
	b = (axis+2) % 3
	
	with pymp.Parallel(npymp) as p:
		for ip in p.range(0, n_p):
			delta_shell = math.fabs(x_p[ip][axis] - center)
			if delta_shell < shell:
				x_s.append(x_p[ip])
				#x_a.append(x_p[ip][a])	
				#x_b.append(x_p[ip][b])	

#	print x_a[1], x_a[111]
	print 'Found %d particles within the shell.' % len(x_s)
	return (x_a, x_b)
'''
