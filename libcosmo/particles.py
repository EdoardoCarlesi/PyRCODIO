from utils import *
import math
import pymp
import numpy as np

def select_particles(x_c, x_p, r, facMpc):
	n_p = len(x_p)
	x_s = []

	print 'Selecting particles within %f kpc/h out of %d' % (r, n_p)

	for ip in range(0, n_p):
		x_p[ip][0] *= facMpc
		x_p[ip][1] *= facMpc
		x_p[ip][2] *= facMpc
		d = vec_distance(x_c, x_p[ip])

		if d < r:
			x_s.append(x_p[ip])

	print 'Found %d particles within the sphere.' % len(x_s)
	return x_s

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
