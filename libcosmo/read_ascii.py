from halo import *
import os
import numpy

def read_ahf(file_name):
	# Open file
	file_ahf = open(file_name, 'r')

	line = file_ahf.readline()
	halos_ahf = []
	count = 0
	
	while line:
    		line = file_ahf.readline()
		line = line.strip()
		column = line.split()
		n_col = len(column)
		
		if n_col > 1:
			# Read halo properties
			idn = long(column[0])
			mass = float(column[3])
			pos = [float(column[5]), float(column[6]), float(column[7])]
			vel = [float(column[8]), float(column[9]), float(column[10])]
			rvir = float(column[11])
			nsub = int(column[2])
			npart = int(column[4])

			# Initialize and append halo to the list
			halo = Halo()
			halo.initialize(idn, mass, pos, vel, rvir, nsub, npart)
			halos_ahf.append(halo)
			count += 1

	n_lines = count
	print "Read %s with a total of %d lines" % (file_name, n_lines)

	return halos_ahf

# Read the LG asciis saved somewhere
def read_lgs(file_name):
	file_txt = open(file_name, 'r')

	line = file_txt.readline()
	lgs_txt = []
	count = 0
	
	lgs = []

	com = [0.0] * 3
	pos = [0.0] * 3
	vel = [0.0] * 3


	while line:
    		line = file_txt.readline()
		line = line.strip()
		column = line.split()
		n_col = len(column)
		
		if n_col > 1:
			lg1 = Halo()
			lg2 = Halo()

			# Read halo properties
			id0 = long(column[0])
			id1 = long(column[1])
			m0 = float(column[2])
			m1 = float(column[3])
			dist = float(column[4])
			vrad = float(column[5])
			nsub0 = int(column[6])
			nsub1 = int(column[7])
			code = column[8]
			com[0] = float(column[9]) * 1000.
			com[1] = float(column[10]) * 1000.
			com[2] = float(column[11]) * 1000.

			# Most of this is initialized to dummy variables
			lg1.initialize(id0, m0, pos, vel, 0.0, nsub0, 0)
			lg2.initialize(id1, m1, pos, vel, 0.0, nsub1, 0)

			lg = LocalGroup(code)
			lg.init_halos(lg1, lg2)
			lg.vrad = vrad
			lg.r = dist
			lg.com = com

			lgs.append(lg)

			del lg1
			del lg2

	return lgs
