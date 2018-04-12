import sys
sys.path.append('../libcosmo/')
from libcosmo.halo import *
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
			angmom = [float(column[21]), float(column[22]), float(column[23])]

			# Initialize and append halo to the list
			halo = Halo()
			halo.initialize(idn, mass, pos, vel, rvir, nsub, npart)
			halo.l = angmom
			halos_ahf.append(halo)
			count += 1

	n_lines = count
	print "Read %s with a total of %d lines" % (file_name, n_lines)

	return halos_ahf

# Reading AHF particle file
def read_particles(file_name):
	ids = []	# List containing all halo IDs & Particle number per halo ID - each list member is a 2 elements array
	parts = []	# Each list member contains an array with all the particle ids
	count_p = 0	# Total number of particles per halo
	count_h = 0	# Total number of haloes in file
	count_l = 0	# Simply count the lines

	this_np = 0
	tot_h = 0	# Double check that the total number of haloes is matched with the counter

	file_part = open(file_name, 'r')

	# First read the header, containing the total number of haloes
	line = file_part.readline()
	line = line.strip()
	column = line.split()
	tot_h = int(column[0])
	
	lines_command = 'wc -l ' + file_name
	out_os = os.popen(lines_command).read()
	(tot_l, fname) = out_os.split()
	tot_l = long(tot_l)
	print 'Reading %s with %ld lines and %d halos. ' % (file_name, tot_l, tot_h)

	while line:
		line = file_part.readline()
		line = line.strip()
		column = line.split()
		
		if count_l > 0 and count_p < this_np:
			this_pid = long(column[0])	# Particle ID
			this_parts.append(this_pid)
			count_p += 1			
			count_l += 1
		else:	
			# All particles have been read in
			if count_p == this_np and count_h > 0:
				this_parts.sort()	# Automatically sort them by ascending order!
				parts.append(this_parts)
			
			# Still reading particle files
			if count_l < tot_l-1:
				this_parts = []
				this_hid = long(column[1])	# Halo ID
				this_np = int(column[0])
				#print 'Line %ld found halo %ld with %d particles' % (count_l, this_hid, this_np)	
				ids.append([this_hid, this_np])
				count_l += 1
				count_h += 1
				count_p = 0			# Reset the particle number
		
	return (ids, parts)


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
