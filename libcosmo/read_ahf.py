from halo import Halo
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
			halo = Halo(idn, mass, pos, vel, rvir, nsub, npart)
			halos_ahf.append(halo)
			count += 1

	n_lines = count
	print "Read %s with a total of %d lines" % (file_name, n_lines)

	return halos_ahf
