# This function creates a mask to be used for the generation of zoom ICs.
# It requires: the lg_center of mass, the name of the initial condition file and z = 0 snapshot, the number of files per IC/snapshot
# It dumps a Lare_NUM.tbl file which then is used to generate a mask file suitable for ginnungagap
# It requires also two scripts as an input: lare_find_sh, which is used to find out the lagrangian region, and lare_gen_sh, which is 
# used to generate the Lare_??_??.dat file necessary for ginnungagap

import os
import subprocess
import numpy as np

def find_lare(lg_com, lare_r, file_ics, file_z0, n_ics, n_z0, file_out, lare_find_sh):
	# The script that looks for the lagrangian region takes in as an input:
	# (1) file_z0 (2) n_file_z0 (3) output file (4) laRe radius (5,6,7) LG center of mass (8) ic file (9) n ic files

	print 'Find lare: %s %f %s %s %d %d %s %s ' % (lg_com, lare_r, file_ics, file_z0, n_ics, n_z0, file_out, lare_find_sh)
	file_tmp = 'out.tmp'
	xLG = lg_com[0]
	yLG = lg_com[1]
	zLG = lg_com[2]

	# Execute the script:
	lare_find_exec = lare_find_sh + " " + file_z0  + " " + str(n_z0) + " "  + file_out + " "  + str(lare_r) + " "  + \
			 str(xLG) + " " + str(yLG) + " " + str(zLG) + " " + file_ics + " " + str(n_ics)
 
	#coords = os.system(lare_find_exec)
	proc = subprocess.Popen(lare_find_exec, stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

	# The output is read as a huge single string. This way we only parse the last 4 characters (x,y,z,r)
	n_lines = len(out)
	n_chars = 27

	lare_properties = out[n_lines-n_chars:n_lines-1]

	# We might not need to convert this into floats - a string works fine for the gen_mask script
	all_lare_x = np.fromstring(lare_properties, dtype=float, sep=' ')
	
	n_lare_x = len(all_lare_x)

	return all_lare_x

def gen_mask(this_run, lare_x, lare_gen_sh, dir_sh):
	print 'Generating mask with %s\n' % lare_gen_sh
	
	lare_gen_exec = lare_gen_sh + " " + this_run + " " + lare_x + " " + dir_sh
	print 'Mask input %s \n' % lare_gen_exec 
	print lare_gen_exec
	os.system(lare_gen_exec)
	







