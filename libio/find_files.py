# These functions find and read AHF files using a series of bash scripts included in the script/ subfolder
import os
import subprocess
import numpy as np

# Read and identify available snapshot numbers
def ahf_snapshots(base_path):
	find_s_sh='scripts/find_s.sh'
	exec_s=find_s_sh+' '+base_path
	proc = subprocess.Popen(exec_s, stdout=subprocess.PIPE, shell=True)
	(f_out, f_err) = proc.communicate()

	file_s = open(f_out.rstrip(), 'r')
	line_s = file_s.readlines()
	ns = len(line_s)
	ss = []

        for js in range(0, ns):
                this_line = line_s[js].rstrip()
		ss.append(str(this_line))
	return ss


# Read and identify redshifts
def ahf_redshifts(base_path):
	find_z_sh='scripts/find_z.sh'
	exec_z=find_z_sh+' '+base_path
	proc = subprocess.Popen(exec_z, stdout=subprocess.PIPE, shell=True)
	(f_out, f_err) = proc.communicate()
	
	print '-'
	print f_out
	print '-'

	file_z = open(f_out.rstrip(), 'r')
	line_z = file_z.readlines()
	nz = len(line_z)
	zs = []

        for iz in range(0, nz):
                this_line = line_z[iz].rstrip()
		zs.append(this_line)
	return zs

