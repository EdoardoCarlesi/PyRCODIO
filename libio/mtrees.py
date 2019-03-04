# Merger tree input - output

import numpy as np
import sys
from .read_ascii import *

def mah_main_file_name(out_path, num_run, i_main):
	main_code = '%02d' % i_main
	file_name = out_path + num_run + '_' + main_code + '_main.txt'
	return file_name

def mah_sub_file_name(out_path, num_run, i_main, i_sub):
	main_code = '%02d' % i_main
	sub_code = '%03d' % i_sub
	file_name = out_path + num_run + '_' + main_code + '_' + sub_code + '_sub.txt'
	return file_name

def n_mah_files(out_path, num_run):
	n_sub_all = []

	lines_main = 'ls ' + out_path + '/' + num_run + '*_main.txt | wc -l'
        n_main = os.popen(lines_main).read()
	n_main = int(n_main.strip())

	for i_main in range(0, n_main):
		main_code = '%02d' % i_main
		lines_sub = 'ls ' + out_path + '/' + num_run + '_' + main_code+'_*_sub.txt | wc -l'
        	n_sub = os.popen(lines_sub).read()
		n_sub = int(n_sub.strip())
		n_sub_all.append(n_sub)

	return (n_main, n_sub_all)

