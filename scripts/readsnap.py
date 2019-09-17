#!/usr/bin/python

from pygadgetreader import *
import sys

file_in=sys.argv[1]

print(file_in)

z   	= readhead(file_in, 'redshift')
box 	= readhead(file_in, 'boxsize')
n_files = readhead(file_in, 'nfiles')
n_tot	= readhead(file_in, 'npartThisFile')

print('z:%5.4f, box:%5.4f, nFiles:%d, nTot:%s' % (z, box, n_files, n_tot))

