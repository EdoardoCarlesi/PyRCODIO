#!/usr/bin/python

#file1='00/Pk-00256-001_snapshot_012.noerr'
#file2='00/Pk-00256-001_snapshot_012.verr'
#file2='01/Pk-00256-001_snapshot_012'

file1='/home/eduardo/CLUES/DATA/ICs/Pk.01_12.deltak.dat'
file2='/home/eduardo/CLUES/DATA/ICs/Pk.09_18.deltak.dat'

f1 = open(file1, 'r')
f2 = open(file2, 'r')

line1 = f1.readline()
line2 = f2.readline()

while line1 and line2:
	try:
		line1 = f1.readline()
		line2 = f2.readline()
		char1 = line1[0]
		char2 = line2[0]
		
		if char1 != '#':
			column1 = line1.strip()
			column1 = line1.split()
	
		if char2 != '#':
			column2 = line2.strip()
			column2 = line2.split()

			k = float(column1[0]) #* 1.e+3 
			k0 = 1.0
			k1 = 1.005
			k2 = 0.995

			print k, float(column2[1])/float(column1[1]), k0, k1, k2

	except:
		print '\n'

