from libio.read_ascii import *
from libcosmo.halo import *
from libcosmo.mtree import *
from libSQL.sqllib import *
from libSQL.mtree import *
import subprocess as sbp
import pandas as pd
import pickle
import os

# How many catalogs do we want to read
endSnaps = 127
totSnaps = 10
iniSnaps = endSnaps - totSnaps
nSnaps = endSnaps - iniSnaps

# File properties
thisDB = '/home/eduardo/CLUES/DATA/HESTIA/8192/trees/hestia_trees.db'
rootPath = '/home/eduardo/CLUES/DATA/HESTIA/8192/catalogs/'
#rootPath = '/home/eduardo/CLUES/DATA/LGF/1024/FIX/'
subPath = '37_11'
outPath = 'output/full_trees/' + subPath + '/'
#baseAHF = 'snapshot_'
baseAHF = 'HESTIA_100Mpc_8192_'+subPath+'.'
suffAHF = '.AHF_halos'

#z_out = rootPath + '../output_z.txt'
z_out = '/home/eduardo/CLUES/DATA/HESTIA/output_hestia.txt'
z_f = open(z_out, 'r')
z_s = z_f.read().split()

# List of halos without direct progenitor
orphHalos = []

# Start reading from the final snapshot
thisNSnap = '%03d' % endSnaps
thisZSnap = z_s[endSnaps]
thisAHF = rootPath + subPath + '/' + baseAHF + thisNSnap + '.z' + thisZSnap + suffAHF

print('Reading AHF halo catalog: ', thisAHF)
descHalos = read_ahf(thisAHF)
idDescHalos = makeHaloDictionary(descHalos)

# Read the merger trees from a SQL database
columnReadPT = 'allNumPart'
columnReadID = 'allHaloIDs'

newSql = SQL_IO(thisDB, totSnaps)
allTrees = []

# Only select certain halos at z=0
#boxCenter = [5.0e+4, 5.0e+4, 5.0e+4]
boxCenter = [4.6e+4, 5.0e+4, 4.7e+4]
radiusMax = 5.0e+3

print("Selecting a subset from a total %d halos..." % (len(descHalos)))

for thisHalo in descHalos:

    # Only look for halos within a given sphere
    if thisHalo.distance(boxCenter) < radiusMax:
        pts = newSql.select_tree(thisHalo.ID, columnReadPT)
        ids = newSql.select_tree(thisHalo.ID, columnReadID)
        thisTree = MergerTree(totSnaps, pts, ids, subPath)
        allTrees.append(thisTree)
try:
	newSql.close()

except:
	'Do nothing'

nTrees = len(allTrees)
allFilesOut = []
allHaloInfo = [[" "] * (nTrees)] * (nSnaps +1)

print("\nDone.\nPrinting %d trees to %s." % (nTrees, outPath))

iTree = 0
for thisTree in allTrees:
    thisID = thisTree.IDs[0]
    thisHalo = descHalos[idDescHalos[thisID]]
    allHaloInfo[0][iTree] = thisHalo.line.rstrip('\n')
    thisName = outPath + 'halo_' + str(thisID) + '_tree.dat'
    thisFile = open(thisName, 'w')
    allFilesOut.append(thisFile)
    print(allHaloInfo[0][iTree], file=thisFile)
    iTree += 1

# Loop on all the snapshots
for jSnap in range(0, nSnaps):
    iSnap = endSnaps - jSnap - 1
    thisNSnap = '%03d' % iSnap
    thisZSnap = z_s[iSnap]
    thisAHF = rootPath + subPath + '/' + baseAHF + thisNSnap + '.z' + thisZSnap + suffAHF

    print('Reading AHF halo catalog: ', thisAHF)
    progHalos = read_ahf(thisAHF)
    idProgHalos = makeHaloDictionary(progHalos)

    iTree = 0
    for thisTree in allTrees:
        thisID = thisTree.IDs[jSnap+1]
        prevID = thisTree.IDs[jSnap]

        # This is an orphan halo 
        if thisID == prevID:
            allHaloInfo[jSnap+1][iTree] = allHaloInfo[jSnap][iTree]

        # If not orphan then copy the new one 
        else:
            thisHalo = progHalos[idProgHalos[thisID]]
            allHaloInfo[jSnap+1][iTree] = thisHalo.line.rstrip('\n')

        iTree += 1

    # Now reset all the indexes
    idProgHalos.clear()
    progHalos.clear()

# Close all files once the reconstructed trees have been written
for thisFile in allFilesOut:
        thisFile.close()
