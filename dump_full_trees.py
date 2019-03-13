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
totSnaps = 126
iniSnaps = endSnaps - totSnaps
nSnaps = endSnaps - iniSnaps

# File properties
rootPath = '/home/eduardo/CLUES/DATA/HESTIA/8192/catalogs/'
#rootPath = '/home/eduardo/CLUES/DATA/LGF/1024/FIX/'
subPath = '17_11'
thisDB = '/home/eduardo/CLUES/DATA/HESTIA/8192/trees/hestia_trees_' + subPath + '.db'
outPath = 'output/full_trees/' + subPath + '/'
#baseAHF = 'snapshot_'
baseAHF = 'HESTIA_100Mpc_8192_'+subPath+'.'
suffAHF = '.AHF_halos'
dwarfList = '/home/eduardo/CLUES/DATA/HESTIA/dwarfs_' + subPath + '.txt'

dwarfIDs = read_dwarf_list(dwarfList)

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

allHalos = []
selectHalos = []

# Select only halos from the dwarf list
for thisID in dwarfIDs:
    thisHalo = descHalos[idDescHalos[thisID]]
    allHalos.append(thisHalo)

# Read the merger trees from a SQL database
columnReadPT = 'allNumPart'
columnReadID = 'allHaloIDs'

# Open the SQL db containing all the merger tree information
newSql = SQL_IO(thisDB, totSnaps)
allTrees = []

# Only select certain halos at z=0 (not implemented)
#boxCenter = [5.0e+4, 5.0e+4, 5.0e+4]
boxCenter = [4.6e+4, 5.0e+4, 4.7e+4]
radiusMax = 5.0e+3

print("Selecting a subset from a total %d halos..." % (len(descHalos)))

# Only use the halos that have been previously selected
for thisHalo in allHalos:
    pts = newSql.select_tree(thisHalo.ID, columnReadPT)
    ids = newSql.select_tree(thisHalo.ID, columnReadID)
    thisTree = MergerTree(totSnaps, pts, ids, subPath)

    if int(thisTree.nPart[0]) > 0:
        allTrees.append(thisTree)
        selectHalos.append(thisHalo)
    else:
        'Tree not found'
#        print("Could not find the tree for: ", thisTree.IDs, thisHalo.info())

# Close the database
try:
    newSql.close()
except:
    'Do nothing'

nTrees = len(allTrees)

# This list of lists will hold all the data to be printed
allHaloInfo = [] 
for iList in range(0, nTrees):
    allHaloInfo.append([])

print("Done.\nPrinting %d trees to %s" % (nTrees, outPath))

# Save some information on orphan halo tracking
orphName = outPath + 'orphan_halo_info.txt'
orphInfo = open(orphName, 'w')

# Find the trees for the selected halos
iTree = 0
for thisTree in allTrees:
    thisID = thisTree.IDs[0]

    if thisID > 0:
        thisHalo = descHalos[idDescHalos[thisID]]
        allHaloInfo[iTree].append(thisHalo.line.rstrip('\n'))
        iTree += 1
    else:
        " "
        #print("Error IDs: ", thisTree.IDs)
        #print("Error IDs: ", thisTree.nPart)
        #allHaloInfo[iTree].append("000000000")

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
        if len(thisTree.IDs) > 1:
            thisID = int(thisTree.IDs[jSnap+1])
            prevID = int(thisTree.IDs[jSnap])
        else:
            thisID = 0

        # Check that the tree history is not truncated
        if thisID == 0:
            'Do not track this halo, which was too small to be worth following'
            #print(thisTree.nPart[jSnap], thisTree.nPart[jSnap+1])

        else: 
            # This is an orphan halo 
            if thisID == prevID:
                #print("Orphan halo ", iTree, thisID, descHalos[iTree].info())

                try:
                    allHaloInfo[iTree].append(allHaloInfo[iTree][jSnap])
                    print(thisZSnap, jSnap, iTree, thisID, selectHalos[iTree].info(), file=orphInfo)
                except:
                    print('Error:', thisZSnap, jSnap, iTree, thisID, selectHalos[iTree].info())

            # If not orphan then copy the new one 
            else:
                try:
                    thisHalo = progHalos[idProgHalos[thisID]]
                    allHaloInfo[iTree].append(thisHalo.line.rstrip('\n'))

                except:
                    'This halo is also too small to be tracked'
                    allHaloInfo[iTree].append(" ")
                    #print("Problem at step %d, tree %d: %s" % (iTree, jSnap, thisHalo.info()))

        iTree += 1

    # Now reset all the indexes
    idProgHalos.clear()
    progHalos.clear()

orphInfo.close()
print("Trees read and reconstructed with halo catalogs, dumping everything to: ", outPath)

# Close all files once the reconstructed trees have been written
for iTree in range(0, nTrees):
    thisID = selectHalos[iTree].ID
    thisName = outPath + 'halo_' + str(thisID) + '_tree.dat'
    thisFile = open(thisName, 'w')
    print(iTree, thisName)

    for iLine in range(0, nSnaps):
        iSnap = endSnaps - iLine
        thisZSnap = z_s[iSnap]
        try:
            lineHaloInfo = allHaloInfo[iTree][iLine]
            print(thisZSnap, lineHaloInfo, file=thisFile)
        except:
            'This line has been cut earlier - probably it was a very small halo'

    thisFile.close()

print("Done.")
