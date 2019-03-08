from libio.read_ascii import *
from libcosmo.halo import *
import os
import pickle

'''
class StrToBytes:
    def __init__(self, fileobj):
        self.fileobj = fileobj
    def read(self, size):
        return self.fileobj.read(size).encode()
    def readline(self, size=-1):
        return self.fileobj.readline(size).encode()
'''

# How many catalogs do we want to read
iniSnaps = 50
endSnaps = 54
totSnaps = 55

# File properties
rootPath = '/home/eduardo/CLUES/DATA/LGF/1024/FIX/'
subPath = '00_00'
baseAHF = 'snapshot_'
suffAHF = '.AHF_halos'

z_out = rootPath + '../output_z.txt'
z_f = open(z_out, 'r')
z_s = z_f.read().split()

# MetroCPP mtree files are loaded with pickle

# AHF halos are (pairwise) stored here
idDescHalos = dict()
idProgHalos = dict()
idOrphHalos = dict()

# List of halos without direct progenitor
orphHalos = []

# Start reading from the final snapshot
thisNSnap = '%03d' % endSnaps
thisZSnap = z_s[endSnaps]
thisAHF = rootPath + subPath + '/' + baseAHF + thisNSnap + '.z' + thisZSnap + suffAHF

# Load the merger tree files from pickle
pklM31 = 'saved/all_trees_lgf_m31.pkl'
pklMW = 'saved/all_trees_lgf_mw.pkl'

fileMW = open(pklMW, 'r')
treeMW = pickle.load(fileMW)
fileM31 = open(pklM31, 'r')
treeM31 = pickle.load(fileM31)

descHalos = read_ahf(thisAHF)
idDescHalos = makeHaloDictionary(descHalos)

print(treeM31[0].info())

# Loop on all the snapshots
for jSnap in range(iniSnaps, endSnaps):
    iSnap = endSnaps - jSnap - 1
    thisNSnap = '%03d' % iSnap
    thisZSnap = z_s[iSnap]
    thisAHF = rootPath + subPath + '/' + baseAHF + thisNSnap + '.z' + thisZSnap + suffAHF

    print('Reading AHF halo catalog: ', thisAHF)
    progHalos = read_ahf(thisAHF)
    idProgHalos = makeHaloDictionary(progHalos)

    '''
    #thisAHF = 
    # Track orphan halos
    '''


    # Now reset all the indexes
    idDescHalos.clear()
    idDescHalos = idProgHalos.copy()
    idProgHalos.clear()
 
    #print('DescHalos: ', len(idDescHalos))
    #print('ProgHalos: ', len(idProgHalos))

