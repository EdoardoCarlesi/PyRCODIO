from libio.read_ascii import *
import pickle
import numpy as np
import os

gridSize = 64

dirBase = '/z/carlesi/CLUES/DATA/512/'
#fileWeb = 'vweb_2048.0000'+str(gridSize)+'.Vweb-ascii'

#dirBase = '/z/carlesi/STORE/SmallScaleVariance/vweb/'
#fileWeb = 'vweb_2048.0000'+str(gridSize)+'.Vweb-ascii'

fileWeb1 = 'vweb_127.000064.Vweb-ascii'
fileWeb2 = 'vweb_lgf_054.000064.Vweb-ascii'

# Do a loop
iSta=0
iEnd=80
gSta=0
gEnd=30

# Main web parameters
boxSize = 100.0
n_half = int(gridSize / 2)

# Subgrid parameters: only store a subcube of sizeSave size
sizeSave=20.0

# The subcube is double in size
nSubGrid= 2 * int(gridSize * sizeSave / boxSize)

subGrid = np.zeros((3, nSubGrid, nSubGrid, nSubGrid))

# Main loop
for iRun in range(iSta, iEnd):
    iRunStr = '%02d' % iRun

    # Sub loop
    for gRun in range(gSta, gEnd):
        gRunStr = '%02d' % gRun

        subRunStr = iRunStr + '_' + gRunStr
        thisFileWeb1 = dirBase + subRunStr + '/' + fileWeb1
        thisFileWeb2 = dirBase + subRunStr + '/' + fileWeb2

        # Check if file exists
        exists1 = os.path.isfile(thisFileWeb1)
        exists2 = os.path.isfile(thisFileWeb2)

        if exists1:
            thisFileWeb = thisFileWeb1

        if exists2:
            thisFileWeb = thisFileWeb2

        # If it does
        if exists1 or exists2:
            print(thisFileWeb, ' found.')
            thisGrid = read_vweb(thisFileWeb, gridSize, boxSize)
            print(thisGrid.evals[:, n_half, n_half, n_half])

            for dx in range(0, nSubGrid):
                ix = gridSize / 2 - nSubGrid / 2 + dx

                for dy in range(0, nSubGrid):
                    iy = gridSize / 2 - nSubGrid / 2 + dy

                    for dz in range(0, nSubGrid):
                        iz = gridSize / 2 - nSubGrid / 2 + dz
                        subGrid[:, dx, dy, dz] = thisGrid.evals[:, ix, iy, iz]

            this_vweb_out = 'saved/lgf_web_' + subRunStr + '.pkl'
            print('Saving to vweb: ', this_vweb_out)
            f_vweb_out = open(this_vweb_out, 'wb')
            pickle.dump(subGrid, f_vweb_out)
