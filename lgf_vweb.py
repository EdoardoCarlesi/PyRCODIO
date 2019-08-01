from libio.read_ascii import *
import pickle
import numpy as np
import os

# Root path for the vweb files
#dirBase='/z/carlesi/CLUES/DATA/LGF/SNAPS/512/'
#dirBase='/work2/eduardo/DATA/512/VWEB/'
#fileWeb='vweb_lgf_054.000064.Vweb-ascii'
fileLGs='saved/lgs_'
fileWeb='saved/lgf_'

fileEV='saved/lgs_evs.pkl'

# The LGs used for the v-web will be attached here
fileSelectedLGs='saved/lgs_select.pkl'

#fileEV='saved/lgs_evs_all.pkl'
#fileSelectedLGs='saved/lgs_select_all.pkl'

# Do a loop
iSta=0
iEnd=10
gSta=0
gEnd=30

# Main web parameters
boxSize=100000.0
gridSize=64

norm = 512.0
ev1 = []; ev2 = []; ev3 = []
lgs = []


# Main loop
for iRun in range(iSta, iEnd):
	iRunStr = '%02d' % iRun

	# Sub loop
	for gRun in range(gSta, gEnd):
		gRunStr = '%02d' % gRun
		
		subRunStr = iRunStr + '_' + gRunStr
		#thisFileWeb = dirBase + subRunStr + '/' + fileWeb
		thisFileWeb = fileWeb + subRunStr + '.pkl'
		thisFileLGs = fileLGs + subRunStr + '.pkl'
		
		# Check if files exist
		exist1 = os.path.isfile(thisFileWeb)
		exist2 = os.path.isfile(thisFileLGs)
		
		# If they do then read the vweb and the lg
		if exist1 and exist2:
                    #thisWeb = read_vweb(thisFileWeb, gridSize, boxSize)
                    f_lg = open(thisFileLGs, 'rb')
                    f_web = open(thisFileWeb, 'rb')
                    thisLG = pickle.load(f_lg)
                        
                    for lg in thisLG:
                        thisCOM = lg.get_com()

                        coord = []
                        for ix in range(0, 3):
                            jx = int(gridSize * thisCOM[ix] / boxSize)
                            coord.append(jx)

                        thisEV = thisWeb.evals[:, coord[0], coord[1], coord[2]]

                        if thisEV[0] < 0.001:
                            norm = 512.0
                        else:
                            norm = 1.0

                        if (thisEV[0] * norm > 0.00  and thisEV[2] < 0.035):
                            print(thisEV[0] * norm, thisEV[1] * norm, thisEV[2] * norm)
                            ev1.append(thisEV[0] * norm)
                            ev2.append(thisEV[1] * norm)
                            ev3.append(thisEV[2] * norm)
                            lgs.append(thisLG)
                                                
evs = [ev1, ev2, ev3]
f_evs = open(fileEV, 'wb')
pickle.dump(evs, f_evs)

f_lgs = open(fileSelectedLGs, 'wb')
pickle.dump(lgs, f_lgs)




