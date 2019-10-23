from libio.read_ascii import *
import pickle
import numpy as np
import os

# Root path for the vweb files
dirBase='/z/carlesi/CLUES/DATA/512/'
#dirBase='/z/carlesi/CLUES/DATA/LGF/SNAPS/512/'
#dirBase='/work2/eduardo/DATA/512/VWEB/'
fileWeb1='vweb_lgf_054.000064.Vweb-ascii'
fileWeb2='vweb_127.000064.Vweb-ascii'
#fileLGs='saved/lgs_'
#'saved/lgs_r_5000.0_mMin4e+11_512_00_00.pkl'
fileLGs='saved/lgs_r_5000.0_mMin4e+11_512_'
#fileWeb='saved/lgf_'

fileEV='saved/lgs_evs_512_5000.pkl'

# The LGs used for the v-web will be attached here
fileSelectedLGs='saved/lgs_select_512_5000.pkl'

#fileEV='saved/lgs_evs_all.pkl'
#fileSelectedLGs='saved/lgs_select_all.pkl'

# Do a loop
iSta=0
iEnd=100
gSta=0
gEnd=30

# Main web parameters
boxSize=100000.0
gridSize=64

norm = 1.e+3
ev1 = []; ev2 = []; ev3 = []
lgs = []


# Main loop
for iRun in range(iSta, iEnd):
	iRunStr = '%02d' % iRun

	# Sub loop
	for gRun in range(gSta, gEnd):
		gRunStr = '%02d' % gRun
		
		subRunStr = iRunStr + '_' + gRunStr
		thisFileWeb1 = dirBase + subRunStr + '/' + fileWeb1
		thisFileWeb2 = dirBase + subRunStr + '/' + fileWeb2
		thisFileLGs = fileLGs + subRunStr + '.pkl'
		
#		print(thisFileWeb)
#		print(thisFileLGs)

		# Check if files exist
		exist1 = os.path.isfile(thisFileWeb1)
		exist2 = os.path.isfile(thisFileWeb2)
		exist3 = os.path.isfile(thisFileLGs)
		
		if exist1:
			thisFileWeb = thisFileWeb1
			exist0 = True

		if exist2:
			thisFileWeb = thisFileWeb2
			exist0 = True

		# If they do then read the vweb and the lg
		if exist0 and exist3:
		    print('Found vWeb file: ', thisFileWeb)
                    f_lg = open(thisFileLGs, 'rb')
                    #f_web = open(thisFileWeb, 'rb')
                    thisLG = pickle.load(f_lg)
                     
                    for lg in thisLG:
			if lg.code != 'EMPTY':
                    	        thisWeb = read_vweb(thisFileWeb, gridSize, boxSize)
 	                        thisCOM = lg.get_com()

        	                coord = []
                	        for ix in range(0, 3):
                        	    jx = int(gridSize * thisCOM[ix] / boxSize)
	                            coord.append(jx)

        	                thisEV = thisWeb.evals[:, coord[0], coord[1], coord[2]]

                	        if abs(thisEV[0]) < 0.001:
                        	    norm = 1.e+3
	                        else:
        	                    norm = 1.0

                	        if (thisEV[0] * norm > -10.00  and thisEV[2] < 10.035):
                        	    print(thisEV[0] * norm, thisEV[1] * norm, thisEV[2] * norm, norm)
	                            ev1.append(thisEV[0] * norm)
        	                    ev2.append(thisEV[1] * norm)
                	            ev3.append(thisEV[2] * norm)
                        	    lgs.append(thisLG)
                                                
evs = [ev1, ev2, ev3]
f_evs = open(fileEV, 'wb')
pickle.dump(evs, f_evs)

f_lgs = open(fileSelectedLGs, 'wb')
pickle.dump(lgs, f_lgs)




