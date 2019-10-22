from libio.read_ascii import *
from libcosmo.halos import *
import os

res='4096'
halo_base='/home/eduardo/CLUES/DATA/' + res + '/'
tree_base='/home/eduardo/CLUES/DATA/trees/' + res + '/'

# 'snapshot_055.0000.z0.000.AHF_halos'
ahf_base = 'snapshot_'
ahf_midd = '.0000.'
ahf_suff = '.AHF_halos'

#'halo_4037937173177511913.ids'
mcpp_base = 'halo_'
mcpp_suff = '.ids'

out_suff = '.allinfo'

run='37_11'
ahf_zs = halo_base + run + '/zs.txt'

# AHF snapshots
sIni = 0
sEnd = 55

gIni = 1
gEnd = 10

f_zs = open(ahf_zs, 'r')
lines = f_zs.readlines()

zs = []
for line in lines:
    zs.append(line.rstrip('\n'))

for g in range(gIni, gEnd):
    gDir = '%02d' % g
    
    halosIDs = dict()

    for s in range(0, sEnd - sIni):
        sInv = sEnd - s
        sSnap = '%03d' % sInv
        zStr = zs[sInv-1]
        file_snap = halo_base + run + '/' + gDir + '/' + ahf_base + sSnap + ahf_midd + zStr + ahf_suff

        if os.path.isfile(file_snap):
            print(s, '. Reading: ', file_snap)

            if s == 0:
                halos0 = read_ahf(file_snap)
            else:
                halos = read_ahf(file_snap)
            
                for halo in halos:
                    halosIDs[str(halo.ID)] = halo
        else:
            print('File: ', file_snap, ' not found.')

    print('Dumping the full mtrees...')

    # Now read the full mtrees
    for halo0 in halos0:
        idStr = str(halo0.ID)
        file_tree = tree_base + run + '/' + gDir + '/' + mcpp_base + idStr + mcpp_suff
        file_out = tree_base + run + '/' + gDir + '/' + mcpp_base + idStr + out_suff
#        print(file_tree)

        if os.path.isfile(file_tree):
            #print('Found: ', file_tree, len(halos0))
            f_tree = open(file_tree, 'r')
            f_out = open(file_out, 'w')

            print(halo0.header_ahf(), file=f_out)
            theseIDs = []; iLine = 0

            for thisID in f_tree.readlines():
                theseIDs.append(thisID.rstrip('\n'))

                if iLine == 0:
                   print(halo0.line.rstrip('\n'), file=f_out)

                if iLine > 0:
                    try:
                        thisHalo = halosIDs[str(thisID.rstrip('\n'))]
                        print(thisHalo.line.rstrip('\n'), file=f_out)

                    except:
                        'This ID has not been found'

                iLine += iLine + 1

            #print('Full tree written to: ', file_out)
            f_out.close()
            f_tree.close()
            
print('Done.')
