import os
import pickle
import numpy as np

resolution = '2048'
base_lgs = 'saved/lgs_' + resolution + '_'
base_sub = 'saved/sub_' + resolution + '_'

iini = 0
iend = 80
gini = 0
gend = 30
sini = 0
send = 10

m_lmc_min = 3e+10
d_lmc_max = 100

m_m33_min = 1.0e+11
d_m33_min = 200
d_m33_max = 450

i_m33 = 0
i_lmc = 0

for i in range(iini, iend):
    i_str = '%02d' % i

    for g in range(gini, gend):
        g_str = '%02d' % g

        for s in range(sini, send):
            s_str = '%02d' % s

            igs_str = i_str + '_' + g_str + '_' + s_str
        
            fn_sub = base_sub + igs_str + '.pkl'
            fn_lgs = base_lgs + igs_str + '.pkl'

            if os.path.exists(fn_sub) and os.path.exists(fn_lgs):
                #print('Found: ', fn_sub, fn_lgs)

                f_lgs = open(fn_lgs, 'rb')
                lgs = pickle.load(f_lgs)
                f_lgs.close()

                f_sub = open(fn_sub, 'rb')
                lg_sat, mw_sat, m31_sat = pickle.load(f_sub)
                f_sub.close()

                lg1 = lgs.LG1
                lg2 = lgs.LG2
                com = lgs.geo_com()

                # MW is the least massive of the pair
                if lg1.m > lg2.m:
                    lg1 = lgs.LG2
                    lg2 = lgs.LG1

                for sub in lg_sat:
                    this_d_lg1  = lg1.distance(sub.x) 
                    this_d_lg2  = lg2.distance(sub.x) 

                    # LMC
                    if this_d_lg1 < d_lmc_max and sub.m > m_lmc_min and sub.m != lg1.m and sub.m != lg2.m:
                        #print('MW : ', (lg1.x[0]-com[0])/1000., (lg1.x[1]-com[1])/1000., (lg1.x[2]-com[2])/1000.)
                        print('LMC: ', (sub.x[0]-com[0])/1000., (sub.x[1]-com[1])/1000., (sub.x[2]-com[2])/1000., this_d_lg1/1000., igs_str, sub.m/1e+11)
                        #print('LMC: ', this_d_lg1, igs_str, sub.info())
                        i_lmc = i_lmc + 1

                    # M33
                    if this_d_lg2 < d_m33_max and sub.m > m_m33_min and this_d_lg2 > d_m33_min and sub.m != lg1.m and sub.m != lg2.m:
                        print('-> M31: ', (lg2.x[0]-com[0])/1000., (lg2.x[1]-com[1])/1000., (lg2.x[2]-com[2])/1000.)
                        #print('-> M33: ', this_d_lg2, igs_str, sub.m/1e+11, sub.x[0] - com[0], sub.x[1] - com[1], sub.x[2] - com[2])
                        print('-> M33: ', this_d_lg2/1000., igs_str, sub.m/1e+11) #, sub.x[0] - com[0], sub.x[1] - com[1], sub.x[2] - com[2])
                        #print('M33: ', this_d_lg2, igs_str, sub.info())
                        i_m33 = i_m33 + 1

print('Found: ', i_lmc, ' LMCs and ', i_m33, ' M33s.')



