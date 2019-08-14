import numpy
import sys
import os
#sys.path.append('../libcosmo/')
#sys.path.append('..')

from libcosmo.grid import *
#from libcosmo.porcodio import *
from libcosmo.halos import Halo


def read_grid(file_name, size, box):
    file_web = open(file_name, 'r')
    grid = Grid(size, box)

    # Read header
    line = file_web.readline()
    tot_n = size * size * size
    index = 0

    while line and index < tot_n:
        line = file_web.readline()
        line = line.strip()
        column = line.split()#; print(column[0])

        # Determine corresponding x, y, z in grid units
        (ix, jy, kz) = grid.reverse_index(index)#; print(ix, jy, kz)

        # Density
        grid.rho[ix, jy, kz] = float(column[0])

        # Velocities
        grid.vel[:, ix, jy, kz] = [float(column[1]), float(column[2]), float(column[3])]

        # Increase line index
        index += 1

    return grid



def read_vweb(file_name, size, box):
    file_web = open(file_name, 'r')
    grid = VWeb(size, box)

    # Read header
    line = file_web.readline()
    tot_n = size * size * size
    index = 0

    while line and index < tot_n:
        line = file_web.readline()
        line = line.strip()
        column = line.split()#; print(column[0])

        # Determine corresponding x, y, z in grid units
        (ix, jy, kz) = grid.reverse_index(index)#; print(ix, jy, kz)

        # Density
        grid.rho[ix, jy, kz] = float(column[0])

        # Velocities
        grid.vel[:, ix, jy, kz] = [float(column[1]), float(column[2]), float(column[3])]

        # Eigenvalues
        grid.evals[0, ix, jy, kz] = float(column[4])
        grid.evals[1, ix, jy, kz] = float(column[5])
        grid.evals[2, ix, jy, kz] = float(column[6])

        # Eigenvectors
        grid.evecs[0, :, ix, jy, kz] = [float(column[7]), float(column[8]), float(column[9])]
        grid.evecs[1, :, ix, jy, kz] = [float(column[10]), float(column[11]), float(column[12])]
        grid.evecs[2, :, ix, jy, kz] = [float(column[13]), float(column[14]), float(column[15])]

        # Increase line index
        index += 1

    return grid


# This function reads haloes within a given mass range from a AHF catalog split into several chunks
def read_ahf_chunks(file_root, file_suffix, n_chunks):
    halos_ahf = []
    count = 0

    for i_chunk in range(0, n_chunks):
        this_chunk = '%04d' % i_chunk
        this_name = file_root + this_chunk + file_suffix
        file_ahf = open(this_name, 'r')

        line = file_ahf.readline()

        while line:
            line = file_ahf.readline()
            line = line.strip()
            column = line.split()
            n_col = len(column)

            if n_col > 1:
                # Read halo properties
                idn = long(column[0])
                mass = float(column[3])
                pos = [float(column[5]), float(column[6]), float(column[7])]
                vel = [float(column[8]), float(column[9]), float(column[10])]
                rvir = float(column[11])
                nsub = int(column[2])
                npart = int(column[4])
                angmom = [float(column[21]), float(column[22]), float(column[23])]
                contam = float(column[38])

                # Initialize and append halo to the list
                halo = Halo()
                halo.initialize(idn, mass, pos, vel, rvir, nsub, npart)
                halo.update_id_index(idn, count)
                #halo.ID = idn
                halo.l = angmom
                halo.contam = contam
                halos_ahf.append(halo)
                count += 1

    n_lines = count
    print("Found a total of %d halos " % (n_lines))

    return halos_ahf



# This function reads haloes within a given mass range from a AHF catalog split into several chunks
def read_ahf_chunks_mass_range(file_root, file_suffix, n_chunks, m_max, m_min):
    halos_ahf = []
    count = 0

    for i_chunk in range(0, n_chunks):
        this_chunk = '%04d' % i_chunk
        this_name = file_root + this_chunk + file_suffix
        file_ahf = open(this_name, 'r')

        line = file_ahf.readline()

        while line:
            line = file_ahf.readline()
            line = line.strip()
            column = line.split()
            n_col = len(column)

            if n_col > 1:
                # Read halo properties
                idn = long(column[0])
                mass = float(column[3])
                pos = [float(column[5]), float(column[6]), float(column[7])]
                vel = [float(column[8]), float(column[9]), float(column[10])]
                rvir = float(column[11])
                nsub = int(column[2])
                npart = int(column[4])
                angmom = [float(column[21]), float(column[22]), float(column[23])]
                contam = float(column[38])

                if mass > m_min and mass < m_max:
                    # Initialize and append halo to the list
                    halo = Halo()
                    halo.initialize(idn, mass, pos, vel, rvir, nsub, npart)
                    halo.update_id_index(idn, count)
                    #halo.ID = idn
                    halo.l = angmom
                    halo.contam = contam
                    halos_ahf.append(halo)
                    count += 1

    n_lines = count
    print("Found a total of %d halos in mass range (%e, %e)" % (n_lines, m_min, m_max))

    return halos_ahf



def read_dwarf_list(file_name):
    file_dw = open(file_name, 'r')
    all_lines = file_dw.readlines() #.rstrip('\n')
    dwarfIDs = []

    iLine = 0
    for line in all_lines:
        if iLine > 0:
            thisLine = line.strip().rstrip('\n').split()
            thisID = thisLine[0]
            dwarfIDs.append(int(thisID))

        iLine += 1

    return dwarfIDs



def read_ahf_mass_range(file_name, m_max, m_min):
        # Open file
    file_ahf = open(file_name, 'r')

    line = file_ahf.readline()
    halos_ahf = []
    count = 0

    while line:
        line = file_ahf.readline()
        line = line.strip()
        column = line.split()
        n_col = len(column)

        if n_col > 1:
            # Read halo properties
            idn = long(column[0])
            mass = float(column[3])
            pos = [float(column[5]), float(column[6]), float(column[7])]
            vel = [float(column[8]), float(column[9]), float(column[10])]
            rvir = float(column[11])
            nsub = int(column[2])
            npart = int(column[4])
            angmom = [float(column[21]), float(column[22]), float(column[23])]
            contam = float(column[38])

            #pos[0] *= 1000. ; pos[1] *= 1000. ; pos[2] *= 1000.

            if mass > m_min and mass < m_max:
                # Initialize and append halo to the list
                halo = Halo()
                halo.initialize(idn, mass, pos, vel, rvir, nsub, npart)
                halo.update_id_index(idn, count)
                #halo.ID = idn
                halo.l = angmom
                halo.contam = contam
                halos_ahf.append(halo)
                count += 1

    n_lines = count
    print("Found a total of %d halos in mass range (%e, %e)" % (n_lines, m_min, m_max))

    return halos_ahf

#ID(1)  hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)        lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)        SurfP(41)       Phi0(42)        cNFW(43)        n_gas(44)       M_gas(45)       lambda_gas(46)  lambdaE_gas(47) Lx_gas(48)      Ly_gas(49)      Lz_gas(50)      b_gas(51)       c_gas(52)       Eax_gas(53)     Eay_gas(54)     Eaz_gas(55)     Ebx_gas(56)     Eby_gas(57)     Ebz_gas(58)     Ecx_gas(59)     Ecy_gas(60)     Ecz_gas(61)     Ekin_gas(62)    Epot_gas(63)    n_star(64)      M_star(65)      lambda_star(66) lambdaE_star(67)        Lx_star(68)     Ly_star(69)     Lz_star(70)     b_star(71)      c_star(72)      Eax_star(73)    Eay_star(74)    Eaz_star(75)    Ebx_star(76)    Eby_star(77)    Ebz_star(78)    Ecx_star(79)    Ecy_star(80)    Ecz_star(81)    Ekin_star(82)   Epot_star(83)   mean_z_gas(84)  mean_z_star(85)

def read_ahf(file_name):
    file_ahf = open(file_name, 'r')

    line = file_ahf.readline()
    halos_ahf = []
    count = 0

    while line:
        full_line = file_ahf.readline()
        line = full_line.strip()
        column = line.split()
        n_col = len(column)

        if n_col > 1 and line[0] != "#":
            # Read halo properties
            idn = int(column[0])
            mass = float(column[3])
            pos = [float(column[5]), float(column[6]), float(column[7])]
            vel = [float(column[8]), float(column[9]), float(column[10])]
            rvir = float(column[11])
            nsub = int(column[2])
            npart = int(column[4])
            vmax = float(column[16])
            angmom = [float(column[21]), float(column[22]), float(column[23])]
            contam = float(column[38])
            ngas = int(column[43])
            mgas = float(column[44])
            nstar = int(column[63])
            mstar = float(column[64])

            #pos[0] *= 1000. ; pos[1] *= 1000. ; pos[2] *= 1000

            # Initialize and append hl to the list
            hl = Halo(); hl.line = full_line
            hl.initialize(idn, mass, pos, vel, rvir, nsub, npart)
            hl.update_id_index(idn, count)
            #hl.ID = idn
            hl.l = angmom
            hl.contam = contam
            hl.vmax = vmax
            hl.m_star = mstar
            hl.m_gas = mgas
            hl.m_dm = mass - mstar - mgas
            hl.nstar = nstar
            hl.ngas = ngas
            halos_ahf.append(hl)
            count += 1

    n_lines = count
    #print "Read %s with a total of %d lines" % (file_name, n_lines)

    return halos_ahf



# Reading AHF particle file
def read_particles_chunks(file_root, file_suff, n_files, n_halos):
    ids = dict()    # List containing all halo IDs & Particle number per halo ID - each list member is a 2 elements array
    parts = []      # Each list member contains an array with all the particle ids
    count_p = 0     # Total number of particles per halo
    count_h = 0     # Total number of haloes in file

    this_np = 0
    tot_h = 0       # Double check that the total number of haloes is matched with the counter

    for i_chunk in range(0, n_files):
        this_chunk = '%04d' % i_chunk
        file_name = file_root + this_chunk + file_suff
        file_part = open(file_name, 'r')
        count_l = 0     # Reset the lines to zero for each file

        # First read the header, containing the total number of haloes
        line = file_part.readline()
        line = line.strip()
        column = line.split()

        lines_command = 'wc -l ' + file_name
        out_os = os.popen(lines_command).read()
        (tot_l, fname) = out_os.split()
        tot_l = long(tot_l)

        print('Reading particles %s with %ld lines and %d halos. ' % (file_name, tot_l, n_halos))

        while line:
            line = file_part.readline()
            line = line.strip()
            column = line.split()

            if count_l > 0 and count_p < this_np:
                this_pid = long(column[0])      # Particle ID
                this_parts.append(this_pid)
                count_p += 1
                count_l += 1
            else:
                # All particles have been read in
                if count_p == this_np and count_h > 0:
                    this_parts.sort()       # Automatically sort them by ascending order!
                    parts.append(this_parts)

                # Still reading particle files
                if count_l < tot_l-1:
                    this_parts = []
                    this_hid = str(column[1])       # Halo ID
                    this_np = int(column[0])
                    this_index = count_h
                    #print 'Line %ld found halo %ld with %d particles' % (count_l, this_hid, this_np)
                    ids.update({this_hid:[this_index, this_np]})

                    count_l += 1
                    count_h += 1
                    count_p = 0                     # Reset the particle number

    print('Expected %d halos, found %d ' % (count_h, n_halos))

    return (ids, parts)




# Reading AHF particle file
def read_particles(file_name):
    ids = dict()    # List containing all halo IDs & Particle number per halo ID - each list member is a 2 elements array
    parts = []      # Each list member contains an array with all the particle ids
    count_p = 0     # Total number of particles per halo
    count_h = 0     # Total number of haloes in file
    count_l = 0     # Simply count the lines

    this_np = 0
    tot_h = 0       # Double check that the total number of haloes is matched with the counter

    file_part = open(file_name, 'r')

    # First read the header, containing the total number of haloes
    line = file_part.readline()
    line = line.strip()
    column = line.split()
    tot_h = int(column[0])

    lines_command = 'wc -l ' + file_name
    out_os = os.popen(lines_command).read()
    (tot_l, fname) = out_os.split()
    tot_l = long(tot_l)
    print('Reading particles %s with %ld lines and %d halos. ' % (file_name, tot_l, tot_h))

    while line:
        line = file_part.readline()
        line = line.strip()
        column = line.split()

        if count_l > 0 and count_p < this_np:
            this_pid = long(column[0])      # Particle ID
            this_parts.append(this_pid)
            count_p += 1
            count_l += 1
        else:
            # All particles have been read in
            if count_p == this_np and count_h > 0:
                this_parts.sort()       # Automatically sort them by ascending order!
                parts.append(this_parts)

            # Still reading particle files
            if count_l < tot_l-1:
                this_parts = []
                this_hid = str(column[1])       # Halo ID
                this_np = int(column[0])
                this_index = count_h
                #print 'Line %ld found halo %ld with %d particles' % (count_l, this_hid, this_np)
                ids.update({this_hid:[this_index, this_np]})

                count_l += 1
                count_h += 1
                count_p = 0                     # Reset the particle number

    return (ids, parts)


# Read the LG asciis saved somewhere
def read_lgs(file_name):
    file_txt = open(file_name, 'r')

    line = file_txt.readline()
    lgs_txt = []
    count = 0

    lgs = []

    com = [0.0] * 3
    pos = [0.0] * 3
    vel = [0.0] * 3


    while line:
        line = file_txt.readline()
        line = line.strip()
        column = line.split()
        n_col = len(column)

        if n_col > 1:
            lg1 = Halo()
            lg2 = Halo()

            # Read halo properties
            id0 = long(column[0])
            id1 = long(column[1])
            m0 = float(column[2])
            m1 = float(column[3])
            dist = float(column[4])
            vrad = float(column[5])
            nsub0 = int(column[6])
            nsub1 = int(column[7])
            code = column[8]
            com[0] = float(column[9]) * 1000.
            com[1] = float(column[10]) * 1000.
            com[2] = float(column[11]) * 1000.

            # Most of this is initialized to dummy variables
            lg1.initialize(id0, m0, pos, vel, 0.0, nsub0, 0)
            lg2.initialize(id1, m1, pos, vel, 0.0, nsub1, 0)

            lg = LocalGroup(code)
            lg.init_halos(lg1, lg2)
            lg.vrad = vrad
            lg.r = dist
            lg.com = com

            lgs.append(lg)

            del lg1
            del lg2

    return lgs
