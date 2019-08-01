import re

class Halo():

    """
    Generic class representing a halo.
    """
    
    properties = dict()

    def __init__(self, halo_id):
        self.properties = dict()
        self._descriptor = "halo_" + str(halo_id)
        self.properties['halo_id'] = halo_id

    def assign_property(self, key, value):
        self.properties[key] = value


def open_(filename, *args):
    """Open a file, determining from the filename whether to use
    gzip decompression"""

    if (filename[-3:] == '.gz'):
        return gzip.open(filename, *args)
    try:
        return open(filename, *args)
    except IOError:
        return gzip.open(filename + ".gz", *args)

def load_ahf_halos(filename):
    f = open_(filename,"rt")
    # get all the property names from the first, commented line
    # remove (#)
    keys = [re.sub('\([0-9]*\)', '', field)
            for field in f.readline().split()]

    halos = []

    # provide translations
    for i, key in enumerate(keys):
        if(key == '#npart'):
            keys[i] = 'npart'
        if(key == '#ID'):
            keys[i] = 'ID'
        if(key == 'a'):
            keys[i] = 'a_axis'
        if(key == 'b'):
            keys[i] = 'b_axis'
        if(key == 'c'):
            keys[i] = 'c_axis'
        if(key == 'Mvir'):
            keys[i] = 'mass'

        # fix for column 0 being a non-column in some versions of the AHF
        # output
        if keys[0] == '#':
            keys = keys[1:]

        print(keys)
        for h, line in enumerate(f):
            values = [float(x) if '.' in x or 'e' in x or 'nan' in x else int(
                x) for x in line.split()]

            halo_id = values[0]
            halo = Halo(halo_id)

            for i, key in enumerate(keys):
                try:
                    halo.assign_property(key, values[i])
                except:
                    'Do nothing'

            halos.append(halo)

        f.close()

        return halos
