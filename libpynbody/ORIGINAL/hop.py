import os.path
import re
import struct

import numpy as np

from . import GrpCatalogue


class HOPCatalogue(GrpCatalogue):
    """A HOP Catalogue as used by Ramses. HOP output files must be in simulation directory, in simulation/hop/ directory
    or specified by fname"""
    def __init__(self, sim, fname=None):
        self._halos = {}

        if fname is None:
            for fname in HOPCatalogue._enumerate_hop_tag_locations_from_sim(sim):
                if os.path.exists(fname):
                    break

            if not os.path.exists(fname):
                raise RuntimeError("Unable to find HOP .tag file in simulation directory")


        sim._create_array('hop_grp', dtype=np.int32)
        sim['hop_grp']=-1
        with open(fname, "rb") as f:
            num_part, = struct.unpack('i', f.read(4))
            if num_part==8:
                # fortran-formatted output
                num_part, num_grps, _, _ = struct.unpack('iiii', f.read(16))
            else:
                # plain binary output
                num_grps, = struct.unpack('i', f.read(4))

            if num_part!=len(sim.dm):
                raise RuntimeError("Mismatching number of particles between snapshot %s and HOP file %s"%(sim.filename, fname))
            sim.dm['hop_grp'] = np.fromfile(f, np.int32, len(sim.dm))
        GrpCatalogue.__init__(self,sim,array="hop_grp")

    @staticmethod
    def _can_load(sim, arr_name='grp'):
        # Hop output must be in output directory or in output_*/hop directory
        exists = any([os.path.exists(fname) for fname in HOPCatalogue._enumerate_hop_tag_locations_from_sim(sim)])
        return exists

    @staticmethod
    def _extract_hop_name_from_sim(sim):
        match = re.search("output_([0-9]*)", sim.filename)
        if match is None:
            raise IOError("Cannot guess the HOP catalogue filename for %s" % sim.filename)
        return "grp%s.tag" % match.group(1)

    @staticmethod
    def _enumerate_hop_tag_locations_from_sim(sim):
        try:
            name = HOPCatalogue._extract_hop_name_from_sim(sim)
        except IOError:
            return []

        s_filename = os.path.abspath(sim.filename)

        return [os.path.join(os.path.dirname(s_filename),name),
                os.path.join(s_filename,name),
                os.path.join(s_filename,'hop',name)]

    def _can_run(self):
        return False