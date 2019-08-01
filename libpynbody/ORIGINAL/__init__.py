"""

halo
====

Implements halo catalogue functions. If you have a supported halo
catalogue on disk or a halo finder installed and correctly configured,
you can access a halo catalogue through f.halos() where f is a
SimSnap.

See the `halo tutorial
<http://pynbody.github.io/pynbody/tutorials/halos.html>`_ for some
examples.

"""

import numpy as np
import weakref
import copy
import logging

from .. import snapshot, util

logger = logging.getLogger("pynbody.halo")

class DummyHalo(object):

    def __init__(self):
        self.properties = {}


class Halo(snapshot.IndexedSubSnap):

    """
    Generic class representing a halo.
    """

    def __init__(self, halo_id, halo_catalogue, *args):
        super(Halo, self).__init__(*args)
        self._halo_catalogue = halo_catalogue
        self._halo_id = halo_id
        self._descriptor = "halo_" + str(halo_id)
        self.properties = copy.copy(self.properties)
        self.properties['halo_id'] = halo_id
        if halo_id in halo_catalogue._halos:
            for key in list(halo_catalogue._halos[halo_id].properties.keys()):
                self.properties[key] = halo_catalogue._halos[halo_id].properties[key]


    def is_subhalo(self, otherhalo):
        """
        Convenience function that calls the corresponding function in
        a halo catalogue.
        """

        return self._halo_catalogue.is_subhalo(self._halo_id, otherhalo._halo_id)


# ----------------------------#
# General HaloCatalogue class #
#-----------------------------#

class HaloCatalogue(object):

    """
    Generic halo catalogue object.
    """

    def __init__(self, sim):
        self._base = weakref.ref(sim)
        self._halos = {}
        self.lazy_off = util.ExecutionControl()

    def calc_item(self, i):
        if i in self._halos:  # and self._halos[i]() is not None :
            if isinstance(self._halos[i],DummyHalo):
                try:
                    return self._get_halo(i)
                except:
                    return self._halos[i]
            else:
                return self._halos[i]
        else:
            h = self._get_halo(i)
            self._halos[i] = h  # weakref.ref(h)
            return h

    def __len__(self):
        return len(self._halos)

    def __iter__(self):
        return self._halo_generator()

    def __getitem__(self, item):
        if isinstance(item, slice):
            for x in self._halo_generator(item.start,item.stop) : pass
            indices = item.indices(len(self._halos))
            res = [self.calc_item(i) for i in range(*indices)]
            return res
        else:
            return self.calc_item(item)

    @property
    def base(self):
        return self._base()

    def _halo_generator(self, i_start=None, i_stop=None) :
        if len(self) == 0 : return
        if i_start is None :
            try :
                self[0]
                i = 0
            except KeyError :
                i = 1
        else :
            i = i_start

        if i_stop is None :
            i_stop = len(self)

        while True:
            try:
                yield self[i]
                i+=1
                if i!=i_stop and len(self[i]) == 0: continue
            except RuntimeError:
                break
            if i == i_stop: return

    def _init_iord_to_fpos(self):
        if not hasattr(self, "_iord_to_fpos"):
            self._iord_to_fpos = np.empty(self.base['iord'].max()+1,dtype=np.int64)
            self._iord_to_fpos[self.base['iord']] = np.arange(len(self.base))

    def is_subhalo(self, childid, parentid):
        """Checks whether the specified 'childid' halo is a subhalo
        of 'parentid' halo.
        """
        if (childid in self._halos[parentid].properties['children']):
            return True
        else:
            return False

    def contains(self, haloid):
        if (haloid in self._halos):
            return True
        else:
            return False

    def __contains__(self, haloid):
        return self.contains(haloid)

    def get_group_array(self):
        """Return an array with an integer for each particle in the simulation
        indicating which halo that particle is associated with. If there are multiple
        levels (i.e. subhalos), the number returned corresponds to the lowest level, i.e.
        the smallest subhalo."""
        raise NotImplementedError

    @staticmethod
    def _can_load(self):
        return False

    @staticmethod
    def _can_run(self):
        return False


class GrpCatalogue(HaloCatalogue):
    """
    A generic catalogue using a .grp file to specify which particles
    belong to which group.
    """
    def __init__(self, sim, array='grp', ignore=None, **kwargs):
        """Construct a GrpCatalogue, extracting halos based on a simulation-wide integer array with their IDs

        *sim* - the SimSnap for which the halos will be constructed
        *array* - the name of the array which should be present, loadable or derivable across the simulation
        *ignore* - a special value indicating "no halo", or None if no such special value is defined
        """
        sim[array] # trigger lazy-loading and/or kick up a fuss if unavailable
        self._halos = {}
        self._array = array
        self._sorted = None
        self._ignore = ignore
        HaloCatalogue.__init__(self,sim)

    def __len__(self):
        if self._ignore is None:
            N = self.base[self._array].max()
        else:
            N = self.base[self._array]
            N = N[N!=self._ignore]
            N = N.max()
        if N<0:
            N=0
        return N

    def precalculate(self):
        """Speed up future operations by precalculating the indices
        for all halos in one operation. This is slow compared to
        getting a single halo, however."""
        self._sorted = np.argsort(
            self.base[self._array], kind='mergesort')  # mergesort for stability
        self._boundaries = util.find_boundaries(
            self.base[self._array][self._sorted])

    def get_group_array(self, family=None):
        if family is not None:
            return self.base[family][self._array]
        else:
            return self.base[self._array]

    def _get_halo_indices(self, i):
        if self.base is None:
            raise RuntimeError("Parent SimSnap has been deleted")

        no_exist = ValueError("Halo %s does not exist" % (str(i)))

        if self._sorted is None:
            # one-off selection
            index = np.where(self.base[self._array] == i)
            return index
        else:
            # pre-calculated
            if i >= len(self._boundaries) or i < 0:
                raise no_exist
            if self._boundaries[i] < 0:
                raise no_exist

            start = self._boundaries[i]
            if start is None:
                raise no_exist

            end = None
            j = i + 1
            while j < len(self._boundaries) and end is None:
                end = self._boundaries[j]
                j += 1

            return self._sorted[start:end]


    def _get_halo(self, i):
        x = Halo(i, self, self.base, self._get_halo_indices(i))
        if len(x) == 0:
            raise ValueError("Halo %s does not exist" % (str(i)))
        x._descriptor = "halo_" + str(i)
        return x


    def load_copy(self, i):
        """Load the a fresh SimSnap with only the particle in halo i"""
        from .. import load
        return load(self.base.filename, take=self._get_halo_indices(i))

    @property
    def base(self):
        return self._base()

    @staticmethod
    def _can_load(sim, arr_name='grp'):
        if (arr_name in sim.loadable_keys()) or (arr_name in list(sim.keys())) :
            return True
        else:
            return False


class AmigaGrpCatalogue(GrpCatalogue):
    def __init__(self, sim, arr_name='amiga.grp',**kwargs):
        GrpCatalogue.__init__(self, sim, arr_name)

    @staticmethod
    def _can_load(sim,arr_name='amiga.grp'):
        return GrpCatalogue._can_load(sim, arr_name)


from pynbody.halo.ahf import AHFCatalogue
from pynbody.halo.hop import HOPCatalogue
from pynbody.halo.legacy import RockstarIntermediateCatalogue
from pynbody.halo.rockstar import RockstarCatalogue
from pynbody.halo.subfind import SubfindCatalogue
from pynbody.halo.subfindhdf import SubFindHDFSubhaloCatalogue, SubFindHDFHaloCatalogue

def _get_halo_classes():
    # AmigaGrpCatalogue MUST be scanned first, because if it exists we probably
    # want to use it, but an AHFCatalogue will probably be on-disk too.
    _halo_classes = [GrpCatalogue, AmigaGrpCatalogue, AHFCatalogue,
                     RockstarCatalogue, SubfindCatalogue, SubFindHDFHaloCatalogue,
                     RockstarIntermediateCatalogue, HOPCatalogue]

    return _halo_classes
