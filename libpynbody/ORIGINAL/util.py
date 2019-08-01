"""

util
====

Various utility routines used internally by pynbody.

"""

import gzip
import struct
import os
import threading
import sys
import time
import functools
import logging
import math
import sys

import numpy as np
import scipy

from .backcompat import fractions
from . import config
from . import units
from .array import SimArray

logger = logging.getLogger('pynbody.util')
from ._util import *


def open_(filename, *args):
    """Open a file, determining from the filename whether to use
    gzip decompression"""

    if (filename[-3:] == '.gz'):
        return gzip.open(filename, *args)
    try:
        return open(filename, *args)
    except IOError:
        return gzip.open(filename + ".gz", *args)


def open_with_size(filename, *args):
    """Open a file for reading, returning also the (decompressed)
    file size"""

    f = open_(filename, *args)
    if isinstance(f, gzip.GzipFile):
        fo = open(f.name, 'rb')
        fo.seek(-4, 2)
        r = fo.read()
        fo.close()
        return f, struct.unpack('<I', r)[0]
    else:
        f.seek(0, os.SEEK_END)
        buflen = f.tell()
        f.seek(4, os.SEEK_SET)
        return f, buflen


def eps_as_simarray(f, eps):
    """Convert th given eps to a SimArray with units of f['pos'] and dtype of f['mass']""" 
    if isinstance(eps, str):
        eps = units.Unit(eps)
    if not isinstance(eps, units.UnitBase):
        eps = eps * f['pos'].units
        logger.info("Considering eps = {}".format(eps))
    eps_value = eps._scale
    eps_unit = eps/eps_value
    eps = SimArray(np.ones(len(f), dtype=f['mass'].dtype) * eps_value, eps_unit)
    return eps


def get_eps(f):
    """The gravitational softening length is determined from (in order of
    preference):
    1. the array f['eps']
    2. f.properties['eps'] (scalar or unit)

    Return a SimArray with correct units and dtype (same dtype as 'mass' array)"""

    try:
        eps = f['eps']
    except KeyError:
        if 'eps' in f.properties:
            eps = eps_as_simarray(f, f.properties['eps'])
        else:
            raise RuntimeError("Cannot retrieve 'eps' from SimSnap")
    return eps


def gcf(a, b):
    while b > 0:
        a, b = b, a % b
    return a


def lcm(a, b):
    return (a * b) // gcf(a, b)


def intersect_slices(s1, s2, array_length=None):
    """Given two python slices s1 and s2, return a new slice which
    will extract the data of an array d which is in both d[s1] and
    d[s2].

    Note that it may not be possible to do this without information on
    the length of the array referred to, hence all slices with
    end-relative indexes are first converted into begin-relative
    indexes. This means that the slice returned may be specific to
    the length specified."""

    assert array_length is not None or \
        (s1.start >= 0 and s2.start >= 0 and s1.stop >= 0 and s2.start >= 0)

    s1_start = s1.start
    s2_start = s2.start
    s1_stop = s1.stop
    s2_stop = s2.stop
    s1_step = s1.step
    s2_step = s2.step

    if s1_step == None:
        s1_step = 1
    if s2_step == None:
        s2_step = 1

    assert s1_step > 0 and s2_step > 0

    if s1_start < 0:
        s1_start = array_length + s1_start
    if s1_start < 0:
        return slice(0, 0)

    if s2_start < 0:
        s2_start = array_length + s2_start
    if s2_start < 0:
        return slice(0, 0)

    if s1_stop < 0:
        s1_stop = array_length + s1_stop
    if s1_stop < 0:
        return slice(0, 0)

    if s2_stop < 0:
        s2_stop = array_length + s2_stop
    if s2_stop < 0:
        return slice(0, 0)

    step = lcm(s1_step, s2_step)

    start = max(s1_start, s2_start)
    stop = min(s1_stop, s2_stop)

    if stop <= start:
        return slice(0, 0)

    s1_offset = start - s1_start
    s2_offset = start - s2_start
    s1_offset_x = int(s1_offset)
    s2_offset_x = int(s2_offset)

    if s1_step == s2_step and s1_offset % s1_step != s2_offset % s1_step:
        # slices are mutually exclusive
        return slice(0, 0)

    # There is surely a more efficient way to do the following, but
    # it eludes me for the moment
    while s1_offset % s1_step != 0 or s2_offset % s2_step != 0:
        start += 1
        s1_offset += 1
        s2_offset += 1
        if s1_offset % s1_step == s1_offset_x % s1_step and s2_offset % s2_step == s2_offset_x % s2_step:
            # slices are mutually exclusive
            return slice(0, 0)

    if step == 1:
        step = None

    return slice(start, stop, step)


def relative_slice(s_relative_to, s):
    """Given a slice s, return a slice s_prime with the property that
    array[s_relative_to][s_prime] == array[s]. Clearly this will
    not be possible for arbitrarily chosen s_relative_to and s, but
    it should be possible for s=intersect_slices(s_relative_to, s_any)
    which is the use case envisioned here (and used by SubSim).
    This code currently does not work with end-relative (i.e. negative)
    start or stop positions."""

    assert (s_relative_to.start >= 0 and s.start >= 0 and s.stop >= 0)

    if s.start == s.stop:
        return slice(0, 0, None)

    s_relative_to_step = s_relative_to.step if s_relative_to.step is not None else 1
    s_step = s.step if s.step is not None else 1

    if (s.start - s_relative_to.start) % s_relative_to_step != 0:
        raise ValueError("Incompatible slices")
    if s_step % s_relative_to_step != 0:
        raise ValueError("Incompatible slices")

    start = (s.start - s_relative_to.start) // s_relative_to_step
    step = s_step // s_relative_to_step
    stop = start + \
        (s_relative_to_step - 1 + s.stop - s.start) // s_relative_to_step

    if step == 1:
        step = None

    return slice(start, stop, step)


def chained_slice(s1, s2):
    """Return a slice s3 with the property that
    ar[s1][s2] == ar[s3] """

    assert (s1.start >= 0 and s2.start >= 0 and s1.stop >= 0 and s2.stop >= 0)
    s1_start = s1.start or 0
    s2_start = s2.start or 0
    s1_step = s1.step or 1
    s2_step = s2.step or 1

    start = s1_start + s2_start * s1_step
    step = s1_step * s2_step
    if s1.stop is None and s2.stop is None:
        stop = None
    elif s1.stop is None:
        stop = start + step * (s2.stop - s2_start) // s2_step
    elif s2.stop is None:
        stop = s1.stop
    else:
        stop_s2 = start + step * (s2.stop - s2_start) // s2_step
        stop_s1 = s1.stop
        stop = stop_s2 if stop_s2 < stop_s1 else stop_s1
    return slice(start, stop, step)


def index_before_slice(s, index):
    """Return an index array new_index with the property that, for a
    slice s (start, stop and step all positive), ar[s][index] ==
    ar[new_index]."""

    start = s.start or 0
    step = s.step or 1

    assert start >= 0
    assert step >= 0
    assert s.stop is None or s.stop >= 0

    new_index = start + index * step
    if s.stop is not None:
        new_index = new_index[np.where(new_index < s.stop)]

    return new_index


def concatenate_indexing(i1, i2):
    """Given either a numpy array or slice for both i1 and i2,
    return either a numpy array or slice i3 with the property that

    ar[i3] == ar[i1][i2].

    As a convenience, if i2 is None, i1 is returned
    """
    if isinstance(i1, tuple) and len(i1) == 1:
        i1 = i1[0]
    if isinstance(i2, tuple) and len(i2) == 1:
        i2 = i2[0]

    if i2 is None:
        return i1
    if isinstance(i1, slice) and isinstance(i2, slice):
        return chained_slice(i1, i2)
    elif isinstance(i1, slice) and isinstance(i2, (np.ndarray, list)):
        return index_before_slice(i1, i2)
    elif isinstance(i1, (np.ndarray, list)) and isinstance(i2, (slice, np.ndarray)):
        return np.asarray(i1)[i2]
    else:
        raise TypeError("Don't know how to chain these index types")


def indexing_length(sl_or_ar):
    """Given either an array or slice, return len(ar[sl_or_ar]) for any
    array ar which is large enough that the slice does not overrun it."""

    if isinstance(sl_or_ar, slice):
        step = sl_or_ar.step or 1
        diff = (sl_or_ar.stop - sl_or_ar.start)
        return diff // step + (diff % step > 0)
    else:
        return len(sl_or_ar)


def arrays_are_same(a1, a2):
    """Returns True if a1 and a2 are numpy views pointing to the exact
    same underlying data; False otherwise."""
    try:
        return a1.__array_interface__['data'] == a2.__array_interface__['data'] \
            and a1.strides == a2.strides
    except AttributeError:
        return False


def set_array_if_not_same(a_store, a_in, index=None):
    """This routine checks whether a_store and a_in ultimately point to the
    same buffer; if not, the contents of a_in are copied into a_store."""
    if index is None:
        index = slice(None)
    if not arrays_are_same(a_store[index], a_in):
        a_store[index] = a_in
        if not hasattr(a_in.units, "_no_unit"):
            a_store.units = a_in.units


def index_of_first(array, find):
    """Returns the index to the first element in array
    which satisfies array[index]>=find. The array must
    be sorted in ascending order."""

    if len(array) == 0:
        return 0

    left = 0
    right = len(array) - 1

    if array[left] >= find:
        return 0

    if array[right] < find:
        return len(array)

    while right - left > 1:
        mid = (left + right) // 2
        if array[mid] >= find:
            right = mid
        else:
            left = mid

    return right


def equipartition(ar, nbins, min=None, max=None):
    """

    Given an array ar, return nbins+1 monotonically increasing bin
    edges such that the number of items in each bin is approximately
    equal.

    """

    a_s = np.sort(ar)

    if max is not None:
        a_s = a_s[a_s <= max]
    if min is not None:
        a_s = a_s[a_s > min]

    return a_s[np.array(np.linspace(0, len(a_s) - 1, nbins + 1), dtype='int')]


def bisect(left, right, f, epsilon=None, eta=0, verbose=False, niter_max=200):
    """

    Finds the value x such that f(x)=0 for a monotonically increasing
    function f, using a binary search.

    The search stops when either the bounding domain is smaller than
    epsilon (by default 10^-7 times the original region) OR a value
    f(x) is found such that |f(x)|<eta (by default eta=0, so this
    criterion is never satisfied).

    """

    if epsilon is None:
        epsilon = (right - left) * 1.e-7

    logger.info("Entering bisection search algorithm")
    for i in range(niter_max):

        if (right - left) < epsilon:
            return (right + left) / 2

        mid = (left + right) / 2
        z = f(mid)

        logger.info("%f %f %f %f" % (left, mid, right, z))

        if (abs(z) < eta):
            return mid
        elif(z < 0):
            left = mid
        else:
            right = mid

    raise ValueError("Bisection algorithm did not converge")


def gauss_jordan(out):
    """A simple Gauss-Jordan matrix inverter. This is provided so that
    matrices of fractions can be inverted (numpy linalg converts
    everything to floats first.)

    Don't use on large matrices -- it's slow!

    Based on public domain code by Jarno Elonen."""

    h, w = out.shape

    assert w > h

    for y in range(0, h):

        maxrow = out[y:, y].argmax() + y

        (out[y], out[maxrow]) = (out[maxrow], out[y].copy())

        if out[y][y] == 0:
            # this will be a problem, see if we can do a row
            # operation to fix it
            for y2 in range(y+1,h):
                if out[y2][y]!=0:
                    out[y]+=out[y2]
                    break

            # no, out of options, must be a singular matrix
            if out[y][y]==0:
                raise np.linalg.linalg.LinAlgError("Singular matrix")

        for y2 in range(y + 1, h):    # Eliminate column y
            c = out[y2][y] / out[y][y]
            out[y2] -= out[y] * c

    for y in range(h - 1, 0 - 1, -1):  # Backsubstitute
        c = out[y][y]
        for y2 in range(0, y):
            for x in range(w - 1, y - 1, -1):
                out[y2][x] -= out[y][x] * out[y2][y] / c
        out[y][y] /= c
        for x in range(h, w):       # Normalize row y
            out[y][x] /= c

    return out


def rational_matrix_inv(matrix):
    """A simple replacement for numpy linalg matrix inverse
    which handles fractions exactly. Not suitable for large
    matrices!"""

    assert len(matrix) == len(matrix[0])
    x = np.ndarray(
        shape=(len(matrix), len(matrix[0]) + len(matrix)), dtype=fractions.Fraction)
    x[:, :] = fractions.Fraction(0)
    for i in range(len(x)):
        x[i, len(x) + i] = fractions.Fraction(1)

    for i in range(len(x)):
        for j in range(len(x)):
            x[i, j] = fractions.Fraction(matrix[i][j])

    return gauss_jordan(x)[:, len(x):]


def random_rotation_matrix():
    """Return a random rotation matrix (Haar measure for 3x3 case), using
    fast algorithm from Graphics Gems III

    (http://tog.acm.org/resources/GraphicsGems/gemsiii/rand_rotation.c)"""

    x = np.random.uniform(size=3)
    theta = x[0]*2*math.pi
    phi = x[1]*2*math.pi
    z = x[2]*2

    r = math.sqrt(z)
    vx = math.sin(phi)*r
    vy = math.cos(phi)*r
    vz = math.sqrt(2.0-z)

    st = math.sin(theta)
    ct = math.cos(theta)

    sx = vx*ct-vy*st
    sy = vx*st+vy*ct

    return np.array([[vx*sx-ct, vx*sy-st, vx*vz],
                     [vy*sx+st, vy*sy-ct, vy*vz],
                     [vz*sx,vz*sy,1.0-z]])


def cutgz(x):
    """Strip the .gz ending off a string"""
    if x[-3:] == '.gz':
        return x[:-3]
    else:
        return x


class ExecutionControl(object):

    def __init__(self):
        self.count = 0
        self.on_exit = None

    def __enter__(self):
        self.count += 1

    def __exit__(self, *excp):
        self.count -= 1
        assert self.count >= 0
        if self.count == 0 and self.on_exit is not None:
            self.on_exit()

    def __bool__(self):
        return self.count > 0

    def __repr__(self):
        return "<ExecutionControl: %s>" % ('True' if self.count > 0 else 'False')


#################################################################
# Code for incomplete gamma function accepting complex arguments
#################################################################

def _gser(a, x, eps=3.e-7, itmax=700):
    """Series representation of the incomplete gamma
    function, based on numerical recipes 3rd ed"""
    if x == 0.0:
        return 0.0
    ap = a
    sum = 1. / a
    delta = sum
    n = 1
    while n <= itmax:
        ap = ap + 1.
        delta = delta * x / ap
        sum = sum + delta
        if (abs(delta) < abs(sum) * eps):
            return (sum * np.exp(-x + a * np.log(x)))
        n = n + 1
    raise RuntimeError("Maximum iterations exceeded in gser")


def _gcf(a, x, eps=3.e-7, itmax=200):
    """Continued fraction representation of the incomplete gamma
    function, based on numerical recipes 3rd ed"""

    gold = 0.
    a0 = 1.
    a1 = x
    b0 = 0.
    b1 = 1.
    fac = 1.
    n = 1
    while n <= itmax:
        an = n
        ana = an - a
        a0 = (a1 + a0 * ana) * fac
        b0 = (b1 + b0 * ana) * fac
        anf = an * fac
        a1 = x * a0 + anf * a1
        b1 = x * b0 + anf * b1
        if (a1 != 0.):
            fac = 1. / a1
            g = b1 * fac
            if (abs((g - gold) / g) < eps):
                return (g * np.exp(-x + a * np.log(x)))
            gold = g
            n = n + 1
    raise RuntimeError("Maximum iterations exceeded in gcf")


def gamma_inc(a, z, eps=3.e-7):
    """Incomplete gamma function accepting complex z, based on algorithm
    given in numerical recipes (3rd ed)"""
    import scipy
    import scipy.special

    if (abs(z) < a + 1.):
        return _gser(a, z, eps)
    else:
        return scipy.special.gamma(a) - _gcf(a, z, eps)


#
# THREAD-SAFE VERSION OF scipy.weave.inline
#

compile_lock = threading.Lock()


def threadsafe_inline(*args, **kwargs):
    """When scipy.weave.inline is called, it may trigger a compile. We
    only want one compilation to be going on at once, otherwise nasty
    race conditions arise. This function wraps scipy.weave.inline to
    be thread-safe."""

    import scipy.weave

    call_frame = sys._getframe().f_back
    if 'local_dict' not in kwargs:
        kwargs['local_dict'] = call_frame.f_locals
    if 'global_dict' not in kwargs:
        kwargs['global_dict'] = call_frame.f_globals

    tid = threading.currentThread().name
    while args[0] not in scipy.weave.inline_tools.function_cache:
        # We need a compilation, so try to acquire the compile lock
        if compile_lock.acquire(False):
            # acquired lock
            try:
                ret = scipy.weave.inline(*args, **kwargs)
            finally:
                compile_lock.release()
            return ret
        else:
            # didn't acquire lock. Wait a while
            time.sleep(1)

    # When we reach this point, we know no compilation will be
    # triggered, so go ahead and call
    return scipy.weave.inline(*args, **kwargs)



##################################################
# fortran reading facilities for ramses and grafic
#################################################

_head_type = np.dtype('i4')

class FortranFile(object):
    """Utilities to help reading fortran files efficiently, using a numpy memmap

    Usage:

    with FortranFile(fname) as f:
        data = f.read_field(dtype, size) # loads a fortran field consisting of size elements of type dtype
        f.skip_fields(num) # skip over <num> fortran fields
        data_as_map = f.read_field_memmapped(dtype, size) # loads the fortran field as a RO memmap straight onto the file
        header = f.read_series(fields_dtype) # load a series of fortran fields defined by the composite numpy dtype
    """

    def __init__(self, filename):
        self._map = np.memmap(filename, mode='r')
        self._offset = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del self._map

    def get_raw_memmapped(self, dtype, length=1):
        start = self._offset
        end = start + length * dtype.itemsize
        result = np.frombuffer(self._map[start:end], dtype=dtype)
        self._offset = end
        return result

    def skip(self, bytes):
        self._offset+=bytes

    def skip_fields(self, n=1):
        for i in range(n):
            alen = self.get_raw_memmapped(_head_type, 1)
            self.skip(alen[0])
            alen2 = self.get_raw_memmapped(_head_type, 1)
            assert alen==alen2

    def read_field(self, dtype, field_length=1):
        return self.get_field_memmapped(dtype, field_length).copy()

    def get_field_memmapped(self, dtype, field_length=1):
        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)

        length = field_length * dtype.itemsize
        alen = self.get_raw_memmapped(_head_type, 1)
        if alen != length:
            raise IOError("Unexpected FORTRAN block length %d!=%d" % (
                alen, length))

        data = self.get_raw_memmapped(dtype, field_length)

        alen = self.get_raw_memmapped(_head_type, 1)
        if alen != length:
            raise IOError("Unexpected FORTRAN block length (tail) %d!=%d" % (
                alen, length))

        return data

    def read_series(self, dtype):
        q = np.empty(1, dtype=dtype)
        for i in range(len(dtype.fields)):
            data = self.read_field(dtype[i], 1)

            # I really don't understand why the following acrobatic should
            # be necessary, but q[0][i] = data[0] doesn't copy arrays properly
            if hasattr(data[0], "__len__"):
                q[0][i][:] = data[0]
            else:
                q[0][i] = data[0]

        return q[0]

def _thread_map(func, *args):

    def r_func(*afunc):
        try:
            this_t = threading.current_thread()
            this_t.ret_value = func(*afunc)
        except Exception as e:
            this_t.ret_excp = e

    threads = []
    for arg_this in zip(*args):
        threads.append(threading.Thread(target=r_func, args=arg_this))
        threads[-1].start()
    rets = []
    excp = None
    for t in threads:
        while t.is_alive():
            # just calling t.join() with no timeout can make it harder to
            # debug deadlocks!
            t.join(1.0)
        if hasattr(t, 'ret_excp'):
            excp = t.ret_excp
        else:
            rets.append(t.ret_value)

    if excp is None:
        return rets
    raise excp  # Note this is a re-raised exception from within a thread


def parallel(p_args=[0],
             threads=config['number_of_threads'], reduce='interleave'):
    """Return a function decorator which makes a function execute in parallel.

    *p_args*: a list of integers specifying which arguments will be
    array-like. These will be sliced up and processed over the number
    of threads specified.

    *reduce*: specifies how to reduce the output, and can take one of
    three values:

    'none': return None

    'interleave': return an array of the size of the input arrays,
     with the elements returned to their correct place

    'sum': sum over the outputs

    """
    def decorator(fn):
        def new_fn(*args):
            threaded_args = []
            for i in range(len(args)):
                if i in p_args:
                    threaded_args.append(
                        [args[i][n::threads] for n in range(threads)])
                else:
                    threaded_args.append([args[i]] * threads)

            ret = _thread_map(fn, *threaded_args)

            if reduce == 'interleave':
                out_len = list(getattr(ret[0], 'shape', [0]))
                out_len[0] = sum(map(len, ret))
                output = np.empty(out_len, dtype=ret[0].dtype)
                for n in range(threads):
                    output[n::threads] = ret[n]
            elif reduce == 'sum':
                output = sum(ret)
            else:
                output = None
            return output
        return functools.wraps(fn)(new_fn)
    return decorator
