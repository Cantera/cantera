"""
  Cantera provides capabilities for simulating problems involving
  chemical kinetics and transport processes.
"""

import types
#ok = 0
from constants import *
from exceptions import *
from gases import *
from set import set
from importFromFile import *

try:
    from Numeric import array, asarray, zeros, ones
except:
    print """

    ERROR: NumPy not found!
    
    Cantera uses the NumPy set of numerical extensions to
    Python. These do not appear to be present on your system, or else
    the 'Numeric' package is not on your Python path. You can get
    NumPy from http://sourceforge.net/projects/numpy, or see
    http://www.python.org/topics/scicomp/numpy.html for more
    information. If you are working on a Windows PC, an installer
    program is available from the sourceforge site that will install
    NumPy for you. The whole process takes less than a minute.

    """
    raise "could not import Numeric"


def writeCSV(f, list):
    """Write list items to file 'f' in comma-separated-value format."""
    for item in list:
        if type(item) == types.StringType:
            f.write(item+', ')
        else:
            f.write(`item`+', ')
    f.write('\n')


def table(keys, values):
    """Create a map with the keys and values specified."""
    x = {}
    pairs = map(None, keys, values)
    for p in pairs:
        k, v = p
        x[k] = v
    return x

def getCanteraError():
    """Return the Cantera error message, if any."""
    import _cantera
    return _cantera.get_Cantera_Error() 

def refCount(a):
    """Return the reference count for an object."""
    import _cantera
    return _cantera.ct_refcnt(a)
