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

ctdata = ''
try:
    import cantera_loc
    ctdata = cantera_loc.CANTERA_DATA
except:
    pass

from Mixture import Mixture

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

def addDirectory(dir):
    import _cantera
    return _cantera.ct_addDirectory(dir)

if ctdata:
    addDirectory(ctdata)
