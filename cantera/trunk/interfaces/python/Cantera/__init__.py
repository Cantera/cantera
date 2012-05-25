"""
  Cantera provides capabilities for simulating problems involving
  chemical kinetics and transport processes.
"""

import types
import _cantera
from num import *
from constants import *
from exceptions import *
from gases import *
from set import set
from importFromFile import *

import os as _os
import sys as _sys

if not os.getenv('PYTHON_CMD'):
    # Setting PYTHON_CMD here avoids issues with .cti -> .xml conversions
    # in cases where the python interpreter isn't in the system path.
    os.environ['PYTHON_CMD'] = _sys.executable

def writeCSV(f, lst):
    """
    Write list items to file *f* in
    comma-separated-value format.  Strings will be written as-is, and
    other types of objects will be converted to strings and then
    written. Each call to writeCSV writes one line of the file.
    """
    for i,item in enumerate(lst):
        if type(item) == types.StringType:
            f.write(item)
        else:
            f.write(repr(item))
        if i != len(lst)-1:
            f.write(',')

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
    return _cantera.get_Cantera_Error()

def refCount(a):
    """Return the reference count for an object."""
    return _cantera.ct_refcnt(a)

def addDirectory(dir):
    """Add a directory to search for Cantera data files."""
    return _cantera.ct_addDirectory(dir)

def writeLogFile(file):
    return _cantera.ct_writelogfile(file)


def reset():
    """Release all cached Cantera data. Equivalent to
    starting a fresh session."""
    _cantera.ct_appdelete()

# workaround for case problems in CVS repository file Mixture.py. On some
# systems it appears as mixture.py, and on others as Mixture.py
try:
    from Mixture import Mixture
except:
    from mixture import Mixture

from num import *
