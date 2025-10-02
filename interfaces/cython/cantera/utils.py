# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import os
import inspect as _inspect

from ._utils import add_data_directory
from .composite import Solution

def import_phases(filename, phase_names):
    """
    Import multiple phases from one file. The phase names should be
    entered as a list of strings.
    """
    return [Solution(filename, p) for p in phase_names]


def add_module_directory():
    """
    Add the directory containing the module from which this function is called
    to the Cantera input file search path.
    """
    add_data_directory(os.path.dirname(os.path.abspath(_inspect.stack()[1][1])))
