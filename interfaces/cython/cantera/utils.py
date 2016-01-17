import os
import inspect

from . import Solution, add_directory


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
    add_directory(os.path.dirname(os.path.abspath(inspect.stack()[1][1])))
