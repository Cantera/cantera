"""
Unit tests for Cantera's Cython-based Python module.

This script gathers all the tests defined 'cantera.test' module, runs them,
and prints a report.
"""
from __future__ import print_function

import sys
import os

cantera_root = os.getcwd().split(os.sep)[:-2]
if sys.version_info[0] == 3:
    sys.path.insert(0, os.sep.join(cantera_root + ['build', 'python3']))
else:
    sys.path.insert(0, os.sep.join(cantera_root + ['build', 'python2']))

from cantera.test.utilities import unittest
import cantera

if __name__ == '__main__':
    print('\n* INFO: using Cantera module found at this location:')
    print('*     ', repr(cantera.__file__), '\n')
    sys.stdout.flush()

    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner(verbosity=2)
    suite = loader.loadTestsFromName('cantera.test')

    results = runner.run(suite)
    sys.exit(len(results.errors) + len(results.failures))
