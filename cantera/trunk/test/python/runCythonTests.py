"""
Unit tests for Cantera's Cython-based Python module.

This script gathers all the tests defined 'cantera.test' module, runs them,
and prints a report. Extra command line arguments can be used to run subsets
of the test suite, e.g.

    python runCythonTests.py thermo kinetics
    python runCythonTests.py onedim reactor
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
import cantera.test

class TestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        unittest.TextTestResult.__init__(self, *args, **kwargs)
        self.outfile = open('python%d-results.txt' % sys.version_info[0], 'w')

    def reformat(self, test_string):
        name, cls = test_string.split()
        cls = cls.replace('(cantera.test.', '').replace(')','')
        return '%s.%s' % (cls, name)

    def addSuccess(self, test):
        self.outfile.write('PASS: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addSuccess(self, test)

    def addFailure(self, test, err):
        self.outfile.write('FAIL: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addFailure(self, test, err)

    def addError(self, test, err):
        self.outfile.write('ERROR: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addFailure(self, test, err)


if __name__ == '__main__':
    print('\n* INFO: using Cantera module found at this location:')
    print('*     ', repr(cantera.__file__), '\n')
    sys.stdout.flush()

    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner(verbosity=2, resultclass=TestResult)
    suite = unittest.TestSuite()
    subsets = []
    for name in dir(cantera.test):
        if name.startswith('test_') and name[5:] in sys.argv:
            subsets.append('cantera.test.' + name)

    if not subsets:
        subsets.append('cantera.test')

    suite = loader.loadTestsFromNames(subsets)

    results = runner.run(suite)
    sys.exit(len(results.errors) + len(results.failures))
