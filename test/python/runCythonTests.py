"""
Unit tests for Cantera's Cython-based Python module.

This script gathers all the tests defined 'cantera.test' module, runs them,
and prints a report. Extra command line arguments can be used to run subsets
of the test suite, e.g.:

all tests from 'test_thermo.py' and 'test_kinetics.py':

    python runCythonTests.py thermo kinetics

all tests from the 'test_reactor.TesTIdealGasReactor' class:

    python runCythonTests.py reactor.TestIdealGasReactor

a single test:

    python runCythonTests.py onedim.TestDiffusionFlame.test_mixture_averaged
"""
from __future__ import print_function

import sys
import os

cantera_root = os.path.relpath(__file__).split(os.sep)[:-1] + ['..', '..']
py_version = 'python3' if sys.version_info[0] == 3 else 'python2'
module_path = os.path.abspath(os.sep.join(cantera_root + ['build', py_version]))

if 'PYTHONPATH' in os.environ:
    os.environ['PYTHONPATH'] = module_path + os.path.pathsep + os.environ['PYTHONPATH']
else:
    os.environ['PYTHONPATH'] = module_path

sys.path.insert(0, module_path)
os.chdir(os.sep.join(cantera_root + ['test', 'work']))

from cantera.test.utilities import unittest
import cantera
import cantera.test

class TestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        unittest.TextTestResult.__init__(self, *args, **kwargs)
        self.outName = 'python%d-results.txt' % sys.version_info[0]
        with open(self.outName, 'w') as f:
            pass # just create an empty output file

    def reformat(self, test_string):
        name, cls = test_string.split()
        cls = cls.replace('(cantera.test.', '').replace(')','')
        return '%s.%s' % (cls, name)

    def addSuccess(self, test):
        with open(self.outName, 'a') as f:
            f.write('PASS: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addSuccess(self, test)

    def addFailure(self, test, err):
        with open(self.outName, 'a') as f:
            f.write('FAIL: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addFailure(self, test, err)

    def addError(self, test, err):
        with open(self.outName, 'a') as f:
            f.write('ERROR: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addFailure(self, test, err)


if __name__ == '__main__':
    print('\n* INFO: using Cantera module found at this location:')
    print('*     ', repr(cantera.__file__), '\n')
    sys.stdout.flush()

    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner(verbosity=2, resultclass=TestResult)
    suite = unittest.TestSuite()
    subsets = []
    for name in sys.argv[1:]:
        subsets.append('cantera.test.test_' + name)

    if not subsets:
        subsets.append('cantera.test')

    suite = loader.loadTestsFromNames(subsets)

    results = runner.run(suite)
    sys.exit(len(results.errors) + len(results.failures))
