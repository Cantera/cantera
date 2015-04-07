"""
Unit tests for the Cantera Python module.

This script gathers all the tests defined in all of the test<foo>.py
files, runs them, and prints a report.
"""
import sys
import os
import unittest

sys.path.insert(0, os.path.abspath('../../interfaces/python'))

import Cantera

Cantera.addDirectory(os.path.join(os.path.split(os.getcwd())[0], 'data'))

if __name__ == '__main__':
    print '\n* INFO: using Cantera module found at this location:'
    print '*     ', repr(Cantera.__file__), '\n'
    sys.stdout.flush()

    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner(verbosity=2)
    suite = loader.loadTestsFromName('testSolution')
    suite.addTests(loader.loadTestsFromName('testKinetics'))
    suite.addTests(loader.loadTestsFromName('testPureFluid'))
    suite.addTests(loader.loadTestsFromName('testEquilibrium'))
    suite.addTests(loader.loadTestsFromName('testReactors'))
    suite.addTests(loader.loadTestsFromName('testConvert'))

    results = runner.run(suite)
    sys.exit(len(results.errors) + len(results.failures))
