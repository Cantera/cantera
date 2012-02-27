"""
Unit tests for the Cantera Python module.

This script gathers all the tests defined in all of the test<foo>.py
files, runs them, and prints a report.
"""
import unittest

import Cantera

if __name__ == '__main__':
    print '\n* INFO: using Cantera module found at this location:'
    print '*     ', repr(Cantera.__file__), '\n'

    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner(verbosity=2)
    suite = loader.loadTestsFromName('testThermo')
    runner.run(suite)
