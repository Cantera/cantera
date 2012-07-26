import os
import unittest
import numpy as np
import ck2cti

import Cantera as ct

class chemkinConverterTest(unittest.TestCase):
    def test_gri30(self):
        if os.path.exists('gri30_test.cti'):
            os.remove('gri30_test.cti')

        ck2cti.convertMech('../../data/inputs/gri30.inp',
                           transportFile='../../data/transport/gri30_tran.dat',
                           outName='gri30_test.cti', quiet=True)

        ref = ct.IdealGasMix('gri30.xml')
        gas = ct.IdealGasMix('gri30_test.cti')

        self.assertEqual(ref.speciesNames(), gas.speciesNames())
        self.assertTrue((ref.reactantStoichCoeffs() == gas.reactantStoichCoeffs()).all())

    def test_missingElement(self):
        if os.path.exists('h2o2_missingElement.cti'):
            os.remove('h2o2_missingElement.cti')

        self.assertRaises(ck2cti.InputParseError,
                          lambda: ck2cti.convertMech('../data/h2o2_missingElement.inp',
                                                    quiet=True))

    def test_missingThermo(self):
        if os.path.exists('h2o2_missingThermo.cti'):
            os.remove('h2o2_missingThermo.cti')

        self.assertRaises(ck2cti.InputParseError,
                          lambda: ck2cti.convertMech('../data/h2o2_missingThermo.inp',
                                                    quiet=True))
