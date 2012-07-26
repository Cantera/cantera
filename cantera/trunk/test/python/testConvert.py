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

