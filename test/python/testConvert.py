import os
import unittest
import numpy as np
import ck2cti

import Cantera as ct

class chemkinConverterTest(unittest.TestCase):
    def checkConversion(self, refFile, testFile):
        ref = ct.IdealGasMix(refFile)
        gas = ct.IdealGasMix(testFile)

        self.assertEqual(ref.elementNames(), gas.elementNames())
        self.assertEqual(ref.speciesNames(), gas.speciesNames())
        self.assertTrue((ref.reactantStoichCoeffs() == gas.reactantStoichCoeffs()).all())

        compositionA = [[ref.nAtoms(i,j) for j in range(ref.nElements())]
                        for i in range(ref.nSpecies())]
        compositionB = [[gas.nAtoms(i,j) for j in range(gas.nElements())]
                        for i in range(gas.nSpecies())]
        self.assertEqual(compositionA, compositionB)

    def test_gri30(self):
        if os.path.exists('gri30_test.cti'):
            os.remove('gri30_test.cti')

        ck2cti.convertMech('../../data/inputs/gri30.inp',
                           transportFile='../../data/transport/gri30_tran.dat',
                           outName='gri30_test.cti', quiet=True)

        self.checkConversion('gri30.xml', 'gri30_test.cti')

    def test_soot(self):
        if os.path.exists('soot_test.cti'):
            os.remove('soot_test.cti')

        ck2cti.convertMech('../data/soot.inp',
                           thermoFile='../data/soot-therm.dat',
                           outName='soot_test.cti', quiet=True)

        self.checkConversion('../data/soot.xml', 'soot_test.cti')

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
