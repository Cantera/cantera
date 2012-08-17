import os
import unittest
import numpy as np
import itertools
import ck2cti

import utilities
import Cantera as ct

class chemkinConverterTest(utilities.CanteraTest):
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

        return ref, gas

    def checkThermo(self, ref, gas, temperatures):
        for T in temperatures:
            ref.set(T=T, P=ct.OneAtm)
            gas.set(T=T, P=ct.OneAtm)
            ref_cp = ref.cp_R()
            gas_cp = gas.cp_R()
            ref_h = ref.enthalpies_RT()
            gas_h = gas.enthalpies_RT()
            ref_s = ref.entropies_R()
            gas_s = gas.entropies_R()
            for i in range(gas.nSpecies()):
                self.assertNear(ref_cp[i], gas_cp[i], 1e-7)
                self.assertNear(ref_h[i], gas_h[i], 1e-7)
                self.assertNear(ref_s[i], gas_s[i], 1e-7)

    def checkKinetics(self, ref, gas, temperatures, pressures):
        for T,P in itertools.product(temperatures, pressures):
            ref.set(T=T, P=P)
            gas.set(T=T, P=P)
            ref_kf = ref.fwdRateConstants()
            ref_kr = ref.revRateConstants()
            gas_kf = gas.fwdRateConstants()
            gas_kr = gas.revRateConstants()
            for i in range(gas.nReactions()):
                self.assertNear(ref_kf[i], gas_kf[i])
                self.assertNear(ref_kr[i], gas_kr[i])


    def test_gri30(self):
        if os.path.exists('gri30_test.cti'):
            os.remove('gri30_test.cti')

        ck2cti.convertMech('../../data/inputs/gri30.inp',
                           transportFile='../../data/transport/gri30_tran.dat',
                           outName='gri30_test.cti', quiet=True)

        ref, gas = self.checkConversion('gri30.xml', 'gri30_test.cti')
        self.checkKinetics(ref, gas, [300, 1500], [5e3, 1e5, 2e6])

    def test_soot(self):
        if os.path.exists('soot_test.cti'):
            os.remove('soot_test.cti')

        ck2cti.convertMech('../data/soot.inp',
                           thermoFile='../data/soot-therm.dat',
                           outName='soot_test.cti', quiet=True)

        ref, gas = self.checkConversion('../data/soot.xml', 'soot_test.cti')
        self.checkThermo(ref, gas, [300, 1100])
        self.checkKinetics(ref, gas, [300, 1100], [5e3, 1e5, 2e6])

    def test_pdep(self):
        if os.path.exists('pdep_test.cti'):
            os.remove('pdep_test.cti')

        ck2cti.convertMech('../data/pdep-test.inp',
                           outName='pdep_test.cti', quiet=True)

        ref, gas = self.checkConversion('../data/pdep-test.xml', 'pdep_test.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])


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

    def test_nasa9(self):
        if os.path.exists('nasa9_test.cti'):
            os.remove('nasa9_test.cti')
        ck2cti.convertMech('../data/nasa9-test.inp',
                           thermoFile='../data/nasa9-test-therm.dat',
                           outName='nasa9_test.cti', quiet=True)

        ref, gas = self.checkConversion('../data/nasa9-test.xml',
                                        'nasa9_test.cti')
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_sri_falloff(self):
        if os.path.exists('sri-falloff.cti'):
            os.remove('sri-falloff.cti')
        ck2cti.convertMech('../data/sri-falloff.inp',
                           thermoFile='../data/dummy-thermo.dat',
                           outName='sri-falloff.cti', quiet=True)

        ref, gas = self.checkConversion('../data/sri-falloff.xml',
                                        'sri-falloff.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_explicit_reverse_rate(self):
        if os.path.exists('explicit-reverse-rate.cti'):
            os.remove('explicit-reverse-rate.cti')
        ck2cti.convertMech('../data/explicit-reverse-rate.inp',
                           thermoFile='../data/dummy-thermo.dat',
                           outName='explicit-reverse-rate.cti', quiet=True)
        ref, gas = self.checkConversion('../data/explicit-reverse-rate.xml',
                                        'explicit-reverse-rate.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

        # Reactions with explicit reverse rate constants are transformed into
        # two irreversible reactions with reactants and products swapped.
        Rr = gas.revRateConstants()
        self.assertEqual(Rr[0], 0.0)
        self.assertEqual(Rr[1], 0.0)
        Rstoich = gas.reactantStoichCoeffs()
        Pstoich = gas.productStoichCoeffs()
        self.assertEqual(list(Rstoich[:,0]), list(Pstoich[:,1]))
        self.assertEqual(list(Rstoich[:,1]), list(Pstoich[:,0]))
