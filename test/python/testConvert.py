import os
import unittest
import numpy as np
import itertools
import ck2cti

import utilities
import Cantera as ct


def convertMech(inputFile, outName=None, **kwargs):
    if os.path.exists(outName):
        os.remove(outName)
    parser = ck2cti.Parser()
    parser.convertMech(inputFile, outName=outName, **kwargs)


class chemkinConverterTest(utilities.CanteraTest):
    def checkConversion(self, refFile, testFile):
        ref = ct.IdealGasMix(refFile)
        gas = ct.IdealGasMix(testFile)

        self.assertEqual(ref.elementNames(), gas.elementNames())
        self.assertEqual(ref.speciesNames(), gas.speciesNames())
        coeffs_ref = ref.reactantStoichCoeffs()
        coeffs_gas = gas.reactantStoichCoeffs()
        self.assertEqual(coeffs_gas.shape, coeffs_ref.shape)
        self.assertTrue((coeffs_gas == coeffs_ref).all())

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
                message = ' for species {0} at T = {1}'.format(i, T)
                self.assertNear(ref_cp[i], gas_cp[i], 1e-7, msg='cp'+message)
                self.assertNear(ref_h[i], gas_h[i], 1e-7, msg='h'+message)
                self.assertNear(ref_s[i], gas_s[i], 1e-7, msg='s'+message)

    def checkKinetics(self, ref, gas, temperatures, pressures, tol=1e-8):
        for T,P in itertools.product(temperatures, pressures):
            ref.set(T=T, P=P)
            gas.set(T=T, P=P)
            ref_kf = ref.fwdRateConstants()
            ref_kr = ref.revRateConstants()
            gas_kf = gas.fwdRateConstants()
            gas_kr = gas.revRateConstants()
            for i in range(gas.nReactions()):
                message = ' for reaction {0} at T = {1}, P = {2}'.format(i, T, P)
                self.assertNear(ref_kf[i], gas_kf[i], rtol=tol, msg='kf '+message)
                self.assertNear(ref_kr[i], gas_kr[i], rtol=tol, msg='kr '+message)

    def test_gri30(self):
        convertMech('../../data/inputs/gri30.inp',
                    transportFile='../../data/transport/gri30_tran.dat',
                    outName='gri30_test.cti', quiet=True)

        ref, gas = self.checkConversion('gri30.xml', 'gri30_test.cti')
        self.checkKinetics(ref, gas, [300, 1500], [5e3, 1e5, 2e6])

    def test_soot(self):
        convertMech('../data/soot.inp',
                    thermoFile='../data/soot-therm.dat',
                    outName='soot_test.cti', quiet=True)

        ref, gas = self.checkConversion('../data/soot.xml', 'soot_test.cti')
        self.checkThermo(ref, gas, [300, 1100])
        self.checkKinetics(ref, gas, [300, 1100], [5e3, 1e5, 2e6])

    def test_pdep(self):
        convertMech('../data/pdep-test.inp',
                    outName='pdep_test.cti', quiet=True)

        ref, gas = self.checkConversion('../data/pdep-test.xml', 'pdep_test.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_missingElement(self):
        self.assertRaises(ck2cti.InputParseError,
                          lambda: convertMech('../data/h2o2_missingElement.inp',
                                              outName='h2o2_missingElement.cti',
                                              quiet=True))

    def test_missingThermo(self):
        self.assertRaises(ck2cti.InputParseError,
                          lambda: convertMech('../data/h2o2_missingThermo.inp',
                                              outName='h2o2_missingThermo.cti',
                                              quiet=True))

    def test_nasa9(self):
        convertMech('../data/nasa9-test.inp',
                    thermoFile='../data/nasa9-test-therm.dat',
                    outName='nasa9_test.cti', quiet=True)

        ref, gas = self.checkConversion('../data/nasa9-test.xml',
                                        'nasa9_test.cti')
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_sri_falloff(self):
        convertMech('../data/sri-falloff.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='sri-falloff.cti', quiet=True)

        ref, gas = self.checkConversion('../data/sri-falloff.xml',
                                        'sri-falloff.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_chemically_activated(self):
        name = 'chemically-activated-reaction'
        convertMech('../data/{0}.inp'.format(name),
                    outName='{0}.cti'.format(name), quiet=True)

        ref, gas = self.checkConversion('../data/{0}.xml'.format(name),
                                        '{0}.cti'.format(name))
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6, 1e7])

    def test_explicit_third_bodies(self):
        convertMech('../data/explicit-third-bodies.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='explicit-third-bodies.cti', quiet=True)

        ref, gas = self.checkConversion('explicit-third-bodies.cti',
                                        '../data/explicit-third-bodies.xml')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_explicit_reverse_rate(self):
        convertMech('../data/explicit-reverse-rate.inp',
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

    def test_explicit_forward_order(self):
        convertMech('../data/explicit-forward-order.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='explicit-forward-order.cti', quiet=True)
        ref, gas = self.checkConversion('../data/explicit-forward-order.xml',
                                        'explicit-forward-order.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_reaction_units(self):
        convertMech('../data/units-default.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='units-default.cti', quiet=True)
        convertMech('../data/units-custom.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='units-custom.cti', quiet=True)

        default, custom = self.checkConversion('units-default.cti',
                                               'units-custom.cti')
        self.checkKinetics(default, custom,
                           [300, 800, 1450, 2800], [5e3, 1e5, 2e6], 1e-7)

    def test_transport_normal(self):
        convertMech('../../data/inputs/h2o2.inp',
                    transportFile='../../data/transport/gri30_tran.dat',
                    outName='h2o2_transport_normal.cti', quiet=True)

        gas = ct.IdealGasMix('h2o2_transport_normal.cti')
        gas.set(X='H2:1.0, O2:1.0', T=300, P=101325)
        self.assertAlmostEqual(gas.thermalConductivity(), 0.07663, 4)

    def test_transport_missing_species(self):
        def convert():
            convertMech('../../data/inputs/h2o2.inp',
                        transportFile='../data/h2o2-missing-species-tran.dat',
                        outName='h2o2_transport_missing_species.cti',
                        quiet=True)

        self.assertRaises(ck2cti.InputParseError, convert)

    def test_transport_duplicate_species(self):
        def convert():
            convertMech('../../data/inputs/h2o2.inp',
                        transportFile='../data/h2o2-duplicate-species-tran.dat',
                        outName='h2o2_transport_duplicate_species.cti',
                        quiet=True)

        # This should fail
        self.assertRaises(ck2cti.InputParseError, convert)

        # This should succeed
        convertMech('../../data/inputs/h2o2.inp',
                    transportFile='../data/h2o2-duplicate-species-tran.dat',
                    outName='h2o2_transport_duplicate_species.cti',
                    quiet=True,
                    permissive=True)

    def test_transport_bad_geometry(self):
        def convert():
            convertMech('../../data/inputs/h2o2.inp',
                        transportFile='../data/h2o2-bad-geometry-tran.dat',
                        outName='h2o2_transport_bad_geometry.cti',
                        quiet=True)

        self.assertRaises(ck2cti.InputParseError, convert)
