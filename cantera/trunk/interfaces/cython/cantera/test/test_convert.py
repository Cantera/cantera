import os
import numpy as np
import itertools

from . import utilities
import cantera as ct
from cantera import ck2cti


def convertMech(inputFile, outName=None, **kwargs):
    if os.path.exists(outName):
        os.remove(outName)
    parser = ck2cti.Parser()
    parser.convertMech(inputFile, outName=outName, **kwargs)


class chemkinConverterTest(utilities.CanteraTest):
    def checkConversion(self, refFile, testFile):
        ref = ct.Solution(refFile)
        gas = ct.Solution(testFile)

        self.assertEqual(ref.element_names, gas.element_names)
        self.assertEqual(ref.species_names, gas.species_names)
        coeffs_ref = ref.reactant_stoich_coeffs()
        coeffs_gas = gas.reactant_stoich_coeffs()
        self.assertEqual(coeffs_gas.shape, coeffs_ref.shape)
        self.assertTrue((coeffs_gas == coeffs_ref).all())

        compositionA = [[ref.n_atoms(i,j) for j in range(ref.n_elements)]
                        for i in range(ref.n_species)]
        compositionB = [[gas.n_atoms(i,j) for j in range(gas.n_elements)]
                        for i in range(gas.n_species)]
        self.assertEqual(compositionA, compositionB)

        return ref, gas

    def checkThermo(self, ref, gas, temperatures):
        for T in temperatures:
            ref.TP = T, ct.one_atm
            gas.TP = T, ct.one_atm
            ref_cp = ref.standard_cp_R
            gas_cp = gas.standard_cp_R
            ref_h = ref.standard_enthalpies_RT
            gas_h = gas.standard_enthalpies_RT
            ref_s = ref.standard_entropies_R
            gas_s = gas.standard_entropies_R
            for i in range(gas.n_species):
                message = ' for species {0} at T = {1}'.format(i, T)
                self.assertNear(ref_cp[i], gas_cp[i], 1e-7, msg='cp'+message)
                self.assertNear(ref_h[i], gas_h[i], 1e-7, msg='h'+message)
                self.assertNear(ref_s[i], gas_s[i], 1e-7, msg='s'+message)

    def checkKinetics(self, ref, gas, temperatures, pressures, tol=1e-8):
        for T,P in itertools.product(temperatures, pressures):
            ref.TP = T, P
            gas.TP = T, P
            ref_kf = ref.forward_rate_constants
            ref_kr = ref.reverse_rate_constants
            gas_kf = gas.forward_rate_constants
            gas_kr = gas.reverse_rate_constants
            for i in range(gas.n_reactions):
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

    def test_species_only(self):
        convertMech(None,
                    thermoFile='../data/dummy-thermo.dat',
                    outName='dummy-thermo.cti', quiet=True)

        cti = "ideal_gas(elements='C H', species='dummy-thermo:R1A R1B P1')"
        gas = ct.Solution(source=cti)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 0)

    def test_missingElement(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../data/h2o2_missingElement.inp',
                        outName='h2o2_missingElement.cti',
                        quiet=True)

    def test_missingThermo(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../data/h2o2_missingThermo.inp',
                        outName='h2o2_missingThermo.cti',
                        quiet=True)

    def test_duplicate_thermo(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../data/duplicate-thermo.inp',
                        outName='duplicate-thermo.cti',
                        quiet=True)

        convertMech('../data/duplicate-thermo.inp',
                    outName='duplicate-thermo.cti',
                    quiet=True, permissive=True)

        gas = ct.Solution('duplicate-thermo.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_duplicate_species(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../data/duplicate-species.inp',
                        outName='duplicate-species.cti',
                        quiet=True)

        convertMech('../data/duplicate-species.inp',
                    outName='duplicate-species.cti',
                    quiet=True, permissive=True)

        gas = ct.Solution('duplicate-species.cti')
        self.assertEqual(gas.species_names, ['foo','bar','baz'])

    def test_pathologicalSpeciesNames(self):
        convertMech('../data/species-names.inp',
                    outName='species-names.cti', quiet=True)
        gas = ct.Solution('species-names.cti')

        self.assertEqual(gas.n_species, 6)
        self.assertEqual(gas.species_name(0), '(Parens)')
        self.assertEqual(gas.species_name(1), '@#$%^-2')
        self.assertEqual(gas.species_name(2), '[xy2]*{.}')
        self.assertEqual(gas.species_name(3), 'plus+')
        self.assertEqual(gas.species_name(4), 'eq=uals')
        self.assertEqual(gas.species_name(5), 'plus')

        self.assertEqual(gas.n_reactions, 6)
        nu = gas.product_stoich_coeffs() - gas.reactant_stoich_coeffs()
        self.assertEqual(list(nu[:,0]), [-1, -1, 2, 0, 0, 0])
        self.assertEqual(list(nu[:,1]), [-2, 3, -1, 0, 0, 0])
        self.assertEqual(list(nu[:,2]), [-1, 0, 0, 1, 0, 0])
        self.assertEqual(list(nu[:,3]), [3, 0, 0, -2, -1, 0])
        self.assertEqual(list(nu[:,4]), [2, 0, 0, -1, 0, -1])
        self.assertEqual(list(nu[:,5]), [1, 0, 0, 1, -1, -1])

    def test_unterminatedSections(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../data/unterminated-sections.inp',
                        outName='unterminated-sections.cti',
                        quiet=True)

        convertMech('../data/unterminated-sections.inp',
                    outName='unterminated-sections.cti',
                    quiet=True, permissive=True)

        gas = ct.Solution('unterminated-sections.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unterminatedSections2(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../data/unterminated-sections2.inp',
                        outName='unterminated-sections2.cti',
                        quiet=True)

        convertMech('../data/unterminated-sections2.inp',
                    outName='unterminated-sections2.cti',
                    quiet=True, permissive=True)

        gas = ct.Solution('unterminated-sections2.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

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
        Rr = gas.reverse_rate_constants
        self.assertEqual(Rr[0], 0.0)
        self.assertEqual(Rr[1], 0.0)
        Rstoich = gas.reactant_stoich_coeffs()
        Pstoich = gas.product_stoich_coeffs()
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

    def test_float_stoich_coeffs(self):
        convertMech('../data/float-stoich.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='float-stoich.cti', quiet=True)
        gas = ct.Solution('float-stoich.cti')

        R = gas.reactant_stoich_coeffs()
        P = gas.product_stoich_coeffs()
        self.assertArrayNear(R[:,0], [0, 1.5, 0.5, 0])
        self.assertArrayNear(P[:,0], [1, 0, 0, 1])
        self.assertArrayNear(R[:,1], [1, 0, 0, 1])
        self.assertArrayNear(P[:,1], [0, 0.33, 1.67, 0])

    def test_transport_normal(self):
        convertMech('../../data/inputs/h2o2.inp',
                    transportFile='../../data/transport/gri30_tran.dat',
                    outName='h2o2_transport_normal.cti', quiet=True)

        gas = ct.Solution('h2o2_transport_normal.cti')
        gas.TPX = 300, 101325, 'H2:1.0, O2:1.0'
        self.assertAlmostEqual(gas.thermal_conductivity, 0.07663, 4)

    def test_transport_embedded(self):
        convertMech('../data/with-transport.inp',
                    outName='with-transport.cti', quiet=True)

        gas = ct.Solution('with-transport.cti')
        gas.X = [0.2, 0.3, 0.5]
        D = gas.mix_diff_coeffs
        for d in D:
            self.assertTrue(d > 0.0)

    def test_transport_missing_species(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../../data/inputs/h2o2.inp',
                        transportFile='../data/h2o2-missing-species-tran.dat',
                        outName='h2o2_transport_missing_species.cti',
                        quiet=True)

    def test_transport_duplicate_species(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../../data/inputs/h2o2.inp',
                        transportFile='../data/h2o2-duplicate-species-tran.dat',
                        outName='h2o2_transport_duplicate_species.cti',
                        quiet=True)

        convertMech('../../data/inputs/h2o2.inp',
                    transportFile='../data/h2o2-duplicate-species-tran.dat',
                    outName='h2o2_transport_duplicate_species.cti',
                    quiet=True,
                    permissive=True)

    def test_transport_bad_geometry(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech('../../data/inputs/h2o2.inp',
                        transportFile='../data/h2o2-bad-geometry-tran.dat',
                        outName='h2o2_transport_bad_geometry.cti',
                        quiet=True)

    def test_empty_reaction_section(self):
        convertMech('../data/h2o2_emptyReactions.inp',
                    outName='h2o2_emptyReactions.cti',
                    quiet=True)

        gas = ct.Solution('h2o2_emptyReactions.cti')
        self.assertEqual(gas.n_species, 9)
        self.assertEqual(gas.n_reactions, 0)

    def test_reaction_comments1(self):
        convertMech('../data/pdep-test.inp',
                    outName='pdep_test.cti', quiet=True)
        text = open('pdep_test.cti').read()
        self.assertIn('Generic mechanism header', text)
        self.assertIn('Single PLOG reaction', text)
        self.assertIn('PLOG with duplicate rates and negative A-factors', text)

    def test_reaction_comments2(self):
        convertMech('../data/explicit-third-bodies.inp',
                    thermoFile='../data/dummy-thermo.dat',
                    outName='explicit_third_bodies.cti', quiet=True)
        text = open('explicit_third_bodies.cti').read()
        self.assertIn('An end of line comment', text)
        self.assertIn('A comment after the last reaction', text)


class CtmlConverterTest(utilities.CanteraTest):
    def test_sofc(self):
        gas_a, anode_bulk, oxide_a = ct.import_phases(
            '../../interfaces/cython/cantera/examples/surface_chemistry/sofc.cti',
            ['gas', 'metal', 'oxide_bulk'])

        self.assertNear(gas_a.P, ct.one_atm)
        self.assertNear(anode_bulk['electron'].X, 1.0)
        self.assertNear(oxide_a.density, 700)

    def test_diamond(self):
        gas, solid = ct.import_phases('diamond.cti', ['gas','diamond'])
        face = ct.Interface('diamond.cti', 'diamond_100', [gas, solid])

        self.assertNear(face.site_density, 3e-8)

    def test_pdep(self):
        gas = ct.Solution('../data/pdep-test.cti')
        self.assertEqual(gas.n_reactions, 6)

    def test_invalid(self):
        try:
            gas = ct.Solution('../data/invalid.cti')
        except RuntimeError as e:
            err = e

        self.assertIn('Multiply-declared species', err.args[0])

    def test_noninteger_atomicity(self):
        gas = ct.Solution('../data/noninteger-atomicity.cti')
        self.assertNear(gas.molecular_weights[gas.species_index('CnHm')],
                        10.65*gas.atomic_weight('C') + 21.8*gas.atomic_weight('H'))
