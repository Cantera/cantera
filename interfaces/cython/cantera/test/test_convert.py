import os
from os.path import join as pjoin
import itertools

from . import utilities
import cantera as ct
from cantera import ck2cti, cti2yaml


def convertMech(inputFile, outName=None, **kwargs):
    if os.path.exists(outName):
        os.remove(outName)
    ck2cti.convertMech(inputFile, outName=outName, **kwargs)


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
        convertMech(pjoin(self.test_data_dir, 'gri30.inp'),
                    transportFile=pjoin(self.test_data_dir, 'gri30_tran.dat'),
                    outName=pjoin(self.test_work_dir, 'gri30_test.cti'), quiet=True)

        ref, gas = self.checkConversion('gri30.xml', 'gri30_test.cti')
        self.checkKinetics(ref, gas, [300, 1500], [5e3, 1e5, 2e6])

    def test_soot(self):
        convertMech(pjoin(self.test_data_dir, 'soot.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'soot-therm.dat'),
                    outName=pjoin(self.test_work_dir, 'soot_test.cti'), quiet=True)

        ref, gas = self.checkConversion('soot.xml', 'soot_test.cti')
        self.checkThermo(ref, gas, [300, 1100])
        self.checkKinetics(ref, gas, [300, 1100], [5e3, 1e5, 2e6])

    def test_pdep(self):
        convertMech(pjoin(self.test_data_dir, 'pdep-test.inp'),
                    outName=pjoin(self.test_work_dir, 'pdep_test.cti'), quiet=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'pdep-test.xml'),
                                        pjoin(self.test_work_dir, 'pdep_test.cti'))
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_species_only(self):
        convertMech(None,
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'dummy-thermo.cti'), quiet=True)

        cti = "ideal_gas(elements='C H', species='dummy-thermo:R1A R1B P1')"
        gas = ct.Solution(source=cti)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 0)

    def test_missingElement(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'Undefined elements'):
            convertMech(pjoin(self.test_data_dir, 'h2o2_missingElement.inp'),
                        outName=pjoin(self.test_work_dir, 'h2o2_missingElement.cti'),
                        quiet=True)

    def test_missingThermo(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'No thermo data'):
            convertMech(pjoin(self.test_data_dir, 'h2o2_missingThermo.inp'),
                        outName=pjoin(self.test_work_dir, 'h2o2_missingThermo.cti'),
                        quiet=True)

    def test_duplicate_thermo(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'additional thermo'):
            convertMech(pjoin(self.test_data_dir, 'duplicate-thermo.inp'),
                        outName=pjoin(self.test_work_dir, 'duplicate-thermo.cti'),
                        quiet=True)

        convertMech(pjoin(self.test_data_dir, 'duplicate-thermo.inp'),
                    outName=pjoin(self.test_work_dir, 'duplicate-thermo.cti'),
                    quiet=True, permissive=True)

        gas = ct.Solution('duplicate-thermo.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_duplicate_species(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'additional declaration'):
            convertMech(pjoin(self.test_data_dir, 'duplicate-species.inp'),
                        outName=pjoin(self.test_work_dir, 'duplicate-species.cti'),
                        quiet=True)

        convertMech(pjoin(self.test_data_dir, 'duplicate-species.inp'),
                    outName=pjoin(self.test_work_dir, 'duplicate-species.cti'),
                    quiet=True, permissive=True)

        gas = ct.Solution('duplicate-species.cti')
        self.assertEqual(gas.species_names, ['foo','bar','baz'])

    def test_pathologicalSpeciesNames(self):
        convertMech(pjoin(self.test_data_dir, 'species-names.inp'),
                    outName=pjoin(self.test_work_dir, 'species-names.cti'), quiet=True)
        gas = ct.Solution('species-names.cti')

        self.assertEqual(gas.n_species, 9)
        self.assertEqual(gas.species_name(0), '(Parens)')
        self.assertEqual(gas.species_name(1), '@#$%^-2')
        self.assertEqual(gas.species_index('co:lons:'), 2)
        self.assertEqual(gas.species_name(3), '[xy2]*{.}')
        self.assertEqual(gas.species_name(4), 'plus+')
        self.assertEqual(gas.species_name(5), 'eq=uals')
        self.assertEqual(gas.species_name(6), 'plus')
        self.assertEqual(gas.species_name(7), 'trans_butene')
        self.assertEqual(gas.species_name(8), 'co')

        self.assertEqual(gas.n_reactions, 12)
        nu = gas.product_stoich_coeffs() - gas.reactant_stoich_coeffs()
        self.assertEqual(list(nu[:,0]), [-1, -1, 0, 2, 0, 0, 0, 0, 0])
        self.assertEqual(list(nu[:,1]), [-2, 3, 0, -1, 0, 0, 0, 0, 0])
        self.assertEqual(list(nu[:,2]), [-1, 0, 0, 0, 1, 0, 0, 0, 0])
        self.assertEqual(list(nu[:,3]), [3, 0, 0, 0, -2, -1, 0, 0, 0])
        self.assertEqual(list(nu[:,4]), [2, 0, 0, 0, -1, 0, -1, 0, 0])
        self.assertEqual(list(nu[:,5]), [1, 0, 0, 0, 1, -1, -1, 0, 0])
        self.assertEqual(list(nu[:,6]), [2, 0, -1, 0, 0, -1, 0, 0, 0])
        self.assertEqual(list(nu[:,7]), [0, 0, 0, 0, -1, 1, 0, 0, 0])
        self.assertEqual(list(nu[:,8]), [0, 0, 0, 0, -1, 1, 0, 0, 0])
        self.assertEqual(list(nu[:,9]), [0, 0, 0, 0, -1, 1, 0, 0, 0])
        self.assertEqual(list(nu[:,10]), [0, 0, -1, 0, 2, 0, 0, -1, 0])
        self.assertEqual(list(nu[:,11]), [0, 0, -1, 0, 2, 0, 0, 0, -1])

    def test_unterminatedSections(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'implicitly ended'):
            convertMech(pjoin(self.test_data_dir, 'unterminated-sections.inp'),
                        outName=pjoin(self.test_work_dir, 'unterminated-sections.cti'),
                        quiet=True)

        convertMech(pjoin(self.test_data_dir, 'unterminated-sections.inp'),
                    outName=pjoin(self.test_work_dir, 'unterminated-sections.cti'),
                    quiet=True, permissive=True)

        gas = ct.Solution('unterminated-sections.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unterminatedSections2(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'implicitly ended'):
            convertMech(pjoin(self.test_data_dir, 'unterminated-sections2.inp'),
                        outName=pjoin(self.test_work_dir, 'unterminated-sections2.cti'),
                        quiet=True)

        convertMech(pjoin(self.test_data_dir, 'unterminated-sections2.inp'),
                    outName=pjoin(self.test_work_dir, 'unterminated-sections2.cti'),
                    quiet=True, permissive=True)

        gas = ct.Solution('unterminated-sections2.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unrecognized_section(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'SPAM'):
            convertMech(pjoin(self.test_data_dir, 'unrecognized-section.inp'),
                        thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                        outName=pjoin(self.test_work_dir, 'unrecognized-section.cti'),
                        quiet=True, permissive=True)

    def test_nasa9(self):
        convertMech(pjoin(self.test_data_dir, 'nasa9-test.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'nasa9-test-therm.dat'),
                    outName=pjoin(self.test_work_dir, 'nasa9_test.cti'), quiet=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'nasa9-test.xml'),
                                        'nasa9_test.cti')
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_nasa9_subset(self):
        convertMech(pjoin(self.test_data_dir, 'nasa9-test-subset.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'nasa9-test-therm.dat'),
                    outName=pjoin(self.test_work_dir, 'nasa9-test-subset.cti'),
                    quiet=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'nasa9-test-subset.xml'),
                                        'nasa9-test-subset.cti')
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_sri_falloff(self):
        convertMech(pjoin(self.test_data_dir, 'sri-falloff.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'sri-falloff.cti'), quiet=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'sri-falloff.xml'),
                                        'sri-falloff.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_chemically_activated(self):
        name = 'chemically-activated-reaction'
        convertMech(pjoin(self.test_data_dir, '{0}.inp'.format(name)),
                    outName=pjoin(self.test_work_dir, '{0}.cti'.format(name)), quiet=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, '{0}.xml'.format(name)),
                                        pjoin(self.test_work_dir, '{0}.cti'.format(name)))
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6, 1e7])

    def test_explicit_third_bodies(self):
        convertMech(pjoin(self.test_data_dir, 'explicit-third-bodies.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'explicit-third-bodies.cti'), quiet=True)

        ref, gas = self.checkConversion('explicit-third-bodies.cti',
                                        pjoin(self.test_data_dir, 'explicit-third-bodies.xml'))
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_explicit_reverse_rate(self):
        convertMech(pjoin(self.test_data_dir, 'explicit-reverse-rate.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'explicit-reverse-rate.cti'), quiet=True)
        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'explicit-reverse-rate.xml'),
                                        'explicit-reverse-rate.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

        # Reactions with explicit reverse rate constants are transformed into
        # two irreversible reactions with reactants and products swapped, unless
        # the explicit reverse rate is zero so only the forward reaction is used.
        Rr = gas.reverse_rate_constants
        self.assertEqual(Rr[0], 0.0)
        self.assertEqual(Rr[1], 0.0)
        self.assertEqual(Rr[2], 0.0)
        self.assertEqual(Rr[3], 0.0)
        self.assertEqual(Rr[4], 0.0)
        Rstoich = gas.reactant_stoich_coeffs()
        Pstoich = gas.product_stoich_coeffs()
        self.assertEqual(list(Rstoich[:, 0]), list(Pstoich[:, 1]))
        self.assertEqual(list(Rstoich[:, 1]), list(Pstoich[:, 0]))
        self.assertEqual(list(Rstoich[:, 2]), list(Pstoich[:, 3]))
        self.assertEqual(list(Rstoich[:, 3]), list(Pstoich[:, 2]))

        self.assertEqual(gas.n_reactions, 5)

    def test_explicit_forward_order(self):
        convertMech(pjoin(self.test_data_dir, 'explicit-forward-order.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'explicit-forward-order.cti'), quiet=True)
        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'explicit-forward-order.xml'),
                                        'explicit-forward-order.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_negative_order(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'Negative reaction order'):
            convertMech(pjoin(self.test_data_dir, 'negative-order.inp'),
                        thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                        outName=pjoin(self.test_work_dir, 'negative-order.cti'), quiet=True)

    def test_negative_order_permissive(self):
        convertMech(pjoin(self.test_data_dir, 'negative-order.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'negative-order.cti'),
                    quiet=True, permissive=True)
        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'explicit-forward-order.xml'),
                                        'explicit-forward-order.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_bad_troe_value(self):
        with self.assertRaises(ValueError):
            convertMech(pjoin(self.test_data_dir, 'bad-troe.inp'),
                        thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                        outName=pjoin(self.test_work_dir, 'bad-troe.cti'), quiet=True)

    def test_invalid_reaction_equation(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'Unparsable'):
            convertMech(pjoin(self.test_data_dir, 'invalid-equation.inp'),
                        thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                        outName=pjoin(self.test_work_dir, 'invalid-equation.cti'), quiet=True)

    def test_reaction_units(self):
        convertMech(pjoin(self.test_data_dir, 'units-default.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'units-default.cti'), quiet=True)
        convertMech(pjoin(self.test_data_dir, 'units-custom.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'units-custom.cti'), quiet=True)

        default, custom = self.checkConversion('units-default.cti',
                                               'units-custom.cti')
        self.checkKinetics(default, custom,
                           [300, 800, 1450, 2800], [5e3, 1e5, 2e6], 1e-7)

    def test_float_stoich_coeffs(self):
        convertMech(pjoin(self.test_data_dir, 'float-stoich.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'float-stoich.cti'), quiet=True)
        gas = ct.Solution('float-stoich.cti')

        R = gas.reactant_stoich_coeffs()
        P = gas.product_stoich_coeffs()
        self.assertArrayNear(R[:,0], [0, 1.5, 0.5, 0])
        self.assertArrayNear(P[:,0], [1, 0, 0, 1])
        self.assertArrayNear(R[:,1], [1, 0, 0, 1])
        self.assertArrayNear(P[:,1], [0, 0.33, 1.67, 0])

    def test_photon(self):
        convertMech(pjoin(self.test_data_dir, 'photo-reaction.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'photo-reaction.cti'), quiet=True,
                    permissive=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'photo-reaction.xml'),
                                        'photo-reaction.cti')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_transport_normal(self):
        convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                    transportFile=pjoin(self.test_data_dir, 'gri30_tran.dat'),
                    outName=pjoin(self.test_work_dir, 'h2o2_transport_normal.cti'), quiet=True)

        gas = ct.Solution('h2o2_transport_normal.cti')
        gas.TPX = 300, 101325, 'H2:1.0, O2:1.0'
        self.assertAlmostEqual(gas.thermal_conductivity, 0.07663, 4)

    def test_transport_embedded(self):
        convertMech(pjoin(self.test_data_dir, 'with-transport.inp'),
                    outName=pjoin(self.test_work_dir, 'with-transport.cti'), quiet=True)

        gas = ct.Solution('with-transport.cti')
        gas.X = [0.2, 0.3, 0.5]
        D = gas.mix_diff_coeffs
        for d in D:
            self.assertTrue(d > 0.0)

    def test_transport_missing_species(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'No transport data'):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-missing-species-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_transport_missing_species.cti'),
                        quiet=True)

    def test_transport_extra_column_entries(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, '572.400'):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-extra-column-entries-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_extra-column-entries-tran.cti'),
                        quiet=True)

    def test_transport_duplicate_species(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'duplicate transport'):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-duplicate-species-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_transport_duplicate_species.cti'),
                        quiet=True)

        convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                    transportFile=pjoin(self.test_data_dir, 'h2o2-duplicate-species-tran.dat'),
                    outName=pjoin(self.test_work_dir, 'h2o2_transport_duplicate_species.cti'),
                    quiet=True,
                    permissive=True)

    def test_transport_bad_geometry(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'geometry flag'):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-bad-geometry-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_transport_bad_geometry.cti'),
                        quiet=True)

    def test_transport_float_geometry(self):
        with self.assertRaisesRegex(ck2cti.InputParseError, 'geometry flag'):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-float-geometry-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_transport_float_geometry.cti'),
                        quiet=True)

    def test_empty_reaction_section(self):
        convertMech(pjoin(self.test_data_dir, 'h2o2_emptyReactions.inp'),
                    outName=pjoin(self.test_work_dir, 'h2o2_emptyReactions.cti'),
                    quiet=True)

        gas = ct.Solution('h2o2_emptyReactions.cti')
        self.assertEqual(gas.n_species, 9)
        self.assertEqual(gas.n_reactions, 0)

    def test_reaction_comments1(self):
        convertMech(pjoin(self.test_data_dir, 'pdep-test.inp'),
                    outName=pjoin(self.test_work_dir, 'pdep_test.cti'), quiet=True)
        with open(pjoin(self.test_work_dir, 'pdep_test.cti')) as f:
            text = f.read()
        self.assertIn('Generic mechanism header', text)
        self.assertIn('Single PLOG reaction', text)
        self.assertIn('Multiple PLOG expressions at the same pressure', text)

    def test_reaction_comments2(self):
        convertMech(pjoin(self.test_data_dir, 'explicit-third-bodies.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'dummy-thermo.dat'),
                    outName=pjoin(self.test_work_dir, 'explicit_third_bodies.cti'), quiet=True)
        with open(pjoin(self.test_work_dir, 'explicit_third_bodies.cti')) as f:
            text = f.read()
        self.assertIn('An end of line comment', text)
        self.assertIn('A comment after the last reaction', text)

    def test_custom_element(self):
        convertMech(pjoin(self.test_data_dir, 'custom-elements.inp'),
                    outName=pjoin(self.test_work_dir, 'custom-elements.cti'), quiet=True)
        gas = ct.Solution('custom-elements.cti')
        self.assertEqual(gas.n_elements, 4)
        self.assertNear(gas.atomic_weight(2), 13.003)
        self.assertEqual(gas.n_atoms('ethane', 'C'), 2)
        self.assertEqual(gas.n_atoms('CC', 'C'), 1)
        self.assertEqual(gas.n_atoms('CC', 'Ci'), 1)

    def test_surface_mech(self):
        convertMech(pjoin(self.test_data_dir, 'surface1-gas.inp'),
                    surfaceFile=pjoin(self.test_data_dir, 'surface1.inp'),
                    outName=pjoin(self.test_work_dir, 'surface1.cti'), quiet=True)

        gas = ct.Solution('surface1.cti', 'gas')
        surf = ct.Interface('surface1.cti', 'PT_SURFACE', [gas])

        self.assertEqual(gas.n_reactions, 11)
        self.assertEqual(surf.n_reactions, 15)
        self.assertEqual(surf.species('O2_Pt').size, 3)

        # Different units for rate constants in each input file
        # 62.1 kJ/gmol = 6.21e7 J/kmol
        self.assertNear(gas.reaction(0).rate.activation_energy, 6.21e7)
        # 67400 J/mol = 6.74e7 J/kmol
        self.assertNear(surf.reaction(1).rate.activation_energy, 6.74e7)

        # Sticking coefficients
        self.assertFalse(surf.reaction(1).is_sticking_coefficient)
        self.assertTrue(surf.reaction(2).is_sticking_coefficient)
        self.assertTrue(surf.reaction(2).use_motz_wise_correction)
        self.assertTrue(surf.reaction(4).is_sticking_coefficient)
        self.assertFalse(surf.reaction(4).use_motz_wise_correction)
        self.assertTrue(surf.reaction(4).duplicate)
        self.assertTrue(surf.reaction(6).use_motz_wise_correction)

        # Coverage dependencies
        covdeps = surf.reaction(1).coverage_deps
        self.assertEqual(len(covdeps), 2)
        self.assertIn('H_Pt', covdeps)
        self.assertEqual(covdeps['OH_Pt'][1], 1.0)
        self.assertNear(covdeps['H_Pt'][2], -6e6) # 6000 J/gmol = 6e6 J/kmol

    def test_surface_mech2(self):
        convertMech(pjoin(self.test_data_dir, 'surface1-gas-noreac.inp'),
                    surfaceFile=pjoin(self.test_data_dir, 'surface1.inp'),
                    outName=pjoin(self.test_work_dir, 'surface1-nogasreac.cti'), quiet=True)

        gas = ct.Solution('surface1-nogasreac.cti', 'gas')
        surf = ct.Interface('surface1-nogasreac.cti', 'PT_SURFACE', [gas])

        self.assertEqual(gas.n_reactions, 0)
        self.assertEqual(surf.n_reactions, 15)


class CtmlConverterTest(utilities.CanteraTest):
    def test_sofc(self):
        gas_a, anode_bulk, oxide_a = ct.import_phases(
            'sofc.cti',
            ['gas', 'metal', 'oxide_bulk'])

        self.assertNear(gas_a.P, ct.one_atm)
        self.assertNear(anode_bulk['electron'].X, 1.0)
        self.assertNear(oxide_a.density_mole, 17.6)

    def test_diamond(self):
        gas, solid = ct.import_phases('diamond.cti', ['gas','diamond'])
        face = ct.Interface('diamond.cti', 'diamond_100', [gas, solid])

        self.assertNear(face.site_density, 3e-8)

    def test_pdep(self):
        gas = ct.Solution('pdep-test.cti')
        self.assertEqual(gas.n_reactions, 7)

    def test_invalid(self):
        try:
            gas = ct.Solution('invalid.cti')
        except ct.CanteraError as e:
            err = e

        self.assertIn('already contains', err.args[0])

    def test_noninteger_atomicity(self):
        gas = ct.Solution('noninteger-atomicity.cti')
        self.assertNear(gas.molecular_weights[gas.species_index('CnHm')],
                        10.65*gas.atomic_weight('C') + 21.8*gas.atomic_weight('H'))

    def test_reaction_orders(self):
        gas = ct.Solution('reaction-orders.cti')
        self.assertEqual(gas.n_reactions, 1)
        R = gas.reaction(0)
        self.assertTrue(R.allow_nonreactant_orders)
        self.assertNear(R.orders.get('OH'), 0.15)
        self.assertTrue(R.allow_negative_orders)
        self.assertNear(R.orders.get('H2'), -0.25)

    def test_long_source_input(self):
        """
        Here we are testing if passing a very long string will result in a
        Solution object. This should result in a temp file creation in most OS's
        """

        gas = ct.Solution(pjoin(self.test_data_dir, 'pdep-test.cti'))

        with open(pjoin(self.test_data_dir, 'pdep-test.cti'), 'r') as f:
            data = f.read()
        data_size_2048kB = data + ' '*2048*1024
        gas2 = ct.Solution(source=data_size_2048kB)

        self.assertEqual(gas.n_reactions, gas2.n_reactions)

    def test_short_source_input(self):
        """
        Here we are testing if passing a short string will result in a Solution
        object. This should not result in a temp file creation in most OS's
        """

        gas = ct.Solution(pjoin(self.test_data_dir, 'pdep-test.cti'))

        with open(pjoin(self.test_data_dir, 'pdep-test.cti'), 'r') as f:
            data = f.read()
        data_size_32kB = data + ' '*18000
        gas2 = ct.Solution(source=data_size_32kB)

        self.assertEqual(gas.n_reactions, gas2.n_reactions)


class cti2yamlTest(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cti2yaml.convert(pjoin(cls.cantera_data, 'gri30.cti'), 'gri30.yaml')

    def checkConversion(self, basename, cls=ct.Solution, ctiphases=(),
                        yamlphases=(), **kwargs):
        ctiPhase = cls(basename + '.cti', phases=ctiphases, **kwargs)
        yamlPhase = cls(basename + '.yaml', phases=yamlphases, **kwargs)

        self.assertEqual(ctiPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctiPhase.species_names, yamlPhase.species_names)
        self.assertEqual(ctiPhase.n_reactions, yamlPhase.n_reactions)
        for C, Y in zip(ctiPhase.species(), yamlPhase.species()):
            self.assertEqual(C.composition, Y.composition)

        for i, (C, Y) in enumerate(zip(ctiPhase.reactions(),
                                       yamlPhase.reactions())):
            self.assertEqual(C.__class__, Y.__class__)
            self.assertEqual(C.reactants, Y.reactants)
            self.assertEqual(C.products, Y.products)
            self.assertEqual(C.duplicate, Y.duplicate)

        for i, sp in zip(range(ctiPhase.n_reactions), ctiPhase.kinetics_species_names):
            self.assertEqual(ctiPhase.reactant_stoich_coeff(sp, i),
                             yamlPhase.reactant_stoich_coeff(sp, i))

        return ctiPhase, yamlPhase

    def checkThermo(self, ctiPhase, yamlPhase, temperatures, tol=1e-7):
        for T in temperatures:
            ctiPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            cp_cti = ctiPhase.partial_molar_cp
            cp_yaml = yamlPhase.partial_molar_cp
            h_cti = ctiPhase.partial_molar_enthalpies
            h_yaml = yamlPhase.partial_molar_enthalpies
            s_cti = ctiPhase.partial_molar_entropies
            s_yaml = yamlPhase.partial_molar_entropies
            self.assertNear(ctiPhase.density, yamlPhase.density)
            for i in range(ctiPhase.n_species):
                message = ' for species {0} at T = {1}'.format(i, T)
                self.assertNear(cp_cti[i], cp_yaml[i], tol, msg='cp'+message)
                self.assertNear(h_cti[i], h_yaml[i], tol, msg='h'+message)
                self.assertNear(s_cti[i], s_yaml[i], tol, msg='s'+message)

    def checkKinetics(self, ctiPhase, yamlPhase, temperatures, pressures, tol=1e-7):
        for T,P in itertools.product(temperatures, pressures):
            ctiPhase.TP = T, P
            yamlPhase.TP = T, P
            kf_cti = ctiPhase.forward_rate_constants
            kr_cti = ctiPhase.reverse_rate_constants
            kf_yaml = yamlPhase.forward_rate_constants
            kr_yaml = yamlPhase.reverse_rate_constants
            for i in range(yamlPhase.n_reactions):
                message = ' for reaction {0} at T = {1}, P = {2}'.format(i, T, P)
                self.assertNear(kf_cti[i], kf_yaml[i], rtol=tol, msg='kf '+message)
                self.assertNear(kr_cti[i], kr_yaml[i], rtol=tol, msg='kr '+message)

    def checkTransport(self, ctiPhase, yamlPhase, temperatures,
                       model='mixture-averaged'):
        ctiPhase.transport_model = model
        yamlPhase.transport_model = model
        for T in temperatures:
            ctiPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            self.assertNear(ctiPhase.viscosity, yamlPhase.viscosity)
            self.assertNear(ctiPhase.thermal_conductivity,
                            yamlPhase.thermal_conductivity)
            Dkm_cti = ctiPhase.mix_diff_coeffs
            Dkm_yaml = yamlPhase.mix_diff_coeffs
            for i in range(ctiPhase.n_species):
                message = 'dkm for species {0} at T = {1}'.format(i, T)
                self.assertNear(Dkm_cti[i], Dkm_yaml[i], msg=message)

    def test_gri30(self):
        ctiPhase, yamlPhase = self.checkConversion('gri30')
        X = {'O2': 0.3, 'H2': 0.1, 'CH4': 0.2, 'CO2': 0.4}
        ctiPhase.X = X
        yamlPhase.X = X
        self.checkThermo(ctiPhase, yamlPhase, [300, 500, 1300, 2000])
        self.checkKinetics(ctiPhase, yamlPhase, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctiPhase, yamlPhase, [298, 1001, 2400])

    def test_pdep(self):
        cti2yaml.convert(pjoin(self.test_data_dir, 'pdep-test.cti'),
                         'pdep-test.yaml')
        ctiPhase, yamlPhase = self.checkConversion('pdep-test')
        self.checkKinetics(ctiPhase, yamlPhase, [300, 1000, 2200],
                           [100, ct.one_atm, 2e5, 2e6, 9.9e6])

    def test_ptcombust(self):
        cti2yaml.convert(pjoin(self.cantera_data, 'ptcombust.cti'),
                         'ptcombust.yaml')
        ctiGas, yamlGas = self.checkConversion('ptcombust')
        ctiSurf, yamlSurf = self.checkConversion('ptcombust', ct.Interface,
            phaseid='Pt_surf', ctiphases=[ctiGas], yamlphases=[yamlGas])

        self.checkKinetics(ctiGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctiSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctiSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    def test_sofc(self):
        cti2yaml.convert(pjoin(self.cantera_data, 'sofc.cti'), 'sofc.yaml')
        ctiGas, yamlGas = self.checkConversion('sofc')
        ctiMetal, yamlMetal = self.checkConversion('sofc', phaseid='metal')
        ctiOxide, yamlOxide = self.checkConversion('sofc', phaseid='oxide_bulk')
        ctiMSurf, yamlMSurf = self.checkConversion('sofc', ct.Interface,
            phaseid='metal_surface', ctiphases=[ctiGas, ctiMetal],
            yamlphases=[yamlGas, yamlMetal])
        ctiOSurf, yamlOSurf = self.checkConversion('sofc', ct.Interface,
            phaseid='oxide_surface', ctiphases=[ctiGas, ctiOxide],
            yamlphases=[yamlGas, yamlOxide])
        cti_tpb, yaml_tpb = self.checkConversion('sofc', ct.Interface,
            phaseid='tpb', ctiphases=[ctiMetal, ctiMSurf, ctiOSurf],
            yamlphases=[yamlMetal, yamlMSurf, yamlOSurf])

        self.checkThermo(ctiMSurf, yamlMSurf, [900, 1000, 1100])
        self.checkThermo(ctiOSurf, yamlOSurf, [900, 1000, 1100])
        ctiMetal.electric_potential = yamlMetal.electric_potential = 2
        self.checkKinetics(cti_tpb, yaml_tpb, [900, 1000, 1100], [1e5])
        ctiMetal.electric_potential = yamlMetal.electric_potential = 4
        self.checkKinetics(cti_tpb, yaml_tpb, [900, 1000, 1100], [1e5])

    def test_liquidvapor(self):
        cti2yaml.convert(pjoin(self.cantera_data, 'liquidvapor.cti'),
                         'liquidvapor.yaml')
        for name in ['water', 'nitrogen', 'methane', 'hydrogen', 'oxygen',
                     'hfc134a', 'carbondioxide', 'heptane']:
            ctiPhase, yamlPhase = self.checkConversion('liquidvapor', phaseid=name)
            self.checkThermo(ctiPhase, yamlPhase,
                             [1.3 * ctiPhase.min_temp, 0.7 * ctiPhase.max_temp])

    def test_Redlich_Kwong_CO2(self):
        cti2yaml.convert(pjoin(self.test_data_dir, 'co2_RK_example.cti'),
                         'co2_RK_example.yaml')
        ctiGas, yamlGas = self.checkConversion('co2_RK_example')
        for P in [1e5, 2e6, 1.3e7]:
            yamlGas.TP = ctiGas.TP = 300, P
            self.checkThermo(ctiGas, yamlGas, [300, 400, 500])

    def test_Redlich_Kwong_ndodecane(self):
        cti2yaml.convert(pjoin(self.cantera_data, 'nDodecane_Reitz.cti'),
                         'nDodecane_Reitz.yaml')
        ctiGas, yamlGas = self.checkConversion('nDodecane_Reitz')
        self.checkThermo(ctiGas, yamlGas, [300, 400, 500])
        self.checkKinetics(ctiGas, yamlGas, [300, 500, 1300], [1e5, 2e6, 1.4e7],
                           1e-6)

    def test_diamond(self):
        cti2yaml.convert(pjoin(self.cantera_data, 'diamond.cti'), 'diamond.yaml')
        ctiGas, yamlGas = self.checkConversion('diamond', phaseid='gas')
        ctiSolid, yamlSolid = self.checkConversion('diamond', phaseid='diamond')
        ctiSurf, yamlSurf = self.checkConversion('diamond',
            ct.Interface, phaseid='diamond_100', ctiphases=[ctiGas, ctiSolid],
            yamlphases=[yamlGas, yamlSolid])
        self.checkThermo(ctiSolid, yamlSolid, [300, 500])
        self.checkThermo(ctiSurf, yamlSurf, [330, 490])
        self.checkKinetics(ctiSurf, yamlSurf, [400, 800], [2e5])

    def test_lithium_ion_battery(self):
        cti2yaml.convert(pjoin(self.cantera_data, 'lithium_ion_battery.cti'),
                         'lithium_ion_battery.yaml')
        name = 'lithium_ion_battery'
        ctiAnode, yamlAnode = self.checkConversion(name, phaseid='anode')
        ctiCathode, yamlCathode = self.checkConversion(name, phaseid='cathode')
        ctiMetal, yamlMetal = self.checkConversion(name, phaseid='electron')
        ctiElyt, yamlElyt = self.checkConversion(name, phaseid='electrolyte')
        ctiAnodeInt, yamlAnodeInt = self.checkConversion(name,
            phaseid='edge_anode_electrolyte',
            ctiphases=[ctiAnode, ctiMetal, ctiElyt],
            yamlphases=[yamlAnode, yamlMetal, yamlElyt])
        ctiCathodeInt, yamlCathodeInt = self.checkConversion(name,
            phaseid='edge_cathode_electrolyte',
            ctiphases=[ctiCathode, ctiMetal, ctiElyt],
            yamlphases=[yamlCathode, yamlMetal, yamlElyt])

        self.checkThermo(ctiAnode, yamlAnode, [300, 330])
        self.checkThermo(ctiCathode, yamlCathode, [300, 330])

        ctiAnode.X = yamlAnode.X = [0.7, 0.3]
        self.checkThermo(ctiAnode, yamlAnode, [300, 330])
        ctiCathode.X = yamlCathode.X = [0.2, 0.8]
        self.checkThermo(ctiCathode, yamlCathode, [300, 330])

        for phase in [ctiAnode, yamlAnode, ctiCathode, yamlCathode, ctiMetal,
                      yamlMetal, ctiElyt, yamlElyt, ctiAnodeInt, yamlAnodeInt,
                      ctiCathodeInt, yamlCathodeInt]:
            phase.TP = 300, 1e5
        ctiMetal.electric_potential = yamlMetal.electric_potential = 0
        ctiElyt.electric_potential = yamlElyt.electric_potential = 1.9
        self.checkKinetics(ctiAnodeInt, yamlAnodeInt, [300], [1e5])

        ctiMetal.electric_potential = yamlMetal.electric_potential = 2.2
        ctiElyt.electric_potential = yamlElyt.electric_potential = 0
        self.checkKinetics(ctiCathodeInt, yamlCathodeInt, [300], [1e5])
