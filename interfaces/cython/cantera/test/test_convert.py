import os
from os.path import join as pjoin
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
        with self.assertRaises(ck2cti.InputParseError):
            convertMech(pjoin(self.test_data_dir, 'h2o2_missingElement.inp'),
                        outName=pjoin(self.test_work_dir, 'h2o2_missingElement.cti'),
                        quiet=True)

    def test_missingThermo(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech(pjoin(self.test_data_dir, 'h2o2_missingThermo.inp'),
                        outName=pjoin(self.test_work_dir, 'h2o2_missingThermo.cti'),
                        quiet=True)

    def test_duplicate_thermo(self):
        with self.assertRaises(ck2cti.InputParseError):
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
        with self.assertRaises(ck2cti.InputParseError):
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

        self.assertEqual(gas.n_species, 7)
        self.assertEqual(gas.species_name(0), '(Parens)')
        self.assertEqual(gas.species_name(1), '@#$%^-2')
        self.assertEqual(gas.species_index('co:lons:'), 2)
        self.assertEqual(gas.species_name(3), '[xy2]*{.}')
        self.assertEqual(gas.species_name(4), 'plus+')
        self.assertEqual(gas.species_name(5), 'eq=uals')
        self.assertEqual(gas.species_name(6), 'plus')

        self.assertEqual(gas.n_reactions, 7)
        nu = gas.product_stoich_coeffs() - gas.reactant_stoich_coeffs()
        self.assertEqual(list(nu[:,0]), [-1, -1, 0, 2, 0, 0, 0])
        self.assertEqual(list(nu[:,1]), [-2, 3, 0, -1, 0, 0, 0])
        self.assertEqual(list(nu[:,2]), [-1, 0, 0, 0, 1, 0, 0])
        self.assertEqual(list(nu[:,3]), [3, 0, 0, 0, -2, -1, 0])
        self.assertEqual(list(nu[:,4]), [2, 0, 0, 0, -1, 0, -1])
        self.assertEqual(list(nu[:,5]), [1, 0, 0, 0, 1, -1, -1])
        self.assertEqual(list(nu[:,6]), [2, 0, -1, 0, 0, -1, 0])

    def test_unterminatedSections(self):
        with self.assertRaises(ck2cti.InputParseError):
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
        with self.assertRaises(ck2cti.InputParseError):
            convertMech(pjoin(self.test_data_dir, 'unterminated-sections2.inp'),
                        outName=pjoin(self.test_work_dir, 'unterminated-sections2.cti'),
                        quiet=True)

        convertMech(pjoin(self.test_data_dir, 'unterminated-sections2.inp'),
                    outName=pjoin(self.test_work_dir, 'unterminated-sections2.cti'),
                    quiet=True, permissive=True)

        gas = ct.Solution('unterminated-sections2.cti')
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_nasa9(self):
        convertMech(pjoin(self.test_data_dir, 'nasa9-test.inp'),
                    thermoFile=pjoin(self.test_data_dir, 'nasa9-test-therm.dat'),
                    outName=pjoin(self.test_work_dir, 'nasa9_test.cti'), quiet=True)

        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'nasa9-test.xml'),
                                        'nasa9_test.cti')
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
        with self.assertRaises(ck2cti.InputParseError):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-missing-species-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_transport_missing_species.cti'),
                        quiet=True)

    def test_transport_extra_column_entries(self):
        with self.assertRaises(ck2cti.InputParseError):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-extra-column-entries-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_extra-column-entries-tran.cti'),
                        quiet=True)

    def test_transport_duplicate_species(self):
        with self.assertRaises(ck2cti.InputParseError):
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
        with self.assertRaises(ck2cti.InputParseError):
            convertMech(pjoin(self.test_data_dir, 'h2o2.inp'),
                        transportFile=pjoin(self.test_data_dir, 'h2o2-bad-geometry-tran.dat'),
                        outName=pjoin(self.test_work_dir, 'h2o2_transport_bad_geometry.cti'),
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
        self.assertIn('PLOG with duplicate rates and negative A-factors', text)

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
        self.assertEqual(surf.n_reactions, 14)

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


class CtmlConverterTest(utilities.CanteraTest):
    def test_sofc(self):
        gas_a, anode_bulk, oxide_a = ct.import_phases(
            'sofc.cti',
            ['gas', 'metal', 'oxide_bulk'])

        self.assertNear(gas_a.P, ct.one_atm)
        self.assertNear(anode_bulk['electron'].X, 1.0)
        self.assertNear(oxide_a.density, 700)

    def test_diamond(self):
        gas, solid = ct.import_phases('diamond.cti', ['gas','diamond'])
        face = ct.Interface('diamond.cti', 'diamond_100', [gas, solid])

        self.assertNear(face.site_density, 3e-8)

    def test_pdep(self):
        gas = ct.Solution('pdep-test.cti')
        self.assertEqual(gas.n_reactions, 6)

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
        R = gas.reaction(0)
        self.assertTrue(R.allow_nonreactant_orders)
        self.assertNear(R.orders.get('OH'), 0.15)
        self.assertTrue(R.allow_negative_orders)
        self.assertNear(R.orders.get('H2'), -0.25)
