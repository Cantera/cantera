import os
from os.path import join as pjoin
import itertools
from pathlib import Path

from . import utilities
import cantera as ct
from cantera import ck2cti, ck2yaml, cti2yaml

class converterTestCommon:
    def convert(self, inputFile, thermo=None, transport=None,
                surface=None, output=None, quiet=True, permissive=None):
        if output is None:
            output = inputFile[:-4]  # strip '.inp'
        if inputFile is not None:
            inputFile = pjoin(self.test_data_dir, inputFile)
        if thermo is not None:
            thermo = pjoin(self.test_data_dir, thermo)
        if transport is not None:
            transport = pjoin(self.test_data_dir, transport)
        if surface is not None:
            surface = pjoin(self.test_data_dir, surface)
        output = pjoin(self.test_work_dir, output + self.ext)
        if os.path.exists(output):
            os.remove(output)
        self._convert(inputFile, thermo=thermo, transport=transport,
            surface=surface, output=output, quiet=quiet, permissive=permissive)

    def checkConversion(self, refFile, testFile):
        ref = ct.Solution(refFile)
        gas = ct.Solution(testFile + self.ext)

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
                self.assertNear(ref_kf[i], gas_kf[i], rtol=tol, msg='kf' + message)
                self.assertNear(ref_kr[i], gas_kr[i], rtol=tol, msg='kr' + message)

    def test_gri30(self):
        self.convert('gri30.inp', thermo='gri30_thermo.dat',
                     transport='gri30_tran.dat', output='gri30_test')

        ref, gas = self.checkConversion('gri30.xml', 'gri30_test')
        self.checkKinetics(ref, gas, [300, 1500], [5e3, 1e5, 2e6])

    def test_soot(self):
        self.convert('soot.inp', thermo='soot-therm.dat', output='soot_test')
        ref, gas = self.checkConversion('soot.xml', 'soot_test')
        self.checkThermo(ref, gas, [300, 1100])
        self.checkKinetics(ref, gas, [300, 1100], [5e3, 1e5, 2e6])

    def test_pdep(self):
        self.convert('pdep-test.inp')
        ref, gas = self.checkConversion(pjoin(self.test_data_dir, 'pdep-test.xml'),
                                        pjoin(self.test_work_dir, 'pdep-test'))
        # Chebyshev coefficients in XML are truncated to 6 digits, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6],
                           tol=2e-4)

    def test_species_only(self):
        self.convert(None, thermo='dummy-thermo.dat', output='dummy-thermo')
        cti = "ideal_gas(elements='C H', species='dummy-thermo:R1A R1B P1')"
        gas = ct.Solution(source=cti)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 0)

    def test_missingThermo(self):
        with self.assertRaisesRegex(self.InputError, 'No thermo data'):
            self.convert('h2o2_missingThermo.inp')

    def test_duplicate_thermo(self):
        with self.assertRaisesRegex(self.InputError, 'additional thermo'):
            self.convert('duplicate-thermo.inp')

        self.convert('duplicate-thermo.inp', permissive=True)

        gas = ct.Solution('duplicate-thermo' + self.ext)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_duplicate_species(self):
        with self.assertRaisesRegex(self.InputError, 'additional declaration'):
            self.convert('duplicate-species.inp')

        self.convert('duplicate-species.inp', permissive=True)

        gas = ct.Solution('duplicate-species' + self.ext)
        self.assertEqual(gas.species_names, ['foo','bar','baz'])

    def test_pathologicalSpeciesNames(self):
        self.convert('species-names.inp')
        gas = ct.Solution('species-names' + self.ext)

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
        with self.assertRaisesRegex(self.InputError, 'implicitly ended'):
            self.convert('unterminated-sections.inp')

        self.convert('unterminated-sections.inp', permissive=True)
        gas = ct.Solution('unterminated-sections' + self.ext)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unterminatedSections2(self):
        with self.assertRaisesRegex(self.InputError, 'implicitly ended'):
            self.convert('unterminated-sections2.inp')

        self.convert('unterminated-sections2.inp', permissive=True)
        gas = ct.Solution('unterminated-sections2' + self.ext)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unrecognized_section(self):
        with self.assertRaisesRegex(self.InputError, 'SPAM'):
            self.convert('unrecognized-section.inp', thermo='dummy-thermo.dat',
                         permissive=True)

    def test_nasa9(self):
        self.convert('nasa9-test.inp',thermo='nasa9-test-therm.dat')
        ref, gas = self.checkConversion('nasa9-test.xml', 'nasa9-test')
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_nasa9_subset(self):
        self.convert('nasa9-test-subset.inp', thermo='nasa9-test-therm.dat')
        ref, gas = self.checkConversion('nasa9-test-subset.xml',
                                        'nasa9-test-subset')
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_sri_falloff(self):
        self.convert('sri-falloff.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('sri-falloff.xml', 'sri-falloff')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_chemically_activated(self):
        self.convert('chemically-activated-reaction.inp')
        ref, gas = self.checkConversion('chemically-activated-reaction.xml',
                                        'chemically-activated-reaction')
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6, 1e7],
                           tol=1e-7)

    def test_explicit_third_bodies(self):
        self.convert('explicit-third-bodies.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('explicit-third-bodies.xml',
                                        'explicit-third-bodies')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_explicit_reverse_rate(self):
        self.convert('explicit-reverse-rate.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('explicit-reverse-rate.xml',
                                        'explicit-reverse-rate')
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
        self.convert('explicit-forward-order.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('explicit-forward-order.xml',
                                        'explicit-forward-order')
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6],
                           tol=2e-7)

    def test_negative_order(self):
        with self.assertRaisesRegex(self.InputError, 'Negative reaction order'):
            self.convert('negative-order.inp', thermo='dummy-thermo.dat')

    def test_negative_order_permissive(self):
        self.convert('negative-order.inp', thermo='dummy-thermo.dat',
            permissive=True)
        ref, gas = self.checkConversion('explicit-forward-order.xml',
                                        'explicit-forward-order')
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6],
                           tol=2e-7)

    def test_bad_troe_value(self):
        with self.assertRaises(ValueError):
            self.convert('bad-troe.inp', thermo='dummy-thermo.dat')

    def test_invalid_reaction_equation(self):
        with self.assertRaisesRegex(self.InputError, 'Unparsable'):
            self.convert('invalid-equation.inp', thermo='dummy-thermo.dat')

    def test_reaction_units(self):
        self.convert('units-default.inp', thermo='dummy-thermo.dat')
        self.convert('units-custom.inp', thermo='dummy-thermo.dat')
        default, custom = self.checkConversion('units-default' + self.ext,
                                               'units-custom')
        self.checkKinetics(default, custom,
                           [300, 800, 1450, 2800], [5e0, 5e3, 1e5, 2e6, 1e8], 1e-7)

    def test_float_stoich_coeffs(self):
        self.convert('float-stoich.inp', thermo='dummy-thermo.dat')
        gas = ct.Solution('float-stoich' + self.ext)

        R = gas.reactant_stoich_coeffs()
        P = gas.product_stoich_coeffs()
        self.assertArrayNear(R[:,0], [0, 1.5, 0.5, 0])
        self.assertArrayNear(P[:,0], [1, 0, 0, 1])
        self.assertArrayNear(R[:,1], [1, 0, 0, 1])
        self.assertArrayNear(P[:,1], [0, 0.33, 1.67, 0])

    def test_photon(self):
        self.convert('photo-reaction.inp', thermo='dummy-thermo.dat',
                     permissive=True)

        ref, gas = self.checkConversion('photo-reaction.xml', 'photo-reaction')
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_transport_normal(self):
        self.convert('h2o2.inp', transport='gri30_tran.dat',
            output='h2o2_transport_normal')

        gas = ct.Solution('h2o2_transport_normal' + self.ext)
        gas.TPX = 300, 101325, 'H2:1.0, O2:1.0'
        self.assertAlmostEqual(gas.thermal_conductivity, 0.07663, 4)

    def test_transport_embedded(self):
        self.convert('with-transport.inp')
        gas = ct.Solution('with-transport' + self.ext)
        gas.X = [0.2, 0.3, 0.5]
        D = gas.mix_diff_coeffs
        for d in D:
            self.assertTrue(d > 0.0)

    def test_transport_missing_species(self):
        with self.assertRaisesRegex(self.InputError, 'No transport data'):
            self.convert('h2o2.inp', transport='h2o2-missing-species-tran.dat',
                output='h2o2_transport_missing_species')

    def test_transport_extra_column_entries(self):
        with self.assertRaisesRegex(self.InputError, '572.400'):
            self.convert('h2o2.inp',
                transport='h2o2-extra-column-entries-tran.dat',
                output='h2o2_extra-column-entries-tran')

    def test_transport_duplicate_species(self):
        with self.assertRaisesRegex(self.InputError, 'duplicate transport'):
            self.convert('h2o2.inp',
                transport='h2o2-duplicate-species-tran.dat',
                output='h2o2_transport_duplicate_species.cti')

        self.convert('h2o2.inp',
            transport='h2o2-duplicate-species-tran.dat',
            output='h2o2_transport_duplicate_species', permissive=True)

    def test_transport_bad_geometry(self):
        with self.assertRaisesRegex(self.InputError, 'geometry flag'):
            self.convert('h2o2.inp',
                transport='h2o2-bad-geometry-tran.dat',
                output='h2o2_transport_bad_geometry')

    def test_transport_float_geometry(self):
        with self.assertRaisesRegex(self.InputError, 'geometry flag'):
            self.convert('h2o2.inp',
                transport='h2o2-float-geometry-tran.dat',
                output='h2o2_transport_float_geometry')

    def test_empty_reaction_section(self):
        self.convert('h2o2_emptyReactions.inp')
        gas = ct.Solution('h2o2_emptyReactions.cti')
        self.assertEqual(gas.n_species, 9)
        self.assertEqual(gas.n_reactions, 0)

    def test_reaction_comments1(self):
        self.convert('pdep-test.inp')
        with open(pjoin(self.test_work_dir, 'pdep-test' + self.ext)) as f:
            text = f.read()
        self.assertIn('Generic mechanism header', text)
        self.assertIn('Single PLOG reaction', text)
        self.assertIn('Multiple PLOG expressions at the same pressure', text)

    def test_reaction_comments2(self):
        self.convert('explicit-third-bodies.inp', thermo='dummy-thermo.dat')
        with open(pjoin(self.test_work_dir, 'explicit-third-bodies' + self.ext)) as f:
            text = f.read()
        self.assertIn('An end of line comment', text)
        self.assertIn('A comment after the last reaction', text)

    def test_custom_element(self):
        self.convert('custom-elements.inp')
        gas = ct.Solution('custom-elements' + self.ext)
        self.assertEqual(gas.n_elements, 4)
        self.assertNear(gas.atomic_weight(2), 13.003)
        self.assertEqual(gas.n_atoms('ethane', 'C'), 2)
        self.assertEqual(gas.n_atoms('CC', 'C'), 1)
        self.assertEqual(gas.n_atoms('CC', 'Ci'), 1)

    def test_surface_mech(self):
        self.convert('surface1-gas.inp', surface='surface1.inp',
                     output='surface1')

        gas = ct.Solution('surface1' + self.ext, 'gas')
        surf = ct.Interface('surface1' + self.ext, 'PT_SURFACE', [gas])

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
        self.convert('surface1-gas-noreac.inp', surface='surface1.inp',
                     output='surface1-nogasreac')

        gas = ct.Solution('surface1-nogasreac' + self.ext, 'gas')
        surf = ct.Interface('surface1-nogasreac' + self.ext, 'PT_SURFACE', [gas])

        self.assertEqual(gas.n_reactions, 0)
        self.assertEqual(surf.n_reactions, 15)


class ck2ctiTest(converterTestCommon, utilities.CanteraTest):
    ext = '.cti'
    InputError = ck2cti.InputParseError

    def _convert(self, inputFile, *, thermo, transport, surface, output, quiet,
                 permissive):
        ck2cti.convertMech(inputFile, thermoFile=thermo,
            transportFile=transport, surfaceFile=surface, outName=output,
            quiet=quiet, permissive=permissive)

    def test_missingElement(self):
        with self.assertRaisesRegex(self.InputError, 'Undefined elements'):
            self.convert('h2o2_missingElement.inp', output='h2o2_missingElement')


class ck2yamlTest(converterTestCommon, utilities.CanteraTest):
    ext = '.yaml'
    InputError = ck2yaml.InputError

    def _convert(self, inputFile, *, thermo, transport, surface, output, quiet,
                 permissive):
        ck2yaml.convert_mech(inputFile, thermo_file=thermo,
            transport_file=transport, surface_file=surface, out_name=output,
            quiet=quiet, permissive=permissive)


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
        cti2yaml.convert(Path(cls.cantera_data).joinpath('gri30.cti'),
                         Path(cls.test_work_dir).joinpath('gri30.yaml'))

    def checkConversion(self, basename, cls=ct.Solution, ctiphases=(),
                        yamlphases=(), **kwargs):
        ctiPhase = cls(basename + '.cti', adjacent=ctiphases, **kwargs)
        yamlPhase = cls(basename + '.yaml', adjacent=yamlphases, **kwargs)

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
        cti2yaml.convert(Path(self.test_data_dir).joinpath('pdep-test.cti'),
                         Path(self.test_work_dir).joinpath('pdep-test.yaml'))
        ctiPhase, yamlPhase = self.checkConversion('pdep-test')
        self.checkKinetics(ctiPhase, yamlPhase, [300, 1000, 2200],
                           [100, ct.one_atm, 2e5, 2e6, 9.9e6])

    def test_ptcombust(self):
        cti2yaml.convert(Path(self.cantera_data).joinpath('ptcombust.cti'),
                         Path(self.test_work_dir).joinpath('ptcombust.yaml'))
        ctiGas, yamlGas = self.checkConversion('ptcombust')
        ctiSurf, yamlSurf = self.checkConversion('ptcombust', ct.Interface,
            name='Pt_surf', ctiphases=[ctiGas], yamlphases=[yamlGas])

        self.checkKinetics(ctiGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctiSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctiSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    def test_ptcombust_motzwise(self):
        cti2yaml.convert(Path(self.test_data_dir).joinpath('ptcombust-motzwise.cti'),
                         Path(self.test_work_dir).joinpath('ptcombust-motzwise.yaml'))
        ctiGas, yamlGas = self.checkConversion('ptcombust-motzwise')
        ctiSurf, yamlSurf = self.checkConversion('ptcombust-motzwise', ct.Interface,
            name='Pt_surf', ctiphases=[ctiGas], yamlphases=[yamlGas])


        self.checkKinetics(ctiGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctiSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctiSurf, yamlSurf, [900], [101325])

    def test_sofc(self):
        cti2yaml.convert(Path(self.cantera_data).joinpath('sofc.cti'),
                         Path(self.test_work_dir).joinpath('sofc.yaml'))
        ctiGas, yamlGas = self.checkConversion('sofc')
        ctiMetal, yamlMetal = self.checkConversion('sofc', name='metal')
        ctiOxide, yamlOxide = self.checkConversion('sofc', name='oxide_bulk')
        ctiMSurf, yamlMSurf = self.checkConversion('sofc', ct.Interface,
            name='metal_surface', ctiphases=[ctiGas, ctiMetal],
            yamlphases=[yamlGas, yamlMetal])
        ctiOSurf, yamlOSurf = self.checkConversion('sofc', ct.Interface,
            name='oxide_surface', ctiphases=[ctiGas, ctiOxide],
            yamlphases=[yamlGas, yamlOxide])
        cti_tpb, yaml_tpb = self.checkConversion('sofc', ct.Interface,
            name='tpb', ctiphases=[ctiMetal, ctiMSurf, ctiOSurf],
            yamlphases=[yamlMetal, yamlMSurf, yamlOSurf])

        self.checkThermo(ctiMSurf, yamlMSurf, [900, 1000, 1100])
        self.checkThermo(ctiOSurf, yamlOSurf, [900, 1000, 1100])
        ctiMetal.electric_potential = yamlMetal.electric_potential = 2
        self.checkKinetics(cti_tpb, yaml_tpb, [900, 1000, 1100], [1e5])
        ctiMetal.electric_potential = yamlMetal.electric_potential = 4
        self.checkKinetics(cti_tpb, yaml_tpb, [900, 1000, 1100], [1e5])

    def test_liquidvapor(self):
        cti2yaml.convert(Path(self.cantera_data).joinpath('liquidvapor.cti'),
                         Path(self.test_work_dir).joinpath('liquidvapor.yaml'))
        for name in ['water', 'nitrogen', 'methane', 'hydrogen', 'oxygen',
                     'hfc134a', 'carbondioxide', 'heptane']:
            ctiPhase, yamlPhase = self.checkConversion('liquidvapor', name=name)
            self.checkThermo(ctiPhase, yamlPhase,
                             [1.3 * ctiPhase.min_temp, 0.7 * ctiPhase.max_temp])

    def test_Redlich_Kwong_CO2(self):
        cti2yaml.convert(Path(self.test_data_dir).joinpath('co2_RK_example.cti'),
                         Path(self.test_work_dir).joinpath('co2_RK_example.yaml'))
        ctiGas, yamlGas = self.checkConversion('co2_RK_example')
        for P in [1e5, 2e6, 1.3e7]:
            yamlGas.TP = ctiGas.TP = 300, P
            self.checkThermo(ctiGas, yamlGas, [300, 400, 500])

    def test_Redlich_Kwong_ndodecane(self):
        cti2yaml.convert(Path(self.cantera_data).joinpath('nDodecane_Reitz.cti'),
                         Path(self.test_work_dir).joinpath('nDodecane_Reitz.yaml'))
        ctiGas, yamlGas = self.checkConversion('nDodecane_Reitz')
        self.checkThermo(ctiGas, yamlGas, [300, 400, 500])
        self.checkKinetics(ctiGas, yamlGas, [300, 500, 1300], [1e5, 2e6, 1.4e7],
                           1e-6)

    def test_diamond(self):
        cti2yaml.convert(Path(self.cantera_data).joinpath('diamond.cti'),
                         Path(self.test_work_dir).joinpath('diamond.yaml'))
        ctiGas, yamlGas = self.checkConversion('diamond', name='gas')
        ctiSolid, yamlSolid = self.checkConversion('diamond', name='diamond')
        ctiSurf, yamlSurf = self.checkConversion('diamond',
            ct.Interface, name='diamond_100', ctiphases=[ctiGas, ctiSolid],
            yamlphases=[yamlGas, yamlSolid])
        self.checkThermo(ctiSolid, yamlSolid, [300, 500])
        self.checkThermo(ctiSurf, yamlSurf, [330, 490])
        self.checkKinetics(ctiSurf, yamlSurf, [400, 800], [2e5])

    def test_lithium_ion_battery(self):
        cti2yaml.convert(Path(self.cantera_data).joinpath('lithium_ion_battery.cti'),
                         Path(self.test_work_dir).joinpath('lithium_ion_battery.yaml'))
        name = 'lithium_ion_battery'
        ctiAnode, yamlAnode = self.checkConversion(name, name='anode')
        ctiCathode, yamlCathode = self.checkConversion(name, name='cathode')
        ctiMetal, yamlMetal = self.checkConversion(name, name='electron')
        ctiElyt, yamlElyt = self.checkConversion(name, name='electrolyte')
        ctiAnodeInt, yamlAnodeInt = self.checkConversion(name,
            name='edge_anode_electrolyte',
            ctiphases=[ctiAnode, ctiMetal, ctiElyt],
            yamlphases=[yamlAnode, yamlMetal, yamlElyt])
        ctiCathodeInt, yamlCathodeInt = self.checkConversion(name,
            name='edge_cathode_electrolyte',
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

    def test_ch4_ion(self):
        cti2yaml.convert(Path(self.test_data_dir).joinpath("ch4_ion.cti"),
                          Path(self.test_work_dir).joinpath("ch4_ion.yaml"))
        ctiGas, yamlGas = self.checkConversion("ch4_ion")
        self.checkThermo(ctiGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctiGas, yamlGas, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctiGas, yamlGas, [298, 1001, 2400])
