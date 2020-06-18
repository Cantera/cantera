import os
from os.path import join as pjoin
import itertools
from pathlib import Path
import logging
import io
import sys

from . import utilities
import cantera as ct
from cantera import ck2cti, ck2yaml, cti2yaml, ctml2yaml

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml


class converterTestCommon:
    def convert(self, inputFile, thermo=None, transport=None,
                surface=None, output=None, extra=None,
                quiet=True, permissive=None):
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
        if extra is not None:
            extra = pjoin(self.test_data_dir, extra)
        output = pjoin(self.test_work_dir, output + self.ext)
        if os.path.exists(output):
            os.remove(output)
        self._convert(inputFile, thermo=thermo, transport=transport,
            surface=surface, output=output, extra=extra,
            quiet=quiet, permissive=permissive)
        return output

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
                self.assertNear(ref_kf[i], gas_kf[i], rtol=tol, msg='kf' + message)
                self.assertNear(ref_kr[i], gas_kr[i], rtol=tol, msg='kr' + message)

    def test_gri30(self):
        output = self.convert('gri30.inp', thermo='gri30_thermo.dat',
                              transport='gri30_tran.dat', output='gri30_test')

        ref, gas = self.checkConversion('gri30.xml', output)
        self.checkKinetics(ref, gas, [300, 1500], [5e3, 1e5, 2e6])

    def test_soot(self):
        output = self.convert('soot.inp', thermo='soot-therm.dat', output='soot_test')
        ref, gas = self.checkConversion('soot.xml', output)
        self.checkThermo(ref, gas, [300, 1100])
        self.checkKinetics(ref, gas, [300, 1100], [5e3, 1e5, 2e6])

    def test_pdep(self):
        output = self.convert('pdep-test.inp')
        ref, gas = self.checkConversion('pdep-test.xml', output)
        # Chebyshev coefficients in XML are truncated to 6 digits, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6],
                           tol=2e-4)

    def test_species_only(self):
        self.convert(None, thermo='dummy-thermo.dat', output='dummy-thermo')
        if self.ext == ".cti":
            cti = "ideal_gas(elements='C H', species='dummy-thermo:R1A R1B P1')"
            gas = ct.Solution(source=cti)
        elif self.ext == ".yaml":
            yaml = ("{phases: [{name: gas, species: "
                    "[{dummy-thermo.yaml/species: [R1A, R1B, P1]}], "
                    "thermo: ideal-gas}]}")
            gas = ct.Solution(yaml=yaml)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 0)

    def test_missingThermo(self):
        with self.assertRaisesRegex(self.InputError, 'No thermo data'):
            self.convert('h2o2_missingThermo.inp')

    def test_duplicate_thermo(self):
        with self.assertRaisesRegex(self.InputError, 'additional thermo'):
            self.convert('duplicate-thermo.inp')

        output = self.convert('duplicate-thermo.inp', permissive=True)

        gas = ct.Solution(output)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_duplicate_species(self):
        with self.assertRaisesRegex(self.InputError, 'additional declaration'):
            self.convert('duplicate-species.inp')

        output = self.convert('duplicate-species.inp', permissive=True)

        gas = ct.Solution(output)
        self.assertEqual(gas.species_names, ['foo','bar','baz'])

    def test_pathologicalSpeciesNames(self):
        output = self.convert('species-names.inp')
        gas = ct.Solution(output)

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

        output = self.convert('unterminated-sections.inp', permissive=True)
        gas = ct.Solution(output)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unterminatedSections2(self):
        with self.assertRaisesRegex(self.InputError, 'implicitly ended'):
            self.convert('unterminated-sections2.inp')

        output = self.convert('unterminated-sections2.inp', permissive=True)
        gas = ct.Solution(output)
        self.assertEqual(gas.n_species, 3)
        self.assertEqual(gas.n_reactions, 2)

    def test_unrecognized_section(self):
        with self.assertRaisesRegex(self.InputError, 'SPAM'):
            self.convert('unrecognized-section.inp', thermo='dummy-thermo.dat',
                         permissive=True)

    def test_nasa9(self):
        output = self.convert('nasa9-test.inp', thermo='nasa9-test-therm.dat')
        ref, gas = self.checkConversion('nasa9-test.xml', output)
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_nasa9_subset(self):
        output = self.convert('nasa9-test-subset.inp', thermo='nasa9-test-therm.dat')
        ref, gas = self.checkConversion('nasa9-test-subset.xml', output)
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_sri_falloff(self):
        output = self.convert('sri-falloff.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('sri-falloff.xml', output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_chemically_activated(self):
        output = self.convert('chemically-activated-reaction.inp')
        ref, gas = self.checkConversion('chemically-activated-reaction.xml',
                                        output)
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6, 1e7],
                           tol=1e-7)

    def test_explicit_third_bodies(self):
        output = self.convert('explicit-third-bodies.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('explicit-third-bodies.xml', output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_explicit_reverse_rate(self):
        output = self.convert('explicit-reverse-rate.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('explicit-reverse-rate.xml', output)
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
        output = self.convert('explicit-forward-order.inp', thermo='dummy-thermo.dat')
        ref, gas = self.checkConversion('explicit-forward-order.xml', output)
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6],
                           tol=2e-7)

    def test_negative_order(self):
        with self.assertRaisesRegex(self.InputError, 'Negative reaction order'):
            self.convert('negative-order.inp', thermo='dummy-thermo.dat')

    def test_negative_order_permissive(self):
        output = self.convert('negative-order.inp', thermo='dummy-thermo.dat',
                              permissive=True)
        ref, gas = self.checkConversion('negative-order.cti', output)
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6],
                           tol=2e-7)

    def test_negative_A_factor(self):
        output = self.convert('negative-rate.inp', thermo='dummy-thermo.dat')
        gas = ct.Solution(output)  # Validate the mechanism
        self.assertLess(gas.reaction(4).rate.pre_exponential_factor, 0)
        self.assertLess(gas.reaction(1).rate.pre_exponential_factor, 0)
        self.assertLess(gas.reaction(2).rate.pre_exponential_factor, 0)
        self.assertLess(gas.forward_rate_constants[5], 0)

    def test_bad_troe_value(self):
        with self.assertRaises(ValueError):
            self.convert('bad-troe.inp', thermo='dummy-thermo.dat')

    def test_invalid_reaction_equation(self):
        with self.assertRaisesRegex(self.InputError, 'Unparsable'):
            self.convert('invalid-equation.inp', thermo='dummy-thermo.dat')

    def test_reaction_units(self):
        out_def = self.convert('units-default.inp', thermo='dummy-thermo.dat')
        out_cus = self.convert('units-custom.inp', thermo='dummy-thermo.dat')
        default, custom = self.checkConversion(out_def, out_cus)
        self.checkKinetics(default, custom,
                           [300, 800, 1450, 2800], [5e0, 5e3, 1e5, 2e6, 1e8], 1e-7)

    def test_float_stoich_coeffs(self):
        output = self.convert('float-stoich.inp', thermo='dummy-thermo.dat')
        gas = ct.Solution(output)

        R = gas.reactant_stoich_coeffs()
        P = gas.product_stoich_coeffs()
        self.assertArrayNear(R[:,0], [0, 1.5, 0.5, 0])
        self.assertArrayNear(P[:,0], [1, 0, 0, 1])
        self.assertArrayNear(R[:,1], [1, 0, 0, 1])
        self.assertArrayNear(P[:,1], [0, 0.33, 1.67, 0])

    def test_photon(self):
        output = self.convert('photo-reaction.inp', thermo='dummy-thermo.dat',
                              permissive=True)

        ref, gas = self.checkConversion('photo-reaction.xml', output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_transport_normal(self):
        output = self.convert('h2o2.inp', transport='gri30_tran.dat',
                              output='h2o2_transport_normal')

        gas = ct.Solution(output)
        gas.TPX = 300, 101325, 'H2:1.0, O2:1.0'
        self.assertAlmostEqual(gas.thermal_conductivity, 0.07663, 4)

    def test_transport_embedded(self):
        output = self.convert('with-transport.inp')
        gas = ct.Solution(output)
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
        output = self.convert('h2o2_emptyReactions.inp')
        gas = ct.Solution(output)
        self.assertEqual(gas.n_species, 9)
        self.assertEqual(gas.n_reactions, 0)

    def test_reaction_comments1(self):
        output = self.convert('pdep-test.inp')
        text = Path(output).read_text()
        self.assertIn('Generic mechanism header', text)
        self.assertIn('Single PLOG reaction', text)
        self.assertIn('Multiple PLOG expressions at the same pressure', text)

    def test_reaction_comments2(self):
        output = self.convert('explicit-third-bodies.inp', thermo='dummy-thermo.dat')
        text = Path(output).read_text()
        self.assertIn('An end of line comment', text)
        self.assertIn('A comment after the last reaction', text)

    def test_custom_element(self):
        output = self.convert('custom-elements.inp')
        gas = ct.Solution(output)
        self.assertEqual(gas.n_elements, 4)
        self.assertNear(gas.atomic_weight(2), 13.003)
        self.assertEqual(gas.n_atoms('ethane', 'C'), 2)
        self.assertEqual(gas.n_atoms('CC', 'C'), 1)
        self.assertEqual(gas.n_atoms('CC', 'Ci'), 1)

    def test_surface_mech(self):
        output = self.convert('surface1-gas.inp', surface='surface1.inp',
                              output='surface1')

        gas = ct.Solution(output, 'gas')
        surf = ct.Interface(output, 'PT_SURFACE', [gas])

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
        output = self.convert('surface1-gas-noreac.inp', surface='surface1.inp',
                              output='surface1-nogasreac')

        gas = ct.Solution(output, 'gas')
        surf = ct.Interface(output, 'PT_SURFACE', [gas])

        self.assertEqual(gas.n_reactions, 0)
        self.assertEqual(surf.n_reactions, 15)

    def test_third_body_plus_falloff_reactions(self):
        self.convert('third_body_plus_falloff_reaction.inp')
        gas = ct.Solution('third_body_plus_falloff_reaction' + self.ext)
        self.assertEqual(gas.n_reactions, 2)

    def test_blank_line_in_header(self):
        self.convert('blank_line_in_header.inp')
        gas = ct.Solution('blank_line_in_header' + self.ext)
        self.assertEqual(gas.n_reactions, 1)


class ck2ctiTest(converterTestCommon, utilities.CanteraTest):
    ext = '.cti'
    InputError = ck2cti.InputParseError

    def _convert(self, inputFile, *, thermo, transport, surface, output, extra,
                 quiet, permissive):
        ck2cti.convertMech(inputFile, thermoFile=thermo,
            transportFile=transport, surfaceFile=surface, outName=output,
            quiet=quiet, permissive=permissive)

    def test_missingElement(self):
        with self.assertRaisesRegex(self.InputError, 'Undefined elements'):
            self.convert('h2o2_missingElement.inp', output='h2o2_missingElement')


class ck2yamlTest(converterTestCommon, utilities.CanteraTest):
    ext = '.yaml'
    InputError = ck2yaml.InputError

    def _convert(self, inputFile, *, thermo, transport, surface, output, extra,
                 quiet, permissive):
        ck2yaml.convert_mech(inputFile, thermo_file=thermo,
            transport_file=transport, surface_file=surface, out_name=output,
            extra_file=extra, quiet=quiet, permissive=permissive)

    def test_extra(self):
        self.convert('gri30.inp', thermo='gri30_thermo.dat',
                     transport='gri30_tran.dat', output='gri30_extra',
                     extra='extra.yaml')

        output = pjoin(self.test_work_dir, 'gri30_extra' + self.ext)
        with open(output, 'rt', encoding="utf-8") as stream:
            yml = yaml.safe_load(stream)

        desc = yml['description'].split('\n')[-1]
        self.assertEqual(desc, 'This is an alternative description.')
        for key in ['foo', 'bar']:
            self.assertIn(key, yml.keys())

    def test_duplicate_reactions(self):
        # Running a test this way instead of using the convertMech function
        # tests the behavior of the ck2yaml.main function and the mechanism
        # validation step.

        # clear global handler created by logging.basicConfig() in ck2cti
        logging.getLogger().handlers.clear()

        # Replace the ck2yaml logger with our own in order to capture the output
        log_stream = io.StringIO()
        logger = logging.getLogger('cantera.ck2yaml')
        original_handler = logger.handlers.pop()
        logformatter = logging.Formatter('%(message)s')
        handler = logging.StreamHandler(log_stream)
        handler.setFormatter(logformatter)
        logger.addHandler(handler)

        with self.assertRaises(SystemExit):
            ck2yaml.main([
                '--input={}/undeclared-duplicate-reactions.inp'.format(self.test_data_dir),
                '--thermo={}/dummy-thermo.dat'.format(self.test_data_dir)])

        # Put the original logger back in place
        logger.handlers.clear()
        logger.addHandler(original_handler)

        message = log_stream.getvalue()
        for token in ('FAILED', 'lines 12 and 14', 'R1A', 'R1B'):
            self.assertIn(token, message)


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
    def checkConversion(self, basename, cls=ct.Solution, ctiphases=(),
                        yamlphases=(), **kwargs):
        ctiPhase = cls(basename + '.cti', adjacent=ctiphases, **kwargs)
        yamlPhase = cls(basename + '.yaml', adjacent=yamlphases, **kwargs)

        self.assertEqual(ctiPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctiPhase.species_names, yamlPhase.species_names)
        self.assertEqual(ctiPhase.n_reactions, yamlPhase.n_reactions)
        for C, Y in zip(ctiPhase.species(), yamlPhase.species()):
            self.assertEqual(C.composition, Y.composition)

        for C, Y in zip(ctiPhase.reactions(), yamlPhase.reactions()):
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
        cti2yaml.convert(Path(self.cantera_data).joinpath('gri30.cti'),
                         Path(self.test_work_dir).joinpath('gri30.yaml'))
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
        output_file = Path(self.test_work_dir).joinpath('liquidvapor.yaml')
        cti2yaml.convert(Path(self.cantera_data).joinpath('liquidvapor.cti'),
                         output_file)
        for name in ['water', 'nitrogen', 'methane', 'hydrogen', 'oxygen',
                     'hfc134a', 'carbondioxide', 'heptane']:
            ctiPhase, yamlPhase = self.checkConversion('liquidvapor', name=name)
            self.checkThermo(ctiPhase, yamlPhase,
                             [1.3 * ctiPhase.min_temp, 0.7 * ctiPhase.max_temp])
        # The output file must be removed after this test because it shadows
        # a file distributed with Cantera that is used in other tests.
        output_file.unlink()

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

class ctml2yamlTest(utilities.CanteraTest):

    def checkConversion(self, basename, cls=ct.Solution, ctmlphases=(),
                        yamlphases=(), **kwargs):
        ctmlPhase = cls(basename + '.xml', adjacent=ctmlphases, **kwargs)
        yamlPhase = cls(basename + '.yaml', adjacent=yamlphases, **kwargs)

        self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
        self.assertEqual(ctmlPhase.n_reactions, yamlPhase.n_reactions)
        for C, Y in zip(ctmlPhase.species(), yamlPhase.species()):
            self.assertEqual(C.composition, Y.composition)

        for C, Y in zip(ctmlPhase.reactions(), yamlPhase.reactions()):
            self.assertEqual(C.__class__, Y.__class__)
            self.assertEqual(C.reactants, Y.reactants)
            self.assertEqual(C.products, Y.products)
            self.assertEqual(C.duplicate, Y.duplicate)

        for i, sp in zip(range(ctmlPhase.n_reactions), ctmlPhase.kinetics_species_names):
            self.assertEqual(ctmlPhase.reactant_stoich_coeff(sp, i),
                             yamlPhase.reactant_stoich_coeff(sp, i))

        return ctmlPhase, yamlPhase

    def checkThermo(self, ctmlPhase, yamlPhase, temperatures, tol=1e-7):
        for T in temperatures:
            ctmlPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            cp_ctml = ctmlPhase.partial_molar_cp
            cp_yaml = yamlPhase.partial_molar_cp
            h_ctml = ctmlPhase.partial_molar_enthalpies
            h_yaml = yamlPhase.partial_molar_enthalpies
            s_ctml = ctmlPhase.partial_molar_entropies
            s_yaml = yamlPhase.partial_molar_entropies
            self.assertNear(ctmlPhase.density, yamlPhase.density)
            for i in range(ctmlPhase.n_species):
                message = ' for species {0} at T = {1}'.format(ctmlPhase.species_names[i], T)
                self.assertNear(cp_ctml[i], cp_yaml[i], tol, msg='cp'+message)
                self.assertNear(h_ctml[i], h_yaml[i], tol, msg='h'+message)
                self.assertNear(s_ctml[i], s_yaml[i], tol, msg='s'+message)

    def checkKinetics(self, ctmlPhase, yamlPhase, temperatures, pressures, tol=1e-7):
        for T,P in itertools.product(temperatures, pressures):
            ctmlPhase.TP = T, P
            yamlPhase.TP = T, P
            kf_ctml = ctmlPhase.forward_rate_constants
            kr_ctml = ctmlPhase.reverse_rate_constants
            kf_yaml = yamlPhase.forward_rate_constants
            kr_yaml = yamlPhase.reverse_rate_constants
            for i in range(yamlPhase.n_reactions):
                message = ' for reaction {0} at T = {1}, P = {2}'.format(i, T, P)
                self.assertNear(kf_ctml[i], kf_yaml[i], rtol=tol, msg='kf '+message)
                self.assertNear(kr_ctml[i], kr_yaml[i], rtol=tol, msg='kr '+message)

    def checkTransport(self, ctmlPhase, yamlPhase, temperatures,
                       model='mixture-averaged'):
        ctmlPhase.transport_model = model
        yamlPhase.transport_model = model
        for T in temperatures:
            ctmlPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            self.assertNear(ctmlPhase.viscosity, yamlPhase.viscosity)
            self.assertNear(ctmlPhase.thermal_conductivity,
                            yamlPhase.thermal_conductivity)
            Dkm_ctml = ctmlPhase.mix_diff_coeffs
            Dkm_yaml = yamlPhase.mix_diff_coeffs
            for i in range(ctmlPhase.n_species):
                message = 'dkm for species {0} at T = {1}'.format(i, T)
                self.assertNear(Dkm_ctml[i], Dkm_yaml[i], msg=message)

    def test_gri30(self):
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath('gri30.xml'),
            Path(self.test_work_dir).joinpath('gri30.yaml'),
        )
        ctmlPhase, yamlPhase = self.checkConversion('gri30')
        X = {'O2': 0.3, 'H2': 0.1, 'CH4': 0.2, 'CO2': 0.4}
        ctmlPhase.X = X
        yamlPhase.X = X
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlPhase, yamlPhase, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctmlPhase, yamlPhase, [298, 1001, 2400])

    def test_pdep(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath('pdep-test.xml'),
            Path(self.test_work_dir).joinpath('pdep-test.yaml'),
        )
        ctmlPhase, yamlPhase = self.checkConversion('pdep-test')
        self.checkKinetics(ctmlPhase, yamlPhase, [300, 1000, 2200],
                           [100, ct.one_atm, 2e5, 2e6, 9.9e6])

    def test_ptcombust(self):
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath('ptcombust.xml'),
            Path(self.test_work_dir).joinpath('ptcombust.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('ptcombust')
        ctmlSurf, yamlSurf = self.checkConversion('ptcombust', ct.Interface,
            name='Pt_surf', ctmlphases=[ctmlGas], yamlphases=[yamlGas])

        self.checkKinetics(ctmlGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctmlSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctmlSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    def test_ptcombust_motzwise(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath('ptcombust-motzwise.xml'),
            Path(self.test_work_dir).joinpath('ptcombust-motzwise.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('ptcombust-motzwise')
        ctmlSurf, yamlSurf = self.checkConversion('ptcombust-motzwise', ct.Interface,
            name='Pt_surf', ctmlphases=[ctmlGas], yamlphases=[yamlGas])

        self.checkKinetics(ctmlGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctmlSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctmlSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    def test_sofc(self):
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath('sofc.xml'),
            Path(self.test_work_dir).joinpath('sofc.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('sofc')
        ctmlMetal, yamlMetal = self.checkConversion('sofc', name='metal')
        ctmlOxide, yamlOxide = self.checkConversion('sofc', name='oxide_bulk')
        ctmlMSurf, yamlMSurf = self.checkConversion('sofc', ct.Interface,
            name='metal_surface', ctmlphases=[ctmlGas, ctmlMetal],
            yamlphases=[yamlGas, yamlMetal])
        ctmlOSurf, yamlOSurf = self.checkConversion('sofc', ct.Interface,
            name='oxide_surface', ctmlphases=[ctmlGas, ctmlOxide],
            yamlphases=[yamlGas, yamlOxide])
        ctml_tpb, yaml_tpb = self.checkConversion('sofc', ct.Interface,
            name='tpb', ctmlphases=[ctmlMetal, ctmlMSurf, ctmlOSurf],
            yamlphases=[yamlMetal, yamlMSurf, yamlOSurf])

        self.checkThermo(ctmlMSurf, yamlMSurf, [900, 1000, 1100])
        self.checkThermo(ctmlOSurf, yamlOSurf, [900, 1000, 1100])
        ctmlMetal.electric_potential = yamlMetal.electric_potential = 2
        self.checkKinetics(ctml_tpb, yaml_tpb, [900, 1000, 1100], [1e5])
        ctmlMetal.electric_potential = yamlMetal.electric_potential = 4
        self.checkKinetics(ctml_tpb, yaml_tpb, [900, 1000, 1100], [1e5])

    def test_liquidvapor(self):
        output_file = Path(self.test_work_dir).joinpath('liquidvapor.yaml')
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath('liquidvapor.xml'), output_file,
        )
        for name in ['water', 'nitrogen', 'methane', 'hydrogen', 'oxygen',
                     'hfc134a', 'carbondioxide', 'heptane']:
            ctmlPhase, yamlPhase = self.checkConversion('liquidvapor', name=name)
            self.checkThermo(ctmlPhase, yamlPhase,
                             [1.3 * ctmlPhase.min_temp, 0.7 * ctmlPhase.max_temp])
        # The output file must be removed after this test because it shadows
        # a file distributed with Cantera that is used in other tests.
        output_file.unlink()

    def test_Redlich_Kwong_CO2(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath('co2_RK_example.xml'),
            Path(self.test_work_dir).joinpath('co2_RK_example.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('co2_RK_example')
        for P in [1e5, 2e6, 1.3e7]:
            yamlGas.TP = ctmlGas.TP = 300, P
            self.checkThermo(ctmlGas, yamlGas, [300, 400, 500])

    def test_Redlich_Kwong_ndodecane(self):
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath('nDodecane_Reitz.xml'),
            Path(self.test_work_dir).joinpath('nDodecane_Reitz.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('nDodecane_Reitz')
        self.checkThermo(ctmlGas, yamlGas, [300, 400, 500])
        self.checkKinetics(ctmlGas, yamlGas, [300, 500, 1300], [1e5, 2e6, 1.4e7],
                           1e-6)

    def test_diamond(self):
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath('diamond.xml'),
            Path(self.test_work_dir).joinpath('diamond.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('diamond', name='gas')
        ctmlSolid, yamlSolid = self.checkConversion('diamond', name='diamond')
        ctmlSurf, yamlSurf = self.checkConversion('diamond',
            ct.Interface, name='diamond_100', ctmlphases=[ctmlGas, ctmlSolid],
            yamlphases=[yamlGas, yamlSolid])
        self.checkThermo(ctmlSolid, yamlSolid, [300, 500])
        self.checkThermo(ctmlSurf, yamlSurf, [330, 490])
        self.checkKinetics(ctmlSurf, yamlSurf, [400, 800], [2e5])

    def test_lithium_ion_battery(self):
        name = 'lithium_ion_battery'
        ctml2yaml.convert(
            Path(self.cantera_data).joinpath(name + ".xml"),
            Path(self.test_work_dir).joinpath(name + ".yaml"),
        )
        ctmlAnode, yamlAnode = self.checkConversion(name, name='anode')
        ctmlCathode, yamlCathode = self.checkConversion(name, name='cathode')
        ctmlMetal, yamlMetal = self.checkConversion(name, name='electron')
        ctmlElyt, yamlElyt = self.checkConversion(name, name='electrolyte')
        ctmlAnodeInt, yamlAnodeInt = self.checkConversion(name,
            name='edge_anode_electrolyte',
            ctmlphases=[ctmlAnode, ctmlMetal, ctmlElyt],
            yamlphases=[yamlAnode, yamlMetal, yamlElyt])
        ctmlCathodeInt, yamlCathodeInt = self.checkConversion(name,
            name='edge_cathode_electrolyte',
            ctmlphases=[ctmlCathode, ctmlMetal, ctmlElyt],
            yamlphases=[yamlCathode, yamlMetal, yamlElyt])

        self.checkThermo(ctmlAnode, yamlAnode, [300, 330])
        self.checkThermo(ctmlCathode, yamlCathode, [300, 330])

        ctmlAnode.X = yamlAnode.X = [0.7, 0.3]
        self.checkThermo(ctmlAnode, yamlAnode, [300, 330])
        ctmlCathode.X = yamlCathode.X = [0.2, 0.8]
        self.checkThermo(ctmlCathode, yamlCathode, [300, 330])

        for phase in [ctmlAnode, yamlAnode, ctmlCathode, yamlCathode, ctmlMetal,
                      yamlMetal, ctmlElyt, yamlElyt, ctmlAnodeInt, yamlAnodeInt,
                      ctmlCathodeInt, yamlCathodeInt]:
            phase.TP = 300, 1e5
        ctmlMetal.electric_potential = yamlMetal.electric_potential = 0
        ctmlElyt.electric_potential = yamlElyt.electric_potential = 1.9
        self.checkKinetics(ctmlAnodeInt, yamlAnodeInt, [300], [1e5])

        ctmlMetal.electric_potential = yamlMetal.electric_potential = 2.2
        ctmlElyt.electric_potential = yamlElyt.electric_potential = 0
        self.checkKinetics(ctmlCathodeInt, yamlCathodeInt, [300], [1e5])

    def test_noxNeg(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath('noxNeg.xml'),
            Path(self.test_work_dir).joinpath('noxNeg.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('noxNeg')
        self.checkThermo(ctmlGas, yamlGas, [300, 1000])
        self.checkKinetics(ctmlGas, yamlGas, [300, 1000], [1e5])

    def test_ch4_ion(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("ch4_ion.xml"),
            Path(self.test_work_dir).joinpath("ch4_ion.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("ch4_ion")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctmlGas, yamlGas, [298, 1001, 2400])

    def test_nasa9(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("nasa9-test.xml"),
            Path(self.test_work_dir).joinpath("nasa9-test.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("nasa9-test")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])

    def test_chemically_activated(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("chemically-activated-reaction.xml"),
            Path(self.test_work_dir).joinpath("chemically-activated-reaction.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("chemically-activated-reaction")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_explicit_forward_order(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("explicit-forward-order.xml"),
            Path(self.test_work_dir).joinpath("explicit-forward-order.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("explicit-forward-order")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_explicit_reverse_rate(self):
        ctml2yaml.convert(Path(self.test_data_dir).joinpath("explicit-reverse-rate.xml"),
                          Path(self.test_work_dir).joinpath("explicit-reverse-rate.yaml"))
        ctmlGas, yamlGas = self.checkConversion("explicit-reverse-rate")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_explicit_third_bodies(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("explicit-third-bodies.xml"),
            Path(self.test_work_dir).joinpath("explicit-third-bodies.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("explicit-third-bodies")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_fractional_stoich_coeffs(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("frac.xml"),
            Path(self.test_work_dir).joinpath("frac.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("frac")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    # @todo Remove after Cantera 2.5 - class FixedChemPotSSTP is deprecated
    def test_fixed_chemical_potential_thermo(self):
        ct.suppress_deprecation_warnings()
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("LiFixed.xml"),
            Path(self.test_work_dir).joinpath("LiFixed.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("LiFixed")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        ct.make_deprecation_warnings_fatal()

    def test_water_IAPWS95_thermo(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("liquid-water.xml"),
            Path(self.test_work_dir).joinpath("liquid-water.yaml"),
        )
        ctmlWater, yamlWater = self.checkConversion("liquid-water")
        self.checkThermo(ctmlWater, yamlWater, [300, 500, 1300, 2000])
        self.assertEqual(ctmlWater.transport_model, yamlWater.transport_model)
        dens = ctmlWater.density
        for T in [298, 1001, 2400]:
            ctmlWater.TD = T, dens
            yamlWater.TD = T, dens
            self.assertNear(ctmlWater.viscosity, yamlWater.viscosity)
            self.assertNear(ctmlWater.thermal_conductivity,
                            yamlWater.thermal_conductivity)

    def test_hmw_nacl_phase(self):
        basename = "HMW_NaCl_sp1977_alt"
        xml_file = Path(self.test_data_dir).joinpath(basename).with_suffix(".xml")
        yaml_file = Path(self.test_data_dir).joinpath(basename).with_suffix(".yaml")
        ctml2yaml.convert(xml_file, yaml_file)

        # Can only be loaded by ThermoPhase due to a bug in TransportFactory
        # ThermoPhase does not have reactions (neither does the input file)
        # so we can't use checkConversion
        ctmlPhase = ct.ThermoPhase(str(xml_file))
        yamlPhase = ct.ThermoPhase(str(yaml_file))
        self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_NaCl_solid_phase(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("NaCl_Solid.xml"),
            Path(self.test_work_dir).joinpath("NaCl_Solid.yaml"),
        )
        ctmlPhase, yamlPhase = self.checkConversion("NaCl_Solid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500, 1300, 2000])

    def test_DH_NaCl_phase(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("debye-huckel-all.xml"),
            Path(self.test_work_dir).joinpath("debye-huckel-all.yaml"),
        )
        for name in [
            "debye-huckel-dilute",
            "debye-huckel-B-dot-ak",
            "debye-huckel-B-dot-a",
            "debye-huckel-pitzer-beta_ij",
            "debye-huckel-beta_ij",
        ]:
            # Can only be loaded by ThermoPhase due to a bug in TransportFactory
            # ThermoPhase does not have reactions (neither does the input file)
            # so we can't use checkConversion
            ctmlPhase = ct.ThermoPhase("debye-huckel-all.xml", name=name)
            yamlPhase = ct.ThermoPhase("debye-huckel-all.yaml", name=name)
            self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
            self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
            self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_Maskell_solid_soln(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("MaskellSolidSolnPhase_valid.xml"),
            Path(self.test_work_dir).joinpath("MaskellSolidSolnPhase_valid.yaml"),
        )

        ctmlPhase, yamlPhase = self.checkConversion("MaskellSolidSolnPhase_valid")
        # Maskell phase doesn't support partial molar properties, so just check density
        for T in [300, 500, 1300, 2000]:
            ctmlPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            self.assertNear(ctmlPhase.density, yamlPhase.density)

    def test_mock_ion(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("mock_ion.xml"),
            Path(self.test_work_dir).joinpath("mock_ion.yaml"),
        )
        ctmlPhase = ct.ThermoPhase("mock_ion.xml")
        yamlPhase = ct.ThermoPhase("mock_ion.yaml")
        # Due to changes in how the species elements are specified, the composition
        # of the species differs from XML to YAML (electrons are used to specify charge
        # in YAML while the charge node is used in XML). Therefore, checkConversion
        # won't work and we have to check a few things manually. There are also no
        # reactions specified for these phases so don't need to do any checks for that.
        self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
        # ions-from-neutral-molecule phase doesn't support partial molar properties,
        # so just check density
        for T in [300, 500, 1300, 2000]:
            ctmlPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            self.assertNear(ctmlPhase.density, yamlPhase.density)

    def test_Redlich_Kister(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("RedlichKisterVPSSTP_valid.xml"),
            Path(self.test_work_dir).joinpath("RedlichKisterVPSSTP_valid.yaml"),
        )

        ctmlPhase, yamlPhase = self.checkConversion("RedlichKisterVPSSTP_valid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_species_names(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath('species-names.xml'),
            Path(self.test_work_dir).joinpath('species-names.yaml'),
        )
        ctmlGas, yamlGas = self.checkConversion('species-names')
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])

    def test_sri_falloff_reaction(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("sri-falloff.xml"),
            Path(self.test_work_dir).joinpath("sri-falloff.yaml"),
        )
        ctmlGas, yamlGas = self.checkConversion("sri-falloff")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_vpss_and_hkft(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("pdss_hkft.xml"),
            Path(self.test_work_dir).joinpath("pdss_hkft.yaml"),
        )
        # @todo Remove "gas" mode test after Cantera 2.5 - "gas" mode of class
        #     IdealSolnGasVPSS is deprecated
        ct.suppress_deprecation_warnings()
        for name in ["vpss_gas_pdss_hkft_phase", "vpss_soln_pdss_hkft_phase"]:
            ctmlPhase = ct.ThermoPhase("pdss_hkft.xml", name=name)
            yamlPhase = ct.ThermoPhase("pdss_hkft.yaml", name=name)
            # Due to changes in how the species elements are specified, the
            # composition of the species differs from XML to YAML (electrons are used
            # to specify charge in YAML while the charge node is used in XML).
            # Therefore, checkConversion won't work and we have to check a few things
            # manually. There are also no reactions specified for these phases so don't
            # need to do any checks for that.
            self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
            self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
            self.checkThermo(ctmlPhase, yamlPhase, [300, 500])
        ct.make_deprecation_warnings_fatal()

    def test_lattice_solid(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("Li7Si3_ls.xml"),
            Path(self.test_work_dir).joinpath("Li7Si3_ls.yaml"),
        )
        # Use ThermoPhase to avoid constructing a default Transport object which
        # throws an error for the LatticeSolidPhase
        basename = "Li7Si3_ls"
        name = "Li7Si3_and_Interstitials(S)"
        ctmlPhase = ct.ThermoPhase(basename + ".xml", name=name)
        yamlPhase = ct.ThermoPhase(basename + ".yaml", name=name)
        self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_margules(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("LiKCl_liquid.xml"),
            Path(self.test_work_dir).joinpath("LiKCl_liquid.yaml"),
        )
        ctmlPhase, yamlPhase = self.checkConversion("LiKCl_liquid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_idealsolidsoln(self):
        with self.assertWarnsRegex(UserWarning, "SolidKinetics type is not implemented"):
            ctml2yaml.convert(
                Path(self.test_data_dir).joinpath("IdealSolidSolnPhaseExample.xml"),
                Path(self.test_work_dir).joinpath("IdealSolidSolnPhaseExample.yaml"),
            )

        # SolidKinetics is not implemented, so can't create a Kinetics class instance.
        basename = "IdealSolidSolnPhaseExample"
        ctmlPhase = ct.ThermoPhase(basename + ".xml")
        yamlPhase = ct.ThermoPhase(basename + ".yaml")

        self.assertEqual(ctmlPhase.element_names, yamlPhase.element_names)
        self.assertEqual(ctmlPhase.species_names, yamlPhase.species_names)
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_idealmolalsoln(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("IdealMolalSolnPhaseExample.xml"),
            Path(self.test_work_dir).joinpath("IdealMolalSolnPhaseExample.yaml"),
        )

        ctmlPhase, yamlPhase = self.checkConversion("IdealMolalSolnPhaseExample")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_transport_models(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("transport_models_test.xml"),
            Path(self.test_work_dir).joinpath("transport_models_test.yaml"),
        )
        for name in ["UnityLewis", "CK_Mix", "CK_Multi", "HighP"]:
            ctmlPhase, yamlPhase = self.checkConversion("transport_models_test", name=name)
            self.checkTransport(ctmlPhase, yamlPhase, [298, 1001, 2500])

    def test_nonreactant_orders(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("reaction-orders.xml"),
            Path(self.test_work_dir).joinpath("reaction-orders.yaml"),
        )

        ctmlPhase, yamlPhase = self.checkConversion("reaction-orders")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])
        self.checkKinetics(ctmlPhase, yamlPhase, [300, 1001, 2500], [1e5, 10e5])

    def test_species_ss_temperature_polynomials(self):
        ctml2yaml.convert(
            Path(self.test_data_dir).joinpath("Li_Liquid.xml"),
            Path(self.test_work_dir).joinpath("Li_Liquid.yaml"),
        )

        ctmlPhase, yamlPhase = self.checkConversion("Li_Liquid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_duplicate_section_ids(self):
        with self.assertWarnsRegex(UserWarning, "Duplicate 'speciesData' id"):
                ctml2yaml.convert(
                    Path(self.test_data_dir).joinpath("duplicate-speciesData-ids.xml"),
                    Path(self.test_work_dir).joinpath("duplicate-speciesData-ids.yaml")
                )
        with self.assertWarnsRegex(UserWarning, "Duplicate 'reactionData' id"):
                ctml2yaml.convert(
                    Path(self.test_data_dir).joinpath("duplicate-reactionData-ids.xml"),
                    Path(self.test_work_dir).joinpath("duplicate-reactionData-ids.yaml")
                )
