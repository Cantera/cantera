from __future__ import annotations

import itertools
from pathlib import Path
import pytest
from pytest import approx

import cantera as ct
from cantera import ck2yaml, cti2yaml, ctml2yaml, yaml2ck, lxcat2yaml
from .utilities import load_yaml


class Testck2yaml:
    @pytest.fixture(autouse=True)
    def inject_fixtures(self, capsys, test_data_path):
        self._capsys = capsys
        self.test_data_path = test_data_path

    def convert(self, inputFile, thermo=None, transport=None,
                surface=None, output=None, quiet=True, extra=None, **kwargs):
        if output is None:
            output = Path(inputFile or thermo).stem  # strip '.inp'
        if inputFile is not None:
            inputFile = self.test_data_path / inputFile
        if thermo is not None:
            thermo = self.test_data_path / thermo
        if transport is not None:
            transport = self.test_data_path / transport
        if surface is not None:
            surface = self.test_data_path / surface
        if extra is not None:
            extra = self.test_data_path / extra
        output = self.test_work_path / (output + "-from-ck.yaml")
        output.unlink(missing_ok=True)
        ck2yaml.convert(inputFile, thermo_file=thermo,
            transport_file=transport, surface_file=surface, out_name=output,
            extra_file=extra, quiet=quiet, **kwargs)
        return output

    def checkConversion(self, refFile, testFile):
        ref = ct.Solution(refFile)
        gas = ct.Solution(testFile)

        assert ref.element_names == gas.element_names
        assert ref.species_names == gas.species_names
        coeffs_ref = ref.reactant_stoich_coeffs
        coeffs_gas = gas.reactant_stoich_coeffs
        assert coeffs_gas.shape == coeffs_ref.shape
        assert (coeffs_gas == coeffs_ref).all()

        compositionA = [[ref.n_atoms(i,j) for j in range(ref.n_elements)]
                        for i in range(ref.n_species)]
        compositionB = [[gas.n_atoms(i,j) for j in range(gas.n_elements)]
                        for i in range(gas.n_species)]
        assert compositionA == compositionB

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
                assert ref_cp[i] == approx(gas_cp[i], rel=1e-7), 'cp' + message
                assert ref_h[i] == approx(gas_h[i], rel=1e-7), 'h' + message
                assert ref_s[i] == approx(gas_s[i], rel=1e-7), 's' + message

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
                assert ref_kf[i] == approx(gas_kf[i], rel=tol), 'kf' + message
                assert ref_kr[i] == approx(gas_kr[i], rel=tol), 'kr' + message

    @pytest.mark.slow_test
    def test_gri30(self):
        output = self.convert('gri30.inp', thermo='gri30_thermo.dat',
                              transport='gri30_tran.dat', output='gri30_test')

        ref, gas = self.checkConversion("gri30.yaml", output)
        self.checkKinetics(ref, gas, [300, 1500], [5e3, 1e5, 2e6])

    def test_soot(self):
        output = self.convert("soot.inp", thermo="soot-therm.dat", output="soot_test")
        ref, gas = self.checkConversion("soot.yaml", output)
        self.checkThermo(ref, gas, [300, 1100])
        self.checkKinetics(ref, gas, [300, 1100], [5e3, 1e5, 2e6])

    def test_pdep(self):
        output = self.convert('pdep-test.inp')
        ref, gas = self.checkConversion('pdep-test.yaml', output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_species_only(self):
        self.convert(None, thermo='dummy-thermo.dat', output='dummy-thermo')
        yaml = ("{phases: [{name: gas, species: "
                "[{dummy-thermo-from-ck.yaml/species: [R1A, R1B, P1]}], "
                "thermo: ideal-gas}]}")
        gas = ct.Solution(yaml=yaml)
        assert gas.n_species == 3
        assert gas.n_reactions == 0

    def test_missing_thermo(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('h2o2_missingThermo.inp')
        captured = self._capsys.readouterr()
        assert "No thermo data found for species 'H2'" in captured.out

    def test_duplicate_thermo(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('duplicate-thermo.inp')
        captured = self._capsys.readouterr()
        assert "Found additional thermo entry for species 'bar'" in captured.out
        assert "2 additional errors about redundant thermo data" in captured.out

        output = self.convert('duplicate-thermo.inp', permissive=True, quiet=False)
        captured = self._capsys.readouterr()
        assert "2 additional warnings about redundant thermo data" in captured.out
        gas = ct.Solution(output)
        assert gas.n_species == 4
        assert gas.n_reactions == 2

    def test_bad_elemental_composition(self):
        with pytest.raises(ck2yaml.InputError,
                           match="Error parsing elemental composition"):
            self.convert(None, thermo='bad-nasa7-composition.dat')

        captured = self._capsys.readouterr()
        assert "Error parsing elemental composition" in captured.out

    def test_bad_nasa7_temperature_ranges(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert(None, thermo='bad-nasa7-Trange.dat')

        captured = self._capsys.readouterr()
        assert "Only one temperature range defined" in captured.out

    def test_bad_nasa7_formatting(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert(None, thermo='bad-nasa7-formatting.dat')

        captured = self._capsys.readouterr()
        assert "Lines could not be parsed as a NASA7 entry" in captured.out
        assert "line 7" in captured.out  # bad entry for O2
        assert "line 15" in captured.out  # bad entry for H2
        assert "line 24" in captured.out  # bad entry for OH

    def test_bad_nasa7_formatting_warn(self):
        self.convert(None, thermo='bad-nasa7-formatting.dat', permissive=True)
        yaml = ("{phases: [{name: gas, species: "
                "[{bad-nasa7-formatting-from-ck.yaml/species: all}], "
                "thermo: ideal-gas}]}")
        gas = ct.Solution(yaml=yaml)
        assert set(gas.species_names) == {"O", "H", "H2O"}

    def test_bad_nasa7_float(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert(None, thermo='bad-nasa7-float.dat')

        captured = self._capsys.readouterr()
        assert "could not convert string to float: '-0.35101023e+ 0'" in captured.out

    def test_bad_nasa9_temperature_ranges(self):
        with pytest.raises(ck2yaml.InputError,
                           match="non-adjacent temperature ranges"):
            self.convert(None, thermo='bad-nasa9-Trange.dat')
        self._capsys.readouterr()

    def test_bad_nasa9_coeffs(self):
        with pytest.raises(ck2yaml.InputError) as info:
            self.convert(None, thermo='bad-nasa9-coeffs.dat')
        err_text = str(info.value)
        captured = self._capsys.readouterr()
        for out in err_text, captured.out:
            assert "Error while reading thermo entry" in out
            assert "\nALCL3" in out
            assert "could not convert string to float" in out

    def test_duplicate_species(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('duplicate-species.inp')
        captured = self._capsys.readouterr()
        assert "multiple declarations for species 'bar'" in captured.out
        assert "Suppressed 2 additional errors about redundant species" in captured.out

        output = self.convert('duplicate-species.inp', permissive=True, quiet=False)
        captured = self._capsys.readouterr()
        assert "Ignoring redundant declaration for species 'bar'" in captured.out
        assert "Suppressed 2 additional warnings" in captured.out

        gas = ct.Solution(output)
        assert gas.species_names == ['foo','bar','baz']

    def test_pathologicalSpeciesNames(self):
        output = self.convert('species-names.inp')
        gas = ct.Solution(output)

        assert gas.n_species == 10
        assert gas.species_name(0) == '(Parens)'
        assert gas.species_name(1) == '@#$%^-2'
        assert gas.species_index('co:lons:') == 2
        assert gas.species_name(3) == '[xy2]*{.}'
        assert gas.species_name(4) == 'plus+'
        assert gas.species_name(5) == 'eq=uals'
        assert gas.species_name(6) == 'plus'
        assert gas.species_name(7) == 'trans_butene'
        assert gas.species_name(8) == 'co'
        assert gas.species_name(9) == "amp&ersand"

        assert gas.n_reactions == 13
        nu = gas.product_stoich_coeffs - gas.reactant_stoich_coeffs
        assert list(nu[:,0]) == [-1, -1, 0, 2, 0, 0, 0, 0, 0, 0]
        assert list(nu[:,1]) == [-2, 3, 0, -1, 0, 0, 0, 0, 0, 0]
        assert list(nu[:,2]) == [-1, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        assert list(nu[:,3]) == [3, 0, 0, 0, -2, -1, 0, 0, 0, 0]
        assert list(nu[:,4]) == [2, 0, 0, 0, -1, 0, -1, 0, 0, 0]
        assert list(nu[:,5]) == [1, 0, 0, 0, 1, -1, -1, 0, 0, 0]
        assert list(nu[:,6]) == [2, 0, -1, 0, 0, -1, 0, 0, 0, 0]
        assert list(nu[:,7]) == [0, 0, 0, 0, -1, 1, 0, 0, 0, 0]
        assert list(nu[:,8]) == [0, 0, 0, 0, -1, 1, 0, 0, 0, 0]
        assert list(nu[:,9]) == [0, 0, 0, 0, -1, 1, 0, 0, 0, 0]
        assert list(nu[:,10]) == [0, 0, -1, 0, 2, 0, 0, -1, 0, 0]
        assert list(nu[:,11]) == [0, 0, -1, 0, 2, 0, 0, 0, -1, 0]
        assert list(nu[:,12]) == [0, 0, 0, 0, 1, 0, 0, 0, 0, -1]

    def test_unterminatedSections(self):
        output = self.convert('unterminated-sections.inp', quiet=False)
        captured = self._capsys.readouterr()
        for section in ["ELEMENTS", "SPECIES", "THERMO", "REACTIONS"]:
            assert f"{section} section implicitly ended" in captured.out
        gas = ct.Solution(output)
        assert gas.n_species == 3
        assert gas.n_reactions == 2

    def test_unterminatedSections2(self):
        output = self.convert('unterminated-sections2.inp', quiet=False)
        captured = self._capsys.readouterr()
        for section in ["ELEMENTS", "SPECIES", "THERMO", "REACTIONS"]:
            assert f"{section} section implicitly ended" in captured.out

        gas = ct.Solution(output)
        assert gas.n_species == 3
        assert gas.n_reactions == 2

    def test_unrecognized_section(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('unrecognized-section.inp', thermo='dummy-thermo.dat',
                         permissive=True)
        captured = self._capsys.readouterr()
        assert "unrecognized keyword 'SPAM'" in captured.out

    def test_nasa9(self):
        output = self.convert("nasa9-test.inp", thermo="nasa9-test-therm.dat")
        ref, gas = self.checkConversion("nasa9-test.yaml", output)
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_nasa9_subset(self):
        output = self.convert("nasa9-test-subset.inp", thermo="nasa9-test-therm.dat")
        ref, gas = self.checkConversion("nasa9-test-subset.yaml", output)
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])

    def test_nasa9_embedded(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert("nasa9-embedded.inp")
        captured = self._capsys.readouterr()
        assert "SPECIES section implicitly ended" in captured.out
        assert "THERMO NASA9 section implicitly ended" in captured.out
        assert "Found additional thermo entry for species 'AR'" in captured.out
        assert "Suppressed 2 additional errors" in captured.out

        output = self.convert("nasa9-embedded.inp", quiet=False, permissive=True)
        captured = self._capsys.readouterr()
        assert "Suppressed 2 additional warnings" in captured.out
        ref, gas = self.checkConversion("nasa9-embedded.yaml", output)
        self.checkThermo(ref, gas, [300, 500, 1200, 5000])
        self.checkKinetics(ref, gas, [300, 1200, 9000], [5e3, 1e5, 2e6])

    def test_sri_falloff(self):
        output = self.convert("sri-falloff.inp", thermo="dummy-thermo.dat")
        ref, gas = self.checkConversion("sri-falloff.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_chemically_activated(self):
        output = self.convert("chemically-activated-reaction.inp")
        ref, gas = self.checkConversion("chemically-activated-reaction.yaml",
                                        output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6, 1e7])

    def test_falloff_no_rate(self):
        with pytest.raises(ck2yaml.InputError, match="implies pressure dependence"):
            self.convert("falloff-error.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()

    def test_explicit_third_bodies(self):
        output = self.convert("explicit-third-bodies.inp", thermo="dummy-thermo.dat")
        ref, gas = self.checkConversion("explicit-third-bodies.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_explicit_reverse_rate(self):
        output = self.convert("explicit-reverse-rate.inp", thermo="dummy-thermo.dat")
        ref, gas = self.checkConversion("explicit-reverse-rate.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

        # Reactions with explicit reverse rate constants are transformed into
        # two irreversible reactions with reactants and products swapped, unless
        # the explicit reverse rate is zero so only the forward reaction is used.
        Rr = gas.reverse_rate_constants
        assert Rr[0] == 0.0
        assert Rr[1] == 0.0
        assert Rr[2] == 0.0
        assert Rr[3] == 0.0
        assert Rr[4] == 0.0
        Rstoich = gas.reactant_stoich_coeffs
        Pstoich = gas.product_stoich_coeffs
        assert list(Rstoich[:, 0]) == list(Pstoich[:, 1])
        assert list(Rstoich[:, 1]) == list(Pstoich[:, 0])
        assert list(Rstoich[:, 2]) == list(Pstoich[:, 3])
        assert list(Rstoich[:, 3]) == list(Pstoich[:, 2])

        assert gas.n_reactions == 5

    def test_explicit_forward_order(self):
        output = self.convert("explicit-forward-order.inp", thermo="dummy-thermo.dat")
        ref, gas = self.checkConversion("explicit-forward-order.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_negative_order(self):
        output = self.convert('negative-order.inp', thermo='dummy-thermo.dat',
                              quiet=False)
        captured = self._capsys.readouterr()
        assert "Negative reaction order" in captured.out
        ref, gas = self.checkConversion("negative-order.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_nonreactant_order(self):
        output = self.convert("nonreactant-order.inp", thermo="dummy-thermo.dat",
                              quiet=False)
        captured = self._capsys.readouterr()
        assert "Non-reactant order" in captured.out
        ref, gas = self.checkConversion("nonreactant-order.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_negative_A_factor(self):
        output = self.convert('negative-rate.inp', thermo='dummy-thermo.dat')
        gas = ct.Solution(output)  # Validate the mechanism
        assert gas.reaction(4).rate.pre_exponential_factor < 0
        assert gas.reaction(1).rate.pre_exponential_factor < 0
        assert gas.reaction(2).rate.pre_exponential_factor < 0
        assert gas.forward_rate_constants[5] < 0

    def test_bad_troe_value(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('bad-troe.inp', thermo='dummy-thermo.dat')
        captured = self._capsys.readouterr()
        assert "Error while reading reaction in bad-troe.inp" in captured.out
        assert "on line 10" in captured.out
        assert "R1A+R1B(+M) <=> H+P1(+M)" in captured.out
        assert "could not convert string to float" in captured.out

    def test_invalid_reaction_equation(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('invalid-equation.inp', thermo='dummy-thermo.dat')
        captured = self._capsys.readouterr()
        assert "Unparsable line: 'R1A <-> R1B" in captured.out

    def test_mismatched_third_body(self):
        with pytest.raises(ck2yaml.InputError, match="Third bodies do not match"):
            self.convert("mismatched-third-body.inp", thermo="gri30_thermo.dat")
        captured = self._capsys.readouterr()
        assert "Third bodies do not match: 'M' and 'AR'" in captured.out

    @pytest.mark.slow_test
    def test_reaction_units(self):
        out_def = self.convert('units-default.inp', thermo='dummy-thermo.dat')
        out_cus = self.convert('units-custom.inp', thermo='dummy-thermo.dat')
        default, custom = self.checkConversion(out_def, out_cus)
        self.checkKinetics(default, custom,
                           [300, 800, 1450, 2800], [5e0, 5e3, 1e5, 2e6, 1e8], 1e-7)

    def test_float_stoich_coeffs(self):
        output = self.convert('float-stoich.inp', thermo='dummy-thermo.dat')
        gas = ct.Solution(output)

        R = gas.reactant_stoich_coeffs
        P = gas.product_stoich_coeffs
        assert R[:,0] == approx([0, 1.5, 0.5, 0])
        assert P[:,0] == approx([1, 0, 0, 1])
        assert R[:,1] == approx([1, 0, 0, 1])
        assert P[:,1] == approx([0, 0.33, 1.67, 0])

    def test_unparsable_reaction(self):
        with pytest.raises(ck2yaml.InputError, match="Unparsable line"):
            self.convert("unparsable-reaction.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()

    def test_reaction_no_reactants(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert("reaction-no-reactants.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()
        assert "No reactant species found in reaction equation" in captured.out

    def test_undefined_species(self):
        with pytest.raises(ck2yaml.InputError, match="Unexpected token"):
            self.convert("undefined-species.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()
        assert "May be due to undeclared species 'XYZ'" in captured.out

    def test_bad_chebyshev_params(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert("bad-chebyshev-params.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()
        assert "Missing TCHEB entry for Chebyshev reaction" in captured.out
        assert "Missing PCHEB entry for Chebyshev reaction" in captured.out
        assert "Expected 3*5 = 15 but got 16." in captured.out

    def test_bad_parameters_multiple_types(self):
        with pytest.raises(ck2yaml.InputError,
                           match="contains parameters for more than one reaction type"):
            self.convert("bad-reaction-params.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()

    def test_photon(self):
        with pytest.raises(ck2yaml.InputError,
                           match="reversible reaction containing a product photon"):
            self.convert("photo-reaction.inp", thermo="dummy-thermo.dat")

        output = self.convert('photo-reaction.inp', thermo='dummy-thermo.dat',
                              permissive=True)
        captured = self._capsys.readouterr()

        ref, gas = self.checkConversion("photo-reaction.yaml", output)
        self.checkKinetics(ref, gas, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_photon_reactant_error(self):
        with pytest.raises(ck2yaml.InputError,
                           match="Reactant photon not supported"):
            self.convert("bad-photo-reaction.inp", thermo="dummy-thermo.dat")
        captured = self._capsys.readouterr()

    def test_transport_normal(self):
        output = self.convert('h2o2.inp', transport='gri30_tran.dat',
                              output='h2o2_transport_normal')

        gas = ct.Solution(output)
        gas.TPX = 300, 101325, 'H2:1.0, O2:1.0'
        assert gas.thermal_conductivity == approx(0.07663, rel=1e-4)

    def test_transport_embedded(self):
        output = self.convert('with-transport.inp')
        gas = ct.Solution(output)
        gas.X = [0.2, 0.3, 0.5]
        D = gas.mix_diff_coeffs
        for d in D:
            assert d > 0.0

    def test_transport_missing_species(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('h2o2.inp', transport='h2o2-missing-species-tran.dat',
                output='h2o2_transport_missing_species')
        captured = self._capsys.readouterr()
        assert "No transport data for species 'H2O2'" in captured.out

    def test_transport_extra_column_entries(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('h2o2.inp',
                transport='h2o2-extra-column-entries-tran.dat',
                output='h2o2_extra-column-entries-tran')
        captured = self._capsys.readouterr()
        assert "572.400" in captured.out
        assert "6 transport parameters were expected, but found 7" in captured.out

    def test_transport_duplicate_species(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('h2o2.inp',
                transport='h2o2-duplicate-species-tran.dat',
                output='h2o2_transport_duplicate_species')
        captured = self._capsys.readouterr()
        assert "Duplicate transport data for species 'H2O'" in captured.out
        assert "2 additional errors about duplicate transport data" in captured.out

        self.convert('h2o2.inp',
            transport='h2o2-duplicate-species-tran.dat',
            output='h2o2_transport_duplicate_species', permissive=True, quiet=False)
        captured = self._capsys.readouterr()
        assert "2 additional warnings about duplicate transport data" in captured.out

    def test_transport_bad_geometry(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('h2o2.inp',
                transport='h2o2-bad-geometry-tran.dat',
                output='h2o2_transport_bad_geometry')
        captured = self._capsys.readouterr()
        assert ("Invalid geometry flag value '3' for species 'HO2'. "
                "Flag value must be 0, 1, or 2.") in captured.out

        with pytest.raises(ck2yaml.InputError):
            self.convert('h2o2.inp',
                transport='h2o2-character-geometry-tran.dat',
                output='h2o2_transport_character_geometry')
        captured = self._capsys.readouterr()
        assert ("Invalid geometry flag 'a' for species 'HO2'. Flag must be an integer."
                in captured.out)

    def test_transport_float_geometry(self):
        # Runs but issues a log message
        self.convert('h2o2.inp',
            transport='h2o2-float-geometry-tran.dat',
            output='h2o2_transport_float_geometry', quiet=False)

        captured = self._capsys.readouterr()
        assert "Incorrect geometry flag" in captured.out
        assert "automatically converted" in captured.out

        # Runs with no logging message
        output = self.convert('h2o2.inp',
            transport='h2o2-float-geometry-tran.dat',
            output='h2o2_transport_float_geometry', permissive=True, quiet=False)
        captured = self._capsys.readouterr()
        assert "Incorrect geometry flag" not in captured.out

        gas = ct.Solution(output)
        assert gas.species("H").transport.geometry == 'atom'
        assert gas.species("H2").transport.geometry == 'linear'
        assert gas.species("H2O").transport.geometry == 'nonlinear'

        with pytest.raises(ck2yaml.InputError, match='Invalid float geometry flag'):
            self.convert('h2o2.inp',
                transport='h2o2-float-arithmetic-error-geometry-tran.dat',
                output='h2o2_transport_float_geometry', permissive=True)
        captured = self._capsys.readouterr()
        assert "0.999999" in captured.out

    def test_empty_reaction_section(self):
        output = self.convert('h2o2_emptyReactions.inp')
        gas = ct.Solution(output)
        assert gas.n_species == 9
        assert gas.n_reactions == 0

    def test_reaction_comments1(self):
        output = self.convert('pdep-test.inp')
        text = output.read_text()
        assert 'Generic mechanism header' in text
        assert 'Single PLOG reaction' in text
        assert 'Multiple PLOG expressions at the same pressure' in text

    def test_reaction_comments2(self):
        output = self.convert('explicit-third-bodies.inp', thermo='dummy-thermo.dat')
        text = output.read_text()
        assert 'An end of line comment' in text
        assert 'A comment after the last reaction' in text

    def test_custom_element(self):
        output = self.convert('custom-elements.inp')
        gas = ct.Solution(output)
        assert gas.n_elements == 4
        assert gas.atomic_weight(2) == approx(13.003)
        assert gas.n_atoms('ethane', 'C') == 2
        assert gas.n_atoms('CC', 'C') == 1
        assert gas.n_atoms('CC', 'Ci') == 1

    def test_surface_mech(self):
        output = self.convert('surface1-gas.inp', surface='surface1.inp',
                              output='surface1')

        surf = ct.Interface(output, 'PT_SURFACE')
        gas = surf.adjacent["gas"]

        assert gas.n_reactions == 11
        assert surf.n_reactions == 15
        assert surf.species('O2_Pt').size == 3

        # Different units for rate constants in each input file
        # 62.1 kJ/gmol = 6.21e7 J/kmol
        assert gas.reaction(0).rate.activation_energy == approx(6.21e7)
        # 67400 J/mol = 6.74e7 J/kmol
        assert surf.reaction(1).rate.activation_energy == approx(6.74e7)

        # Sticking coefficients
        assert surf.reaction(4).duplicate
        assert not isinstance(surf.reaction(1).rate, ct.StickingArrheniusRate)
        assert isinstance(surf.reaction(2).rate, ct.StickingArrheniusRate)
        assert surf.reaction(2).rate.motz_wise_correction
        assert isinstance(surf.reaction(4).rate, ct.StickingArrheniusRate)
        assert not surf.reaction(4).rate.motz_wise_correction
        assert surf.reaction(6).rate.motz_wise_correction

        # Coverage dependencies
        covdeps = surf.reaction(1).rate.coverage_dependencies
        assert len(covdeps) == 2
        assert "H_Pt" in covdeps
        assert covdeps["OH_Pt"]["m"] == 1.0
        assert covdeps["H_Pt"]["E"] == approx(-6e6) # 6000 J/gmol = 6e6 J/kmol

    def test_surface_mech2(self):
        output = self.convert('surface1-gas-noreac.inp', surface='surface1.inp',
                              output='surface1-nogasreac')

        gas = ct.Solution(output, 'gas')
        surf = ct.Interface(output, 'PT_SURFACE', [gas])

        assert gas.n_reactions == 0
        assert surf.n_reactions == 15

        covdeps = surf.reaction(1).rate.coverage_dependencies
        assert "H_Pt" in covdeps
        assert covdeps["OH_Pt"]["m"] == 1.0
        assert covdeps["H_Pt"]["E"] == approx(-6e6)

    def test_surf_bad_files(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('surface1.inp')

        captured = self._capsys.readouterr()
        assert "must be specified using the '--surface' option" in captured.out

        with pytest.raises(SystemExit):
            ck2yaml.main([
                f"--surface={self.test_data_path}/surface1.inp",
                f"--thermo={self.test_data_path}/dummy-thermo.dat"])
        captured = self._capsys.readouterr()
        assert "Cannot specify a surface mechanism without a gas phase" in captured.out

    def test_surface_mech3(self):
        # This tests the case where the thermo data for both the gas and surface are
        # combined in a file separate from the gas and surface definitions.

        output = self.convert('surface2-gas.inp', thermo='surface2-thermo.dat',
                              surface='surface2.inp', output='surface2')
        captured = self._capsys.readouterr()
        assert "SITE section implicitly ended" in captured.out
        surf = ct.Interface(output, 'PT_SURFACE')

        assert surf.n_species == 6
        assert surf.n_reactions ==  15
        assert surf.reaction(4).duplicate is True

    def test_missing_site_density(self):
        with pytest.raises(ck2yaml.InputError, match="no site density"):
            self.convert("surface1-gas.inp", surface="missing-site-density.inp")
        captured = self._capsys.readouterr()

    def test_bad_reaction_option(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert("surface1-gas.inp", surface="surface1-bad-option.inp")
        captured = self._capsys.readouterr()
        assert "Unrecognized token 'XYZ' on REACTIONS line" in captured.out

    def test_third_body_plus_falloff_reactions(self):
        output = self.convert("third_body_plus_falloff_reaction.inp")
        gas = ct.Solution(output)
        assert gas.n_reactions == 2

    def test_blank_line_in_header(self):
        output = self.convert("blank_line_in_header.inp")
        gas = ct.Solution(output)
        assert gas.n_reactions == 1

    def test_missing_input_files(self):
        with pytest.raises(IOError, match="Missing input file"):
            self.convert('nonexistent-file-813.inp')

        with pytest.raises(SystemExit):
            ck2yaml.main([
                f"--input=nonexistent-file-813.inp"])
        captured = self._capsys.readouterr()
        assert "Missing input file: 'nonexistent-file-813.inp'" in captured.out

    def test_extra(self):
        output = self.convert("h2o2.inp", output="h2o2_extra", extra="extra.yaml")
        yml = load_yaml(output)

        desc = yml['description'].split('\n')[-1]
        assert desc == 'This is an alternative description.'
        for key in ['foo', 'bar']:
            assert key in yml.keys()

    def test_extra_reserved(self):
        with pytest.raises(ck2yaml.InputError,
                           match="must not redefine reserved field.+'species'"):
            self.convert("h2o2.inp", output="h2o2_extra1", extra="extra-reserved.yaml")
        captured = self._capsys.readouterr()

    def test_extra_description_str(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert("h2o2.inp", output="h2o2_extra1",
                         extra="extra-bad-description.yaml")
        captured = self._capsys.readouterr()
        assert "alternate description" in captured.out
        assert "needs to be a string" in captured.out

    def test_sri_zero(self):
        # This test tests it can convert the SRI parameters when D or E equal to 0
        output = self.convert('sri_convert_test.txt')
        mech = load_yaml(output)
        D = mech['reactions'][0]['SRI']['D']
        E = mech['reactions'][0]['SRI']['E']
        assert D == 0
        assert E == 0

    def test_duplicate_reactions(self):
        # Running a test this way instead of using the convertMech function
        # tests the behavior of the ck2yaml.main function and the mechanism
        # validation step.

        with pytest.raises(SystemExit):
            ck2yaml.main([
                f"--input={self.test_data_path}/undeclared-duplicate-reactions.inp",
                f"--thermo={self.test_data_path}/dummy-thermo.dat",
                f"--output={self.test_work_path}/undeclared-duplicate-reactions.yaml"])

        captured = self._capsys.readouterr()
        for token in ('FAILED', 'Line 12', 'Line 14', 'R1A', 'R1B'):
            assert token in captured.out

    def test_multiple_duplicate_reactions(self):
        with pytest.raises(SystemExit):
            ck2yaml.main([
                f"--input={self.test_data_path}/undeclared-duplicate-reactions2.inp",
                f"--thermo={self.test_data_path}/dummy-thermo.dat",
                f"--output={self.test_work_path}/undeclared-duplicate-reactions2.yaml"])

        captured = self._capsys.readouterr()
        for token in ('FAILED', 'Line 12', 'Line 14', 'Line 11', 'Line 15',
                    'R1A + R1B', 'R3 + H'):
            assert token in captured.out

    def test_single_Tint(self):
        output = self.convert(None, thermo="thermo_single_Tint.dat",
                              output="thermo_single_Tint",
                              single_intermediate_temperature=True)
        mech = load_yaml(output)

        # Al(cr)
        thermo = mech["species"][0]["thermo"]
        assert thermo["temperature-ranges"] == [200.0, 933.61]
        assert len(thermo["data"]) == 1
        assert thermo["data"][0][0] == 1.01040191

        # AlBr3(L)
        thermo = mech["species"][1]["thermo"]
        assert thermo["temperature-ranges"] == [370.6, 5000.0]
        assert len(thermo["data"]) == 1
        assert thermo["data"][0][0] == 15.02975

        # AlF3(b)
        thermo = mech["species"][2]["thermo"]
        assert thermo["temperature-ranges"] == [728.0, 1000.0, 2523.0]
        assert len(thermo["data"]) == 2
        assert thermo["data"][1][0] == 10.41947

        # AlF3(L)
        thermo = mech["species"][3]["thermo"]
        assert thermo["temperature-ranges"] == [2523.0, 5000.0]
        assert len(thermo["data"]) == 1
        assert thermo["data"][0][0] == 15.096679

    def test_error_for_big_element_number(self):
        with pytest.raises(ck2yaml.InputError):
            self.convert('big_element_num_err.inp')

        captured = self._capsys.readouterr()
        assert "no more than 3 digits." in captured.out


class Testyaml2ck:
    """Test yaml2ck by converting to CK then back to YAML to read with Cantera."""
    ext: str = "-from-yaml2ck.yaml"

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, test_data_path):
        self.test_data_path = test_data_path

    def _convert_to_ck(
        self,
        input_file: Path,
        phase_name: str = "",
        output: tuple[str, str, str] | tuple = (),
    ) -> tuple[Path | None, Path | None, Path | None]:
        mechanism_path: Path | str
        if not output:
            stem = Path(input_file).stem  # strip '.inp'
            mechanism_path = self.test_work_path / (stem + "-from-yaml.ck")
            thermo_path = transport_path = None
        else:
            if len(output) != 3:
                raise ValueError(
                    "convert_to_ck output must be a tuple of length three "
                    "containing the mechanism, thermo, and transport file names."
                )
            mechanism_path, thermo_path, transport_path = output

        mech, thermo, transport = yaml2ck.convert(
            input_file,
            phase_name=phase_name,
            mechanism_path=mechanism_path,
            thermo_path=thermo_path,
            transport_path=transport_path,
            overwrite=True,
            sort_elements=None,
            sort_species=None
        )

        return mech, thermo, transport

    def convert(
        self,
        input_file: Path,
        phase_name: str = "",
        mech: str | Path | None = None,
        thermo: str | Path | None = None,
        transport: str | Path | None = None,
        permissive: bool = False,
    ) -> str:
        if mech is not None:
            mech, thermo, transport = self._convert_to_ck(
                input_file,
                phase_name,
                (mech, thermo, transport),
            )
        else:
            mech, thermo, transport = self._convert_to_ck(input_file, phase_name)

        output = self.test_work_path / (Path(input_file).stem + self.ext)
        ck2yaml.convert(
            mech,
            thermo_file=thermo,
            transport_file=transport,
            out_name=output,
            quiet=True,
            permissive=permissive,
        )
        return mech

    def check_conversion(self, basename, cls=ct.Solution, **kwargs):
        # The round-trip YAML->CK->YAML will always have the single phase name 'gas'
        # even if the input YAML phase has a different name
        if "name" in kwargs:
            phase_name = kwargs.pop("name")
        else:
            phase_name = ""
        ckname = self.test_work_path / (basename.stem + self.ext)
        ck_phase = cls(ckname, **kwargs)
        yaml_phase = cls(basename, phase_name, **kwargs)

        assert set(ck_phase.element_names) == set(yaml_phase.element_names)
        assert set(ck_phase.species_names) == set(yaml_phase.species_names)

        yamlSpecies = [yaml_phase.species(s) for s in ck_phase.species_names]
        for C, Y in zip(ck_phase.species(), yamlSpecies):
            assert C.composition == Y.composition

        assert ck_phase.n_reactions == yaml_phase.n_reactions
        for C, Y in zip(ck_phase.reactions(), yaml_phase.reactions()):
            assert C.__class__ == Y.__class__
            assert C.reactants == Y.reactants
            assert C.products == Y.products
            assert C.duplicate == Y.duplicate

        for i, sp in zip(range(ck_phase.n_reactions), ck_phase.kinetics_species_names):
            assert ck_phase.reactant_stoich_coeff(sp, i) == yaml_phase.reactant_stoich_coeff(sp, i)

        return ck_phase, yaml_phase

    def check_thermo(self, ck_phase, yaml_phase, temperatures, tol=1e-7):
        yaml_idx = {ck_phase.species_index(s): yaml_phase.species_index(s) for s in ck_phase.species_names}

        for T in temperatures:
            ck_phase.TP = T, ct.one_atm
            yaml_phase.TP = T, ct.one_atm
            cp_ck = ck_phase.partial_molar_cp
            cp_yaml = yaml_phase.partial_molar_cp
            h_ck = ck_phase.partial_molar_enthalpies
            h_yaml = yaml_phase.partial_molar_enthalpies
            s_ck = ck_phase.partial_molar_entropies
            s_yaml = yaml_phase.partial_molar_entropies
            assert ck_phase.density == approx(yaml_phase.density)
            for i in range(ck_phase.n_species):
                message = ' for species {0} at T = {1}'.format(i, T)
                assert cp_ck[i] == approx(cp_yaml[yaml_idx[i]], rel=tol), 'cp' + message
                assert h_ck[i] == approx(h_yaml[yaml_idx[i]], rel=tol), 'h' + message
                assert s_ck[i] == approx(s_yaml[yaml_idx[i]], rel=tol), 's' + message

    def check_kinetics(self, ck_phase, yaml_phase, temperatures, pressures, tol=1e-7):
        for T, P in itertools.product(temperatures, pressures):
            ck_phase.TP = T, P
            yaml_phase.TP = T, P
            kf_ck = ck_phase.forward_rate_constants
            kr_ck = ck_phase.reverse_rate_constants
            kf_yaml = yaml_phase.forward_rate_constants
            kr_yaml = yaml_phase.reverse_rate_constants
            for i in range(yaml_phase.n_reactions):
                message = f"for reaction {i+1}: {yaml_phase.reaction(i)} at T = {T}, P = {P}"
                assert kf_ck[i] == approx(kf_yaml[i], rel=tol), 'kf ' + message
                assert kr_ck[i] == approx(kr_yaml[i], rel=tol), 'kr ' + message

    def check_transport(self, ck_phase, yaml_phase, temperatures, model="mixture-averaged"):
        yaml_idx = {ck_phase.species_index(s): yaml_phase.species_index(s) for s in ck_phase.species_names}
        ck_phase.transport_model = model
        yaml_phase.transport_model = model
        for T in temperatures:
            ck_phase.TP = T, ct.one_atm
            yaml_phase.TP = T, ct.one_atm
            assert ck_phase.viscosity == approx(yaml_phase.viscosity)
            assert ck_phase.thermal_conductivity == approx(yaml_phase.thermal_conductivity)
            Dkm_ck = ck_phase.mix_diff_coeffs
            Dkm_yaml = yaml_phase.mix_diff_coeffs
            for i in range(ck_phase.n_species):
                message = 'dkm for species {0} at T = {1}'.format(i, T)
                assert Dkm_ck[i] == approx(Dkm_yaml[yaml_idx[i]]), message

    @pytest.mark.slow_test
    def test_gri30(self, cantera_data_path):
        input_file = cantera_data_path / "gri30.yaml"
        self.convert(input_file)
        X = {'O2': 0.3, 'H2': 0.1, 'CH4': 0.2, 'CO2': 0.4}
        ck_phase, yaml_phase = self.check_conversion(input_file)
        ck_phase.X = X
        yaml_phase.X = X
        self.check_thermo(ck_phase, yaml_phase, [300, 500, 1300, 2000])
        self.check_kinetics(ck_phase, yaml_phase, [900, 1800], [2e5, 20e5])
        self.check_transport(ck_phase, yaml_phase, [298, 1001, 2400])

    def test_nonreactant_orders(self):
        input_file = self.test_data_path / "reaction-orders.yaml"
        self.convert(input_file, permissive=True)
        ck_phase, yaml_phase = self.check_conversion(input_file)
        self.check_thermo(ck_phase, yaml_phase, [300, 500])
        self.check_kinetics(ck_phase, yaml_phase, [300, 1001, 2500], [1e5, 10e5])

    def test_phase_id(self, cantera_data_path):
        input_file = cantera_data_path / "nDodecane_Reitz.yaml"
        self.convert(input_file, "nDodecane_IG")
        ck_phase, yaml_phase = self.check_conversion(input_file, name="nDodecane_IG")
        ck_phase.X = "h2:1"
        yaml_phase.X = "h2:1"
        self.check_kinetics(
            ck_phase, yaml_phase, [300, 800, 1450, 2800], [5e3, 1e5, 2e6], tol=4e-6
        )

    def test_third_body_reactions(self):
        input_file = self.test_data_path / "explicit-third-bodies.yaml"
        mech = self.convert(input_file)
        with open(mech) as fid:
            lines = fid.readlines()
        for i, line in enumerate(lines):
            if line.startswith("R1A + R1B"):
                next = lines[i + 1]
                assert next.startswith("LOW") or next.strip() == "DUPLICATE"
        ck_phase, yaml_phase = self.check_conversion(input_file)
        self.check_kinetics(
            ck_phase, yaml_phase, [300, 800, 1450, 2800], [5e3, 1e5, 2e6]
        )

    def test_undeclared_third_body(self):
        input_file = self.test_data_path / "undeclared-third-body.yaml"
        mech = self.convert(input_file)
        with open(mech) as fid:
            lines = fid.readlines()
        for i, line in enumerate(lines):
            if line.startswith("H + O2 + M <=> HO2 + M"):
                next = lines[i + 1]
                assert "N2" not in next
                assert all(sp in next for sp in ["AR", "C2H6", "CO", "CO2", "H2O", "O2"])
        ck_phase, yaml_phase = self.check_conversion(input_file)
        self.check_kinetics(ck_phase, yaml_phase, [900, 1800], [2e5, 20e5])

    def test_pdep(self):
        input_file = self.test_data_path / "pdep-test.yaml"
        self.convert(input_file)
        ck_phase, yaml_phase = self.check_conversion(input_file)
        # Chebyshev coefficients in XML are truncated to 6 digits, limiting accuracy
        self.check_kinetics(ck_phase, yaml_phase, [300, 1000, 2200],
                            [100, ct.one_atm, 2e5, 2e6, 9.9e6], tol=2e-4)

    def test_sri_falloff(self):
        input_file = self.test_data_path / "sri-falloff.yaml"
        self.convert(input_file)
        ck_phase, yaml_phase = self.check_conversion(input_file)
        self.check_kinetics(ck_phase, yaml_phase, [300, 800, 1450, 2800], [5e3, 1e5, 2e6])

    def test_chemically_activated(self):
        input_file = self.test_data_path / "chemically-activated-reaction.yaml"
        self.convert(input_file)
        ck_phase, yaml_phase = self.check_conversion(input_file)
        # pre-exponential factor in XML is truncated to 7 sig figs, limiting accuracy
        self.check_kinetics(
            ck_phase, yaml_phase, [300, 800, 1450, 2800], [5e3, 1e5, 2e6, 1e7], tol=1e-7
        )

    def test_yaml_2_ck_reactions(self):
        input_file = self.test_data_path / "yaml-ck-reactions.yaml"
        self.convert(input_file)
        ck_phase, yaml_phase = self.check_conversion(input_file)
        X = {'O2': 0.3, 'H': 0.1, 'H2': 0.2, 'AR': 0.4}
        ck_phase.X = X
        yaml_phase.X = X
        self.check_thermo(ck_phase, yaml_phase, [300, 500, 1300, 2000])
        self.check_kinetics(ck_phase, yaml_phase, [900, 1800], [2e5, 20e5], tol=2e-7)
        self.check_transport(ck_phase, yaml_phase, [298, 1001, 2400])

    def test_phase_no_reactions(self):
        input_file = self.test_data_path / "yaml-ck-reactions.yaml"
        yaml_phase = ct.Solution(input_file, name="no-reactions")
        assert yaml_phase.n_reactions == 0

        ck_file = self.test_work_path / 'no-reactions.ck'
        ck_file.unlink(missing_ok=True)
        yaml_phase.write_chemkin(ck_file, quiet=True)
        assert ck_file.exists()

    def test_write_chemkin(self):
        # test alternative converter
        yaml_phase = ct.Solution('h2o2.yaml')
        ck_file = self.test_work_path / 'test.ck'
        ck_file.unlink(missing_ok=True)
        yaml_phase.write_chemkin(ck_file, quiet=True)
        yaml_phase.write_chemkin(
            ck_file, sort_species='alphabetical', overwrite=True, quiet=True)
        assert ck_file.exists()

        yaml_file = self.test_work_path / 'test.yaml'
        yaml_file.unlink(missing_ok=True)
        ck2yaml.convert(ck_file, out_name=yaml_file, quiet=True)
        assert yaml_file.exists()
        ck_phase = ct.Solution(yaml_file)

        X = {'O2': 0.3, 'H': 0.1, 'H2': 0.2, 'AR': 0.4}
        ck_phase.X = X
        yaml_phase.X = X
        self.check_thermo(ck_phase, yaml_phase, [300, 500, 1300, 2000])
        self.check_kinetics(ck_phase, yaml_phase, [900, 1800], [2e5, 20e5], tol=2e-7)
        self.check_transport(ck_phase, yaml_phase, [298, 1001, 2400])

    def test_write_notes(self):
        input_file = self.test_data_path / 'species-names.yaml'
        yaml_phase = ct.Solution(input_file)
        assert yaml_phase.species("eq=uals").input_data["thermo"]["note"] == 120521
        assert yaml_phase.species("plus").input_data["thermo"]["note"] == 12.05

        ck_file = self.test_work_path / 'species-names.ck'
        ck_file.unlink(missing_ok=True)
        yaml_phase.write_chemkin(ck_file, quiet=True)

        yaml_file = self.test_work_path / 'species-names.yaml'
        yaml_file.unlink(missing_ok=True)
        ck2yaml.convert(ck_file, out_name=yaml_file, quiet=True)
        assert yaml_file.exists()

        ck_phase = ct.Solution(yaml_file)
        assert ck_phase.species("eq=uals").input_data["thermo"]["note"] == "120521"
        assert ck_phase.species("plus").input_data["thermo"]["note"] == "12.05"

    def test_multi_line_notes(self):
        input_file = self.test_data_path / 'multi-line-notes.yaml'
        yaml_phase = ct.Solution(input_file)
        assert yaml_phase.species("H2").input_data["thermo"]["note"] == "Line 1\nLine 2"

        ck_file = self.test_work_path / 'multi-line-notes.ck'
        ck_file.unlink(missing_ok=True)
        yaml_phase.write_chemkin(ck_file, quiet=True)

        yaml_file = self.test_work_path / 'multi-line-notes.yaml'
        yaml_file.unlink(missing_ok=True)
        ck2yaml.convert(ck_file, out_name=yaml_file, quiet=True)
        assert yaml_file.exists()

        ck_phase = ct.Solution(yaml_file)
        assert ck_phase.species("H2").input_data["thermo"]["note"] == "Line 1\nLine 2"


class Testcti2yaml:

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, test_data_path):
        self.test_data_path = test_data_path

    def convert(self, basename, src_dir=None, encoding=None):
        if src_dir is None:
            src_dir = self.test_data_path

        cti2yaml.convert(
            filename=Path(src_dir) / f"{basename}.cti",
            output_name=self.test_work_path / f"{basename}-from-cti.yaml",
            encoding=encoding,
        )

    def checkConversion(self, basename, cls=ct.Solution, ctiphases=(),
                        yamlphases=(), **kwargs):
        ctiPhase = cls(f"{basename}-from-cti.yaml", adjacent=ctiphases, **kwargs)
        yamlPhase = cls(f"{basename}.yaml", adjacent=yamlphases, **kwargs)

        assert ctiPhase.element_names == yamlPhase.element_names
        assert ctiPhase.species_names == yamlPhase.species_names
        assert ctiPhase.n_reactions == yamlPhase.n_reactions
        for C, Y in zip(ctiPhase.species(), yamlPhase.species()):
            assert C.composition == Y.composition

        for C, Y in zip(ctiPhase.reactions(), yamlPhase.reactions()):
            assert C.__class__ == Y.__class__
            assert C.reactants == Y.reactants
            assert C.products == Y.products
            assert C.duplicate == Y.duplicate

        for i, sp in zip(range(ctiPhase.n_reactions), ctiPhase.kinetics_species_names):
            assert ctiPhase.reactant_stoich_coeff(sp, i) == yamlPhase.reactant_stoich_coeff(sp, i)

        return ctiPhase, yamlPhase

    def checkThermo(self, ctiPhase, yamlPhase, temperatures, tol=1e-7, check_cp=True):
        for T in temperatures:
            ctiPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            if check_cp:
                cp_cti = ctiPhase.partial_molar_cp
                cp_yaml = yamlPhase.partial_molar_cp
            else:
                with pytest.raises(NotImplementedError):
                    yamlPhase.partial_molar_cp
            h_cti = ctiPhase.partial_molar_enthalpies
            h_yaml = yamlPhase.partial_molar_enthalpies
            s_cti = ctiPhase.partial_molar_entropies
            s_yaml = yamlPhase.partial_molar_entropies
            assert ctiPhase.density == approx(yamlPhase.density)
            for i in range(ctiPhase.n_species):
                message = ' for species {0} at T = {1}'.format(i, T)
                if check_cp:
                    cp_cti[i] == approx(cp_yaml[i], rel=tol), 'cp' + message
                assert h_cti[i] == approx(h_yaml[i], rel=tol), 'h' + message
                assert s_cti[i] == approx(s_yaml[i], rel=tol), 's' + message

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
                assert kf_cti[i] == approx(kf_yaml[i], rel=tol), 'kf ' + message
                assert kr_cti[i] == approx(kr_yaml[i], rel=tol), 'kr ' + message

    def checkTransport(self, ctiPhase, yamlPhase, temperatures,
                       model='mixture-averaged'):
        ctiPhase.transport_model = model
        yamlPhase.transport_model = model
        for T in temperatures:
            ctiPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            assert ctiPhase.viscosity == approx(yamlPhase.viscosity)
            assert ctiPhase.thermal_conductivity == approx(yamlPhase.thermal_conductivity)
            Dkm_cti = ctiPhase.mix_diff_coeffs
            Dkm_yaml = yamlPhase.mix_diff_coeffs
            for i in range(ctiPhase.n_species):
                message = 'dkm for species {0} at T = {1}'.format(i, T)
                assert Dkm_cti[i] == approx(Dkm_yaml[i]), message

    @pytest.mark.slow_test
    def test_gri30(self):
        self.convert("gri30")
        ctiPhase, yamlPhase = self.checkConversion('gri30')
        X = {'O2': 0.3, 'H2': 0.1, 'CH4': 0.2, 'CO2': 0.4}
        ctiPhase.X = X
        yamlPhase.X = X
        self.checkThermo(ctiPhase, yamlPhase, [300, 500, 1300, 2000])
        self.checkKinetics(ctiPhase, yamlPhase, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctiPhase, yamlPhase, [298, 1001, 2400])

    def test_pdep(self):
        self.convert("pdep-test")
        ctiPhase, yamlPhase = self.checkConversion('pdep-test')
        # Agreement limited by low precision used by ck2cti for Chebyshev coeffs
        self.checkKinetics(ctiPhase, yamlPhase, [300, 1000, 2200],
                           [100, ct.one_atm, 2e5, 2e6, 9.9e6], tol=2e-4)

    def test_ptcombust(self):
        self.convert("ptcombust")
        ctiSurf, yamlSurf = self.checkConversion("ptcombust", ct.Interface,
            name="Pt_surf")
        yamlGas = yamlSurf.adjacent["gas"]
        ctiGas = ctiSurf.adjacent["gas"]

        self.checkKinetics(ctiGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctiSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctiSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    @pytest.mark.slow_test
    def test_ptcombust_motzwise(self):
        self.convert("ptcombust-motzwise")
        ctiSurf, yamlSurf = self.checkConversion("ptcombust-motzwise", ct.Interface,
            name="Pt_surf")
        yamlGas = yamlSurf.adjacent["gas"]
        ctiGas = ctiSurf.adjacent["gas"]

        self.checkKinetics(ctiGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctiSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctiSurf, yamlSurf, [900], [101325])

    def test_sofc(self):
        self.convert("sofc")
        cti_tpb, yaml_tpb = self.checkConversion("sofc", ct.Interface, name="tpb")
        ctiMetal, ctiMSurf, ctiOSurf = cti_tpb.adjacent.values()
        yamlMetal, yamlMSurf, yamlOSurf = yaml_tpb.adjacent.values()

        assert "oxide_bulk" in ctiOSurf.adjacent
        assert "gas" in ctiOSurf.adjacent

        self.checkThermo(ctiMSurf, yamlMSurf, [900, 1000, 1100])
        self.checkThermo(ctiOSurf, yamlOSurf, [900, 1000, 1100])
        ctiMetal.electric_potential = yamlMetal.electric_potential = 2
        self.checkKinetics(cti_tpb, yaml_tpb, [900, 1000, 1100], [1e5])
        ctiMetal.electric_potential = yamlMetal.electric_potential = 4
        self.checkKinetics(cti_tpb, yaml_tpb, [900, 1000, 1100], [1e5])

    @pytest.mark.slow_test
    def test_liquidvapor(self):
        self.convert("liquidvapor")
        for name in ["water", "nitrogen", "methane", "hydrogen", "oxygen", "heptane"]:
            ctiPhase, yamlPhase = self.checkConversion("liquidvapor", name=name)
            self.checkThermo(ctiPhase, yamlPhase,
                             [1.3 * ctiPhase.min_temp, 0.7 * ctiPhase.max_temp])

    def test_Redlich_Kwong_CO2(self):
        self.convert("co2_RK_example")
        ctiGas, yamlGas = self.checkConversion('co2_RK_example')
        for P in [1e5, 2e6, 1.3e7]:
            yamlGas.TP = ctiGas.TP = 300, P
            self.checkThermo(ctiGas, yamlGas, [300, 400, 500], check_cp=False)

    def test_diamond(self):
        self.convert("diamond")
        ctiSurf, yamlSurf = self.checkConversion("diamond", ct.Interface,
            name="diamond_100")
        ctiSolid = ctiSurf.adjacent["diamond"]
        yamlSolid = yamlSurf.adjacent["diamond"]
        self.checkThermo(ctiSolid, yamlSolid, [300, 500])
        self.checkThermo(ctiSurf, yamlSurf, [330, 490])
        self.checkKinetics(ctiSurf, yamlSurf, [400, 800], [2e5])

    def test_lithium_ion_battery(self):
        name = 'lithium_ion_battery'
        self.convert(name, encoding="utf-8")
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
        self.convert("ch4_ion")
        ctiGas, yamlGas = self.checkConversion("ch4_ion")
        self.checkThermo(ctiGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctiGas, yamlGas, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctiGas, yamlGas, [298, 1001, 2400])

    def test_description(self):
        self.convert("haca2")
        ctiGas, yamlGas = self.checkConversion("haca2")
        assert ctiGas.input_header["description"].startswith("HACA Mechanism")
        assert yamlGas.input_header["description"].startswith("HACA Mechanism")

    def test_nonreactant_orders(self):
        self.convert("reaction-orders")
        ctiGas, yamlGas = self.checkConversion("reaction-orders")
        assert ctiGas.input_header["description"].startswith("Input file to test")
        self.checkThermo(ctiGas, yamlGas, [300, 500])
        self.checkKinetics(ctiGas, yamlGas, [300, 1001, 2500], [1e5, 10e5])


class Testctml2yaml:

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, test_data_path):
        self.test_data_path = test_data_path

    def convert(self, basename, src_dir=None):
        if src_dir is None:
            src_dir = self.test_data_path

        ctml2yaml.convert(
            Path(src_dir) / f"{basename}.xml",
            self.test_work_path / f"{basename}-from-xml.yaml",
        )

    def checkConversion(self, basename, cls=ct.Solution, ctmlphases=(),
                        yamlphases=(), **kwargs):
        ctmlPhase = cls(f"{basename}-from-xml.yaml", adjacent=ctmlphases, **kwargs)
        yamlPhase = cls(f"{basename}.yaml", adjacent=yamlphases, **kwargs)

        assert ctmlPhase.element_names == yamlPhase.element_names
        assert ctmlPhase.species_names == yamlPhase.species_names
        assert ctmlPhase.n_reactions == yamlPhase.n_reactions
        for C, Y in zip(ctmlPhase.species(), yamlPhase.species()):
            assert C.composition == Y.composition

        for C, Y in zip(ctmlPhase.reactions(), yamlPhase.reactions()):
            assert C.__class__ == Y.__class__
            assert C.reactants == Y.reactants
            assert C.products == Y.products
            assert C.duplicate == Y.duplicate

        for i, sp in zip(range(ctmlPhase.n_reactions), ctmlPhase.kinetics_species_names):
            assert (ctmlPhase.reactant_stoich_coeff(sp, i)
                    == yamlPhase.reactant_stoich_coeff(sp, i))

        return ctmlPhase, yamlPhase

    def checkThermo(self, ctmlPhase, yamlPhase, temperatures, pressure=ct.one_atm,
                    tol=1e-7, check_cp=True):
        for T in temperatures:
            ctmlPhase.TP = T, pressure
            yamlPhase.TP = T, pressure
            if check_cp:
                cp_ctml = ctmlPhase.partial_molar_cp
                cp_yaml = yamlPhase.partial_molar_cp
            else:
                with pytest.raises(NotImplementedError):
                    yamlPhase.partial_molar_cp
            h_ctml = ctmlPhase.partial_molar_enthalpies
            h_yaml = yamlPhase.partial_molar_enthalpies
            s_ctml = ctmlPhase.partial_molar_entropies
            s_yaml = yamlPhase.partial_molar_entropies
            assert ctmlPhase.density == approx(yamlPhase.density)
            for i in range(ctmlPhase.n_species):
                message = ' for species {0} at T = {1}'.format(ctmlPhase.species_names[i], T)
                if check_cp:
                    assert cp_ctml[i] == approx(cp_yaml[i], rel=tol), 'cp' + message
                assert h_ctml[i] == approx(h_yaml[i], rel=tol), 'h' + message
                assert s_ctml[i] == approx(s_yaml[i], rel=tol), 's' + message

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
                assert kf_ctml[i] == approx(kf_yaml[i], rel=tol), 'kf ' + message
                assert kr_ctml[i] == approx(kr_yaml[i], rel=tol), 'kr ' + message

    def checkTransport(self, ctmlPhase, yamlPhase, temperatures,
                       model='mixture-averaged'):
        ctmlPhase.transport_model = model
        yamlPhase.transport_model = model
        for T in temperatures:
            ctmlPhase.TP = T, ct.one_atm
            yamlPhase.TP = T, ct.one_atm
            assert ctmlPhase.viscosity == approx(yamlPhase.viscosity)
            assert ctmlPhase.thermal_conductivity ==  approx(yamlPhase.thermal_conductivity)
            Dkm_ctml = ctmlPhase.mix_diff_coeffs
            Dkm_yaml = yamlPhase.mix_diff_coeffs
            for i in range(ctmlPhase.n_species):
                message = 'dkm for species {0} at T = {1}'.format(i, T)
                assert Dkm_ctml[i] == approx(Dkm_yaml[i]), message

    @pytest.mark.slow_test
    def test_gri30(self):
        self.convert("gri30")
        ctmlPhase, yamlPhase = self.checkConversion('gri30')
        X = {'O2': 0.3, 'H2': 0.1, 'CH4': 0.2, 'CO2': 0.4}
        ctmlPhase.X = X
        yamlPhase.X = X
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlPhase, yamlPhase, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctmlPhase, yamlPhase, [298, 1001, 2400])

    def test_pdep(self):
        self.convert("pdep-test")
        ctmlPhase, yamlPhase = self.checkConversion('pdep-test')
        # Chebyshev coefficients in XML are truncated to 6 digits, limiting accuracy
        self.checkKinetics(ctmlPhase, yamlPhase, [300, 1000, 2200],
                           [100, ct.one_atm, 2e5, 2e6, 9.9e6], tol=2e-4)

    def test_ptcombust(self):
        self.convert("ptcombust")
        ctmlSurf, yamlSurf = self.checkConversion("ptcombust", ct.Interface,
            name="Pt_surf")
        ctmlGas = ctmlSurf.adjacent["gas"]
        yamlGas = yamlSurf.adjacent["gas"]

        self.checkKinetics(ctmlGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctmlSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctmlSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    def test_ptcombust_motzwise(self):
        self.convert("ptcombust-motzwise")
        ctmlGas, yamlGas = self.checkConversion('ptcombust-motzwise')
        ctmlSurf, yamlSurf = self.checkConversion('ptcombust-motzwise', ct.Interface,
            name='Pt_surf', ctmlphases=[ctmlGas], yamlphases=[yamlGas])

        self.checkKinetics(ctmlGas, yamlGas, [500, 1200], [1e4, 3e5])
        self.checkThermo(ctmlSurf, yamlSurf, [400, 800, 1600])
        self.checkKinetics(ctmlSurf, yamlSurf, [500, 1200], [1e4, 3e5])

    def test_sofc(self):
        self.convert("sofc")
        ctml_tpb, yaml_tpb = self.checkConversion("sofc", ct.Interface, name="tpb")
        ctmlMetal, ctmlMSurf, ctmlOSurf = ctml_tpb.adjacent.values()
        yamlMetal, yamlMSurf, yamlOSurf = yaml_tpb.adjacent.values()

        assert "oxide_bulk" in ctmlOSurf.adjacent
        assert "gas" in  ctmlOSurf.adjacent

        self.checkThermo(ctmlMSurf, yamlMSurf, [900, 1000, 1100])
        self.checkThermo(ctmlOSurf, yamlOSurf, [900, 1000, 1100])
        ctmlMetal.electric_potential = yamlMetal.electric_potential = 2
        self.checkKinetics(ctml_tpb, yaml_tpb, [900, 1000, 1100], [1e5])
        ctmlMetal.electric_potential = yamlMetal.electric_potential = 4
        self.checkKinetics(ctml_tpb, yaml_tpb, [900, 1000, 1100], [1e5])

    def test_liquidvapor(self):
        self.convert("liquidvapor")
        for name in ["water", "nitrogen", "methane", "hydrogen", "oxygen", "heptane"]:
            ctmlPhase, yamlPhase = self.checkConversion("liquidvapor", name=name)
            self.checkThermo(ctmlPhase, yamlPhase,
                             [1.3 * ctmlPhase.min_temp, 0.7 * ctmlPhase.max_temp])

    def test_Redlich_Kwong_CO2(self):
        self.convert("co2_RK_example")
        ctmlGas, yamlGas = self.checkConversion('co2_RK_example')
        for P in [1e5, 2e6, 1.3e7]:
            yamlGas.TP = ctmlGas.TP = 300, P
            self.checkThermo(ctmlGas, yamlGas, [300, 400, 500], check_cp=False)

    def test_diamond(self):
        self.convert("diamond")
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
        self.convert(name)
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
        self.convert("noxNeg")
        ctmlGas, yamlGas = self.checkConversion('noxNeg')
        self.checkThermo(ctmlGas, yamlGas, [300, 1000])
        self.checkKinetics(ctmlGas, yamlGas, [300, 1000], [1e5])

    def test_ch4_ion(self):
        self.convert("ch4_ion")
        ctmlGas, yamlGas = self.checkConversion("ch4_ion")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])
        self.checkTransport(ctmlGas, yamlGas, [298, 1001, 2400])

    def test_nasa9(self):
        self.convert("nasa9-test")
        ctmlGas, yamlGas = self.checkConversion("nasa9-test")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])

    def test_chemically_activated(self):
        self.convert("chemically-activated-reaction")
        ctmlGas, yamlGas = self.checkConversion("chemically-activated-reaction")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_explicit_forward_order(self):
        self.convert("explicit-forward-order")
        ctmlGas, yamlGas = self.checkConversion("explicit-forward-order")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        # Accuracy limited by precision of ck2cti
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5], tol=2e-7)

    def test_explicit_reverse_rate(self):
        self.convert("explicit-reverse-rate")
        ctmlGas, yamlGas = self.checkConversion("explicit-reverse-rate")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_explicit_third_bodies(self):
        self.convert("explicit-third-bodies")
        ctmlGas, yamlGas = self.checkConversion("explicit-third-bodies")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_fractional_stoich_coeffs(self):
        self.convert("frac")
        ctmlGas, yamlGas = self.checkConversion("frac")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_water_IAPWS95_thermo(self):
        self.convert("liquid-water")
        ctmlWater, yamlWater = self.checkConversion("liquid-water")
        self.checkThermo(ctmlWater, yamlWater, [300, 500, 1300, 2000], pressure=22064000.0)
        assert ctmlWater.transport_model == yamlWater.transport_model
        ctmlWater.TP = yamlWater.TP = 300, 22064000.0
        dens = ctmlWater.density
        for T in [298, 1001, 2400]:
            ctmlWater.TD = T, dens
            yamlWater.TD = T, dens
            assert ctmlWater.viscosity == approx(yamlWater.viscosity)
            assert ctmlWater.thermal_conductivity == approx(yamlWater.thermal_conductivity)

    def test_hmw_nacl_phase(self):
        basename = "HMW_NaCl_sp1977_alt"
        self.convert(basename)
        ctmlPhase, yamlPhase = self.checkConversion(basename)
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_NaCl_solid_phase(self):
        self.convert("NaCl_Solid")
        ctmlPhase, yamlPhase = self.checkConversion("NaCl_Solid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500, 1300, 2000])

    def test_DH_NaCl_phase(self):
        self.convert("debye-huckel-all")
        for name in [
            "debye-huckel-dilute",
            "debye-huckel-B-dot-ak",
            "debye-huckel-B-dot-a",
            "debye-huckel-pitzer-beta_ij",
            "debye-huckel-beta_ij",
        ]:
            ctmlPhase, yamlPhase = self.checkConversion("debye-huckel-all", name=name)
            self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_Redlich_Kister(self):
        self.convert("RedlichKisterVPSSTP_valid")
        ctmlPhase, yamlPhase = self.checkConversion("RedlichKisterVPSSTP_valid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_species_names(self):
        self.convert("species-names")
        ctmlGas, yamlGas = self.checkConversion('species-names')
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])

    def test_sri_falloff_reaction(self):
        self.convert("sri-falloff")
        ctmlGas, yamlGas = self.checkConversion("sri-falloff")
        self.checkThermo(ctmlGas, yamlGas, [300, 500, 1300, 2000])
        self.checkKinetics(ctmlGas, yamlGas, [900, 1800], [2e5, 20e5])

    def test_vpss_and_hkft(self):
        self.convert("pdss_hkft")
        ctmlPhase, yamlPhase = self.checkConversion("pdss_hkft")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_lattice_solid(self):
        self.convert("Li7Si3_ls")
        ctmlPhase, yamlPhase = self.checkConversion("Li7Si3_ls",
                                                    name="Li7Si3_and_Interstitials(S)")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_margules(self):
        self.convert("LiKCl_liquid")
        ctmlPhase, yamlPhase = self.checkConversion("LiKCl_liquid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_idealsolidsoln(self):
        with pytest.warns(UserWarning, match="SolidKinetics type is not implemented"):
            self.convert("IdealSolidSolnPhaseExample")

        # SolidKinetics is not implemented, so can't create a Kinetics class instance.
        basename = "IdealSolidSolnPhaseExample"
        ctmlPhase = ct.ThermoPhase(basename + "-from-xml.yaml")
        yamlPhase = ct.ThermoPhase(basename + ".yaml")

        assert ctmlPhase.element_names == yamlPhase.element_names
        assert ctmlPhase.species_names == yamlPhase.species_names
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_idealmolalsoln(self):
        self.convert("IdealMolalSolnPhaseExample")
        ctmlPhase, yamlPhase = self.checkConversion("IdealMolalSolnPhaseExample")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_transport_models(self):
        self.convert("transport_models_test")
        for name in ["UnityLewis", "CK_Mix", "CK_Multi", "HighP"]:
            ctmlPhase, yamlPhase = self.checkConversion("transport_models_test", name=name)
            self.checkTransport(ctmlPhase, yamlPhase, [298, 1001, 2500])

    def test_nonreactant_orders(self):
        self.convert("reaction-orders")
        ctmlPhase, yamlPhase = self.checkConversion("reaction-orders")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])
        self.checkKinetics(ctmlPhase, yamlPhase, [300, 1001, 2500], [1e5, 10e5])

    def test_species_ss_temperature_polynomials(self):
        self.convert("Li_Liquid")
        ctmlPhase, yamlPhase = self.checkConversion("Li_Liquid")
        self.checkThermo(ctmlPhase, yamlPhase, [300, 500])

    def test_duplicate_section_ids(self):
        with pytest.warns(UserWarning, match="Duplicate 'speciesData' id"):
            self.convert("duplicate-speciesData-ids")
        with pytest.warns(UserWarning, match="Duplicate 'reactionData' id"):
            self.convert("duplicate-reactionData-ids")

class Testlxcat2yaml:

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, test_data_path):
        self.test_data_path = test_data_path

    def convert(self, inputFile=None, database=None, mechFile=None, phase=None,
                insert=True, output=None):
        if inputFile is not None:
            inputFile = self.test_data_path / inputFile
        if mechFile is not None:
            mechFile = self.test_data_path / mechFile
        if output is None:
            output = Path(inputFile).stem  # strip '.xml'
        # output to work dir
        output = self.test_work_path / output

        lxcat2yaml.convert(inputFile, database, mechFile, phase, insert, output)
        return output

    def test_mechanism_with_lxcat(self):
        # get Solution from the mechanism file
        phase = "isotropic-electron-energy-plasma"
        mechFile = "lxcat-test-convert.yaml"
        gas1 = ct.Solution(self.test_data_path / mechFile,
                           phase=phase, transport_model=None)

        # get a stand-alone collisions
        standAloneFile = "stand-alone-lxcat.yaml"
        self.convert(inputFile='lxcat-test-convert.xml', database="test",
                     mechFile=mechFile, insert=False,
                     output=standAloneFile)

        # add collisions to the reaction list
        rxn_list = ct.Reaction.list_from_file(self.test_work_path / standAloneFile,
                                       gas1, section="collisions")
        for R in rxn_list:
            gas1.add_reaction(R)

        # get Solution from the output file
        output = "output-lxcat.yaml"
        self.convert(inputFile="lxcat-test-convert.xml", database="test",
                     mechFile=mechFile, insert=True,
                     output=output)
        gas2 = ct.Solution(self.test_work_path / output,
                           phase=phase, transport_model=None)

        # check number of reactions
        assert gas1.n_reactions == gas2.n_reactions == 4
        for i in range(1, gas1.n_reactions):
            assert (gas1.reaction(i).rate.energy_levels
                    == approx(gas2.reaction(i).rate.energy_levels))
            assert (gas1.reaction(i).rate.cross_sections
                    == approx(gas2.reaction(i).rate.cross_sections))

    def test_stand_alone_lxcat(self):
        outfile = "stand-alone-lxcat-without-mech.yaml"
        self.convert(inputFile='lxcat-test-convert.xml',
                     database="test", insert=False,
                     output=outfile)

        # get Solution from the mechanism file
        phase = "isotropic-electron-energy-plasma"
        mechFile = "lxcat-test-convert.yaml"
        gas = ct.Solution(self.test_data_path / mechFile,
                           phase=phase, transport_model=None)

        # add collisions to the reaction list
        rxn_list = ct.Reaction.list_from_file(self.test_work_path / outfile,
                                              gas, section="collisions")

        # verify the data
        assert len(rxn_list) == 2
        assert rxn_list[0].equation == "O2 + e => O2(Total-Ionization)+ + 2 e"
        assert rxn_list[0].reaction_type == "electron-collision-plasma"
        assert rxn_list[0].rate.energy_levels == approx([15., 20.])
        assert rxn_list[0].rate.cross_sections == approx([0.0, 5.5e-22])

        assert rxn_list[1].equation == "O2 + e => O2-"
        assert rxn_list[1].reaction_type == "electron-collision-plasma"
        assert rxn_list[1].rate.energy_levels == approx([0.0, 1.0])
        assert rxn_list[1].rate.cross_sections == approx([0.0, 1.0e-22])
