import gc
import numpy as np
import pytest
from pytest import approx

import cantera as ct


# Workaround to support both Numpy 1.x and 2.4.0+
# TODO: Replace when dropping Numpy 1.x support
trapezoid = getattr(np, "trapezoid", None) or np.trapz

@pytest.fixture(scope='function')
def setup_thermo_phase_tests(request):
    request.cls.phase = ct.Solution('h2o2.yaml', transport_model=None)

@pytest.mark.usefixtures('setup_thermo_phase_tests')
class TestThermoPhase:

    def test_source(self):
        assert self.phase.source == 'h2o2.yaml'

    def test_input_header(self):
        extra = self.phase.input_header
        assert extra["description"].startswith("Hydrogen-Oxygen submechanism")
        assert extra["cantera-version"] == "2.5.0"
        assert extra["generator"] == "ck2yaml"

    def test_missing_phases_key(self):
        yaml = '''
        species:
        - name: H
          composition: {H: 1}
          thermo:
            model: NASA7
            temperature-ranges: [200.0, 1000.0, 6000.0]
            data:
            - [2.5, 0.0, 0.0, 0.0, 0.0, 2.547366e+04, -0.44668285]
            - [2.5, 0.0, 0.0, 0.0, 0.0, 2.547366e+04, -0.44668285]
            note: L6/94
        '''
        with pytest.raises(ct.CanteraError, match="Key 'phases' not found"):
            _ = ct.Solution(yaml=yaml)

    def test_deprecated_phase(self):
        yaml = """
        phases:
        - name: phasename
          thermo: ideal-gas
          species: [{h2o2.yaml/species: all}]
          deprecated: This phase is deprecated because I said so.
        """
        with pytest.raises(ct.CanteraError, match="(?s)phasename.*said so"):
            ct.Solution(yaml=yaml)

    def test_base_attributes(self):
        assert isinstance(self.phase.name, str)
        assert isinstance(self.phase.phase_of_matter, str)
        assert isinstance(self.phase.thermo_model, str)
        assert isinstance(self.phase.kinetics_model, str)
        assert isinstance(self.phase.transport_model, str)
        assert isinstance(self.phase.composite, tuple)
        assert len(self.phase.composite) == 3
        assert self.phase.composite == (self.phase.thermo_model,
                                        self.phase.kinetics_model,
                                        self.phase.transport_model)
        self.phase.name = 'spam'
        assert self.phase.name == 'spam'
        with pytest.raises(AttributeError):
            self.phase.type = 'eggs'

    def test_phases(self):
        assert self.phase.n_phases == 1
        assert self.phase.phase_of_matter == "gas"

    def test_states(self):
        assert self.phase._native_state == ('T', 'D', 'Y')
        assert 'TPY' in self.phase._full_states.values()
        assert 'TD' in self.phase._partial_states.values()
        assert self.phase._native_mode == "TDY"

    def test_species(self):
        assert self.phase.n_species == 10
        for i, name in enumerate(self.phase.species_names):
            assert name == self.phase.species_name(i)
            assert i == self.phase.species_index(name)
            assert i == self.phase.species_index(i)
        with pytest.raises(ct.CanteraError,
                           match='IndexError thrown by Phase::checkSpeciesIndex'):
            self.phase.species(self.phase.n_species)

    def test_elements(self):
        assert self.phase.n_elements == 4
        for i, symbol in enumerate(self.phase.element_names):
            assert symbol == self.phase.element_name(i)
            assert i == self.phase.element_index(symbol)
            assert i == self.phase.element_index(i)

    def test_n_atoms(self):
        data = [(1, 'O', 'O'), (2, 'O', 'O2'), (1, b'H', b'OH'),
                (2, 'H', 'H2O'), (2, 'O', 'H2O2'), (1, 'Ar', 'AR'),
                (0, 'O', 'H'), (0, 'H', 'AR'), (0, 'Ar', 'HO2')]
        for (n, elem, species) in data:
            assert self.phase.n_atoms(species, elem) == n
            mElem = self.phase.element_index(elem)
            kSpec = self.phase.species_index(species)
            assert self.phase.n_atoms(kSpec, mElem) == n

        with pytest.raises(ct.CanteraError, match="Species 'C' not found"):
            self.phase.n_atoms('C', 'H2')
        with pytest.raises(ct.CanteraError, match="Element 'CH4' not found"):
            self.phase.n_atoms('H', 'CH4')

    def test_elemental_mass_fraction(self):
        self.phase.Y = 'H2O:0.5, O2:0.5'
        Zo = self.phase.elemental_mass_fraction('O')
        Zh = self.phase.elemental_mass_fraction('H')
        Zar = self.phase.elemental_mass_fraction('Ar')

        mO = self.phase.element_index('O')
        assert Zo == self.phase.elemental_mass_fraction(mO)
        assert Zo == approx(0.5 + 0.5 * (15.999 / 18.015))
        assert Zh == approx(0.5 * (2.016 / 18.015))
        assert Zar == 0.0

        with pytest.raises(ct.CanteraError, match="Element 'C' not found"):
            self.phase.elemental_mass_fraction('C')
        with pytest.raises(ct.CanteraError, match="outside valid range"):
            self.phase.elemental_mass_fraction(5)

    def test_elemental_mole_fraction(self):
        self.phase.X = 'H2O:0.5, O2:0.5'
        Zo = self.phase.elemental_mole_fraction('O')
        Zh = self.phase.elemental_mole_fraction('H')
        Zar = self.phase.elemental_mole_fraction('Ar')

        mO = self.phase.element_index('O')
        assert Zo == self.phase.elemental_mole_fraction(mO)
        assert Zo == approx((0.5 + 1) / (0.5 * 3 + 0.5 * 2))
        assert Zh == approx((2 * 0.5) / (0.5 * 3 + 0.5 * 2))
        assert Zar == 0.0

        with pytest.raises(ct.CanteraError, match="Element 'C' not found"):
            self.phase.elemental_mole_fraction('C')
        with pytest.raises(ct.CanteraError, match="outside valid range"):
            self.phase.elemental_mole_fraction(5)

    def test_elemental_mass_mole_fraction(self):
        # expected relationship between elemental mass and mole fractions
        comps = ['H2O:0.5, O2:0.5', 'H2:0.1, O2:0.4, H2O2:0.3, AR:0.2',
                 'O2:0.1, H2:0.9']
        for comp in comps:
            self.phase.X = comp

            denom = sum(self.phase.elemental_mole_fraction(i)
                        * self.phase.atomic_weight(i)
                        for i in range(self.phase.n_elements))

            for i in range(self.phase.n_elements):
                assert self.phase.elemental_mass_fraction(i) == approx(
                       self.phase.elemental_mole_fraction(i)
                       * self.phase.atomic_weight(i) / denom)

    def test_weights(self):
        atomic_weights = self.phase.atomic_weights
        molecular_weights = self.phase.molecular_weights
        assert self.phase.n_elements == len(atomic_weights)
        assert self.phase.n_species == len(molecular_weights)

        for i, mw in enumerate(molecular_weights):
            test_weight = 0.0
            for j, aw in enumerate(atomic_weights):
                test_weight += aw * self.phase.n_atoms(i, j)
            assert test_weight == approx(mw)

    def test_charges(self):
        gas = ct.Solution('ch4_ion.yaml')
        charges = gas.charges
        test = {'E': -1., 'N2': 0., 'H3O+': 1.}
        for species, charge in test.items():
            assert species in gas.species_names
            index = gas.species_index(species)
            assert charges[index] == charge

    def test_setComposition(self):
        X = np.zeros(self.phase.n_species)
        X[2] = 1.0
        self.phase.X = X
        Y = self.phase.Y
        assert list(X) == list(Y)

    def test_setComposition_singleton(self):
        X = np.zeros((1, self.phase.n_species, 1))
        X[0, 2, 0] = 1.0
        self.phase.X = X
        Y = self.phase.Y
        assert list(X[0, :, 0]) == list(Y)

    def test_setCompositionString(self):
        self.phase.X = 'h2:1.0, o2:1.0'
        X = self.phase.X
        assert X[0] == approx(0.5)
        assert X[3] == approx(0.5)

        with pytest.raises(ct.CanteraError, match="Species 'CO2' not found"):
            self.phase.X = 'H2:1.0, CO2:1.5'

    def test_setCompositionStringBad(self):
        X0 = self.phase.X
        with pytest.raises(ct.CanteraError, match='Trouble processing'):
            self.phase.X = 'H2:1.0, O2:asdf'
        assert X0 == approx(self.phase.X)

        with pytest.raises(ct.CanteraError, match='Trouble processing'):
            self.phase.X = 'H2:1e-x4'
        assert X0 == approx(self.phase.X)

        with pytest.raises(ct.CanteraError, match='decimal point in exponent'):
            self.phase.X = 'H2:1e-1.4'
        assert X0 == approx(self.phase.X)

        with pytest.raises(ct.CanteraError, match='Duplicate key'):
            self.phase.X = 'H2:0.5, O2:1.0, H2:0.1'
        assert X0 == approx(self.phase.X)

    def test_setCompositionDict(self):
        self.phase.X = {b'H2': 1.0, b'O2': 3.0}
        X = self.phase.X
        assert X[0] == approx(0.25)
        assert X[3] == approx(0.75)

        self.phase.Y = {'H2': 1.0, 'O2': 3.0}
        Y = self.phase.Y
        assert Y[0] == approx(0.25)
        assert Y[3] == approx(0.75)

    def test_getCompositionDict(self):
        self.phase.X = 'oh:1e-9, O2:0.4, AR:0.6'
        assert len(self.phase.mole_fraction_dict(1e-7)) == 2
        assert len(self.phase.mole_fraction_dict()) == 3

        self.phase.Y = 'O2:0.4, AR:0.6'
        Y1 = self.phase.mass_fraction_dict()
        assert Y1['O2'] == approx(0.4)
        assert Y1['AR'] == approx(0.6)

    def test_setCompositionNoNorm(self):
        X = np.zeros(self.phase.n_species)
        X[2] = 1.0
        X[0] = 0.01
        self.phase.set_unnormalized_mole_fractions(X)
        assert self.phase.X == approx(X)
        assert sum(X) == approx(1.01)

        Y = np.zeros(self.phase.n_species)
        Y[2] = 1.0
        Y[0] = 0.01
        self.phase.set_unnormalized_mass_fractions(Y)
        assert self.phase.Y == approx(Y)
        assert sum(Y) == approx(1.01)

    def test_setCompositionNoNormBad(self):
        X = np.zeros(self.phase.n_species - 1)
        with pytest.raises(ValueError, match='incorrect length'):
            self.phase.set_unnormalized_mole_fractions(X)

        with pytest.raises(ValueError, match='incorrect length'):
            self.phase.set_unnormalized_mass_fractions([1, 2, 3])

    def test_setCompositionDict_bad1(self):
        with pytest.raises(ct.CanteraError, match="Species 'HCl' not found"):
            self.phase.X = {'H2': 1.0, 'HCl': 3.0}

    def test_setCompositionDict_bad2(self):
        with pytest.raises(TypeError):
            self.phase.Y = {'H2': 1.0, 'O2': 'xx'}

    def test_setCompositionSlice(self):
        self.phase['h2', 'o2'].X = 0.1, 0.9
        X = self.phase.X
        assert X[0] == approx(0.1)
        assert X[3] == approx(0.9)

    def test_setCompositionSingleSpecies(self):
        yaml_def = """
            phases:
            - name: gas
              species: [{h2o2.yaml/species: [AR]}]
              thermo: ideal-gas
        """
        gas = ct.Solution(yaml=yaml_def)
        gas.X = [1]
        gas.Y = np.array([[1.001]])
        assert gas.Y[0] == 1.0

    def test_setCompositionSlice_bad(self):
        X0 = self.phase.X
        with pytest.raises(ValueError, match='incorrect length'):
            self.phase['H2', 'O2'].Y = [0.1, 0.2, 0.3]
        assert self.phase.X == approx(X0)

    def test_setCompositionEmpty_bad(self):
        X0 = self.phase.X
        with pytest.raises(ValueError, match='incorrect length'):
            self.phase.Y = np.array([])
        assert self.phase.X == approx(X0)

    @pytest.mark.slow_test
    def test_set_equivalence_ratio_stoichiometric(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        for fuel in ('C2H6', 'H2:0.7, CO:0.3', 'NH3:0.4, CH3OH:0.6'):
            for oxidizer in ('O2:1.0, N2:3.76', 'H2O2:1.0'):
                gas.set_equivalence_ratio(1.0, fuel, oxidizer)
                gas.equilibrate('TP')
                # Almost everything should end up as CO2, H2O and N2
                assert sum(gas['H2O','CO2','N2'].X) > 0.999999

    def test_set_equivalence_ratio_lean(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        excess = 0
        for phi in np.linspace(0.9, 0, 5):
            gas.set_equivalence_ratio(phi, 'CH4:0.8, CH3OH:0.2', 'O2:1.0, N2:3.76')
            gas.equilibrate('TP')
            assert gas['O2'].X[0] > excess
            excess = gas['O2'].X[0]
        assert sum(gas['O2','N2'].X) == approx(1.0)

    def test_set_equivalence_ratio_sulfur(self):
        sulfur_species = [k for k in ct.Species.list_from_file("nasa_gas.yaml")
                          if k.name in ("SO", "SO2")]
        gas = ct.Solution(thermo="ideal-gas",
                          species=ct.Species.list_from_file("gri30.yaml") + sulfur_species)
        fuel = "CH3:0.5, SO:0.25, OH:0.125, N2:0.125"
        ox = "O2:0.5, SO2:0.25, CO2:0.125, CH:0.125"

        def test_sulfur_results(gas, fuel, ox, basis):
            gas.set_equivalence_ratio(2.0, fuel, ox, basis)
            Z = gas.mixture_fraction(fuel, ox, basis)
            assert gas.stoich_air_fuel_ratio(fuel, ox, basis)/((1.0-Z)/Z) == approx(2.0)
            gas.set_mixture_fraction(Z, fuel, ox, basis)
            assert gas['SO2'].X[0] == approx(31.0/212.0)
            assert gas['O2'].X[0] == approx(31.0/106.0)
            assert gas['SO'].X[0] == approx(11.0/106.0)
            assert gas['CO2'].X[0] == approx(31.0/424.0)
            assert gas['CH3'].X[0] == approx(11.0/53.0)
            assert gas['N2'].X[0] == approx(11.0/212.0)
            assert gas['CH'].X[0] == approx(31.0/424.0)
            assert gas['OH'].X[0] == approx(11.0/212.0)
            assert gas.equivalence_ratio(fuel, ox, basis) == approx(2.0)

        test_sulfur_results(gas, fuel, ox, 'mole')

        gas.TPX = None, None, fuel
        fuel = gas.Y
        gas.TPX = None, None, ox
        ox = gas.Y
        test_sulfur_results(gas, fuel, ox, 'mass')

    def test_equivalence_ratio(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        for phi in np.linspace(0.5, 2.0, 5):
            gas.set_equivalence_ratio(phi, 'CH4:0.8, CH3OH:0.2', 'O2:1.0, N2:3.76')
            assert phi == approx(gas.equivalence_ratio('CH4:0.8, CH3OH:0.2', 'O2:1.0, N2:3.76'))
        # Check sulfur species
        sulfur_species = [k for k in ct.Species.list_from_file("nasa_gas.yaml")
                          if k.name in ("SO", "SO2")]
        gas = ct.Solution(thermo="ideal-gas", kinetics="gas",
                          species=ct.Species.list_from_file("gri30.yaml") + sulfur_species)
        for phi in np.linspace(0.5, 2.0, 5):
            gas.set_equivalence_ratio(phi, 'CH3:0.5, SO:0.25, OH:0.125, N2:0.125', 'O2:0.5, SO2:0.25, CO2:0.125')
            assert phi == approx(gas.equivalence_ratio('CH3:0.5, SO:0.25, OH:0.125, N2:0.125', 'O2:0.5, SO2:0.25, CO2:0.125'))
        gas.X = 'CH4:1' # pure fuel
        assert gas.equivalence_ratio() == np.inf

    def test_get_set_equivalence_ratio_functions(self):
        fuel = "CH4:0.2,O2:0.02,N2:0.1,CO:0.05,CO2:0.02"
        ox = "O2:0.21,N2:0.79,CO:0.04,CH4:0.01,CO2:0.03"

        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 300, 1e5, fuel
        Y_Cf = gas.elemental_mass_fraction("C")
        Y_Of = gas.elemental_mass_fraction("O")
        gas.TPX = 300, 1e5, ox
        Y_Co = gas.elemental_mass_fraction("C")
        Y_Oo = gas.elemental_mass_fraction("O")

        def test_equil_results(gas, fuel, ox, Y_Cf, Y_Of, Y_Co, Y_Oo, basis):
            gas.TP = 300, 1e5
            gas.set_equivalence_ratio(1.3, fuel, ox, basis)
            T = gas.T

            # set mixture to burnt state to make sure that equivalence ratio and
            # mixture fraction are independent of reaction progress
            gas.equilibrate("HP")

            phi = gas.equivalence_ratio(fuel, ox, basis)
            phi_loc = gas.equivalence_ratio()
            mf = gas.mixture_fraction(fuel, ox, basis)
            mf_C = gas.mixture_fraction(fuel, ox, basis, element="C")
            mf_O = gas.mixture_fraction(fuel, ox, basis, element="O")
            l = gas.stoich_air_fuel_ratio(fuel, ox, basis)

            gas.set_mixture_fraction(mf, fuel,ox, basis)
            phi2 = gas.equivalence_ratio(fuel, ox, basis)

            assert phi == approx(1.3)
            assert phi2 == approx(1.3)
            assert phi_loc == approx(1.1726068608195617)
            assert mf == approx(0.13415725911057605)
            assert mf_C == approx((gas.elemental_mass_fraction("C")-Y_Co)/(Y_Cf-Y_Co))
            assert mf_O == approx((gas.elemental_mass_fraction("O")-Y_Oo)/(Y_Of-Y_Oo))
            assert l == approx(8.3901204498353561)
            assert gas.P == approx(1e5)
            assert T == approx(300.0)

        test_equil_results(gas, fuel, ox, Y_Cf, Y_Of, Y_Co, Y_Oo, 'mole')

        # do the same for mass-based functions

        gas.TPX = None, None, fuel
        fuel = gas.Y
        Y_Cf = gas.elemental_mass_fraction("C")
        Y_Of = gas.elemental_mass_fraction("O")
        gas.TPX = None, None, ox
        ox = gas.Y
        Y_Co = gas.elemental_mass_fraction("C")
        Y_Oo = gas.elemental_mass_fraction("O")

        test_equil_results(gas, fuel, ox, Y_Cf, Y_Of, Y_Co, Y_Oo, 'mass')

    def test_equivalence_ratio_simple_dilution(self):
        gas = ct.Solution("gri30.yaml")

        phi = 2
        inv_afr = 2 * phi # inverse molar AFR for H2/O2
        X = "H2:4,O2:1,CO:3,CO2:4,N2:5,CH4:6"
        T, P = 300, 1e5
        gas.TPX = T, P, X
        original_X = gas.X
        assert gas.equivalence_ratio(include_species=["H2", "O2"]) == approx(2)
        assert gas.equivalence_ratio("H2", "O2",
                                     include_species=["H2", "O2"]) == approx(2)
        assert gas.T == approx(T)
        assert gas.P == approx(P)
        assert gas.X == approx(original_X)

        def test_simple_dilution(fraction, basis):
            if isinstance(fraction, str):
                fraction_dict = {fraction[:fraction.find(":")]:
                                 float(fraction[fraction.find(":")+1:])}
            else:
                fraction_dict = fraction

            fraction_type  = list(fraction_dict.keys())[0]
            fraction_value = float(list(fraction_dict.values())[0])

            M_H2 = gas.molecular_weights[gas.species_index("H2")]
            M_O2 = gas.molecular_weights[gas.species_index("O2")]

            gas.TP = T, P
            gas.set_equivalence_ratio(phi, "H2", "O2", fraction=fraction,
                                      diluent="CO2", basis=basis)
            if basis == "mole" and fraction_type == "diluent":
                assert gas["H2"].X[0] == approx((1 - fraction_value)
                                                * inv_afr / (inv_afr + 1))
                assert gas["O2"].X[0] == approx((1 - fraction_value) / (inv_afr + 1))
                assert gas["CO2"].X[0] == approx(fraction_value)
            elif basis == "mass" and fraction_type == "diluent":
                assert gas["H2"].Y[0] / gas["O2"].Y[0] == approx(inv_afr * M_H2 / M_O2)
                assert gas["CO2"].Y[0] == approx(fraction_value)
            elif basis == "mole" and fraction_type == "fuel":
                assert gas["H2"].X[0] == approx(fraction_value)
                assert gas["O2"].X[0] == approx(fraction_value / inv_afr)
                assert gas["CO2"].X[0] == approx(1 - fraction_value * (1 + 1 / inv_afr))
            elif basis == "mass" and fraction_type == "fuel":
                assert gas["H2"].Y[0] == approx(fraction_value)
                assert gas["H2"].Y[0] / gas["O2"].Y[0] == approx(inv_afr * M_H2 / M_O2)
            elif basis == "mole" and fraction_type == "oxidizer":
                assert gas["H2"].X[0] == approx(fraction_value * inv_afr)
                assert gas["O2"].X[0] == approx(fraction_value)
                assert gas["CO2"].X[0] == approx(1 - fraction_value * (1 + inv_afr))
            elif basis == "mass" and fraction_type == "oxidizer":
                assert gas["O2"].Y[0] == approx(fraction_value)
                assert gas["H2"].Y[0] / gas["O2"].Y[0] == approx(inv_afr * M_H2 / M_O2)

            Y = gas.Y
            assert phi == approx(gas.equivalence_ratio("H2", "O2", include_species=["H2", "O2"], basis=basis))
            assert phi == approx(gas.equivalence_ratio(include_species=["H2", "O2"], basis=basis))
            assert Y == approx(gas.Y)
            assert gas.T == approx(T)
            assert gas.P == approx(P)

        # brute force all possible input combinations
        test_simple_dilution("diluent:0.3", "mole")
        test_simple_dilution({"diluent": 0.3}, "mole")
        test_simple_dilution("diluent:0.3", "mass")
        test_simple_dilution({"diluent": 0.3}, "mass")

        test_simple_dilution("fuel:0.1", "mole")
        test_simple_dilution({"fuel": 0.1}, "mole")
        test_simple_dilution("fuel:0.1", "mass")
        test_simple_dilution({"fuel": 0.1}, "mass")

        test_simple_dilution("oxidizer:0.1", "mole")
        test_simple_dilution({"oxidizer": 0.1}, "mole")
        test_simple_dilution("oxidizer:0.1", "mass")
        test_simple_dilution({"oxidizer": 0.1}, "mass")

    def test_equivalence_ratio_arbitrary_dilution(self):
        fuel = "CH4:1,O2:0.01,N2:0.1,CO:0.05,CO2:0.02"
        oxidizer = "O2:0.8,N2:0.2,CO:0.01,CH4:0.005,CO2:0.03"
        diluent = "CO2:1,CO:0.025,N2:0.3,CH4:0.07,O2:0.09"
        gas = ct.Solution("gri30.yaml")

        gas.X = fuel
        X_fuel = gas.X
        gas.Y = fuel
        Y_fuel = gas.Y
        gas.X = oxidizer
        X_oxidizer = gas.X
        gas.Y = oxidizer
        Y_oxidizer = gas.Y
        gas.X = diluent
        X_diluent = gas.X
        gas.Y = diluent
        Y_diluent = gas.Y

        phi = 2
        fraction = 0.6
        gas.set_equivalence_ratio(phi, fuel, oxidizer)
        X_Mix = gas.X
        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"diluent": fraction},
                                  diluent=diluent)
        X_expected = X_Mix * (1 - fraction) + fraction * X_diluent
        assert gas.X == approx(X_expected)

        gas.set_equivalence_ratio(phi, fuel, oxidizer, basis="mass")
        Y_Mix = gas.Y
        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"diluent": fraction},
                                  diluent=diluent, basis="mass")
        Y_expected = Y_Mix * (1 - fraction) + fraction * Y_diluent
        assert gas.Y == approx(Y_expected)

        phi = 0.8
        fraction = 0.05
        gas.set_equivalence_ratio(phi, fuel, oxidizer, basis="mass")
        AFR = gas.stoich_air_fuel_ratio(fuel, oxidizer, basis="mass") / phi
        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"fuel": fraction},
                                  diluent=diluent, basis="mass")
        Y_expected = fraction * (Y_fuel + AFR * Y_oxidizer) \
                     + (1 - fraction * (1 + AFR)) * Y_diluent
        assert gas.Y == approx(Y_expected)

        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"oxidizer": fraction},
                                  diluent=diluent, basis="mass")
        Y_expected = fraction * (Y_fuel / AFR + Y_oxidizer) \
                     + (1 - fraction * (1 + 1 / AFR)) * Y_diluent
        assert gas.Y == approx(Y_expected)

        gas.X = fuel
        M_fuel = gas.mean_molecular_weight
        gas.X = oxidizer
        M_oxidizer = gas.mean_molecular_weight

        gas.set_equivalence_ratio(phi, fuel, oxidizer)
        AFR = M_fuel / M_oxidizer * gas.stoich_air_fuel_ratio(fuel, oxidizer) / phi

        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"fuel": fraction},
                                  diluent=diluent)
        X_expected = fraction * (X_fuel + AFR * X_oxidizer) \
                     + (1 - fraction * (1 + AFR)) * X_diluent
        assert gas.X == approx(X_expected)

        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"oxidizer": fraction},
                                  diluent=diluent)
        X_expected = fraction * (X_fuel / AFR + X_oxidizer) \
                     + (1 - fraction * (1 + 1 / AFR)) * X_diluent
        assert gas.X == approx(X_expected)

    def test_full_report(self):
        report = self.phase.report(threshold=0.0)
        assert self.phase.name in report
        assert 'temperature' in report
        assert 'minor' not in report
        for name in self.phase.species_names:
            assert name in report

    def test_default_report(self):
        self.phase.X = 'H2:0.1, O2:0.9, HO2:1e-10, H2O2:1e-20'
        report = self.phase.report()
        assert 'minor' in report
        for name in (' H2 ', ' O2 ', ' HO2 '):
            assert name in report
        for name in (' H2O2 ', ' OH ', ' AR '):
            assert name not in report

    def test_name(self):
        self.phase.name = 'something'
        assert self.phase.name == 'something'
        assert 'something' in self.phase.report()

    def test_badLength(self):
        X = np.zeros(5)
        with pytest.raises(ValueError, match='incorrect length'):
            self.phase.X = X
        with pytest.raises(ValueError, match='incorrect length'):
            self.phase.Y = X

    def test_mass_basis(self):
        assert self.phase.basis == 'mass'
        assert self.phase.density_mass == approx(self.phase.density)
        assert self.phase.enthalpy_mass == approx(self.phase.h)
        assert self.phase.entropy_mass == approx(self.phase.s)
        assert self.phase.int_energy_mass == approx(self.phase.u)
        assert self.phase.volume_mass == approx(self.phase.v)
        assert self.phase.cv_mass == approx(self.phase.cv)
        assert self.phase.cp_mass == approx(self.phase.cp)

    def test_molar_basis(self):
        self.phase.basis = 'molar'
        assert self.phase.basis == 'molar'
        assert self.phase.density_mole == approx(self.phase.density)
        assert self.phase.enthalpy_mole == approx(self.phase.h)
        assert self.phase.entropy_mole == approx(self.phase.s)
        assert self.phase.int_energy_mole == approx(self.phase.u)
        assert self.phase.volume_mole == approx(self.phase.v)
        assert self.phase.cv_mole == approx(self.phase.cv)
        assert self.phase.cp_mole == approx(self.phase.cp)

    def check_setters(self, T1, rho1, Y1):
        T0, rho0, Y0 = self.phase.TDY
        self.phase.TDY = T1, rho1, Y1
        X1 = self.phase.X
        P1 = self.phase.P
        h1 = self.phase.h
        s1 = self.phase.s
        u1 = self.phase.u
        v1 = self.phase.v

        def check_state(T, rho, Y):
            assert self.phase.T == approx(T)
            assert self.phase.Te == approx(T)
            assert self.phase.density == approx(rho)
            assert self.phase.Y == approx(Y)

        self.phase.TDY = T0, rho0, Y0
        self.phase.TPY = T1, P1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.UVY = u1, v1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.HPY = h1, P1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SPY = s1, P1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.TPX = T1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.UVX = u1, v1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.HPX = h1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SPX = s1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SVX = s1, v1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SVY = s1, v1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.DPX = rho1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.DPY = rho1, P1, Y1
        check_state(T1, rho1, Y1)

    def test_setState_mass(self):
        self.check_setters(T1 = 500.0, rho1 = 1.5,
                           Y1 = [0.1, 0.0, 0.0, 0.1, 0.4, 0.2, 0.0, 0.0, 0.2, 0.0])

    def test_setState_mole(self):
        self.phase.basis = 'molar'
        self.check_setters(T1 = 750.0, rho1 = 0.02,
                           Y1 = [0.2, 0.1, 0.0, 0.3, 0.1, 0.0, 0.0, 0.2, 0.1, 0.0])

    def test_setters_hold_constant(self):
        props = ('T','P','s','h','u','v','X','Y')
        pairs = [('TP', 'T', 'P'), ('SP', 's', 'P'),
                 ('UV', 'u', 'v')]

        self.phase.TDX = 1000, 1.5, 'H2O:0.1, O2:0.95, AR:3.0'
        values = {}
        for p in props:
            values[p] = getattr(self.phase, p)

        for pair, first, second in pairs:
            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            first_val = getattr(self.phase, first)
            second_val = getattr(self.phase, second)

            setattr(self.phase, pair, (values[first], None))
            assert getattr(self.phase, first) == approx(values[first])
            assert getattr(self.phase, second) == approx(second_val)

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair, (None, values[second]))
            assert getattr(self.phase, first) == approx(first_val)
            assert getattr(self.phase, second) == approx(values[second])

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'X', (None, None, values['X']))
            assert getattr(self.phase, first) == approx(first_val)
            assert getattr(self.phase, second) == approx(second_val)

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'Y', (None, None, values['Y']))
            assert getattr(self.phase, first) == approx(first_val)
            assert getattr(self.phase, second) == approx(second_val)

    def test_setter_errors(self):
        with pytest.raises(TypeError):
            self.phase.TD = 400

        with pytest.raises(AssertionError, match='incorrect number'):
            self.phase.TP = 300, 101325, 'CH4:1.0'

        with pytest.raises(AssertionError, match='incorrect number'):
            self.phase.HPY = 1.2e6, 101325

        with pytest.raises(AssertionError, match='incorrect number'):
            self.phase.UVX = -4e5, 4.4, 'H2:1.0', -1

    def test_invalid_property(self):
        x = self.phase
        with pytest.raises(AttributeError):
            x.foobar = 300
        with pytest.raises(AttributeError):
            x.foobar

    def check_getters(self):
        T,D,X = self.phase.TDX
        assert T == approx(self.phase.T)
        assert D == approx(self.phase.density)
        assert X == approx(self.phase.X)

        T,D,Y = self.phase.TDY
        assert T == approx(self.phase.T)
        assert D == approx(self.phase.density)
        assert Y == approx(self.phase.Y)

        T,D = self.phase.TD
        assert T == approx(self.phase.T)
        assert D == approx(self.phase.density)

        T,P,X = self.phase.TPX
        assert T == approx(self.phase.T)
        assert P == approx(self.phase.P)
        assert X == approx(self.phase.X)

        T,P,Y = self.phase.TPY
        assert T == approx(self.phase.T)
        assert P == approx(self.phase.P)
        assert Y == approx(self.phase.Y)

        T,P = self.phase.TP
        assert T == approx(self.phase.T)
        assert P == approx(self.phase.P)

        H,P,X = self.phase.HPX
        assert H == approx(self.phase.h)
        assert P == approx(self.phase.P)
        assert X == approx(self.phase.X)

        H,P,Y = self.phase.HPY
        assert H == approx(self.phase.h)
        assert P == approx(self.phase.P)
        assert Y == approx(self.phase.Y)

        H,P = self.phase.HP
        assert H == approx(self.phase.h)
        assert P == approx(self.phase.P)

        U,V,X = self.phase.UVX
        assert U == approx(self.phase.u)
        assert V == approx(self.phase.v)
        assert X == approx(self.phase.X)

        U,V,Y = self.phase.UVY
        assert U == approx(self.phase.u)
        assert V == approx(self.phase.v)
        assert Y == approx(self.phase.Y)

        U,V = self.phase.UV
        assert U == approx(self.phase.u)
        assert V == approx(self.phase.v)

        S,P,X = self.phase.SPX
        assert S == approx(self.phase.s)
        assert P == approx(self.phase.P)
        assert X == approx(self.phase.X)

        S,P,Y = self.phase.SPY
        assert S == approx(self.phase.s)
        assert P == approx(self.phase.P)
        assert Y == approx(self.phase.Y)

        S,P = self.phase.SP
        assert S == approx(self.phase.s)
        assert P == approx(self.phase.P)

        S,V,X = self.phase.SVX
        assert S == approx(self.phase.s)
        assert V == approx(self.phase.v)
        assert X == approx(self.phase.X)

        S,V,Y = self.phase.SVY
        assert S == approx(self.phase.s)
        assert V == approx(self.phase.v)
        assert Y == approx(self.phase.Y)

        S,V = self.phase.SV
        assert S == approx(self.phase.s)
        assert V == approx(self.phase.v)

        D,P,X = self.phase.DPX
        assert D == approx(self.phase.density)
        assert P == approx(self.phase.P)
        assert X == approx(self.phase.X)

        D,P,Y = self.phase.DPY
        assert D == approx(self.phase.density)
        assert P == approx(self.phase.P)
        assert Y == approx(self.phase.Y)

        D,P = self.phase.DP
        assert D == approx(self.phase.density)
        assert P == approx(self.phase.P)

        Te = self.phase.Te
        assert Te == approx(self.phase.Te)

    def test_getState_mass(self):
        self.phase.TDY = 350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.check_getters()

    def test_getState_mole(self):
        self.phase.basis = 'molar'
        self.phase.TDX = 350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'
        self.check_getters()

    def test_getState(self):
        assert self.phase.P == approx(ct.one_atm)
        assert self.phase.T == approx(300)

    def test_partial_molar(self):
        self.phase.TDY = 350.0, 0.6, 'H2:0.1, H2O2:0.1, AR:0.8'
        assert sum(self.phase.partial_molar_enthalpies * self.phase.X) == approx(
               self.phase.enthalpy_mole)

        assert sum(self.phase.partial_molar_entropies * self.phase.X) == approx(
               self.phase.entropy_mole)

        assert sum(self.phase.partial_molar_int_energies * self.phase.X) == approx(
               self.phase.int_energy_mole)

        assert sum(self.phase.chemical_potentials * self.phase.X) == approx(
               self.phase.gibbs_mole)

        assert sum(self.phase.partial_molar_cp * self.phase.X) == approx(
               self.phase.cp_mole)

    def test_nondimensional(self):
        self.phase.TDY = 850.0, 0.2, 'H2:0.1, H2O:0.6, AR:0.3'
        H = (sum(self.phase.standard_enthalpies_RT * self.phase.X) *
             ct.gas_constant * self.phase.T)
        assert H == approx(self.phase.enthalpy_mole)

        U = (sum(self.phase.standard_int_energies_RT * self.phase.X) *
             ct.gas_constant * self.phase.T)
        assert U == approx(self.phase.int_energy_mole)

        cp = sum(self.phase.standard_cp_R * self.phase.X) * ct.gas_constant
        assert cp == approx(self.phase.cp_mole)

    def test_activities(self):
        self.phase.TDY = 850.0, 0.2, 'H2:0.1, H2O:0.6, AR:0.3'
        assert self.phase.X == approx(self.phase.activities)

        assert self.phase.activity_coefficients == approx(np.ones(self.phase.n_species))

    def test_isothermal_compressibility(self):
        assert self.phase.isothermal_compressibility == approx(1.0/self.phase.P)

    def test_thermal_expansion_coeff(self):
        assert self.phase.thermal_expansion_coeff == approx(1.0/self.phase.T)

    def test_ref_info(self):
        assert self.phase.reference_pressure == approx(ct.one_atm)
        assert self.phase.min_temp == approx(300.0)
        assert self.phase.max_temp == approx(3500.0)

    def test_uncopyable(self):
        import copy
        with pytest.raises(NotImplementedError):
            copy.copy(self.phase)

    def test_add_species(self):
        ref = ct.Solution('gri30.yaml', transport_model=None)
        n_orig = self.phase.n_species
        self.phase.add_species(ref.species('CO2'))
        self.phase.add_species(ref.species('CO'))

        assert self.phase.n_species == n_orig + 2
        assert 'CO2' in self.phase.species_names
        assert 'CO' in self.phase.species_names

        state = 400, 2e5, 'H2:0.7, CO2:0.2, CO:0.1'
        ref.TPY = state
        self.phase.TPY = state
        assert self.phase.enthalpy_mass == approx(ref.enthalpy_mass)
        assert self.phase.entropy_mole == approx(ref.entropy_mole)
        assert ref[self.phase.species_names].partial_molar_cp == approx(
               self.phase.partial_molar_cp)

    def test_add_species_disabled(self):
        ref = ct.Solution('gri30.yaml', transport_model=None)

        self.phase.transport_model = "unity-Lewis-number"
        reactor = ct.IdealGasReactor(self.phase, clone=False)
        with pytest.raises(ct.CanteraError, match='Cannot add species'):
            reactor.phase.add_species(ref.species('CH4'))
        del reactor
        gc.collect()
        self.phase.add_species(ref.species('CH4'))

        flame = ct.FreeFlame(self.phase, width=0.1)
        with pytest.raises(ct.CanteraError, match='Cannot add species'):
            self.phase.add_species(ref.species('CO'))
        del flame
        gc.collect()
        self.phase.add_species(ref.species('CO'))

        mix = ct.Mixture([(self.phase, 2.0)])
        with pytest.raises(ct.CanteraError, match='Cannot add species'):
            self.phase.add_species(ref.species('CH2O'))
        del mix
        gc.collect()
        self.phase.add_species(ref.species('CH2O'))

    def test_add_species_duplicate(self):
        species = self.phase.species('H2O2')
        with pytest.raises(ct.CanteraError, match='already contains'):
            self.phase.add_species(species)

    def test_auxiliary_data(self):
        """
        This test checks that the values of `aAlpha_mix` and `b_mix` should, when
        used in the Peng-Robinson equation of state expression, yield the same
        pressure as the pressure that was set for the thermodynamic state.
        """
        gas = ct.Solution('co2_PR_example.yaml', 'CO2-PR')
        gas.TPX = 300, 101325, 'H2:1.0'
        params = gas.auxiliary_data

        # Get the Peng-Robinson equation of state parameters
        aAlpha_mix = params['aAlpha_mix']
        b_mix = params['b_mix']

        molar_volume = 1.0 / gas.density_mole
        # Evaluate the pressure using the Peng-Robinson equation of state
        # parameters
        evaluated_pressure = (
            ct.gas_constant * gas.T / (molar_volume - b_mix) -
            aAlpha_mix / (molar_volume ** 2 + 2 * b_mix * molar_volume - b_mix ** 2)
        )

        assert abs(gas.P - evaluated_pressure) < 1e-6

    def _partial_lookup_peng_robinson_phase(self):
        """
        Build a mixed-parameter Peng-Robinson phase:
        - CH4 / O2 use critical-property lookup
        - all other species use ideal-like PR placeholders (a=b=0)
        """
        species = ct.Species.list_from_file("co2_PR_example.yaml")
        phase_species = []
        for sp in species:
            data = dict(sp.input_data)
            if data["name"] in {"CH4", "O2"}:
                data.pop("equation-of-state", None)
            else:
                data["equation-of-state"] = {
                    "model": "Peng-Robinson",
                    "a": 0.0,
                    "b": 0.0,
                    "acentric-factor": 0.0,
                }
            phase_species.append(ct.Species.from_dict(data))

        return ct.Solution(
            thermo="Peng-Robinson",
            kinetics="bulk",
            species=phase_species,
            reactions=[],
            name="partial-lookup-pr",
        )

    @pytest.mark.parametrize(
        "temperature,pressure,composition",
        [
            (
                198.5,
                6.0e7,
                [
                    0.13124237144841314,
                    0.10978856515950584,
                    0.13745471373026613,
                    0.23651613973352531,
                    0.25457002638236204,
                    0.06814754884267295,
                    0.062280634703254616,
                ],
            ),
            (220.0, 6.0e7, [0.2, 0.08, 0.12, 0.2, 0.25, 0.1, 0.05]),
            (260.0, 6.0e7, [0.1, 0.2, 0.1, 0.2, 0.2, 0.1, 0.1]),
            (240.0, 4.0e7, [0.12, 0.14, 0.14, 0.2, 0.22, 0.08, 0.1]),
            (300.0, 2.0e7, [0.2, 0.1, 0.15, 0.2, 0.2, 0.1, 0.05]),
        ],
    )
    def test_partial_lookup_pr_state_setting(self, temperature, pressure, composition):
        gas = self._partial_lookup_peng_robinson_phase()
        gas.TPX = temperature, pressure, np.array(composition)
        assert np.isfinite(gas.density)
        assert gas.density > 0.0
        assert gas.P == approx(pressure)

    @pytest.mark.parametrize("p", [2.0e7, 4.0e7, 6.0e7])
    def test_partial_lookup_pr_high_pressure_equilibrate(self, p):
        gas = self._partial_lookup_peng_robinson_phase()

        fuel_temp = 255.0
        oxidizer_temp = 142.0

        gas.TPY = fuel_temp, p, "CH4:1.0"
        yin_f = gas.Y
        gas.TPY = oxidizer_temp, p, "O2:1.0"
        yin_o = gas.Y
        zst = 1.0 / (1.0 + gas.stoich_air_fuel_ratio(yin_f, yin_o, "mass"))
        yst = zst * yin_f + (1.0 - zst) * yin_o

        tbar = 0.5 * (fuel_temp + oxidizer_temp)
        gas.TPY = tbar, p, yst
        gas.equilibrate("HP")

        assert np.isfinite(gas.T)
        assert gas.T > tbar
        assert gas.P == approx(p)


@pytest.fixture(scope='function')
def setup_thermo_tests(request):
    request.cls.gas = ct.ThermoPhase("h2o2.yaml")
    request.cls.gas.TPX = 450, 2e5, 'H2:1.0, O2:0.4, AR:3, H2O:0.1'

@pytest.mark.usefixtures("setup_thermo_tests")
class TestThermo:

    def test_setSV_lowT(self):
        """
        Set state in terms of (s,v) when the end temperature is below the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        s1, v1 = self.gas.SV
        self.gas.SV = s1, 3 * v1

        assert self.gas.s == approx(s1)
        assert self.gas.v == approx(3 * v1)
        assert self.gas.T < self.gas.min_temp

    def test_setSV_low_invalid(self):
        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        self.gas.SV = 4600, None
        with pytest.raises(ct.CanteraError):
            self.gas.SV = -1000, None

    def test_setSV_highT(self):
        """
        Set state in terms of (s,v) when the end temperature is above the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 2900, 1e5, 'H2:1.0, O2:0.4, AR:3'
        s1, v1 = self.gas.SV
        self.gas.SV = s1, 0.3 * v1

        assert self.gas.s == approx(s1)
        assert self.gas.v == approx(0.3 * v1)
        assert self.gas.T > self.gas.max_temp

    def test_setHP_lowT(self):
        """
        Set state in terms of (s,v) when the end temperature is below the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        deltaH = 1.25e5
        h1, p1 = self.gas.HP
        self.gas.HP = h1 - deltaH, None

        assert self.gas.h == approx(h1 - deltaH)
        assert self.gas.P == approx(p1)
        assert self.gas.T < self.gas.min_temp

    def test_setHP_low_invalid(self):
        """
        Set state in terms of (h,p) when the enthalpy would imply a negative
        temperature
        """

        self.gas.TPX = 300, 101325, 'H2:1.0'
        with pytest.raises(ct.CanteraError):
            self.gas.HP = -4e6, 101325

    def test_setHP_highT(self):
        """
        Set state in terms of (s,v) when the end temperature is above the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 2800, 1e5, 'H2:1.0, O2:0.4, AR:3'
        deltaH = 8.25e5
        h1, p1 = self.gas.HP
        self.gas.HP = h1 + deltaH, None

        assert self.gas.h == approx(h1 + deltaH)
        assert self.gas.P == approx(p1)
        assert self.gas.T > self.gas.max_temp

    def test_volume(self):
        """
        This phase should follow the ideal gas law
        """
        g = self.gas
        assert g.P == approx(g.density_mole * ct.gas_constant * g.T)

        assert g.P / g.density == approx(
               ct.gas_constant / g.mean_molecular_weight * g.T)

        assert g.density == approx(1.0 / g.volume_mass)

    def test_energy(self):
        g = self.gas
        mmw = g.mean_molecular_weight
        assert g.enthalpy_mass == approx(g.enthalpy_mole / mmw)
        assert g.int_energy_mass == approx(g.int_energy_mole / mmw)
        assert g.gibbs_mass == approx(g.gibbs_mole / mmw)
        assert g.entropy_mass == approx(g.entropy_mole / mmw)

        assert g.cv_mass == approx(g.cv_mole / mmw)
        assert g.cp_mass == approx(g.cp_mole / mmw)
        assert g.cv_mole + ct.gas_constant == approx(g.cp_mole)

    def test_nondimensional(self):
        g = self.gas
        R = ct.gas_constant

        assert np.dot(g.standard_cp_R, g.X) == approx(g.cp_mole / R)
        assert np.dot(g.standard_enthalpies_RT, g.X) == approx(
               g.enthalpy_mole / (R*g.T))

        Smix_R = - np.dot(g.X, np.log(g.X+1e-20))
        assert np.dot(g.standard_entropies_R, g.X) + Smix_R == approx(
               g.entropy_mole / R)
        assert np.dot(g.standard_gibbs_RT, g.X) - Smix_R == approx(
               g.gibbs_mole / (R*g.T))


@pytest.fixture(scope='function')
def setup_interface_tests(request):
    request.cls.interface = ct.Interface("diamond.yaml", "diamond_100")

@pytest.mark.usefixtures("setup_interface_tests")
class TestInterfacePhase:

    def test_properties(self):
        self.interface.site_density = 100
        assert self.interface.site_density == approx(100)

    def test_coverages_array(self):
        C = np.zeros(self.interface.n_species)
        C[1] = 0.25
        C[3] = 0.125
        C[4] = 0.125
        self.interface.coverages = C
        C = self.interface.coverages
        # should now be normalized
        assert C[1] == approx(0.5)
        assert C[3] == approx(0.25)
        assert C[4] == approx(0.25)
        assert sum(C) == approx(1.0)

    def test_mole_fractions(self):
        self.interface.X = 'c6HM:0.3, c6H*:0.7'
        assert sum(self.interface.concentrations) == approx(self.interface.site_density)

    def test_coverages_string(self):
        self.interface.coverages = 'c6HM:0.2, c6H*:0.8'
        C = self.interface.coverages
        assert C[self.interface.species_index('c6HM')] == approx(0.2)
        assert C[self.interface.species_index('c6H*')] == approx(0.8)

    def test_coverages_dict(self):
        self.interface.coverages = {'c6**':1.0, 'c6*M':3.0}
        C = self.interface.coverages
        assert C[self.interface.species_index('c6**')] == approx(0.25)
        assert C[self.interface.species_index('c6*M')] == approx(0.75)


class TestInterfacePhase2:

    @pytest.fixture(scope='function')
    def surf(self):
        """Fixture to initialize the Cantera interface phase."""
        return ct.Interface("surface-phases.yaml", "Pt-multi-sites")

    """ Test special cases of interface phases """
    def test_multi_site_species(self, surf):
        # O2(s) consumes two surface sites
        surf.coverages = {"Pt(s)": 0.5, "H(s)": 0.1, "O2(s)": 0.4}
        X = surf.mole_fraction_dict()
        moles = 0.5 + 0.1 + 0.4 / 2
        assert np.isclose(X["Pt(s)"], 0.5 / moles)
        assert np.isclose(X["H(s)"], 0.1 / moles)
        assert np.isclose(X["O2(s)"], 0.2 / moles)

    def test_multi_site_unnormalized(self, surf):
        theta = [0.5, 0.1, 0.1, 0.32]
        surf.set_unnormalized_coverages(theta)
        assert surf.coverages == approx(theta)

        X_unnormalized = surf.X
        assert np.isclose(sum(theta), sum(surf.X))

        surf.coverages = theta
        X_normalized = surf.X
        assert X_normalized * sum(theta) == approx(X_unnormalized)

    def test_multi_site_density(self, surf):
        surf.coverages = {"O2(s)": 1.0}  # O2(s) covers two sites
        assert sum(surf.concentrations) == approx(0.5 * surf.site_density)

        surf.coverages = {"Pt(s)": 1.0}  # Pt(s) covers one site
        assert sum(surf.concentrations) == approx(surf.site_density)

    def test_finite_difference_derivative(self, surf):
        theta0 = np.array([0.5, 0.1, 0.1, 0.3])
        surf.coverages = theta0
        C0 = surf.concentrations
        dtheta = 0.01

        for k in range(surf.n_species):
            theta = theta0.copy()
            theta[k] += dtheta
            surf.set_unnormalized_coverages(theta)
            C = surf.concentrations
            dC_dtheta = (C - C0) / dtheta
            assert dC_dtheta[k] == approx(surf.site_density / surf.species(k).size)
            for j in range(surf.n_species):
                if j != k:
                    assert abs(dC_dtheta[j]) < 1e-8 * surf.site_density


class TestPlasmaPhase:
    @pytest.fixture(scope='function')
    def phase(self):
        phase = ct.Solution('oxygen-plasma.yaml', 'isotropic-electron-energy-plasma',
                           transport_model=None)
        phase.isotropic_shape_factor = 1.0
        return phase

    @property
    def collision_data(self):
        return {
            "equation": "O2 + E => E + O2",
            "type": "electron-collision-plasma",
            "energy-levels": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
            "cross-sections": [0.0, 3.83e-20, 4.47e-20, 4.79e-20, 5.07e-20, 5.31e-20,
                               5.49e-20, 5.64e-20, 5.77e-20, 5.87e-20, 5.97e-20],
            "duplicate": True,
        }

    def test_converting_electron_energy_to_temperature(self, phase):
        phase.mean_electron_energy = 1.0
        Te = 2.0 / 3.0 * ct.electron_charge / ct.boltzmann
        assert phase.Te == approx(Te)

    def test_converting_electron_temperature_to_energy(self, phase):
        phase.Te = 10000
        energy = phase.Te * 3.0 / 2.0 / ct.electron_charge * ct.boltzmann
        assert phase.mean_electron_energy == approx(energy)

    def test_set_get_electron_energy_levels(self, phase):
        levels = np.linspace(0.01, 10, num=9)
        phase.electron_energy_levels = levels
        assert levels == approx(phase.electron_energy_levels)

    def test_isotropic_velocity_electron_energy_distribution(self, phase):
        levels = np.linspace(0.01, 10, num=9)
        phase.electron_energy_levels = levels
        phase.Te = 2e5
        mean_electron_energy = 3.0 / 2.0 * (phase.Te * ct.gas_constant /
                               (ct.avogadro * ct.electron_charge))
        assert mean_electron_energy == approx(phase.mean_electron_energy)

    def test_discretized_electron_energy_distribution(self, phase):
        levels = np.array([0.0, 1.0, 10.0])
        dist = np.array([0.0, 0.9, 0.01])
        phase.normalize_electron_energy_distribution_enabled = False
        phase.quadrature_method = "trapezoidal"
        phase.set_discretized_electron_energy_distribution(levels, dist)
        assert levels == approx(phase.electron_energy_levels)
        assert dist == approx(phase.electron_energy_distribution)
        mean_energy = 2.0 / 5.0 * trapezoid(dist, np.power(levels, 5./2.))
        assert phase.mean_electron_energy == approx(mean_energy, rel=1e-4)
        electron_temp = 2.0 / 3.0 * (phase.mean_electron_energy *
                        ct.avogadro * ct.electron_charge / ct.gas_constant)
        assert phase.Te == approx(electron_temp)

    def test_electron_thermodynamic_properties(self, phase):
        assert phase.standard_gibbs_RT[0] == approx(
               phase.standard_enthalpies_RT[0] - phase.standard_entropies_R[0])

    def test_add_multiple_electron_species(self, phase):
        electron = ct.Species('Electron', 'E:1')
        electron.thermo = ct.ConstantCp(100, 200, 101325, coeffs=(300, 1, 1, 1))
        with pytest.raises(ct.CanteraError,
                           match='Only one electron species is allowed'):
            phase.add_species(electron)

    def test_elastic_power_loss_low_T(self, phase):
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        assert phase.elastic_power_loss == approx(6846332332)

    def test_elastic_power_loss_high_T(self, phase):
        # when T is as high as Te the energy loss rate becomes small
        phase.TPX = 4000, ct.one_atm, "O2:1, E:1e-5"
        assert phase.elastic_power_loss == approx(2865540)

    def test_elastic_power_loss_replace_rate(self, phase):
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        rate = ct.ReactionRate.from_dict(self.collision_data)
        phase.reaction(1).rate = rate
        assert phase.elastic_power_loss == approx(11765800095)

    def test_elastic_power_loss_add_reaction(self, phase):
        phase2 = ct.Solution(thermo="plasma", kinetics="bulk",
                             species=phase.species(), reactions=[])
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        phase.add_reaction(ct.Reaction.from_dict(self.collision_data, phase))
        assert phase.elastic_power_loss == approx(18612132428)

    def test_elastic_power_loss_change_levels(self, phase):
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        phase.electron_energy_levels = np.linspace(0,10,101)
        assert phase.elastic_power_loss == approx(113058853)

    def test_elastic_power_loss_change_dist(self, phase):
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        levels = np.array([0.0, 0.1, 1.0, 9.0, 10.0])
        dist = np.array([0.0, 0.2, 0.7, 0.01, 0.01])
        # set the electron energy levels first to test if
        # set_discretized_electron_energy_distribution triggers
        # updating the interpolated cross sections
        phase.electron_energy_levels = np.array([0.0, 0.1, 1.0, 7.0, 10.0])
        phase.normalize_electron_energy_distribution_enabled = False
        phase.set_discretized_electron_energy_distribution(levels, dist)
        assert phase.elastic_power_loss == approx(7568518396)

    def test_elastic_power_loss_change_mean_electron_energy(self, phase):
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        phase.mean_electron_energy = 2.0
        assert phase.elastic_power_loss == approx(5826212349)

    def test_elastic_power_loss_change_shape_factor(self, phase):
        phase.TPX = 1000, ct.one_atm, "O2:1, E:1e-5"
        phase.isotropic_shape_factor = 1.1
        assert phase.elastic_power_loss == approx(7408711810)

    def test_eedf_solver(self):

        phase = ct.Solution('air-plasma.yaml')
        phase.TPX = 300., 101325., 'N2:0.79, O2:0.21, N2+:1E-10, Electron:1E-10'
        phase.reduced_electric_field = 200.0 * 1e-21  # Reduced electric field [V.m^2]
        phase.update_electron_energy_distribution()

        grid = phase.electron_energy_levels
        eedf = phase.electron_energy_distribution

        reference_grid = np.logspace(-1, np.log10(60))

        reference_eedf = np.array([
            9.1027381e-02, 9.1026393e-02, 9.1025267e-02, 9.1023985e-02, 9.1022523e-02,
            9.1020858e-02, 9.1015025e-02, 9.1006713e-02, 9.0997242e-02, 9.0986450e-02,
            9.0974154e-02, 9.0954654e-02, 9.0923885e-02, 9.0888824e-02, 9.0842837e-02,
            9.0775447e-02, 9.0695937e-02, 9.0578309e-02, 9.0398980e-02, 9.0118320e-02,
            8.9293838e-02, 8.7498617e-02, 8.3767419e-02, 7.5765714e-02, 6.4856820e-02,
            5.5592157e-02, 4.9309310e-02, 4.5268611e-02, 4.2261381e-02, 3.9440745e-02,
            3.6437762e-02, 3.3181527e-02, 2.9616717e-02, 2.5795007e-02, 2.1676205e-02,
            1.7347058e-02, 1.3022044e-02, 8.9705614e-03, 5.5251937e-03, 3.1894295e-03,
            1.7301525e-03, 8.4647152e-04, 3.6030983e-04, 1.2894755e-04, 3.7416645e-05,
            8.4693678e-06, 1.4299900e-06, 1.7026957e-07, 1.3992350e-08, 1.5340110e-09
        ])

        interp = np.interp(reference_grid, grid, eedf, left=0.0, right=0.0)

        mask = reference_eedf > 1e-8
        rel_error = np.abs(interp[mask] - reference_eedf[mask]) / reference_eedf[mask]

        assert max(rel_error) < 0.01

        l2_norm = np.linalg.norm(interp - reference_eedf)
        assert l2_norm < 1e-3


class TestImport:
    """
    Tests the various ways of creating a Solution object
    """
    def check(self, gas, phase, T, P, nSpec, nElem):
        assert gas.name == phase
        assert gas.T == approx(T)
        assert gas.P == approx(P)
        assert gas.n_species == nSpec
        assert gas.n_elements == nElem

    def test_import_from_species(self):
        gas1 = ct.Solution('h2o2.yaml', transport_model=None)
        gas1.TPX = 350, 101325, 'H2:0.3, O2:0.7'
        gas1.equilibrate('HP')

        species = ct.Species.list_from_file("h2o2.yaml")
        gas2 = ct.ThermoPhase(thermo='ideal-gas', species=species)
        gas2.TPX = 350, 101325, 'H2:0.3, O2:0.7'
        gas2.equilibrate('HP')
        assert gas1.n_elements == gas2.n_elements
        assert gas1.species_names == gas2.species_names
        assert gas1.T == approx(gas2.T)
        assert gas1.X == approx(gas2.X)

    def test_yaml_ideal_gas_simple(self):
        gas = ct.ThermoPhase('ideal-gas.yaml', 'simple')
        self.check(gas, 'simple', 500, 10 * ct.one_atm, 3, 2)

    def test_yaml_ideal_gas_remote_species(self):
        gas = ct.ThermoPhase('ideal-gas.yaml', 'species-remote')
        self.check(gas, 'species-remote', 300, ct.one_atm, 4, 2)

    def test_yaml_duplicate(self):
        with pytest.raises(ct.CanteraError, match='duplicate'):
            gas = ct.ThermoPhase('ideal-gas.yaml', 'duplicate-species')


@pytest.fixture(scope='function')
def setup_species_tests(request):
    request.cls.gas = ct.Solution('h2o2.yaml', transport_model=None)

@pytest.mark.usefixtures("setup_species_tests")
class TestSpecies:

    def test_standalone(self):
        s = ct.Species('CH4', {'C':1, 'H':4})

        assert s.name == 'CH4'
        c = s.composition
        assert len(c) == 2
        assert c['C'] == 1
        assert c['H'] == 4
        assert s.molecular_weight == approx(16.043)

    def test_molecular_weight_with_unstable_element(self):
        s = ct.Species("CPo", {"C": 1, "Po": 1})
        with pytest.raises(ct.CanteraError, match="has no stable isotopes"):
            s.molecular_weight

    def test_molecular_weight_with_electron(self):
        yaml = """
            name: Li+[elyt]
            composition: {Li: 1, E: -1}
        """
        s = ct.Species.from_yaml(yaml)
        assert np.isclose(s.molecular_weight, 6.939451420091127)

    def test_custom_element_weight_is_set_on_species(self):
        gas = ct.Solution("ideal-gas.yaml", "element-override")
        # Check that the molecular weight stored in the phase definition is the same
        # as the one on the Species instance
        assert gas["AR"].molecular_weights[0] == approx(gas.species("AR").molecular_weight)
        # Check that the custom value is actually used
        assert gas.species("AR").molecular_weight == approx(36.0)

    def test_species_can_be_added_to_phase(self):
        s = ct.Species.from_dict({
            "name": "AR2",
            "composition": {"Ar": 2},
            "thermo": {"model": "constant-cp", "h0": 100}
        })
        # Access the molecular weight to make sure it's been computed by the Species
        assert s.molecular_weight == approx(39.95 * 2)
        # This should not cause a warning because the Ar element definition in
        # self.gas is the same as the default
        self.gas.add_species(s)
        assert self.gas["AR2"].molecular_weights[0] == approx(
               self.gas.species("AR2").molecular_weight)

    def test_species_warns_when_changing_molecular_weight(self):
        s = ct.Species.from_dict({
            "name": "AR2",
            "composition": {"Ar": 2},
            "thermo": {"model": "constant-cp", "h0": 100}
        })
        # Access the molecular weight to make sure it's been computed by the Species
        assert s.molecular_weight == approx(39.95 * 2)
        gas = ct.Solution("ideal-gas.yaml", "element-override")
        # The warning here is because the weight of the Argon element has been changed in
        # the phase definition, but the molecular weight of the species has already been
        # computed, so loading the phase definition changes the molecular weight on the
        # species.
        with pytest.warns(UserWarning, match="Molecular weight.*changing"):
            gas.add_species(s)

    def test_species_can_be_added_to_phase_custom_element(self):
        s = ct.Species.from_dict({
            "name": "AR2",
            "composition": {"Ar": 2},
            "thermo": {"model": "constant-cp", "h0": 100}
        })
        # DO NOT access the molecular weight on the Species instance before adding it
        # to the phase to make sure the weight has not been computed by the Species
        gas = ct.Solution("ideal-gas.yaml", "element-override")
        gas.add_species(s)
        assert  gas["AR2"].molecular_weights[0] == approx(
                gas.species("AR2").molecular_weight)
        assert s.molecular_weight == approx(gas.species("AR2").molecular_weight)

    def test_defaults(self):
        s = ct.Species('H2')
        assert s.size == 1.0
        assert s.charge == 0.0

        assert s.thermo is None
        assert s.transport is None

    def test_index_accessor(self):
        for k in range(self.gas.n_species):
            s = self.gas.species(k)
            assert s.name == self.gas.species_name(k)

            for m,n in s.composition.items():
                assert n == self.gas.n_atoms(k,m)

    def test_species_noargs(self):
        for k,s in enumerate(self.gas.species()):
            assert s.name == self.gas.species_name(k)

    def test_name_accessor(self):
        for name in self.gas.species_names:
            s = self.gas.species(name)
            assert s.name == name

    def test_listfromFile_yaml(self):
        S = ct.Species.list_from_file("h2o2.yaml")
        assert {sp.name for sp in S} == set(self.gas.species_names)

    def test_list_from_yaml(self):
        yaml = '''
        - name: H2O
          composition: {H: 2, O: 1}
          thermo: {model: constant-cp, h0: 100}
        - name: HO2
          composition: {H: 1, O: 2}
          thermo: {model: constant-cp, h0: 200}
        '''
        species = ct.Species.list_from_yaml(yaml)
        assert species[0].name == 'H2O'
        assert species[1].composition == {'H': 1, 'O': 2}
        assert species[0].thermo.h(300) == approx(100)

    def test_list_from_yaml_section(self, test_data_path):
        species = ct.Species.list_from_yaml(
            (test_data_path / "ideal-gas.yaml").read_text(),
            'species')

        assert species[0].name == 'O2'
        assert species[1].composition == {'N': 1, 'O': 1}

    def test_from_yaml(self):
        yaml = """
        name: H2O
        composition: {H: 2, O: 1}
        thermo: {model: constant-cp, h0: 100}
        """
        species = ct.Species.from_yaml(yaml)
        assert species.name == 'H2O'
        assert species.composition == {'H': 2, 'O': 1}
        assert species.thermo.h(300) == approx(100)

    def test_from_dict(self):
        data = {
            "name": "H2O",
            "composition": {"H": 2, "O": 1},
            "thermo": {"model": "constant-cp", "h0": 100},
        }
        species = ct.Species.from_dict(data)
        assert species.name == 'H2O'
        assert species.composition == {'H': 2, 'O': 1}
        assert species.thermo.h(300) == approx(100)

    def test_modify_thermo(self):
        S = {sp.name: sp for sp in ct.Species.list_from_file("h2o2.yaml")}
        self.gas.TPX = 400, 2*ct.one_atm, 'H2:1.0'
        g0 = self.gas.gibbs_mole

        self.gas.TPX = None, None, 'O2:1.0'
        assert g0 != approx(self.gas.gibbs_mole)
        # Replace O2 thermo with the data from H2
        S['O2'].thermo = S['H2'].thermo
        self.gas.modify_species(self.gas.species_index('O2'), S['O2'])
        assert g0 == approx(self.gas.gibbs_mole)

    def test_modify_thermo_invalid(self):
        S = {sp.name: sp for sp in ct.Species.list_from_file("h2o2.yaml")}

        orig = S['H2']
        thermo = orig.thermo
        copy = ct.Species('foobar', orig.composition)
        copy.thermo = thermo
        with pytest.raises(ct.CanteraError, match='modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

        copy = ct.Species('H2', {'H': 3})
        copy.thermo = thermo
        with pytest.raises(ct.CanteraError, match='modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

        copy = ct.Species('H2', orig.composition)
        copy.thermo = ct.ConstantCp(thermo.min_temp, thermo.max_temp,
            thermo.reference_pressure, [300, 123, 456, 789])
        with pytest.raises(ct.CanteraError, match='modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

        copy = ct.Species('H2', orig.composition)
        copy.thermo = ct.NasaPoly2(thermo.min_temp+200, thermo.max_temp,
            thermo.reference_pressure, thermo.coeffs)
        with pytest.raises(ct.CanteraError, match='modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

    def test_alias(self):
        self.gas.add_species_alias('H2', 'hydrogen')
        assert self.gas.species_index('hydrogen') == 0
        self.gas.X = 'hydrogen:.5, O2:.5'
        assert self.gas.X[0] == approx(0.5)
        with pytest.raises(ct.CanteraError, match='Invalid alias'):
            self.gas.add_species_alias('H2', 'O2')
        with pytest.raises(ct.CanteraError, match='Unable to add alias'):
            self.gas.add_species_alias('spam', 'eggs')

    def test_isomers(self):
        gas = ct.Solution('nDodecane_Reitz.yaml')
        iso = gas.find_isomers({'C':4, 'H':9, 'O':2})
        assert len(iso) == 2
        iso = gas.find_isomers('C:4, H:9, O:2')
        assert len(iso) == 2
        iso = gas.find_isomers({'C':7, 'H':15})
        assert len(iso) == 1
        iso = gas.find_isomers({'C':7, 'H':16})
        assert len(iso) == 0

@pytest.fixture(scope='function')
def setup_species_thermo_tests(request):
    request.cls.gas = ct.Solution('h2o2.yaml', transport_model=None)
    request.cls.gas.X = 'H2O:1.0'

@pytest.mark.usefixtures("setup_species_thermo_tests")
class TestSpeciesThermo:
    h2o_coeffs = [
        1000.0, 3.03399249E+00, 2.17691804E-03, -1.64072518E-07,
        -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00,
        4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09,
        1.77197817E-12, -3.02937267E+04, -8.49032208E-01
    ]

    def test_create(self):
        st = ct.NasaPoly2(300, 3500, 101325, self.h2o_coeffs)

        for T in [300, 500, 900, 1200, 2000]:
            self.gas.TP = T, 101325
            assert st.cp(T) == approx(self.gas.cp_mole)
            assert st.h(T) == approx(self.gas.enthalpy_mole)
            assert st.s(T) == approx(self.gas.entropy_mole)

    def test_invalid(self):
        with pytest.raises(ValueError, match='incorrect length'):
            # not enough coefficients
            st = ct.NasaPoly2(300, 3500, 101325,
                              [1000.0, 3.03399249E+00, 2.17691804E-03])

    def test_wrap(self):
        st = self.gas.species('H2O').thermo

        assert isinstance(st, ct.NasaPoly2)

        for T in [300, 500, 900, 1200, 2000]:
            self.gas.TP = T, 101325
            assert st.cp(T) == approx(self.gas.cp_mole)
            assert st.h(T) == approx(self.gas.enthalpy_mole)
            assert st.s(T) == approx(self.gas.entropy_mole)

    def test_coeffs(self):
        st = ct.NasaPoly2(300, 3500, 101325, self.h2o_coeffs)
        assert st.min_temp == 300
        assert st.max_temp == 3500
        assert st.reference_pressure == 101325
        assert self.h2o_coeffs == approx(st.coeffs)
        assert st.n_coeffs == len(st.coeffs)
        assert st._check_n_coeffs(st.n_coeffs)

    def test_nasa9_load(self):
        gas = ct.Solution("airNASA9.yaml")
        st = gas.species(3).thermo
        assert isinstance(st, ct.Nasa9PolyMultiTempRegion)
        assert st.n_coeffs == len(st.coeffs)
        assert st._check_n_coeffs(st.n_coeffs)

    def test_nasa9_create(self):
        gas = ct.Solution("airNASA9.yaml")
        st = gas.species(3).thermo
        t_min = st.min_temp
        t_max = st.max_temp
        p_ref = st.reference_pressure
        coeffs = st.coeffs
        st2 = ct.Nasa9PolyMultiTempRegion(t_min, t_max, p_ref, coeffs)
        assert isinstance(st2, ct.Nasa9PolyMultiTempRegion)
        assert st.min_temp == t_min
        assert st.max_temp == t_max
        assert st.reference_pressure == p_ref
        for T in range(300, 20000, 1000):
            assert st.cp(T) == approx(st2.cp(T))
            assert st.h(T) == approx(st2.h(T))
            assert st.s(T) == approx(st2.s(T))

    def test_shomate_load(self):
        sol = ct.Solution('thermo-models.yaml', 'molten-salt-Margules')
        st = sol.species(0).thermo
        assert isinstance(st, ct.ShomatePoly2)
        assert st.n_coeffs == len(st.coeffs)
        assert st._check_n_coeffs(st.n_coeffs)

    def test_shomate_create(self):
        sol = ct.Solution('thermo-models.yaml', 'molten-salt-Margules')
        st = sol.species(0).thermo
        t_min = st.min_temp
        t_max = st.max_temp
        p_ref = st.reference_pressure
        coeffs = st.coeffs
        st2 = ct.ShomatePoly2(t_min, t_max, p_ref, coeffs)
        assert isinstance(st2, ct.ShomatePoly2)
        assert st.min_temp == t_min
        assert st.max_temp == t_max
        assert st.reference_pressure == p_ref
        for T in [300, 500, 700, 900]:
            assert st.cp(T) == approx(st2.cp(T))
            assert st.h(T) == approx(st2.h(T))
            assert st.s(T) == approx(st2.s(T))

    def test_piecewise_gibbs_load(self):
        sol = ct.Solution('thermo-models.yaml', 'HMW-NaCl-electrolyte')
        st = sol.species(1).thermo
        assert isinstance(st, ct.Mu0Poly)
        assert st.n_coeffs == len(st.coeffs)
        assert st._check_n_coeffs(st.n_coeffs)

    def test_piecewise_gibbs_create1(self):
        # use OH- ion data from test/thermo/phaseConstructors.cpp
        h298 = -230.015e6
        T1 = 298.15
        mu1 = -91.50963
        T2 = 333.15
        mu2 = -85
        pref = 101325
        coeffs = [2, h298, T1, mu1*ct.gas_constant*T1, T2, mu2*ct.gas_constant*T2]
        st2 = ct.Mu0Poly(200, 3500, pref, coeffs)
        assert isinstance(st2, ct.Mu0Poly)
        assert st2.n_coeffs == len(coeffs)
        assert st2.n_coeffs == len(st2.coeffs)

    def test_piecewise_gibbs_create2(self):
        sol = ct.Solution('thermo-models.yaml', 'HMW-NaCl-electrolyte')
        st = sol.species(1).thermo
        t_min = st.min_temp
        t_max = st.max_temp
        p_ref = st.reference_pressure
        coeffs = st.coeffs
        st2 = ct.Mu0Poly(t_min, t_max, p_ref, coeffs)
        assert isinstance(st2, ct.Mu0Poly)
        assert st.min_temp == t_min
        assert st.max_temp == t_max
        assert st.reference_pressure == p_ref
        for T in [300, 500, 700, 900]:
            assert st.cp(T) == approx(st2.cp(T))
            assert st.h(T) == approx(st2.h(T))
            assert st.s(T) == approx(st2.s(T))

def test_null_species_thermo():
    # Some species don't use a SpeciesThermoInterp object
    liq = ct.Solution('debye-huckel-all.yaml', 'debye-huckel-B-dot-ak-IAPWS')
    assert liq.species('H2O(L)').thermo is None

@pytest.fixture(scope='class')
def setup_quantity_tests(request):
    request.cls.gas = ct.Solution('gri30.yaml', transport_model=None)

@pytest.fixture(scope='function')
def setup_quantity_tests_data(request, setup_quantity_tests):
    request.cls.gas.TPX = 300, 101325, 'O2:1.0, N2:3.76'
    request.cls.gas.basis = 'mass'

@pytest.mark.usefixtures("setup_quantity_tests_data")
class TestQuantity:

    def test_mass_moles(self):
        q1 = ct.Quantity(self.gas, mass=5)
        assert q1.mass == approx(5)
        assert q1.moles == approx(5 / q1.mean_molecular_weight)

        q1.mass = 7
        assert q1.moles == approx(7 / q1.mean_molecular_weight)

        q1.moles = 9
        assert q1.moles == approx(9)
        assert q1.mass == approx(9 * q1.mean_molecular_weight)

    def test_extensive(self):
        q1 = ct.Quantity(self.gas, mass=5)
        assert q1.mass == approx(5)

        assert q1.volume * q1.density == approx(q1.mass)
        assert q1.V * q1.density == approx(q1.mass)
        assert q1.int_energy == approx(q1.moles * q1.int_energy_mole)
        assert q1.enthalpy == approx(q1.moles * q1.enthalpy_mole)
        assert q1.entropy == approx(q1.moles * q1.entropy_mole)
        assert q1.gibbs == approx(q1.moles * q1.gibbs_mole)
        assert q1.int_energy == approx(q1.U)
        assert q1.enthalpy == approx(q1.H)
        assert q1.entropy == approx(q1.S)
        assert q1.gibbs == approx(q1.G)

    def test_set_equivalence_ratio(self):
        q1 = ct.Quantity(self.gas, mass=3)
        T1, P1 = q1.TP
        q1.set_equivalence_ratio(2.0, 'CH4:1.0', 'O2:1.0, N2:3.76')
        assert q1.T == approx(T1)
        assert q1.P == approx(P1)
        assert q1.X[q1.species_index('CH4')] == approx(1.0 / (1 + 4.76))

    def test_set_mixture_fraction(self):
        q1 = ct.Quantity(self.gas, mass=3)
        T1, P1 = q1.TP
        q1.set_mixture_fraction(1.0, 'CH3OH:1.0', 'O2:1.0, N2:3.76')
        assert q1.T == approx(T1)
        assert q1.P == approx(P1)
        assert q1.mass == 3
        assert q1.X[q1.species_index('CH3OH')] == approx(1.0)

    def test_basis(self):
        q1 = ct.Quantity(self.gas, mass=5)
        T1, P1 = q1.TP
        h1 = q1.h  # J/kg

        q1.basis = 'molar'
        assert q1.T == approx(T1)
        assert q1.P == approx(P1)
        assert q1.h == approx(h1 * q1.mean_molecular_weight)

    def test_multiply(self):
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = q1 * 2.5
        assert q1.mass * 2.5 == approx(q2.mass)
        assert q1.moles * 2.5 == approx(q2.moles)
        assert q1.entropy * 2.5 == approx(q2.entropy)
        assert q1.X == approx(q2.X)

    def test_multiply_HP(self):
        self.gas.TPX = 500, 101325, 'CH4:1.0, O2:0.4'
        q1 = ct.Quantity(self.gas, mass=2, constant='HP')
        q2 = ct.Quantity(self.gas, mass=1, constant='HP')
        q2.equilibrate('HP')
        q3 = 0.2 * q1 + q2 * 0.4
        assert q1.P == approx(q3.P)
        assert q1.enthalpy_mass == approx(q3.enthalpy_mass)
        assert q2.enthalpy_mass == approx(q3.enthalpy_mass)

    def test_iadd(self):
        q0 = ct.Quantity(self.gas, mass=5)
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = ct.Quantity(self.gas, mass=5)
        q2.TPX = 500, 101325, 'CH4:1.0'

        q1 += q2
        assert q0.mass + q2.mass == approx(q1.mass)
        # addition is at constant UV
        assert q0.U + q2.U == approx(q1.U)
        assert q0.V + q2.V == approx(q1.V)
        assert q0.X*q0.moles + q2.X*q2.moles == approx(q1.X*q1.moles)

    def test_add(self):
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = ct.Quantity(self.gas, mass=5)
        q2.TPX = 500, 101325, 'CH4:1.0'

        q3 = q1 + q2
        assert q1.mass + q2.mass == approx(q3.mass)
        # addition is at constant UV
        assert q1.U + q2.U == approx(q3.U)
        assert q1.V + q2.V == approx(q3.V)
        assert q1.X*q1.moles + q2.X*q2.moles == approx(q3.X*q3.moles)

    def test_add_molar(self):
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = ct.Quantity(self.gas, mass=5)
        q2.TPX = 900, 101325, 'CH4:1.0'

        q1.basis = 'molar'
        q2.basis = 'molar'

        # addition at constant UV
        q3 = q1 + q2
        assert q1.mass + q2.mass == approx(q3.mass)

        assert q1.U + q2.U == approx(q3.U)
        assert q1.V + q2.V == approx(q3.V)

        # addition at constant HP
        q1.constant = q2.constant = 'HP'
        q4 = q1 + q2
        assert q1.mass + q2.mass == approx(q4.mass)

        assert q1.H + q2.H == approx(q4.H)
        assert q4.P == approx(q1.P)
        assert q1.X*q1.moles + q2.X*q2.moles == approx(q4.X*q4.moles)

    def test_add_errors(self):
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = ct.Quantity(self.gas, mass=5)
        q1.constant = q2.constant = 'HP'
        q2.TP = q1.T, 1.2 * q1.P
        with pytest.raises(ValueError, match="pressure is not equal"):
            q1 + q2

    def test_equilibrate(self):
        self.gas.TPX = 300, 101325, 'CH4:1.0, O2:0.2, N2:1.0'
        q1 = ct.Quantity(self.gas)
        self.gas.equilibrate('HP')
        T2 = self.gas.T

        assert q1.T == approx(300)
        q1.equilibrate('HP')
        assert q1.T == approx(T2)

    def test_invalid_setter(self):
        q1 = ct.Quantity(self.gas, mass =3)
        with pytest.raises(AttributeError):
            q1.HPQ = self.gas.H, self.gas.P, 1

        with pytest.raises(AttributeError):
            q1.set_unnormalized_mass_fractions(np.ones(q1.n_species))

        with pytest.raises(AttributeError):
            q1.set_unnormalized_mole_fractions(np.ones(q1.n_species))

    def test_incompatible(self):
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        q1 = ct.Quantity(self.gas)
        q2 = ct.Quantity(gas2)

        with pytest.raises(ValueError, match='different phase definitions'):
            q1+q2

    def test_disables_add_species(self):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        q1 = ct.Quantity(gas)

        species_x = ct.Species("X", {"H": 3})
        species_x.thermo = ct.ConstantCp(200, 5000, ct.one_atm, coeffs=(0,0,0,0))
        N = gas.n_species
        with pytest.raises(ct.CanteraError, match="is being used"):
            gas.add_species(species_x)
        assert gas.n_species == N

        # Adding species works again after the Solution is no longer in use
        del q1
        gas.add_species(species_x)
        assert gas.n_species == N + 1


class TestMisc:
    def test_stringify_bad(self):
        gas = ct.Solution("h2o2.yaml")
        with pytest.raises(AttributeError):
            gas.add_species_alias("H2", 3)

    def test_case_sensitive_names(self):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        assert not gas.case_sensitive_species_names
        assert gas.species_index('h2') == 0
        gas.X = 'h2:.5, o2:.5'
        assert gas.X[0] == approx(0.5)
        gas.Y = 'h2:.5, o2:.5'
        assert gas.Y[0] == approx(0.5)

        gas.case_sensitive_species_names = True
        assert gas.case_sensitive_species_names
        with pytest.raises(ct.CanteraError, match="Species 'h2' not found"):
            gas.species_index('h2')
        with pytest.raises(ct.CanteraError, match="Species 'h2' not found"):
            gas.X = 'h2:1.0, o2:1.0'
        with pytest.raises(ct.CanteraError, match="Species 'h2' not found"):
            gas.Y = 'h2:1.0, o2:1.0'

        gas_yaml = """
            phases:
            - name: gas
              thermo: ideal-gas
              elements: [S, C, Cs]
              species: [{nasa_gas.yaml/species: all}]
              skip-undeclared-elements: true
        """
        ct.suppress_thermo_warnings(True)
        gas = ct.Solution(yaml=gas_yaml)
        with pytest.raises(ct.CanteraError, match='is not unique'):
            gas.species_index('cs')
        gas.case_sensitive_species_names = True
        with pytest.raises(ct.CanteraError, match="Species 'cs' not found"):
            gas.species_index('cs')

    def test_unstable_element_in_phase(self):
        gas_yaml = """
            phases:
            - name: gas
              thermo: ideal-gas
              elements: [Po]
            species:
            - name: Po
              composition: {Po: 1}
        """
        with pytest.raises(ct.CanteraError, match="has no stable isotopes"):
            ct.Solution(yaml=gas_yaml)


@pytest.fixture(scope='class')
def setup_element_tests(request):
    request.cls.ar_sym = ct.Element('Ar')
    request.cls.ar_name = ct.Element('argon')
    request.cls.ar_num = ct.Element(18)

@pytest.mark.usefixtures("setup_element_tests")
class TestElement:

    def test_element_multiple_possibilities(self):
        # Carbon starts with Ca, the symbol for calcium.
        carbon = ct.Element('Carbon')
        assert carbon.name == 'carbon'
        assert carbon.symbol == 'C'

    def test_element_weight(self):
        assert self.ar_sym.weight == approx(39.95)
        assert self.ar_name.weight == approx(39.95)
        assert self.ar_num.weight == approx(39.95)

    def test_element_symbol(self):
        assert self.ar_sym.symbol == 'Ar'
        assert self.ar_name.symbol == 'Ar'
        assert self.ar_num.symbol == 'Ar'

    def test_element_name(self):
        assert self.ar_sym.name == 'argon'
        assert self.ar_name.name == 'argon'
        assert self.ar_num.name == 'argon'

    def test_element_atomic_number(self):
        assert self.ar_sym.atomic_number == 18
        assert self.ar_name.atomic_number == 18
        assert self.ar_num.atomic_number == 18

    def test_element_name_not_present(self):
        with pytest.raises(ct.CanteraError, match='element not found'):
            ct.Element('I am not an element')

    def test_element_atomic_number_small(self):
        with pytest.raises(ct.CanteraError, match='IndexError'):
            ct.Element(0)

    def test_element_atomic_number_big(self):
        num_elements = ct.Element.num_elements_defined
        with pytest.raises(ct.CanteraError, match='IndexError'):
            ct.Element(num_elements + 1)

    def test_element_no_weight(self):
        with pytest.raises(ct.CanteraError, match='no stable isotopes'):
            ct.Element('Tc')

    def test_element_bad_input(self):
        with pytest.raises(TypeError, match='input argument to Element'):
            ct.Element(1.2345)

    def test_get_isotope(self):
        d_sym = ct.Element('D')
        assert d_sym.atomic_number == 1
        assert d_sym.weight == approx(2.0141017781)
        assert d_sym.name == 'deuterium'
        assert d_sym.symbol == 'D'

        d_name = ct.Element('deuterium')
        assert d_name.atomic_number == 1
        assert d_name.weight == approx(2.0141017781)
        assert d_name.name == 'deuterium'
        assert d_name.symbol == 'D'

    def test_elements_lists(self):
        syms = ct.Element.element_symbols
        names = ct.Element.element_names
        num_elements = ct.Element.num_elements_defined
        assert len(syms) == num_elements
        assert len(names) == num_elements


@pytest.fixture(scope='class')
def setup_solution_array_tests(request):
    request.cls.gas = ct.Solution('h2o2.yaml')

@pytest.mark.usefixtures("setup_solution_array_tests")
class TestSolutionArray:

    def test_passthrough(self):
        states = ct.SolutionArray(self.gas, 3)
        assert states.n_species == self.gas.n_species
        assert states.reaction(10).equation == self.gas.reaction(10).equation

    def test_meta(self):
        meta = {'foo': 'bar', 'spam': 'eggs'}
        states = ct.SolutionArray(self.gas, 3, meta=meta)
        assert states.meta['foo'] == 'bar'
        assert states.meta['spam'] == 'eggs'

    def test_get_state(self):
        states = ct.SolutionArray(self.gas, 4)
        H, P = states.HP
        assert H.shape == (4,)
        assert P.shape == (4,)

        S, P, Y = states.SPY
        assert S.shape == (4,)
        assert P.shape == (4,)
        assert Y.shape == (4, self.gas.n_species)

    def test_idealgas_getters(self):
        N = 11
        states = ct.SolutionArray(self.gas, N)
        getters = 'TDPUVHSXY'  # omit getters that contain Q

        # obtain setters from thermo objects
        all_getters = [k for k in dir(self.gas)
                       if not set(k) - set(getters) and len(k)>1]

        # ensure that getters do not raise attribute errors
        for g in all_getters:
            out = getattr(states, g)
            assert len(out) == len(g)

    def test_properties_onedim(self):
        N = 11
        states = ct.SolutionArray(self.gas, N)
        T = np.linspace(300, 2200, N)
        P = np.logspace(3, 8, N)
        X = 'H2:0.5, O2:0.4, AR:0.1, H2O2:0.01, OH:0.001'
        states.TPX = T, P, X

        assert states.T == approx(T)
        assert states.P == approx(P)

        h = states.enthalpy_mass
        ropr = states.reverse_rates_of_progress
        Dkm = states.mix_diff_coeffs
        for i in range(N):
            self.gas.TPX = T[i], P[i], X
            assert self.gas.enthalpy_mass == approx(h[i])
            assert self.gas.reverse_rates_of_progress == approx(ropr[i])
            assert self.gas.mix_diff_coeffs == approx(Dkm[i])

    def test_properties_ndim(self):
        states = ct.SolutionArray(self.gas, (2,3,5))

        T = np.linspace(300, 2200, 5)
        P = np.logspace(3, 8, 2)[:,np.newaxis, np.newaxis]
        X = np.random.random((3,1,self.gas.n_species))
        states.TPX = T, P, X

        TT, PP = states.TP

        h = states.enthalpy_mass
        ropr = states.reverse_rates_of_progress
        Dkm = states.mix_diff_coeffs

        assert h.shape == (2,3,5)
        assert ropr.shape == (2,3,5,self.gas.n_reactions)
        assert Dkm.shape == (2,3,5,self.gas.n_species)

        for i,j,k in np.ndindex(TT.shape):
            self.gas.TPX = T[k], P[i][0][0], X[j]
            assert self.gas.enthalpy_mass == approx(h[i,j,k])
            assert self.gas.reverse_rates_of_progress == approx(ropr[i,j,k])
            assert self.gas.mix_diff_coeffs == approx(Dkm[i,j,k])

    def test_array_properties_exist(self):
        grid_shape = (7, 3)
        states = ct.SolutionArray(self.gas, grid_shape)

        skip = {
            # Skipped because they are complicated (conversion not implemented)
            "forward_rates_of_progress_ddX", "net_rates_of_progress_ddX",
            "reverse_rates_of_progress_ddX", "state", "net_rates_of_progress_ddCi",
            "forward_rates_of_progress_ddCi", "reverse_rates_of_progress_ddCi"
        }
        skip.update(ct.SolutionArray._passthrough)

        for attr in dir(self.gas):
            if attr.startswith("_") or attr in skip:
                continue

            try:
                soln_value = getattr(self.gas, attr)
            except (ct.CanteraError, ct.ThermoModelMethodError, NotImplementedError):
                continue

            if not isinstance(soln_value, (float, np.ndarray)):
                continue

            assert hasattr(states, attr), attr
            array_value = getattr(states, attr)
            assert array_value.shape == grid_shape + np.asarray(soln_value).shape, attr

    def test_slicing_onedim(self):
        states = ct.SolutionArray(self.gas, 5)
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        T0 = states.T
        H0 = states.enthalpy_mass

        # Verify that original object is updated when slices change
        state = states[1]
        state.TD = 300, 0.5
        assert states.T[0] == approx(500)
        assert states.T[1] == approx(300)
        assert states.P[2] == approx(2e5)
        assert states.density[1] == approx(0.5)

        # Verify that the slices are updated when the original object changes
        states.TD = 900, None
        assert state.T == approx(900)
        assert states.density[1] == approx(0.5)

    def test_slicing_ndim(self):
        states = ct.SolutionArray(self.gas, (2,5))
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        T0 = states.T
        H0 = states.enthalpy_mass

        # Verify that original object is updated when slices change
        row2 = states[1]
        row2.TD = 300, 0.5
        T = states.T
        D = states.density
        assert T[0] == approx(T0[0])
        assert T[1] == approx(300*np.ones(5))
        assert D[1] == approx(0.5*np.ones(5))

        col3 = states[:,2]
        col3.TD = 400, 2.5
        T = states.T
        D = states.density
        assert T[:,2] == approx(400*np.ones(2))
        assert D[:,2] == approx(2.5*np.ones(2))

        # Verify that the slices are updated when the original object changes
        states.TP = 900, None
        assert col3.T == approx(900*np.ones(2))
        assert row2.T == approx(900*np.ones(5))

    def test_extra_create_by_dict(self):
        extra = {"grid": np.arange(10), "velocity": np.random.rand(10)}
        states = ct.SolutionArray(self.gas, 10, extra=extra)
        keys = states.extra
        assert keys[0] == 'grid'
        assert states.grid == approx(np.arange(10))

    def test_extra_no_shape(self):
        # The shape of the value for "prop" here is (), which is falsey
        # and causes the use of np.full()
        states = ct.SolutionArray(self.gas, 3, extra={"prop": 1})
        assert states.prop.shape == (3,)
        assert states.prop == approx(np.array((1, 1, 1)))

        # Check a multidimensional SolutionArray
        states = ct.SolutionArray(self.gas, (2, 2), extra={"prop": 3})
        assert states.prop.shape == (2, 2)
        assert states.prop == approx(np.array(((3, 3), (3, 3))))

    def test_extra_not_empty(self):
        """Test that a non-empty SolutionArray raises a ValueError if
           initial values for properties are not supplied.
        """
        with pytest.raises(ValueError, match="Initial values for extra components"):
            ct.SolutionArray(self.gas, 3, extra=["prop"])
        with pytest.raises(ValueError, match="Initial values for extra components"):
            ct.SolutionArray(self.gas, 3, extra=np.array(["prop", "prop2"]))

    def test_extra_create_multidim(self):
        # requires matching first dimensions
        extra_val = [[1, 2, 3] for i in range(5)]
        states = ct.SolutionArray(self.gas, 5, extra={"prop": extra_val})
        assert states.prop.shape == (5, 3,)
        states = ct.SolutionArray(self.gas, 5, extra={"prop": np.array(extra_val)})
        assert states.prop.shape == (5, 3,)
        states = ct.SolutionArray(self.gas, (3, 4,), extra={"prop": np.ones((3, 4, 5,))})
        assert states.prop.shape == (3, 4, 5,)
        states = ct.SolutionArray(self.gas, 1, extra={"prop": extra_val[0]})
        assert states.prop.shape == (1, 3,)
        states = ct.SolutionArray(self.gas, 1, extra={"prop": [2]})
        assert states.prop.shape == (1,)
        with pytest.raises(ValueError, match="Unable to map"):
            ct.SolutionArray(self.gas, (3, 3), extra={"prop": np.arange(3)})

    def test_extra_create_by_iterable(self):
        states = ct.SolutionArray(self.gas, extra=("prop1"))
        assert states.prop1.shape == (0,)

        # An integer is not an iterable, and only bare strings are
        # turned into iterables
        with pytest.raises(ValueError, match="Extra components"):
            ct.SolutionArray(self.gas, extra=2)

    def test_extra_not_string(self):
        with pytest.raises(TypeError, match="is not a string"):
            ct.SolutionArray(self.gas, extra=[1])

    def test_extra_no_objects(self):
        with pytest.raises(ValueError, match="not supported"):
            prop = np.array([0, [1, 2], (3, 4)], dtype=object)
            states = ct.SolutionArray(self.gas, 3, extra={"prop": prop})

    def test_extra_reserved_names(self):
        with pytest.raises(ValueError, match="name is already used"):
            ct.SolutionArray(self.gas, extra=["creation_rates"])

        with pytest.raises(ValueError, match="name is already used"):
            ct.SolutionArray(self.gas, extra={"creation_rates": 0})

    def test_extra_create_by_string(self):
        states = ct.SolutionArray(self.gas, extra="prop")
        assert states.prop.shape == (0,)

    def test_extra_setattr(self):
        states = ct.SolutionArray(self.gas, 7, extra={'prop': range(7)})
        states.prop = 0
        assert states.prop == approx(np.zeros((7,)))
        mod_array = np.linspace(0, 10, 7).astype(np.int64)
        states.prop = mod_array
        assert states.prop == approx(mod_array)
        with pytest.raises(ValueError, match="Incompatible shapes"):
            states.prop = [1, 2]

    def test_assign_to_slice(self):
        # assign to slices does not work after Cantera 3.0
        # as Python and C++ representations do not reference shared memory
        states = ct.SolutionArray(self.gas, 7, extra={'prop': range(7)})
        array = np.arange(7)
        assert states.prop == approx(array)
        states[1].prop = -5
        assert states.prop[1] == -5
        with pytest.raises(ValueError, match="read-only"):
            states.prop[1] = -10
        # assign to multi-dimensional extra
        extra_val = [[1, 2, 3] for i in range(5)]
        states = ct.SolutionArray(self.gas, 5, extra={"prop": extra_val})
        states[1].prop = [-1, -1, -1]
        assert (states[1].prop == np.array([-1, -1, -1])).all()
        with pytest.raises(ValueError, match="read-only"):
            states.prop[:, 1] = -2

    def test_extra_create_by_ndarray(self):
        properties_array = np.array(["prop1", "prop2", "prop3"])
        states = ct.SolutionArray(self.gas, shape=(0,), extra=properties_array)
        assert states.prop1.shape == (0,)
        assert states.prop2.shape == (0,)
        assert states.prop3.shape == (0,)
        # Ensure that a 2-dimensional array is flattened
        properties_array = np.array((["prop1"], ["prop2"]))
        states = ct.SolutionArray(self.gas, extra=properties_array)
        assert states.prop1.shape == (0,)
        assert states.prop2.shape == (0,)

    def test_append(self):
        states = ct.SolutionArray(self.gas, 5)
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        assert states.cp_mass.shape == (5,)

        states.append(T=1100, P=3e5, X='AR:1.0')
        assert states.cp_mass.shape == (6,)
        assert states.P[-1] == approx(3e5)
        assert states.T[-1] == approx(1100)

        self.gas.TPX = 1200, 5e5, 'O2:0.3, AR:0.7'
        states.append(self.gas.state)
        assert states.cp_mass.shape == (7,)
        assert states.P[-1] == approx(5e5)
        assert states.X[-1, self.gas.species_index('AR')] == approx(0.7)

        self.gas.TPX = 300, 1e4, 'O2:0.5, AR:0.5'
        HPY = self.gas.HPY
        self.gas.TPX = 1200, 5e5, 'O2:0.3, AR:0.7'  # to make sure it gets changed
        states.append(HPY=HPY)
        assert states.cp_mass.shape == (8,)
        assert states.P[-1] == approx(1e4)
        assert states.T[-1] == approx(300)

    def test_append_with_extra(self):
        states = ct.SolutionArray(self.gas, 5, extra={"prop": "value"})
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        assert states.shape == (5,)
        states.append(T=1100, P=3e5, X="AR:1.0", prop="value2")
        assert states.prop[-1] == "value2"
        assert states.prop.shape == (6,)
        states.append(T=1100, P=3e5, X="AR:1.0", prop=100)
        # NumPy converts to the existing type of the array
        assert states.prop[-1] == "100"
        assert states.prop.shape == (7,)
        # two-dimensional input array
        states = ct.SolutionArray(self.gas, 1, extra={"prop": [1, 2, 3]})
        states.append(T=1100, P=3e5, X="AR:1.0", prop=[4, 5, 6])
        assert states.shape == (2,)

    def test_append_failures(self):
        states = ct.SolutionArray(self.gas, 5, extra={"prop": "value"})
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        assert states.shape == (5,)

        with pytest.raises(TypeError, match="Missing keyword arguments for extra"):
            states.append(T=1100, P=3e5, X="AR:1.0")
        # Failing to append a state shouldn't change the size
        assert states.shape == (5,)

        with pytest.raises(KeyError, match="does not specify"):
            # I is not a valid property
            states.append(TPI=(1100, 3e5, "AR:1.0"), prop="value2")
        # Failing to append a state shouldn't change the size
        assert states.shape == (5,)

        with pytest.raises(KeyError, match="is not a valid"):
            # I is not a valid property
            states.append(T=1100, P=3e5, I="AR:1.0", prop="value2")
        # Failing to append a state shouldn't change the size
        assert states.shape == (5,)

        with pytest.raises(ct.CanteraError, match="incompatible value"):
            # prop has incompatible dimensions
            states.append(T=1100, P=3e5, X="AR:1.0", prop=[1, 2, 3])
        # Failing to append a state shouldn't change the size
        assert states.shape == (5,)

        states = ct.SolutionArray(self.gas, 1, extra={"prop": [1, 2, 3]})
        with pytest.raises(ct.CanteraError, match="incompatible value"):
            # prop has incorrect type (no implicit conversion; changed from Cantera 2.6)
            states.append(T=1100, P=3e5, X="AR:1.0", prop=['a', 'b', 'c'])
        # Failing to append a state shouldn't change the size
        assert states.shape == (1,)

        states = ct.SolutionArray(self.gas, 1, extra={"prop": [1, 2, 3]})
        with pytest.raises(ct.CanteraError, match="does not match"):
            # prop has incorrect dimensions
            states.append(T=1100, P=3e5, X="AR:1.0", prop=[1, 2])
        # Failing to append a state shouldn't change the size
        assert states.shape == (1,)

    def test_purefluid(self):
        water = ct.Water()
        states = ct.SolutionArray(water, 5)
        states.TQ = 400, np.linspace(0, 1, 5)

        P = states.P
        for i in range(1, 5):
            assert P[0] == approx(P[i])

        states.TP = np.linspace(400, 500, 5), 101325
        assert states.Q.squeeze() == approx(np.ones(5))

    def test_phase_of_matter(self):
        water = ct.Water()
        states = ct.SolutionArray(water, 5)
        T = [300, 500, water.critical_temperature*2, 300]
        P = [101325, 101325, 101325, water.critical_pressure*2]
        states[:4].TP = T, P
        states[4].TQ = 300, .4
        pom = ['liquid', 'gas', 'supercritical', 'supercritical', 'liquid-gas-mix']
        assert list(states.phase_of_matter) == pom

    def test_purefluid_getters(self):
        N = 11
        water = ct.Water()
        states = ct.SolutionArray(water, N)
        getters = 'TDPUVHSQ'  # omit getters that contain X or Y

        # obtain setters from thermo objects
        all_getters = [k for k in dir(water)
                       if not set(k) - set(getters) and len(k)>1]

        # ensure that getters do not raise attribute errors
        for g in all_getters:
            out = getattr(states, g)
            assert len(out) == len(g)

    def test_purefluid_yaml_state(self):
        yaml = """
        phases:
        - name: water
          thermo: pure-fluid
          species: [{{liquidvapor.yaml/species: [H2O]}}]
          state: {{T: {T}, P: 101325, Q: {Q}}}
          pure-fluid-name: water
        """
        w = ct.PureFluid(yaml=yaml.format(T=373.177233, Q=0.5))
        assert w.Q == approx(0.5)

        with pytest.raises(ct.CanteraError, match="setState"):
            ct.PureFluid(yaml=yaml.format(T=373, Q=0.5))

        w = ct.PureFluid(yaml=yaml.format(T=370, Q=0.0))
        assert w.P == approx(101325)

        with pytest.raises(ct.CanteraError, match="setState"):
            ct.PureFluid(yaml=yaml.format(T=370, Q=1.0))

    def test_sort(self):
        np.random.seed(0)
        t = np.random.random(101)
        T = np.linspace(300., 1000., 101)
        P = ct.one_atm * (1. + 10.*np.random.random(101))

        states = ct.SolutionArray(self.gas, 101, extra={'t': t})
        states.TP = T, P

        states.sort('t')
        assert (states.t[1:] - states.t[:-1] > 0).all()
        assert not (states.T[1:] - states.T[:-1] > 0).all()
        assert not np.allclose(states.P, P)

        states.sort('T')
        assert not (states.t[1:] - states.t[:-1] > 0).all()
        assert (states.T[1:] - states.T[:-1] > 0).all()
        assert states.P == approx(P)

        states.sort('T', reverse=True)
        assert (states.T[1:] - states.T[:-1] < 0).all()

    def test_set_equivalence_ratio(self):
        states = ct.SolutionArray(self.gas, 8)
        phi = np.linspace(0.5, 2, 8)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        # The mole fraction arrays need to be squeezed here to reduce their
        # dimensionality from a 2-d column array to a vector for comparison
        # with phi.
        states.set_equivalence_ratio(phi, fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(1)
        assert comp == approx(phi)
        states.set_equivalence_ratio(phi[0], fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(1)
        assert comp == approx(np.full_like(phi, phi[0]))
        states.set_equivalence_ratio(phi.tolist(), fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(1)
        assert comp == approx(phi)

    def test_set_equivalence_ratio_wrong_shape_raises(self):
        states = ct.SolutionArray(self.gas, 8)
        phi = np.linspace(0.5, 2, 7)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        with pytest.raises(ValueError, match=r"shape mismatch"):
            states.set_equivalence_ratio(phi, fuel, oxidizer)

    def test_set_equivalence_ratio_2d(self):
        states = ct.SolutionArray(self.gas, (2, 4))
        phi = np.linspace(0.5, 2, 8).reshape((2, 4))
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        # The mole fraction arrays need to be squeezed here to reduce their
        # dimensionality from a 2-d column array to a vector for comparison
        # with phi.
        states.set_equivalence_ratio(phi, fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(2)
        assert comp == approx(phi)

    def test_set_mixture_fraction(self):
        states = ct.SolutionArray(self.gas, 8)
        mixture_fraction = np.linspace(0.5, 1, 8)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        # The mass fraction arrays need to be squeezed here to reduce their
        # dimensionality from a 2-d column array to a vector for comparison
        # with mixture_fraction.
        states.set_mixture_fraction(mixture_fraction, fuel, oxidizer)
        assert states("H2").Y.squeeze(1) == approx(mixture_fraction)
        states.set_mixture_fraction(mixture_fraction[0], fuel, oxidizer)
        assert states("H2").Y.squeeze(1) == approx(np.full_like(mixture_fraction,
                                                                mixture_fraction[0]),)
        states.set_mixture_fraction(mixture_fraction.tolist(), fuel, oxidizer)
        assert states("H2").Y.squeeze(1) == approx(mixture_fraction)

    def test_set_mixture_fraction_wrong_shape_raises(self):
        states = ct.SolutionArray(self.gas, 8)
        mixture_fraction = np.linspace(0.5, 1, 7)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        with pytest.raises(ValueError, match=r"shape mismatch"):
            states.set_mixture_fraction(mixture_fraction, fuel, oxidizer)

    def test_set_mixture_fraction_2D(self):
        states = ct.SolutionArray(self.gas, (2, 4))
        mixture_fraction = np.linspace(0.5, 1, 8).reshape((2, 4))
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        states.set_mixture_fraction(mixture_fraction, fuel, oxidizer)
        # The mass fraction array needs to be squeezed here to reduce its
        # dimensionality from a 3-d array to a 2-d array for comparison
        # with mixture_fraction.
        assert states("H2").Y.squeeze(2) == approx(mixture_fraction)

    def test_species_slicing(self):
        states = ct.SolutionArray(self.gas, (2,5))
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        states.equilibrate('HP')
        assert states('H2').X.squeeze() == approx(
               states.X[...,self.gas.species_index('H2')])

        kk = (self.gas.species_index('OH'), self.gas.species_index('O'))
        assert states('OH','O').partial_molar_cp == approx(
               states.partial_molar_cp[...,kk])

    def test_slice_SolutionArray(self):
        soln = ct.SolutionArray(self.gas, 10)
        arr = soln[2:9:3]
        assert len(arr.T) == 3

    def test_zero_length_slice_SolutionArray(self):
        states = ct.SolutionArray(self.gas, 4)
        arr1 = states[3:3]
        assert len(arr1.T) == 0
        assert arr1.X.shape == (0,10)
        assert arr1.n_reactions == 29

        states.TP = [100,300,900,323.23], ct.one_atm
        arr2 = states[slice(0)]
        assert len(arr2.T) == 0
