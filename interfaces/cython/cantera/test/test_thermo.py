from collections import OrderedDict
import numpy as np
import gc

import cantera as ct
from . import utilities
from .utilities import allow_deprecated
import pytest


class TestThermoPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.yaml', transport_model=None)

    def test_source(self):
        self.assertEqual(self.phase.source, 'h2o2.yaml')

    def test_input_header(self):
        extra = self.phase.input_header
        self.assertTrue(extra["description"].startswith("Hydrogen-Oxygen submechanism"))
        self.assertEqual(extra["cantera-version"], "2.5.0")
        self.assertEqual(extra["generator"], "ck2yaml")

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
        with self.assertRaisesRegex(ct.CanteraError, "Key 'phases' not found"):
            _ = ct.Solution(yaml=yaml)

    def test_deprecated_phase(self):
        yaml = """
        phases:
        - name: phasename
          thermo: ideal-gas
          species: [{h2o2.yaml/species: all}]
          deprecated: This phase is deprecated because I said so.
        """
        with self.assertRaisesRegex(ct.CanteraError, "(?s)phasename.*said so"):
            ct.Solution(yaml=yaml)

    def test_base_attributes(self):
        self.assertIsInstance(self.phase.name, str)
        self.assertIsInstance(self.phase.phase_of_matter, str)
        self.assertIsInstance(self.phase.thermo_model, str)
        self.assertIsInstance(self.phase.kinetics_model, str)
        self.assertIsInstance(self.phase.transport_model, str)
        self.assertIsInstance(self.phase.composite, tuple)
        self.assertEqual(len(self.phase.composite), 3)
        self.assertEqual(self.phase.composite,
                         (self.phase.thermo_model,
                          self.phase.kinetics_model,
                          self.phase.transport_model))
        self.phase.name = 'spam'
        self.assertEqual(self.phase.name, 'spam')
        with self.assertRaises(AttributeError):
            self.phase.type = 'eggs'

    def test_phases(self):
        self.assertEqual(self.phase.n_phases, 1)
        self.assertEqual(self.phase.phase_of_matter, "gas")

    def test_states(self):
        self.assertEqual(self.phase._native_state, ('T', 'D', 'Y'))
        self.assertIn('TPY', self.phase._full_states.values())
        self.assertIn('TD', self.phase._partial_states.values())

    def test_species(self):
        self.assertEqual(self.phase.n_species, 10)
        for i,name in enumerate(self.phase.species_names):
            self.assertEqual(name, self.phase.species_name(i))
            self.assertEqual(i, self.phase.species_index(name))
            self.assertEqual(i, self.phase.species_index(i))
        with self.assertRaisesRegex(ct.CanteraError, 'IndexError thrown by Phase::checkSpeciesIndex'):
            self.phase.species(self.phase.n_species)

    def test_elements(self):
        self.assertEqual(self.phase.n_elements, 4)
        for i,symbol in enumerate(self.phase.element_names):
            self.assertEqual(symbol, self.phase.element_name(i))
            self.assertEqual(i, self.phase.element_index(symbol))
            self.assertEqual(i, self.phase.element_index(i))

    def test_n_atoms(self):
        data = [(1, 'O', 'O'), (2, 'O', 'O2'), (1, b'H', b'OH'),
                (2, 'H', 'H2O'), (2, 'O', 'H2O2'), (1, 'Ar', 'AR'),
                (0, 'O', 'H'), (0, 'H', 'AR'), (0, 'Ar', 'HO2')]
        for (n, elem, species) in data:
            self.assertEqual(self.phase.n_atoms(species, elem), n)
            mElem = self.phase.element_index(elem)
            kSpec = self.phase.species_index(species)
            self.assertEqual(self.phase.n_atoms(kSpec, mElem), n)

        with self.assertRaisesRegex(ValueError, 'No such species'):
            self.phase.n_atoms('C', 'H2')
        with self.assertRaisesRegex(ValueError, 'No such element'):
            self.phase.n_atoms('H', 'CH4')

    def test_elemental_mass_fraction(self):
        self.phase.Y = 'H2O:0.5, O2:0.5'
        Zo = self.phase.elemental_mass_fraction('O')
        Zh = self.phase.elemental_mass_fraction('H')
        Zar = self.phase.elemental_mass_fraction('Ar')

        mO = self.phase.element_index('O')
        self.assertEqual(Zo, self.phase.elemental_mass_fraction(mO))
        self.assertNear(Zo, 0.5 + 0.5 * (15.999 / 18.015))
        self.assertNear(Zh, 0.5 * (2.016 / 18.015))
        self.assertEqual(Zar, 0.0)

        with self.assertRaisesRegex(ValueError, 'No such element'):
            self.phase.elemental_mass_fraction('C')
        with self.assertRaisesRegex(ValueError, 'No such element'):
            self.phase.elemental_mass_fraction(5)

    def test_elemental_mole_fraction(self):
        self.phase.X = 'H2O:0.5, O2:0.5'
        Zo = self.phase.elemental_mole_fraction('O')
        Zh = self.phase.elemental_mole_fraction('H')
        Zar = self.phase.elemental_mole_fraction('Ar')

        mO = self.phase.element_index('O')
        self.assertEqual(Zo, self.phase.elemental_mole_fraction(mO))
        self.assertNear(Zo, (0.5 + 1) / (0.5*3 + 0.5*2))
        self.assertNear(Zh, (2*0.5) / (0.5*3 + 0.5*2))
        self.assertEqual(Zar, 0.0)

        with self.assertRaisesRegex(ValueError, 'No such element'):
            self.phase.elemental_mole_fraction('C')
        with self.assertRaisesRegex(ValueError, 'No such element'):
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
                self.assertNear(self.phase.elemental_mass_fraction(i),
                                self.phase.elemental_mole_fraction(i)
                                * self.phase.atomic_weight(i) / denom)

    def test_weights(self):
        atomic_weights = self.phase.atomic_weights
        molecular_weights = self.phase.molecular_weights
        self.assertEqual(self.phase.n_elements, len(atomic_weights))
        self.assertEqual(self.phase.n_species, len(molecular_weights))

        for i,mw in enumerate(molecular_weights):
            test_weight = 0.0
            for j,aw in enumerate(atomic_weights):
                test_weight += aw * self.phase.n_atoms(i,j)
            self.assertNear(test_weight, mw)

    def test_charges(self):
        gas = ct.Solution('ch4_ion.yaml')
        charges = gas.charges
        test = {'E': -1., 'N2': 0., 'H3O+': 1.}
        for species, charge in test.items():
            self.assertIn(species, gas.species_names)
            index = gas.species_index(species)
            self.assertEqual(charges[index], charge)

    def test_setComposition(self):
        X = np.zeros(self.phase.n_species)
        X[2] = 1.0
        self.phase.X = X
        Y = self.phase.Y

        self.assertEqual(list(X), list(Y))

    def test_setComposition_singleton(self):
        X = np.zeros((1,self.phase.n_species,1))
        X[0,2,0] = 1.0
        self.phase.X = X
        Y = self.phase.Y

        self.assertEqual(list(X[0,:,0]), list(Y))

    def test_setCompositionString(self):
        self.phase.X = 'h2:1.0, o2:1.0'
        X = self.phase.X
        self.assertNear(X[0], 0.5)
        self.assertNear(X[3], 0.5)

        with self.assertRaisesRegex(ct.CanteraError, 'Unknown species'):
            self.phase.X = 'H2:1.0, CO2:1.5'

    def test_setCompositionStringBad(self):
        X0 = self.phase.X
        with self.assertRaisesRegex(ct.CanteraError, 'Trouble processing'):
            self.phase.X = 'H2:1.0, O2:asdf'
        self.assertArrayNear(X0, self.phase.X)

        with self.assertRaisesRegex(ct.CanteraError, 'Trouble processing'):
            self.phase.X = 'H2:1e-x4'
        self.assertArrayNear(X0, self.phase.X)

        with self.assertRaisesRegex(ct.CanteraError, 'decimal point in exponent'):
            self.phase.X = 'H2:1e-1.4'
        self.assertArrayNear(X0, self.phase.X)

        with self.assertRaisesRegex(ct.CanteraError, 'Duplicate key'):
            self.phase.X = 'H2:0.5, O2:1.0, H2:0.1'
        self.assertArrayNear(X0, self.phase.X)

    def test_setCompositionDict(self):
        self.phase.X = {b'H2':1.0, b'O2':3.0}
        X = self.phase.X
        self.assertNear(X[0], 0.25)
        self.assertNear(X[3], 0.75)

        self.phase.Y = {'H2':1.0, 'O2':3.0}
        Y = self.phase.Y
        self.assertNear(Y[0], 0.25)
        self.assertNear(Y[3], 0.75)

    def test_getCompositionDict(self):
        self.phase.X = 'oh:1e-9, O2:0.4, AR:0.6'
        self.assertEqual(len(self.phase.mole_fraction_dict(1e-7)), 2)
        self.assertEqual(len(self.phase.mole_fraction_dict()), 3)

        self.phase.Y = 'O2:0.4, AR:0.6'
        Y1 = self.phase.mass_fraction_dict()
        self.assertNear(Y1['O2'], 0.4)
        self.assertNear(Y1['AR'], 0.6)

    def test_setCompositionNoNorm(self):
        X = np.zeros(self.phase.n_species)
        X[2] = 1.0
        X[0] = 0.01
        self.phase.set_unnormalized_mole_fractions(X)
        self.assertArrayNear(self.phase.X, X)
        self.assertNear(sum(X), 1.01)

        Y = np.zeros(self.phase.n_species)
        Y[2] = 1.0
        Y[0] = 0.01
        self.phase.set_unnormalized_mass_fractions(Y)
        self.assertArrayNear(self.phase.Y, Y)
        self.assertNear(sum(Y), 1.01)

    def test_setCompositionNoNormBad(self):
        X = np.zeros(self.phase.n_species - 1)
        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            self.phase.set_unnormalized_mole_fractions(X)

        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            self.phase.set_unnormalized_mass_fractions([1,2,3])

    def test_setCompositionDict_bad1(self):
        with self.assertRaisesRegex(ct.CanteraError, 'Unknown species'):
            self.phase.X = {'H2':1.0, 'HCl':3.0}

    def test_setCompositionDict_bad2(self):
        with self.assertRaises(TypeError):
            self.phase.Y = {'H2':1.0, 'O2':'xx'}

    def test_setCompositionSlice(self):
        self.phase['h2', 'o2'].X = 0.1, 0.9
        X = self.phase.X
        self.assertNear(X[0], 0.1)
        self.assertNear(X[3], 0.9)

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
        self.assertEqual(gas.Y[0], 1.0)

    def test_setCompositionSlice_bad(self):
        X0 = self.phase.X
        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            self.phase['H2','O2'].Y = [0.1, 0.2, 0.3]
        self.assertArrayNear(self.phase.X, X0)

    def test_setCompositionEmpty_bad(self):
        X0 = self.phase.X
        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            self.phase.Y = np.array([])
        self.assertArrayNear(self.phase.X, X0)

    @utilities.slow_test
    def test_set_equivalence_ratio_stoichiometric(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        for fuel in ('C2H6', 'H2:0.7, CO:0.3', 'NH3:0.4, CH3OH:0.6'):
            for oxidizer in ('O2:1.0, N2:3.76', 'H2O2:1.0'):
                gas.set_equivalence_ratio(1.0, fuel, oxidizer)
                gas.equilibrate('TP')
                # Almost everything should end up as CO2, H2O and N2
                self.assertGreater(sum(gas['H2O','CO2','N2'].X), 0.999999)

    def test_set_equivalence_ratio_lean(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        excess = 0
        for phi in np.linspace(0.9, 0, 5):
            gas.set_equivalence_ratio(phi, 'CH4:0.8, CH3OH:0.2', 'O2:1.0, N2:3.76')
            gas.equilibrate('TP')
            self.assertGreater(gas['O2'].X[0], excess)
            excess = gas['O2'].X[0]
        self.assertNear(sum(gas['O2','N2'].X), 1.0)

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
            self.assertNear(gas.stoich_air_fuel_ratio(fuel, ox, basis)/((1.0-Z)/Z),  2.0)
            gas.set_mixture_fraction(Z, fuel, ox, basis)
            self.assertNear(gas['SO2'].X[0], 31.0/212.0)
            self.assertNear(gas['O2'].X[0],  31.0/106.0)
            self.assertNear(gas['SO'].X[0],  11.0/106.0)
            self.assertNear(gas['CO2'].X[0], 31.0/424.0)
            self.assertNear(gas['CH3'].X[0], 11.0/53.0)
            self.assertNear(gas['N2'].X[0],  11.0/212.0)
            self.assertNear(gas['CH'].X[0],  31.0/424.0)
            self.assertNear(gas['OH'].X[0],  11.0/212.0)
            self.assertNear(gas.equivalence_ratio(fuel, ox, basis),  2.0)

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
            self.assertNear(phi, gas.equivalence_ratio('CH4:0.8, CH3OH:0.2', 'O2:1.0, N2:3.76'))
        # Check sulfur species
        sulfur_species = [k for k in ct.Species.list_from_file("nasa_gas.yaml")
                          if k.name in ("SO", "SO2")]
        gas = ct.Solution(thermo="ideal-gas", kinetics="gas",
                          species=ct.Species.list_from_file("gri30.yaml") + sulfur_species)
        for phi in np.linspace(0.5, 2.0, 5):
            gas.set_equivalence_ratio(phi, 'CH3:0.5, SO:0.25, OH:0.125, N2:0.125', 'O2:0.5, SO2:0.25, CO2:0.125')
            self.assertNear(phi, gas.equivalence_ratio('CH3:0.5, SO:0.25, OH:0.125, N2:0.125', 'O2:0.5, SO2:0.25, CO2:0.125'))
        gas.X = 'CH4:1' # pure fuel
        self.assertEqual(gas.equivalence_ratio(), np.inf)

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

            self.assertNear(phi, 1.3)
            self.assertNear(phi2, 1.3)
            self.assertNear(phi_loc, 1.1726068608195617)
            self.assertNear(mf, 0.13415725911057605)
            self.assertNear(mf_C, (gas.elemental_mass_fraction("C")-Y_Co)/(Y_Cf-Y_Co))
            self.assertNear(mf_O, (gas.elemental_mass_fraction("O")-Y_Oo)/(Y_Of-Y_Oo))
            self.assertNear(l, 8.3901204498353561)
            self.assertNear(gas.P, 1e5)
            self.assertNear(T, 300.0)

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
        self.assertNear(gas.equivalence_ratio(include_species=["H2", "O2"]), 2)
        self.assertNear(gas.equivalence_ratio("H2", "O2",
                                              include_species=["H2", "O2"]), 2)
        self.assertNear(gas.T, T)
        self.assertNear(gas.P, P)
        self.assertArrayNear(gas.X, original_X)

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
                self.assertNear(gas["H2"].X[0], (1 - fraction_value)
                                * inv_afr / (inv_afr + 1))
                self.assertNear(gas["O2"].X[0], (1 - fraction_value) / (inv_afr + 1))
                self.assertNear(gas["CO2"].X[0], fraction_value)
            elif basis == "mass" and fraction_type == "diluent":
                self.assertNear(gas["H2"].Y[0] / gas["O2"].Y[0], inv_afr * M_H2 / M_O2)
                self.assertNear(gas["CO2"].Y[0], fraction_value)
            elif basis == "mole" and fraction_type == "fuel":
                self.assertNear(gas["H2"].X[0], fraction_value)
                self.assertNear(gas["O2"].X[0], fraction_value / inv_afr)
                self.assertNear(gas["CO2"].X[0], 1 - fraction_value * (1 + 1 / inv_afr))
            elif basis == "mass" and fraction_type == "fuel":
                self.assertNear(gas["H2"].Y[0], fraction_value)
                self.assertNear(gas["H2"].Y[0] / gas["O2"].Y[0], inv_afr * M_H2 / M_O2)
            elif basis == "mole" and fraction_type == "oxidizer":
                self.assertNear(gas["H2"].X[0], fraction_value * inv_afr)
                self.assertNear(gas["O2"].X[0], fraction_value)
                self.assertNear(gas["CO2"].X[0], 1 - fraction_value * (1 + inv_afr))
            elif basis == "mass" and fraction_type == "oxidizer":
                self.assertNear(gas["O2"].Y[0], fraction_value)
                self.assertNear(gas["H2"].Y[0] / gas["O2"].Y[0], inv_afr * M_H2 / M_O2)

            Y = gas.Y
            self.assertNear(gas.equivalence_ratio("H2", "O2",
                                                  include_species=["H2", "O2"],
                                                  basis=basis), phi)
            self.assertNear(gas.equivalence_ratio(include_species=["H2", "O2"],
                                                  basis=basis), phi)
            self.assertArrayNear(Y, gas.Y)
            self.assertNear(gas.T, T)
            self.assertNear(gas.P, P)

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
        self.assertArrayNear(gas.X, X_expected)

        gas.set_equivalence_ratio(phi, fuel, oxidizer, basis="mass")
        Y_Mix = gas.Y
        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"diluent": fraction},
                                  diluent=diluent, basis="mass")
        Y_expected = Y_Mix * (1 - fraction) + fraction * Y_diluent
        self.assertArrayNear(gas.Y, Y_expected)

        phi = 0.8
        fraction = 0.05
        gas.set_equivalence_ratio(phi, fuel, oxidizer, basis="mass")
        AFR = gas.stoich_air_fuel_ratio(fuel, oxidizer, basis="mass") / phi
        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"fuel": fraction},
                                  diluent=diluent, basis="mass")
        Y_expected = fraction * (Y_fuel + AFR * Y_oxidizer) \
                     + (1 - fraction * (1 + AFR)) * Y_diluent
        self.assertArrayNear(gas.Y, Y_expected)

        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"oxidizer": fraction},
                                  diluent=diluent, basis="mass")
        Y_expected = fraction * (Y_fuel / AFR + Y_oxidizer) \
                     + (1 - fraction * (1 + 1 / AFR)) * Y_diluent
        self.assertArrayNear(gas.Y, Y_expected)

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
        self.assertArrayNear(gas.X, X_expected)

        gas.set_equivalence_ratio(phi, fuel, oxidizer, fraction={"oxidizer": fraction},
                                  diluent=diluent)
        X_expected = fraction * (X_fuel / AFR + X_oxidizer) \
                     + (1 - fraction * (1 + 1 / AFR)) * X_diluent
        self.assertArrayNear(gas.X, X_expected)

    def test_full_report(self):
        report = self.phase.report(threshold=0.0)
        self.assertIn(self.phase.name, report)
        self.assertIn('temperature', report)
        self.assertNotIn('minor', report)
        for name in self.phase.species_names:
            self.assertIn(name, report)

    def test_default_report(self):
        self.phase.X = 'H2:0.1, O2:0.9, HO2:1e-10, H2O2:1e-20'
        report = self.phase.report()
        self.assertIn('minor', report)
        for name in (' H2 ', ' O2 ', ' HO2 '):
            self.assertIn(name, report)
        for name in (' H2O2 ', ' OH ', ' AR '):
            self.assertNotIn(name, report)

    def test_name(self):
        self.phase.name = 'something'
        self.assertEqual(self.phase.name, 'something')
        self.assertIn('something', self.phase.report())

    def test_badLength(self):
        X = np.zeros(5)
        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            self.phase.X = X
        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            self.phase.Y = X

    def test_mass_basis(self):
        self.assertEqual(self.phase.basis, 'mass')
        self.assertNear(self.phase.density_mass, self.phase.density)
        self.assertNear(self.phase.enthalpy_mass, self.phase.h)
        self.assertNear(self.phase.entropy_mass, self.phase.s)
        self.assertNear(self.phase.int_energy_mass, self.phase.u)
        self.assertNear(self.phase.volume_mass, self.phase.v)
        self.assertNear(self.phase.cv_mass, self.phase.cv)
        self.assertNear(self.phase.cp_mass, self.phase.cp)

    def test_molar_basis(self):
        self.phase.basis = 'molar'
        self.assertEqual(self.phase.basis, 'molar')
        self.assertNear(self.phase.density_mole, self.phase.density)
        self.assertNear(self.phase.enthalpy_mole, self.phase.h)
        self.assertNear(self.phase.entropy_mole, self.phase.s)
        self.assertNear(self.phase.int_energy_mole, self.phase.u)
        self.assertNear(self.phase.volume_mole, self.phase.v)
        self.assertNear(self.phase.cv_mole, self.phase.cv)
        self.assertNear(self.phase.cp_mole, self.phase.cp)

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
            self.assertNear(self.phase.T, T)
            self.assertNear(self.phase.Te, T)
            self.assertNear(self.phase.density, rho)
            self.assertArrayNear(self.phase.Y, Y)

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
            self.assertNear(getattr(self.phase, first), values[first])
            self.assertNear(getattr(self.phase, second), second_val)

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair, (None, values[second]))
            self.assertNear(getattr(self.phase, first), first_val)
            self.assertNear(getattr(self.phase, second), values[second])

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'X', (None, None, values['X']))
            self.assertNear(getattr(self.phase, first), first_val)
            self.assertNear(getattr(self.phase, second), second_val)

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'Y', (None, None, values['Y']))
            self.assertNear(getattr(self.phase, first), first_val)
            self.assertNear(getattr(self.phase, second), second_val)

    def test_setter_errors(self):
        with self.assertRaises(TypeError):
            self.phase.TD = 400

        with self.assertRaisesRegex(AssertionError, 'incorrect number'):
            self.phase.TP = 300, 101325, 'CH4:1.0'

        with self.assertRaisesRegex(AssertionError, 'incorrect number'):
            self.phase.HPY = 1.2e6, 101325

        with self.assertRaisesRegex(AssertionError, 'incorrect number'):
            self.phase.UVX = -4e5, 4.4, 'H2:1.0', -1

    def test_invalid_property(self):
        x = self.phase
        with self.assertRaises(AttributeError):
            x.foobar = 300
        with self.assertRaises(AttributeError):
            x.foobar

    def check_getters(self):
        T,D,X = self.phase.TDX
        self.assertNear(T, self.phase.T)
        self.assertNear(D, self.phase.density)
        self.assertArrayNear(X, self.phase.X)

        T,D,Y = self.phase.TDY
        self.assertNear(T, self.phase.T)
        self.assertNear(D, self.phase.density)
        self.assertArrayNear(Y, self.phase.Y)

        T,D = self.phase.TD
        self.assertNear(T, self.phase.T)
        self.assertNear(D, self.phase.density)

        T,P,X = self.phase.TPX
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        T,P,Y = self.phase.TPY
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        T,P = self.phase.TP
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)

        H,P,X = self.phase.HPX
        self.assertNear(H, self.phase.h)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        H,P,Y = self.phase.HPY
        self.assertNear(H, self.phase.h)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        H,P = self.phase.HP
        self.assertNear(H, self.phase.h)
        self.assertNear(P, self.phase.P)

        U,V,X = self.phase.UVX
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(X, self.phase.X)

        U,V,Y = self.phase.UVY
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(Y, self.phase.Y)

        U,V = self.phase.UV
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)

        S,P,X = self.phase.SPX
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        S,P,Y = self.phase.SPY
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        S,P = self.phase.SP
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)

        S,V,X = self.phase.SVX
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(X, self.phase.X)

        S,V,Y = self.phase.SVY
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(Y, self.phase.Y)

        S,V = self.phase.SV
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)

        D,P,X = self.phase.DPX
        self.assertNear(D, self.phase.density)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        D,P,Y = self.phase.DPY
        self.assertNear(D, self.phase.density)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        D,P = self.phase.DP
        self.assertNear(D, self.phase.density)
        self.assertNear(P, self.phase.P)

        Te = self.phase.Te
        self.assertNear(Te, self.phase.Te)

    def test_getState_mass(self):
        self.phase.TDY = 350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.check_getters()

    def test_getState_mole(self):
        self.phase.basis = 'molar'
        self.phase.TDX = 350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'
        self.check_getters()

    def test_getState(self):
        self.assertNear(self.phase.P, ct.one_atm)
        self.assertNear(self.phase.T, 300)

    def test_partial_molar(self):
        self.phase.TDY = 350.0, 0.6, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.assertNear(sum(self.phase.partial_molar_enthalpies * self.phase.X),
                        self.phase.enthalpy_mole)

        self.assertNear(sum(self.phase.partial_molar_entropies * self.phase.X),
                        self.phase.entropy_mole)

        self.assertNear(sum(self.phase.partial_molar_int_energies * self.phase.X),
                        self.phase.int_energy_mole)

        self.assertNear(sum(self.phase.chemical_potentials * self.phase.X),
                        self.phase.gibbs_mole)

        self.assertNear(sum(self.phase.partial_molar_cp * self.phase.X),
                        self.phase.cp_mole)

    def test_nondimensional(self):
        self.phase.TDY = 850.0, 0.2, 'H2:0.1, H2O:0.6, AR:0.3'
        H = (sum(self.phase.standard_enthalpies_RT * self.phase.X) *
             ct.gas_constant * self.phase.T)
        self.assertNear(H, self.phase.enthalpy_mole)

        U = (sum(self.phase.standard_int_energies_RT * self.phase.X) *
             ct.gas_constant * self.phase.T)
        self.assertNear(U, self.phase.int_energy_mole)

        cp = sum(self.phase.standard_cp_R * self.phase.X) * ct.gas_constant
        self.assertNear(cp, self.phase.cp_mole)

    def test_activities(self):
        self.phase.TDY = 850.0, 0.2, 'H2:0.1, H2O:0.6, AR:0.3'
        self.assertArrayNear(self.phase.X, self.phase.activities)

        self.assertArrayNear(self.phase.activity_coefficients,
                             np.ones(self.phase.n_species))

    def test_isothermal_compressibility(self):
        self.assertNear(self.phase.isothermal_compressibility, 1.0/self.phase.P)

    def test_thermal_expansion_coeff(self):
        self.assertNear(self.phase.thermal_expansion_coeff, 1.0/self.phase.T)

    def test_ref_info(self):
        self.assertNear(self.phase.reference_pressure, ct.one_atm)
        self.assertNear(self.phase.min_temp, 300.0)
        self.assertNear(self.phase.max_temp, 3500.0)

    def test_uncopyable(self):
        import copy
        with self.assertRaises(NotImplementedError):
            copy.copy(self.phase)

    def test_add_species(self):
        ref = ct.Solution('gri30.yaml', transport_model=None)
        n_orig = self.phase.n_species
        self.phase.add_species(ref.species('CO2'))
        self.phase.add_species(ref.species('CO'))

        self.assertEqual(self.phase.n_species, n_orig + 2)
        self.assertIn('CO2', self.phase.species_names)
        self.assertIn('CO', self.phase.species_names)

        state = 400, 2e5, 'H2:0.7, CO2:0.2, CO:0.1'
        ref.TPY = state
        self.phase.TPY = state
        self.assertNear(self.phase.enthalpy_mass, ref.enthalpy_mass)
        self.assertNear(self.phase.entropy_mole, ref.entropy_mole)
        self.assertArrayNear(ref[self.phase.species_names].partial_molar_cp,
                             self.phase.partial_molar_cp)

    def test_add_species_disabled(self):
        ref = ct.Solution('gri30.yaml', transport_model=None)

        self.phase.transport_model = "unity-Lewis-number"
        reactor = ct.IdealGasReactor(self.phase)
        with self.assertRaisesRegex(ct.CanteraError, 'Cannot add species'):
            self.phase.add_species(ref.species('CH4'))
        del reactor
        gc.collect()
        self.phase.add_species(ref.species('CH4'))

        flame = ct.FreeFlame(self.phase, width=0.1)
        with self.assertRaisesRegex(ct.CanteraError, 'Cannot add species'):
            self.phase.add_species(ref.species('CO'))
        del flame
        gc.collect()
        self.phase.add_species(ref.species('CO'))

        mix = ct.Mixture([(self.phase, 2.0)])
        with self.assertRaisesRegex(ct.CanteraError, 'Cannot add species'):
            self.phase.add_species(ref.species('CH2O'))
        del mix
        gc.collect()
        self.phase.add_species(ref.species('CH2O'))

    def test_add_species_duplicate(self):
        species = self.phase.species('H2O2')
        with self.assertRaisesRegex(ct.CanteraError, 'already contains'):
            self.phase.add_species(species)


class TestThermo(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.ThermoPhase("h2o2.yaml")
        self.gas.TPX = 450, 2e5, 'H2:1.0, O2:0.4, AR:3, H2O:0.1'

    def test_setSV_lowT(self):
        """
        Set state in terms of (s,v) when the end temperature is below the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        s1, v1 = self.gas.SV
        self.gas.SV = s1, 3 * v1

        self.assertNear(self.gas.s, s1)
        self.assertNear(self.gas.v, 3 * v1)
        self.assertTrue(self.gas.T < self.gas.min_temp)

    def test_setSV_low_invalid(self):
        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        self.gas.SV = 4600, None
        with self.assertRaises(ct.CanteraError):
            self.gas.SV = -1000, None

    def test_setSV_highT(self):
        """
        Set state in terms of (s,v) when the end temperature is above the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 2900, 1e5, 'H2:1.0, O2:0.4, AR:3'
        s1, v1 = self.gas.SV
        self.gas.SV = s1, 0.3 * v1

        self.assertNear(self.gas.s, s1)
        self.assertNear(self.gas.v, 0.3 * v1)
        self.assertTrue(self.gas.T > self.gas.max_temp)

    def test_setHP_lowT(self):
        """
        Set state in terms of (s,v) when the end temperature is below the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        deltaH = 1.25e5
        h1, p1 = self.gas.HP
        self.gas.HP = h1 - deltaH, None

        self.assertNear(self.gas.h, h1 - deltaH)
        self.assertNear(self.gas.P, p1)
        self.assertTrue(self.gas.T < self.gas.min_temp)

    def test_setHP_low_invalid(self):
        """
        Set state in terms of (h,p) when the enthalpy would imply a negative
        temperature
        """

        self.gas.TPX = 300, 101325, 'H2:1.0'
        with self.assertRaises(ct.CanteraError):
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

        self.assertNear(self.gas.h, h1 + deltaH)
        self.assertNear(self.gas.P, p1)
        self.assertTrue(self.gas.T > self.gas.max_temp)

    def test_volume(self):
        """
        This phase should follow the ideal gas law
        """
        g = self.gas
        self.assertNear(g.P, g.density_mole * ct.gas_constant * g.T)

        self.assertNear(
            g.P / g.density,
            ct.gas_constant / g.mean_molecular_weight * g.T)

        self.assertNear(g.density, 1.0 / g.volume_mass)

    def test_energy(self):
        g = self.gas
        mmw = g.mean_molecular_weight
        self.assertNear(g.enthalpy_mass, g.enthalpy_mole / mmw)
        self.assertNear(g.int_energy_mass, g.int_energy_mole / mmw)
        self.assertNear(g.gibbs_mass, g.gibbs_mole / mmw)
        self.assertNear(g.entropy_mass, g.entropy_mole / mmw)

        self.assertNear(g.cv_mass, g.cv_mole / mmw)
        self.assertNear(g.cp_mass, g.cp_mole / mmw)
        self.assertNear(g.cv_mole + ct.gas_constant, g.cp_mole)

    def test_nondimensional(self):
        g = self.gas
        R = ct.gas_constant

        self.assertNear(np.dot(g.standard_cp_R, g.X), g.cp_mole / R)
        self.assertNear(np.dot(g.standard_enthalpies_RT, g.X),
                        g.enthalpy_mole / (R*g.T))

        Smix_R = - np.dot(g.X, np.log(g.X+1e-20))
        self.assertNear(np.dot(g.standard_entropies_R, g.X) + Smix_R,
                        g.entropy_mole / R)
        self.assertNear(np.dot(g.standard_gibbs_RT, g.X) - Smix_R,
                        g.gibbs_mole / (R*g.T))


class TestInterfacePhase(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution("diamond.yaml", "gas")
        self.solid = ct.Solution("diamond.yaml", "diamond")
        self.interface = ct.Interface("diamond.yaml", "diamond_100",
                                      (self.gas, self.solid))

    def test_properties(self):
        self.interface.site_density = 100
        self.assertNear(self.interface.site_density, 100)

    def test_coverages_array(self):
        C = np.zeros(self.interface.n_species)
        C[1] = 0.25
        C[3] = 0.125
        C[4] = 0.125
        self.interface.coverages = C
        C = self.interface.coverages
        # should now be normalized
        self.assertNear(C[1], 0.5)
        self.assertNear(C[3], 0.25)
        self.assertNear(C[4], 0.25)
        self.assertNear(sum(C), 1.0)

    def test_coverages_string(self):
        self.interface.coverages = 'c6HM:0.2, c6H*:0.8'
        C = self.interface.coverages
        self.assertNear(C[self.interface.species_index('c6HM')], 0.2)
        self.assertNear(C[self.interface.species_index('c6H*')], 0.8)

    def test_coverages_dict(self):
        self.interface.coverages = {'c6**':1.0, 'c6*M':3.0}
        C = self.interface.coverages
        self.assertNear(C[self.interface.species_index('c6**')], 0.25)
        self.assertNear(C[self.interface.species_index('c6*M')], 0.75)


class TestPlasmaPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('oxygen-plasma.yaml',
                                 'isotropic-electron-energy-plasma',
                                 transport_model=None)

    def test_converting_electron_energy_to_temperature(self):
        self.phase.mean_electron_energy = 1.0
        Te = 2.0 / 3.0 * ct.electron_charge / ct.boltzmann
        self.assertNear(self.phase.Te, Te)

    def test_converting_electron_temperature_to_energy(self):
        self.phase.Te = 10000
        energy = self.phase.Te * 3.0 / 2.0 / ct.electron_charge * ct.boltzmann
        self.assertNear(self.phase.mean_electron_energy, energy)

    def test_set_get_electron_energy_levels(self):
        levels = np.linspace(0.01, 10, num=9)
        self.phase.electron_energy_levels = levels
        self.assertArrayNear(levels, self.phase.electron_energy_levels)

    def test_isotropic_velocity_electron_energy_distribution(self):
        levels = np.linspace(0.01, 10, num=9)
        self.phase.electron_energy_levels = levels
        self.phase.Te = 2e5
        mean_electron_energy = 3.0 / 2.0 * (self.phase.Te * ct.gas_constant /
                               (ct.avogadro * ct.electron_charge))
        self.assertNear(mean_electron_energy, self.phase.mean_electron_energy)

    def test_discretized_electron_energy_distribution(self):
        levels = np.array([0.0, 1.0, 10.0])
        dist = np.array([0.0, 0.9, 0.01])
        self.phase.normalize_electron_energy_distribution_enabled = False
        self.phase.quadrature_method = "trapezoidal"
        self.phase.set_discretized_electron_energy_distribution(levels, dist)
        self.assertArrayNear(levels, self.phase.electron_energy_levels)
        self.assertArrayNear(dist, self.phase.electron_energy_distribution)
        mean_energy = 2.0 / 5.0 * np.trapz(dist, np.power(levels, 5./2.))
        self.assertNear(self.phase.mean_electron_energy, mean_energy, 1e-4)
        electron_temp = 2.0 / 3.0 * (self.phase.mean_electron_energy *
                        ct.avogadro * ct.electron_charge / ct.gas_constant)
        self.assertNear(self.phase.Te, electron_temp)

    def test_electron_thermodynamic_properties(self):
        self.assertNear(self.phase.standard_gibbs_RT[0],
                        self.phase.standard_enthalpies_RT[0] -
                        self.phase.standard_entropies_R[0])

    def test_add_multiple_electron_species(self):
        electron = ct.Species('Electron', 'E:1')
        electron.thermo = ct.ConstantCp(100, 200, 101325, coeffs=(300, 1, 1, 1))
        with self.assertRaisesRegex(ct.CanteraError,
                                    'Only one electron species is allowed'):
            self.phase.add_species(electron)


class ImportTest(utilities.CanteraTest):
    """
    Test the various ways of creating a Solution object
    """
    def check(self, gas, phase, T, P, nSpec, nElem):
        self.assertEqual(gas.name, phase)
        self.assertNear(gas.T, T)
        self.assertNear(gas.P, P)
        self.assertEqual(gas.n_species, nSpec)
        self.assertEqual(gas.n_elements, nElem)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_import_phase_cti(self):
        gas1 = ct.Solution('air-no-reactions.cti', 'air')
        self.check(gas1, 'air', 300, 101325, 8, 3)

        gas2 = ct.Solution('air-no-reactions.cti', 'notair')
        self.check(gas2, 'notair', 900, 5*101325, 7, 2)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_import_phase_cti2(self):
        # This should import the first phase, that is, 'air'
        gas = ct.Solution('air-no-reactions.cti')
        self.check(gas, 'air', 300, 101325, 8, 3)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_import_phase_xml(self):
        gas1 = ct.Solution('air-no-reactions.xml', 'air')
        self.check(gas1, 'air', 300, 101325, 8, 3)

        gas2 = ct.Solution('air-no-reactions.xml', 'notair')
        self.check(gas2, 'notair', 900, 5*101325, 7, 2)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_import_phase_cti_text(self):
        cti_def = """
ideal_gas(name='spam', elements='O H',
          species='gri30: all',
          options='skip_undeclared_elements',
          initial_state=state(temperature=350, pressure=2e6))
"""
        gas = ct.Solution(source=cti_def)
        self.check(gas, 'spam', 350, 2e6, 8, 2)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_import_phase_xml_text(self):
        xml_def = """
<?xml version="1.0"?>
<ctml>
  <validate reactions="yes" species="yes"/>
  <phase dim="3" id="spam">
    <elementArray datasrc="elements.xml">O</elementArray>
    <speciesArray datasrc="gri30.xml#species_data">all
      <skip element="undeclared"/>
    </speciesArray>
    <state>
      <temperature units="K">350.0</temperature>
      <pressure units="Pa">2000000.0</pressure>
    </state>
    <thermo model="IdealGas"/>
    <kinetics model="GasKinetics"/>
    <transport model="None"/>
  </phase>
</ctml>"""
        gas = ct.Solution(source=xml_def)
        self.check(gas, 'spam', 350, 2e6, 2, 1)

    def test_import_from_species(self):
        gas1 = ct.Solution('h2o2.yaml', transport_model=None)
        gas1.TPX = 350, 101325, 'H2:0.3, O2:0.7'
        gas1.equilibrate('HP')

        species = ct.Species.list_from_file("h2o2.yaml")
        gas2 = ct.ThermoPhase(thermo='IdealGas', species=species)
        gas2.TPX = 350, 101325, 'H2:0.3, O2:0.7'
        gas2.equilibrate('HP')
        self.assertEqual(gas1.n_elements, gas2.n_elements)
        self.assertEqual(gas1.species_names, gas2.species_names)
        self.assertNear(gas1.T, gas2.T)
        self.assertArrayNear(gas1.X, gas2.X)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_checkReactionBalance(self):
        with self.assertRaisesRegex(ct.CanteraError, 'reaction is unbalanced'):
            ct.Solution('h2o2_unbalancedReaction.xml')

    def test_yaml_ideal_gas_simple(self):
        gas = ct.ThermoPhase('ideal-gas.yaml', 'simple')
        self.check(gas, 'simple', 500, 10 * ct.one_atm, 3, 2)

    def test_yaml_ideal_gas_remote_species(self):
        gas = ct.ThermoPhase('ideal-gas.yaml', 'species-remote')
        self.check(gas, 'species-remote', 300, ct.one_atm, 4, 2)

    def test_yaml_duplicate(self):
        with self.assertRaisesRegex(ct.CanteraError, 'duplicate'):
            gas = ct.ThermoPhase('ideal-gas.yaml', 'duplicate-species')


class TestSpecies(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('h2o2.yaml', transport_model=None)

    def test_standalone(self):
        s = ct.Species('CH4', {'C':1, 'H':4})

        self.assertEqual(s.name, 'CH4')
        c = s.composition
        self.assertEqual(len(c), 2)
        self.assertEqual(c['C'], 1)
        self.assertEqual(c['H'], 4)

    def test_defaults(self):
        s = ct.Species('H2')
        self.assertEqual(s.size, 1.0)
        self.assertEqual(s.charge, 0.0)

        self.assertIsNone(s.thermo)
        self.assertIsNone(s.transport)

    def test_index_accessor(self):
        for k in range(self.gas.n_species):
            s = self.gas.species(k)
            self.assertEqual(s.name, self.gas.species_name(k))

            for m,n in s.composition.items():
                self.assertEqual(n, self.gas.n_atoms(k,m))

    def test_species_noargs(self):
        for k,s in enumerate(self.gas.species()):
            self.assertEqual(s.name, self.gas.species_name(k))

    def test_name_accessor(self):
        for name in self.gas.species_names:
            s = self.gas.species(name)
            self.assertEqual(s.name, name)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_fromCti(self):
        h2_cti = """
            species(
                name="H2",
                atoms="H:2",
                thermo=(
                    NASA([200.00, 1000.00],
                         [2.344331120E+00, 7.980520750E-03, -1.947815100E-05,
                          2.015720940E-08, -7.376117610E-12, -9.179351730E+02,
                          6.830102380E-01]),
                    NASA([1000.00, 3500.00],
                         [3.337279200E+00, -4.940247310E-05, 4.994567780E-07,
                         -1.795663940E-10, 2.002553760E-14, -9.501589220E+02,
                         -3.205023310E+00])
                ),
                transport=gas_transport(geom="linear",
                                        diam=2.92,
                                        well_depth=38.00,
                                        polar=0.79,
                                        rot_relax=280.00),
                note = "TPIS78"
            )"""
        s1 = self.gas.species('H2')
        s2 = ct.Species.fromCti(h2_cti)
        self.assertEqual(s2.name, 'H2')
        self.assertEqual(s1.composition, s2.composition)
        self.assertEqual(s1.thermo.cp(350), s2.thermo.cp(350))

    @pytest.mark.usefixtures("allow_deprecated")
    def test_fromXml(self):
        import xml.etree.ElementTree as ET
        root = ET.parse(self.cantera_data_path / "h2o2.xml").getroot()
        h2_node = root.find('.//species[@name="H2"]')
        h2_string = ET.tostring(h2_node)

        s1 = self.gas.species('H2')
        s2 = ct.Species.fromXml(h2_string)

        self.assertEqual(s2.name, 'H2')
        self.assertEqual(s1.composition, s2.composition)
        self.assertEqual(s1.thermo.cp(350), s2.thermo.cp(350))

    @pytest.mark.usefixtures("allow_deprecated")
    def test_listFromFile_cti(self):
        S = ct.Species.listFromFile('h2o2.cti')
        self.assertEqual(S[3].name, self.gas.species_name(3))

    @pytest.mark.usefixtures("allow_deprecated")
    def test_listFromFile_xml(self):
        S = ct.Species.listFromFile('h2o2.xml')
        self.assertEqual(S[3].name, self.gas.species_name(3))

    def test_listfromFile_yaml(self):
        S = ct.Species.list_from_file("h2o2.yaml")
        self.assertEqual({sp.name for sp in S}, set(self.gas.species_names))

    @pytest.mark.usefixtures("allow_deprecated")
    def test_listFromCti(self):
        S = ct.Species.listFromCti((self.cantera_data_path / "h2o2.cti").read_text())
        self.assertEqual(S[3].name, self.gas.species_name(3))

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
        self.assertEqual(species[0].name, 'H2O')
        self.assertEqual(species[1].composition, {'H': 1, 'O': 2})
        self.assertNear(species[0].thermo.h(300), 100)

    def test_list_from_yaml_section(self):
        species = ct.Species.list_from_yaml(
            (self.test_data_path / "ideal-gas.yaml").read_text(),
            'species')

        self.assertEqual(species[0].name, 'O2')
        self.assertEqual(species[1].composition, {'N': 1, 'O': 1})

    def test_from_yaml(self):
        yaml = """
        name: H2O
        composition: {H: 2, O: 1}
        thermo: {model: constant-cp, h0: 100}
        """
        species = ct.Species.from_yaml(yaml)
        self.assertEqual(species.name, 'H2O')
        self.assertEqual(species.composition, {'H': 2, 'O': 1})
        self.assertNear(species.thermo.h(300), 100)

    def test_from_dict(self):
        data = {
            "name": "H2O",
            "composition": {"H": 2, "O": 1},
            "thermo": {"model": "constant-cp", "h0": 100},
        }
        species = ct.Species.from_dict(data)
        self.assertEqual(species.name, 'H2O')
        self.assertEqual(species.composition, {'H': 2, 'O': 1})
        self.assertNear(species.thermo.h(300), 100)

    @pytest.mark.usefixtures("allow_deprecated")
    def test_listFromXml(self):
        S = ct.Species.listFromXml((self.cantera_data_path / "h2o2.xml").read_text())
        self.assertEqual(S[4].name, self.gas.species_name(4))

    def test_modify_thermo(self):
        S = {sp.name: sp for sp in ct.Species.list_from_file("h2o2.yaml")}
        self.gas.TPX = 400, 2*ct.one_atm, 'H2:1.0'
        g0 = self.gas.gibbs_mole

        self.gas.TPX = None, None, 'O2:1.0'
        self.assertNotAlmostEqual(g0, self.gas.gibbs_mole)
        # Replace O2 thermo with the data from H2
        S['O2'].thermo = S['H2'].thermo
        self.gas.modify_species(self.gas.species_index('O2'), S['O2'])
        self.assertNear(g0, self.gas.gibbs_mole)

    def test_modify_thermo_invalid(self):
        S = {sp.name: sp for sp in ct.Species.list_from_file("h2o2.yaml")}

        orig = S['H2']
        thermo = orig.thermo
        copy = ct.Species('foobar', orig.composition)
        copy.thermo = thermo
        with self.assertRaisesRegex(ct.CanteraError, 'modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

        copy = ct.Species('H2', {'H': 3})
        copy.thermo = thermo
        with self.assertRaisesRegex(ct.CanteraError, 'modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

        copy = ct.Species('H2', orig.composition)
        copy.thermo = ct.ConstantCp(thermo.min_temp, thermo.max_temp,
            thermo.reference_pressure, [300, 123, 456, 789])
        with self.assertRaisesRegex(ct.CanteraError, 'modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

        copy = ct.Species('H2', orig.composition)
        copy.thermo = ct.NasaPoly2(thermo.min_temp+200, thermo.max_temp,
            thermo.reference_pressure, thermo.coeffs)
        with self.assertRaisesRegex(ct.CanteraError, 'modifySpecies'):
            self.gas.modify_species(self.gas.species_index('H2'), copy)

    def test_alias(self):
        self.gas.add_species_alias('H2', 'hydrogen')
        self.assertEqual(self.gas.species_index('hydrogen'), 0)
        self.gas.X = 'hydrogen:.5, O2:.5'
        self.assertNear(self.gas.X[0], 0.5)
        with self.assertRaisesRegex(ct.CanteraError, 'Invalid alias'):
            self.gas.add_species_alias('H2', 'O2')
        with self.assertRaisesRegex(ct.CanteraError, 'Unable to add alias'):
            self.gas.add_species_alias('spam', 'eggs')

    def test_isomers(self):
        gas = ct.Solution('nDodecane_Reitz.yaml')
        iso = gas.find_isomers({'C':4, 'H':9, 'O':2})
        self.assertEqual(len(iso), 2)
        iso = gas.find_isomers('C:4, H:9, O:2')
        self.assertEqual(len(iso), 2)
        iso = gas.find_isomers({'C':7, 'H':15})
        self.assertEqual(len(iso), 1)
        iso = gas.find_isomers({'C':7, 'H':16})
        self.assertEqual(len(iso), 0)


class TestSpeciesThermo(utilities.CanteraTest):
    h2o_coeffs = [
        1000.0, 3.03399249E+00, 2.17691804E-03, -1.64072518E-07,
        -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00,
        4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09,
        1.77197817E-12, -3.02937267E+04, -8.49032208E-01
    ]
    def setUp(self):
        self.gas = ct.Solution('h2o2.yaml', transport_model=None)
        self.gas.X = 'H2O:1.0'

    def test_create(self):
        st = ct.NasaPoly2(300, 3500, 101325, self.h2o_coeffs)

        for T in [300, 500, 900, 1200, 2000]:
            self.gas.TP = T, 101325
            self.assertNear(st.cp(T), self.gas.cp_mole)
            self.assertNear(st.h(T), self.gas.enthalpy_mole)
            self.assertNear(st.s(T), self.gas.entropy_mole)

    def test_invalid(self):
        with self.assertRaisesRegex(ValueError, 'incorrect length'):
            # not enough coefficients
            st = ct.NasaPoly2(300, 3500, 101325,
                              [1000.0, 3.03399249E+00, 2.17691804E-03])

    def test_wrap(self):
        st = self.gas.species('H2O').thermo

        self.assertIsInstance(st, ct.NasaPoly2)

        for T in [300, 500, 900, 1200, 2000]:
            self.gas.TP = T, 101325
            self.assertNear(st.cp(T), self.gas.cp_mole)
            self.assertNear(st.h(T), self.gas.enthalpy_mole)
            self.assertNear(st.s(T), self.gas.entropy_mole)

    def test_coeffs(self):
        st = ct.NasaPoly2(300, 3500, 101325, self.h2o_coeffs)
        self.assertEqual(st.min_temp, 300)
        self.assertEqual(st.max_temp, 3500)
        self.assertEqual(st.reference_pressure, 101325)
        self.assertArrayNear(self.h2o_coeffs, st.coeffs)
        self.assertEqual(st.n_coeffs, len(st.coeffs))
        self.assertTrue(st._check_n_coeffs(st.n_coeffs))

    def test_nasa9_load(self):
        gas = ct.Solution("airNASA9.yaml")
        st = gas.species(3).thermo
        self.assertIsInstance(st, ct.Nasa9PolyMultiTempRegion)
        self.assertEqual(st.n_coeffs, len(st.coeffs))
        self.assertTrue(st._check_n_coeffs(st.n_coeffs))

    def test_nasa9_create(self):
        gas = ct.Solution("airNASA9.yaml")
        st = gas.species(3).thermo
        t_min = st.min_temp
        t_max = st.max_temp
        p_ref = st.reference_pressure
        coeffs = st.coeffs
        st2 = ct.Nasa9PolyMultiTempRegion(t_min, t_max, p_ref, coeffs)
        self.assertIsInstance(st2, ct.Nasa9PolyMultiTempRegion)
        self.assertEqual(st.min_temp, t_min)
        self.assertEqual(st.max_temp, t_max)
        self.assertEqual(st.reference_pressure, p_ref)
        for T in range(300, 20000, 1000):
            self.assertNear(st.cp(T), st2.cp(T))
            self.assertNear(st.h(T), st2.h(T))
            self.assertNear(st.s(T), st2.s(T))

    def test_shomate_load(self):
        sol = ct.Solution('thermo-models.yaml', 'molten-salt-Margules')
        st = sol.species(0).thermo
        self.assertIsInstance(st, ct.ShomatePoly2)
        self.assertEqual(st.n_coeffs, len(st.coeffs))
        self.assertTrue(st._check_n_coeffs(st.n_coeffs))

    def test_shomate_create(self):
        sol = ct.Solution('thermo-models.yaml', 'molten-salt-Margules')
        st = sol.species(0).thermo
        t_min = st.min_temp
        t_max = st.max_temp
        p_ref = st.reference_pressure
        coeffs = st.coeffs
        st2 = ct.ShomatePoly2(t_min, t_max, p_ref, coeffs)
        self.assertIsInstance(st2, ct.ShomatePoly2)
        self.assertEqual(st.min_temp, t_min)
        self.assertEqual(st.max_temp, t_max)
        self.assertEqual(st.reference_pressure, p_ref)
        for T in [300, 500, 700, 900]:
            self.assertNear(st.cp(T), st2.cp(T))
            self.assertNear(st.h(T), st2.h(T))
            self.assertNear(st.s(T), st2.s(T))

    def test_piecewise_gibbs_load(self):
        sol = ct.Solution('thermo-models.yaml', 'HMW-NaCl-electrolyte')
        st = sol.species(1).thermo
        self.assertIsInstance(st, ct.Mu0Poly)
        self.assertEqual(st.n_coeffs, len(st.coeffs))
        self.assertTrue(st._check_n_coeffs(st.n_coeffs))

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
        self.assertIsInstance(st2, ct.Mu0Poly)
        self.assertEqual(st2.n_coeffs, len(coeffs))
        self.assertEqual(st2.n_coeffs, len(st2.coeffs))

    def test_piecewise_gibbs_create2(self):
        sol = ct.Solution('thermo-models.yaml', 'HMW-NaCl-electrolyte')
        st = sol.species(1).thermo
        t_min = st.min_temp
        t_max = st.max_temp
        p_ref = st.reference_pressure
        coeffs = st.coeffs
        st2 = ct.Mu0Poly(t_min, t_max, p_ref, coeffs)
        self.assertIsInstance(st2, ct.Mu0Poly)
        self.assertEqual(st.min_temp, t_min)
        self.assertEqual(st.max_temp, t_max)
        self.assertEqual(st.reference_pressure, p_ref)
        for T in [300, 500, 700, 900]:
            self.assertNear(st.cp(T), st2.cp(T))
            self.assertNear(st.h(T), st2.h(T))
            self.assertNear(st.s(T), st2.s(T))


class TestQuantity(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('gri30.yaml', transport_model=None)

    def setUp(self):
        self.gas.TPX = 300, 101325, 'O2:1.0, N2:3.76'

    def test_mass_moles(self):
        q1 = ct.Quantity(self.gas, mass=5)
        self.assertNear(q1.mass, 5)
        self.assertNear(q1.moles, 5 / q1.mean_molecular_weight)

        q1.mass = 7
        self.assertNear(q1.moles, 7 / q1.mean_molecular_weight)

        q1.moles = 9
        self.assertNear(q1.moles, 9)
        self.assertNear(q1.mass, 9 * q1.mean_molecular_weight)

    def test_extensive(self):
        q1 = ct.Quantity(self.gas, mass=5)
        self.assertNear(q1.mass, 5)

        self.assertNear(q1.volume * q1.density, q1.mass)
        self.assertNear(q1.V * q1.density, q1.mass)
        self.assertNear(q1.int_energy, q1.moles * q1.int_energy_mole)
        self.assertNear(q1.enthalpy, q1.moles * q1.enthalpy_mole)
        self.assertNear(q1.entropy, q1.moles * q1.entropy_mole)
        self.assertNear(q1.gibbs, q1.moles * q1.gibbs_mole)
        self.assertNear(q1.int_energy, q1.U)
        self.assertNear(q1.enthalpy, q1.H)
        self.assertNear(q1.entropy, q1.S)
        self.assertNear(q1.gibbs, q1.G)

    def test_multiply(self):
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = q1 * 2.5
        self.assertNear(q1.mass * 2.5, q2.mass)
        self.assertNear(q1.moles * 2.5, q2.moles)
        self.assertNear(q1.entropy * 2.5, q2.entropy)
        self.assertArrayNear(q1.X, q2.X)

    def test_multiply_HP(self):
        self.gas.TPX = 500, 101325, 'CH4:1.0, O2:0.4'
        q1 = ct.Quantity(self.gas, mass=2, constant='HP')
        q2 = ct.Quantity(self.gas, mass=1, constant='HP')
        q2.equilibrate('HP')
        q3 = 0.2 * q1 + q2 * 0.4
        self.assertNear(q1.P, q3.P)
        self.assertNear(q1.enthalpy_mass, q3.enthalpy_mass)
        self.assertNear(q2.enthalpy_mass, q3.enthalpy_mass)

    def test_iadd(self):
        q0 = ct.Quantity(self.gas, mass=5)
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = ct.Quantity(self.gas, mass=5)
        q2.TPX = 500, 101325, 'CH4:1.0'

        q1 += q2
        self.assertNear(q0.mass + q2.mass, q1.mass)
        # addition is at constant UV
        self.assertNear(q0.U + q2.U, q1.U)
        self.assertNear(q0.V + q2.V, q1.V)
        self.assertArrayNear(q0.X*q0.moles + q2.X*q2.moles, q1.X*q1.moles)

    def test_add(self):
        q1 = ct.Quantity(self.gas, mass=5)
        q2 = ct.Quantity(self.gas, mass=5)
        q2.TPX = 500, 101325, 'CH4:1.0'

        q3 = q1 + q2
        self.assertNear(q1.mass + q2.mass, q3.mass)
        # addition is at constant UV
        self.assertNear(q1.U + q2.U, q3.U)
        self.assertNear(q1.V + q2.V, q3.V)
        self.assertArrayNear(q1.X*q1.moles + q2.X*q2.moles, q3.X*q3.moles)

    def test_equilibrate(self):
        self.gas.TPX = 300, 101325, 'CH4:1.0, O2:0.2, N2:1.0'
        q1 = ct.Quantity(self.gas)
        self.gas.equilibrate('HP')
        T2 = self.gas.T

        self.assertNear(q1.T, 300)
        q1.equilibrate('HP')
        self.assertNear(q1.T, T2)

    def test_invalid_setter(self):
        q1 = ct.Quantity(self.gas, mass =3)
        with self.assertRaises(AttributeError):
            q1.HPQ = self.gas.H, self.gas.P, 1

    def test_incompatible(self):
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        q1 = ct.Quantity(self.gas)
        q2 = ct.Quantity(gas2)

        with self.assertRaisesRegex(ValueError, 'different phase definitions'):
            q1+q2


class TestMisc(utilities.CanteraTest):
    def test_stringify_bad(self):
        with self.assertRaises(AttributeError):
            ct.Solution(3)

    def test_case_sensitive_names(self):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        self.assertFalse(gas.case_sensitive_species_names)
        self.assertEqual(gas.species_index('h2'), 0)
        gas.X = 'h2:.5, o2:.5'
        self.assertNear(gas.X[0], 0.5)
        gas.Y = 'h2:.5, o2:.5'
        self.assertNear(gas.Y[0], 0.5)

        gas.case_sensitive_species_names = True
        self.assertTrue(gas.case_sensitive_species_names)
        with self.assertRaises(ValueError):
            gas.species_index('h2')
        with self.assertRaisesRegex(ct.CanteraError, 'Unknown species'):
            gas.X = 'h2:1.0, o2:1.0'
        with self.assertRaisesRegex(ct.CanteraError, 'Unknown species'):
            gas.Y = 'h2:1.0, o2:1.0'

        gas_cti = """ideal_gas(
            name="gas",
            elements=" S C Cs ",
            species=" nasa: all ",
            options=["skip_undeclared_elements"],
            initial_state=state(temperature=300, pressure=(1, "bar"))
        )"""
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
        with self.assertRaisesRegex(ct.CanteraError, 'is not unique'):
            gas.species_index('cs')
        gas.case_sensitive_species_names = True
        with self.assertRaises(ValueError):
            gas.species_index('cs')


class TestElement(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.ar_sym = ct.Element('Ar')
        cls.ar_name = ct.Element('argon')
        cls.ar_num = ct.Element(18)

    def test_element_multiple_possibilities(self):
        carbon = ct.Element('Carbon')
        self.assertEqual(carbon.name, 'carbon')
        self.assertEqual(carbon.symbol, 'C')

    def test_element_weight(self):
        self.assertNear(self.ar_sym.weight, 39.95)
        self.assertNear(self.ar_name.weight, 39.95)
        self.assertNear(self.ar_num.weight, 39.95)

    def test_element_symbol(self):
        self.assertEqual(self.ar_sym.symbol, 'Ar')
        self.assertEqual(self.ar_name.symbol, 'Ar')
        self.assertEqual(self.ar_num.symbol, 'Ar')

    def test_element_name(self):
        self.assertEqual(self.ar_sym.name, 'argon')
        self.assertEqual(self.ar_name.name, 'argon')
        self.assertEqual(self.ar_num.name, 'argon')

    def test_element_atomic_number(self):
        self.assertEqual(self.ar_sym.atomic_number, 18)
        self.assertEqual(self.ar_name.atomic_number, 18)
        self.assertEqual(self.ar_num.atomic_number, 18)

    def test_element_name_not_present(self):
        with self.assertRaisesRegex(ct.CanteraError, 'element not found'):
            ct.Element('I am not an element')

    def test_element_atomic_number_small(self):
        with self.assertRaisesRegex(ct.CanteraError, 'IndexError'):
            ct.Element(0)

    def test_element_atomic_number_big(self):
        num_elements = ct.Element.num_elements_defined
        with self.assertRaisesRegex(ct.CanteraError, 'IndexError'):
            ct.Element(num_elements + 1)

    def test_element_no_weight(self):
        with self.assertRaisesRegex(ct.CanteraError, 'no stable isotopes'):
            ct.Element('Tc')

    def test_element_bad_input(self):
        with self.assertRaisesRegex(TypeError, 'input argument to Element'):
            ct.Element(1.2345)

    def test_get_isotope(self):
        d_sym = ct.Element('D')
        self.assertEqual(d_sym.atomic_number, 1)
        self.assertNear(d_sym.weight, 2.0141017781)
        self.assertEqual(d_sym.name, 'deuterium')
        self.assertEqual(d_sym.symbol, 'D')

        d_name = ct.Element('deuterium')
        self.assertEqual(d_name.atomic_number, 1)
        self.assertNear(d_name.weight, 2.0141017781)
        self.assertEqual(d_name.name, 'deuterium')
        self.assertEqual(d_name.symbol, 'D')

    def test_elements_lists(self):
        syms = ct.Element.element_symbols
        names = ct.Element.element_names
        num_elements = ct.Element.num_elements_defined
        self.assertEqual(len(syms), num_elements)
        self.assertEqual(len(names), num_elements)


class TestSolutionArray(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('h2o2.yaml')

    def test_passthrough(self):
        states = ct.SolutionArray(self.gas, 3)
        self.assertEqual(states.n_species, self.gas.n_species)
        self.assertEqual(states.reaction(10).equation,
                         self.gas.reaction(10).equation)

    def test_meta(self):
        meta = {'foo': 'bar', 'spam': 'eggs'}
        states = ct.SolutionArray(self.gas, 3, meta=meta)
        self.assertEqual(states.meta['foo'], 'bar')
        self.assertEqual(states.meta['spam'], 'eggs')

    def test_get_state(self):
        states = ct.SolutionArray(self.gas, 4)
        H, P = states.HP
        self.assertEqual(H.shape, (4,))
        self.assertEqual(P.shape, (4,))

        S, P, Y = states.SPY
        self.assertEqual(S.shape, (4,))
        self.assertEqual(P.shape, (4,))
        self.assertEqual(Y.shape, (4, self.gas.n_species))

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
            self.assertEqual(len(out), len(g))

    def test_properties_onedim(self):
        N = 11
        states = ct.SolutionArray(self.gas, N)
        T = np.linspace(300, 2200, N)
        P = np.logspace(3, 8, N)
        X = 'H2:0.5, O2:0.4, AR:0.1, H2O2:0.01, OH:0.001'
        states.TPX = T, P, X

        self.assertArrayNear(states.T, T)
        self.assertArrayNear(states.P, P)

        h = states.enthalpy_mass
        ropr = states.reverse_rates_of_progress
        Dkm = states.mix_diff_coeffs
        for i in range(N):
            self.gas.TPX = T[i], P[i], X
            self.assertNear(self.gas.enthalpy_mass, h[i])
            self.assertArrayNear(self.gas.reverse_rates_of_progress, ropr[i])
            self.assertArrayNear(self.gas.mix_diff_coeffs, Dkm[i])

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

        self.assertEqual(h.shape, (2,3,5))
        self.assertEqual(ropr.shape, (2,3,5,self.gas.n_reactions))
        self.assertEqual(Dkm.shape, (2,3,5,self.gas.n_species))

        for i,j,k in np.ndindex(TT.shape):
            self.gas.TPX = T[k], P[i], X[j]
            self.assertNear(self.gas.enthalpy_mass, h[i,j,k])
            self.assertArrayNear(self.gas.reverse_rates_of_progress, ropr[i,j,k])
            self.assertArrayNear(self.gas.mix_diff_coeffs, Dkm[i,j,k])

    def test_slicing_onedim(self):
        states = ct.SolutionArray(self.gas, 5)
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        T0 = states.T
        H0 = states.enthalpy_mass

        # Verify that original object is updated when slices change
        state = states[1]
        state.TD = 300, 0.5
        self.assertNear(states.T[0], 500)
        self.assertNear(states.T[1], 300)
        self.assertNear(states.P[2], 2e5)
        self.assertNear(states.density[1], 0.5)

        # Verify that the slices are updated when the original object changes
        states.TD = 900, None
        self.assertNear(state.T, 900)
        self.assertNear(states.density[1], 0.5)

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
        self.assertArrayNear(T[0], T0[0])
        self.assertArrayNear(T[1], 300*np.ones(5))
        self.assertArrayNear(D[1], 0.5*np.ones(5))

        col3 = states[:,2]
        col3.TD = 400, 2.5
        T = states.T
        D = states.density
        self.assertArrayNear(T[:,2], 400*np.ones(2))
        self.assertArrayNear(D[:,2], 2.5*np.ones(2))

        # Verify that the slices are updated when the original object changes
        states.TP = 900, None
        self.assertArrayNear(col3.T, 900*np.ones(2))
        self.assertArrayNear(row2.T, 900*np.ones(5))

    def test_extra_create_by_dict(self):
        extra = OrderedDict([('grid', np.arange(10)),
                             ('velocity', np.random.rand(10))])
        states = ct.SolutionArray(self.gas, 10, extra=extra)
        keys = list(states._extra.keys())
        self.assertEqual(keys[0], 'grid')
        self.assertArrayNear(states.grid, np.arange(10))

    def test_extra_no_shape(self):
        # The shape of the value for "prop" here is (), which is falsey
        # and causes the use of np.full()
        states = ct.SolutionArray(self.gas, 3, extra={"prop": 1})
        self.assertEqual(states.prop.shape, (3,))
        self.assertArrayNear(states.prop, np.array((1, 1, 1)))

        # Check a multidimensional SolutionArray
        states = ct.SolutionArray(self.gas, (2, 2), extra={"prop": 3})
        self.assertEqual(states.prop.shape, (2, 2))
        self.assertArrayNear(states.prop, np.array(((3, 3), (3, 3))))

    def test_extra_not_empty(self):
        """Test that a non-empty SolutionArray raises a ValueError if
           initial values for properties are not supplied.
        """
        with self.assertRaisesRegex(ValueError, "Initial values for extra properties"):
            ct.SolutionArray(self.gas, 3, extra=["prop"])
        with self.assertRaisesRegex(ValueError, "Initial values for extra properties"):
            ct.SolutionArray(self.gas, 3, extra=np.array(["prop", "prop2"]))

    def test_extra_create_multidim(self):
        # requires matching first dimensions
        extra_val = [[1, 2, 3] for i in range(5)]
        states = ct.SolutionArray(self.gas, 5, extra={"prop": extra_val})
        self.assertEqual(states.prop.shape, (5, 3,))
        states = ct.SolutionArray(self.gas, 5, extra={"prop": np.array(extra_val)})
        self.assertEqual(states.prop.shape, (5, 3,))
        states = ct.SolutionArray(self.gas, (3, 4,), extra={"prop": np.ones((3, 4, 5,))})
        self.assertEqual(states.prop.shape, (3, 4, 5,))
        states = ct.SolutionArray(self.gas, 1, extra={"prop": extra_val[0]})
        self.assertEqual(states.prop.shape, (1, 3,))
        states = ct.SolutionArray(self.gas, 1, extra={"prop": [2]})
        self.assertEqual(states.prop.shape, (1,))
        with self.assertRaisesRegex(ValueError, "Unable to map"):
            ct.SolutionArray(self.gas, (3, 3), extra={"prop": np.arange(3)})

    def test_extra_create_by_iterable(self):
        states = ct.SolutionArray(self.gas, extra=("prop1"))
        self.assertEqual(states.prop1.shape, (0,))

        # An integer is not an iterable, and only bare strings are
        # turned into iterables
        with self.assertRaisesRegex(ValueError, "Extra properties"):
            ct.SolutionArray(self.gas, extra=2)

    def test_extra_not_string(self):
        with self.assertRaisesRegex(TypeError, "is not a string"):
            ct.SolutionArray(self.gas, extra=[1])

    def test_extra_no_objects(self):
        with self.assertRaisesRegex(ValueError, "not supported"):
            prop = np.array([0, [1, 2], (3, 4)], dtype=object)
            states = ct.SolutionArray(self.gas, 3, extra={"prop": prop})

    def test_extra_reserved_names(self):
        with self.assertRaisesRegex(ValueError, "name is already used"):
            ct.SolutionArray(self.gas, extra=["creation_rates"])

        with self.assertRaisesRegex(ValueError, "name is already used"):
            ct.SolutionArray(self.gas, extra={"creation_rates": 0})

    def test_extra_create_by_string(self):
        states = ct.SolutionArray(self.gas, extra="prop")
        self.assertEqual(states.prop.shape, (0,))

    def test_extra_setattr(self):
        states = ct.SolutionArray(self.gas, 7, extra={'prop': range(7)})
        states.prop = 0
        self.assertArrayNear(states.prop, np.zeros((7,)))
        mod_array = np.linspace(0, 10, 7)
        states.prop = mod_array
        self.assertArrayNear(states.prop, mod_array)
        with self.assertRaisesRegex(ValueError, "Incompatible shapes"):
            states.prop = [1, 2]

    def test_assign_to_slice(self):
        states = ct.SolutionArray(self.gas, 7, extra={'prop': range(7)})
        array = np.arange(7)
        self.assertArrayNear(states.prop, array)
        states.prop[1] = -5
        states.prop[3:5] = [0, 1]
        array_mod = np.array([0, -5, 2, 0, 1, 5, 6])
        self.assertArrayNear(states.prop, array_mod)
        # assign to multi-dimensional extra
        extra_val = [[1, 2, 3] for i in range(5)]
        states = ct.SolutionArray(self.gas, 5, extra={"prop": extra_val})
        states.prop[:, 1] = -1
        self.assertArrayNear(states.prop[:, 1], -1 * np.ones((5,)))
        states.prop[2, :] = -2
        self.assertArrayNear(states.prop[2, :], -2 * np.ones((3,)))

    def test_extra_create_by_ndarray(self):
        properties_array = np.array(["prop1", "prop2", "prop3"])
        states = ct.SolutionArray(self.gas, shape=(0,), extra=properties_array)
        self.assertEqual(states.prop1.shape, (0,))
        self.assertEqual(states.prop2.shape, (0,))
        self.assertEqual(states.prop3.shape, (0,))
        # Ensure that a 2-dimensional array is flattened
        properties_array = np.array((["prop1"], ["prop2"]))
        states = ct.SolutionArray(self.gas, extra=properties_array)
        self.assertEqual(states.prop1.shape, (0,))
        self.assertEqual(states.prop2.shape, (0,))

    def test_append(self):
        states = ct.SolutionArray(self.gas, 5)
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        self.assertEqual(states.cp_mass.shape, (5,))

        states.append(T=1100, P=3e5, X='AR:1.0')
        self.assertEqual(states.cp_mass.shape, (6,))
        self.assertNear(states.P[-1], 3e5)
        self.assertNear(states.T[-1], 1100)

        self.gas.TPX = 1200, 5e5, 'O2:0.3, AR:0.7'
        states.append(self.gas.state)
        self.assertEqual(states.cp_mass.shape, (7,))
        self.assertNear(states.P[-1], 5e5)
        self.assertNear(states.X[-1, self.gas.species_index('AR')], 0.7)

        self.gas.TPX = 300, 1e4, 'O2:0.5, AR:0.5'
        HPY = self.gas.HPY
        self.gas.TPX = 1200, 5e5, 'O2:0.3, AR:0.7'  # to make sure it gets changed
        states.append(HPY=HPY)
        self.assertEqual(states.cp_mass.shape, (8,))
        self.assertNear(states.P[-1], 1e4)
        self.assertNear(states.T[-1], 300)

    def test_append_with_extra(self):
        states = ct.SolutionArray(self.gas, 5, extra={"prop": "value"})
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        self.assertEqual(states._shape, (5,))
        states.append(T=1100, P=3e5, X="AR:1.0", prop="value2")
        self.assertEqual(states.prop[-1], "value2")
        self.assertEqual(states.prop.shape, (6,))
        states.append(T=1100, P=3e5, X="AR:1.0", prop=100)
        # NumPy converts to the existing type of the array
        self.assertEqual(states.prop[-1], "100")
        self.assertEqual(states.prop.shape, (7,))
        # two-dimensional input array
        states = ct.SolutionArray(self.gas, 1, extra={"prop": [1, 2, 3]})
        states.append(T=1100, P=3e5, X="AR:1.0", prop=['a', 'b', 'c'])
        self.assertEqual(states._shape, (2,))

    def test_append_failures(self):
        states = ct.SolutionArray(self.gas, 5, extra={"prop": "value"})
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        self.assertEqual(states._shape, (5,))

        with self.assertRaisesRegex(TypeError, "Missing keyword arguments for extra"):
            states.append(T=1100, P=3e5, X="AR:1.0")
        # Failing to append a state shouldn't change the size
        self.assertEqual(states._shape, (5,))

        with self.assertRaisesRegex(KeyError, "does not specify"):
            # I is not a valid property
            states.append(TPI=(1100, 3e5, "AR:1.0"), prop="value2")
        # Failing to append a state shouldn't change the size
        self.assertEqual(states._shape, (5,))

        with self.assertRaisesRegex(KeyError, "is not a valid"):
            # I is not a valid property
            states.append(T=1100, P=3e5, I="AR:1.0", prop="value2")
        # Failing to append a state shouldn't change the size
        self.assertEqual(states._shape, (5,))

        with self.assertRaisesRegex(ValueError, "incompatible value"):
            # prop has incompatible dimensions
            states.append(T=1100, P=3e5, X="AR:1.0", prop=[1, 2, 3])
        # Failing to append a state shouldn't change the size
        self.assertEqual(states._shape, (5,))

        states = ct.SolutionArray(self.gas, 1, extra={"prop": [1, 2, 3]})
        with self.assertRaisesRegex(ValueError, "does not match"):
            # prop has incorrect dimensions
            states.append(T=1100, P=3e5, X="AR:1.0", prop=['a', 'b'])
        # Failing to append a state shouldn't change the size
        self.assertEqual(states._shape, (1,))

    def test_purefluid(self):
        water = ct.Water()
        states = ct.SolutionArray(water, 5)
        states.TQ = 400, np.linspace(0, 1, 5)

        P = states.P
        for i in range(1, 5):
            self.assertNear(P[0], P[i])

        states.TP = np.linspace(400, 500, 5), 101325
        self.assertArrayNear(states.Q.squeeze(), np.ones(5))

    def test_phase_of_matter(self):
        water = ct.Water()
        states = ct.SolutionArray(water, 5)
        T = [300, 500, water.critical_temperature*2, 300]
        P = [101325, 101325, 101325, water.critical_pressure*2]
        states[:4].TP = T, P
        states[4].TQ = 300, .4
        pom = ['liquid', 'gas', 'supercritical', 'supercritical', 'liquid-gas-mix']
        self.assertEqual(list(states.phase_of_matter), pom)

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
            self.assertEqual(len(out), len(g))

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
        self.assertNear(w.Q, 0.5)

        with self.assertRaisesRegex(ct.CanteraError, "setState"):
            ct.PureFluid(yaml=yaml.format(T=373, Q=0.5))

        w = ct.PureFluid(yaml=yaml.format(T=370, Q=0.0))
        self.assertNear(w.P, 101325)

        with self.assertRaisesRegex(ct.CanteraError, "setState"):
            ct.PureFluid(yaml=yaml.format(T=370, Q=1.0))

    def test_sort(self):
        np.random.seed(0)
        t = np.random.random(101)
        T = np.linspace(300., 1000., 101)
        P = ct.one_atm * (1. + 10.*np.random.random(101))

        states = ct.SolutionArray(self.gas, 101, extra={'t': t})
        states.TP = T, P

        states.sort('t')
        self.assertTrue((states.t[1:] - states.t[:-1] > 0).all())
        self.assertFalse((states.T[1:] - states.T[:-1] > 0).all())
        self.assertFalse(np.allclose(states.P, P))

        states.sort('T')
        self.assertFalse((states.t[1:] - states.t[:-1] > 0).all())
        self.assertTrue((states.T[1:] - states.T[:-1] > 0).all())
        self.assertArrayNear(states.P, P)

        states.sort('T', reverse=True)
        self.assertTrue((states.T[1:] - states.T[:-1] < 0).all())

    def test_set_equivalence_ratio(self):
        states = ct.SolutionArray(self.gas, 8)
        phi = np.linspace(0.5, 2, 8)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        # The mole fraction arrays need to be squeezed here to reduce their
        # dimensionality from a 2-d column array to a vector for comparison
        # with phi.
        states.set_equivalence_ratio(phi, fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(1)
        self.assertArrayNear(comp, phi)
        states.set_equivalence_ratio(phi[0], fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(1)
        self.assertArrayNear(comp, np.full_like(phi, phi[0]))
        states.set_equivalence_ratio(phi.tolist(), fuel, oxidizer)
        comp = (states("H2").X / (2 * states("O2").X)).squeeze(1)
        self.assertArrayNear(comp, phi)

    def test_set_equivalence_ratio_wrong_shape_raises(self):
        states = ct.SolutionArray(self.gas, 8)
        phi = np.linspace(0.5, 2, 7)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        with self.assertRaisesRegex(ValueError, r"shape mismatch"):
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
        self.assertArrayNear(comp, phi)

    def test_set_mixture_fraction(self):
        states = ct.SolutionArray(self.gas, 8)
        mixture_fraction = np.linspace(0.5, 1, 8)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        # The mass fraction arrays need to be squeezed here to reduce their
        # dimensionality from a 2-d column array to a vector for comparison
        # with mixture_fraction.
        states.set_mixture_fraction(mixture_fraction, fuel, oxidizer)
        self.assertArrayNear(states("H2").Y.squeeze(1), mixture_fraction)
        states.set_mixture_fraction(mixture_fraction[0], fuel, oxidizer)
        self.assertArrayNear(
            states("H2").Y.squeeze(1),
            np.full_like(mixture_fraction, mixture_fraction[0]),
        )
        states.set_mixture_fraction(mixture_fraction.tolist(), fuel, oxidizer)
        self.assertArrayNear(states("H2").Y.squeeze(1), mixture_fraction)

    def test_set_mixture_fraction_wrong_shape_raises(self):
        states = ct.SolutionArray(self.gas, 8)
        mixture_fraction = np.linspace(0.5, 1, 7)
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        with self.assertRaisesRegex(ValueError, r"shape mismatch"):
            states.set_mixture_fraction(mixture_fraction, fuel, oxidizer)

    def test_set_mixture_fraction_2D(self):
        states = ct.SolutionArray(self.gas, (2, 4))
        mixture_fraction = np.linspace(0.5, 1, 8).reshape((2, 4))
        fuel, oxidizer = "H2:1.0", "O2:1.0"
        states.set_mixture_fraction(mixture_fraction, fuel, oxidizer)
        # The mass fraction array needs to be squeezed here to reduce its
        # dimensionality from a 3-d array to a 2-d array for comparison
        # with mixture_fraction.
        self.assertArrayNear(states("H2").Y.squeeze(2), mixture_fraction)

    def test_species_slicing(self):
        states = ct.SolutionArray(self.gas, (2,5))
        states.TPX = np.linspace(500, 1000, 5), 2e5, 'H2:0.5, O2:0.4'
        states.equilibrate('HP')
        self.assertArrayNear(states('H2').X.squeeze(),
                             states.X[...,self.gas.species_index('H2')])

        kk = (self.gas.species_index('OH'), self.gas.species_index('O'))
        self.assertArrayNear(states('OH','O').partial_molar_cp,
                             states.partial_molar_cp[...,kk])

    def test_slice_SolutionArray(self):
        soln = ct.SolutionArray(self.gas, 10)
        arr = soln[2:9:3]
        self.assertEqual(len(arr.T), 3)

    def test_zero_length_slice_SolutionArray(self):
        states = ct.SolutionArray(self.gas, 4)
        arr1 = states[3:3]
        self.assertEqual(len(arr1.T), 0)
        self.assertEqual(arr1.X.shape, (0,10))
        self.assertEqual(arr1.n_reactions, 29)

        states.TP = [100,300,900,323.23], ct.one_atm
        arr2 = states[slice(0)]
        self.assertEqual(len(arr2.T), 0)
