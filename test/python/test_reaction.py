from math import exp
from pathlib import Path
import sys
import textwrap
import gc

import cantera as ct
import numpy as np
from . import utilities
from .utilities import has_temperature_derivative_warnings
import pytest


class TestImplicitThirdBody(utilities.CanteraTest):
    # tests for three-body reactions with specified collision partner

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution("gri30.yaml")

    def test_implicit_three_body(self):
        # check equivalency of auto-detected and explicit specification
        yaml1 = """
            equation: H + 2 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            """
        rxn1 = ct.Reaction.from_yaml(yaml1, self.gas)
        self.assertEqual(rxn1.reaction_type, "three-body-Arrhenius")
        self.assertEqual(rxn1.third_body.default_efficiency, 0.)
        self.assertEqual(rxn1.third_body.efficiencies, {"O2": 1})

        yaml2 = """
            equation: H + O2 + M <=> HO2 + M
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            default-efficiency: 0
            efficiencies: {O2: 1.0}
            """
        rxn2 = ct.Reaction.from_yaml(yaml2, self.gas)
        self.assertEqual(rxn1.third_body.efficiencies, rxn2.third_body.efficiencies)
        self.assertEqual(rxn1.third_body.default_efficiency, rxn2.third_body.default_efficiency)

    def test_duplicate(self):
        # @todo simplify this test
        #     duplicates are currently only checked for import from file
        gas1 = ct.Solution(thermo="ideal-gas", kinetics="gas",
                           species=self.gas.species(), reactions=[])

        yaml1 = """
            equation: H + O2 + H2O <=> HO2 + H2O
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            """
        rxn1 = ct.Reaction.from_yaml(yaml1, gas1)

        yaml2 = """
            equation: H + O2 + M <=> HO2 + M
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            default-efficiency: 0
            efficiencies: {H2O: 1}
            """
        rxn2 = ct.Reaction.from_yaml(yaml2, gas1)

        self.assertEqual(rxn1.reaction_type, rxn2.reaction_type)
        self.assertEqual(rxn1.reactants, rxn2.reactants)
        self.assertEqual(rxn1.products, rxn2.products)
        self.assertEqual(rxn1.third_body.efficiencies, rxn2.third_body.efficiencies)
        self.assertEqual(rxn1.third_body.default_efficiency, rxn2.third_body.default_efficiency)

        gas1.add_reaction(rxn1)
        gas1.add_reaction(rxn2)

        fname = "duplicate.yaml"
        gas1.write_yaml(fname)

        with self.assertRaisesRegex(Exception, "Undeclared duplicate reactions"):
            ct.Solution(fname)

        Path(fname).unlink()

    def test_short_serialization(self):
        # check that serialized output is compact
        yaml = """
            equation: H + O2 + H2O <=> HO2 + H2O
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        input_data = rxn.input_data

        self.assertNotIn("type", input_data)
        self.assertNotIn("default-efficiency", input_data)
        self.assertNotIn("efficiencies", input_data)

    def test_non_integer_stoich(self):
        # check that non-integer coefficients prevent automatic conversion
        yaml = """
            equation: 2 H + 1.5 O2 <=> H2O + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "Arrhenius")

    def test_not_three_body(self):
        # check that insufficient reactants prevent automatic conversion
        yaml = """
            equation: HCNO + H <=> H + HNCO  # Reaction 270
            rate-constant: {A: 2.1e+15, b: -0.69, Ea: 2850.0}
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "Arrhenius")

    def test_user_override(self):
        # check that type specification prevents automatic conversion
        yaml = """
            equation: H + 2 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            type: elementary
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "Arrhenius")
        assert "type" in rxn.input_data
        assert rxn.input_data["type"] == "elementary"


class ReactionRateTests:
    # test suite for reaction rate expressions

    _cls = None # reaction rate object to be tested
    _type = None # name of reaction rate
    _index = None # index of reaction in "kineticsfromscratch.yaml"
    _parts = {}
    _input = None # input parameters (dict corresponding to YAML)
    _yaml = None # yaml string specifying parameters

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.soln = ct.Solution("kineticsfromscratch.yaml")

    def setUp(self):
        self.soln.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5, H2O2:1e-7"
        self.soln.TP = 900, 2 * ct.one_atm

    def finalize(self, rate):
        # perform additional setup after construction (whenever applicable)
        return rate

    def from_parts(self):
        # create reaction rate object from parts
        return self.finalize(self._cls(**self._parts))

    def from_input(self, input=None):
        # create reaction rate object from input_data
        if input is None:
            input = self._input
        else:
            self.assertIsInstance(input, dict)
        return self.finalize(self._cls(input_data=input))

    def from_yaml(self):
        # create reaction rate object from yaml
        return self.finalize(ct.ReactionRate.from_yaml(self._yaml))

    def from_dict(self, input=None):
        # create reaction rate object from dictionary
        if input is None:
            input = self.from_yaml().input_data
        else:
            self.assertIsInstance(input, dict)
        return self.finalize(ct.ReactionRate.from_dict(input))

    def eval(self, rate):
        # evaluate rate expression
        return rate(self.soln.T)

    def check_rate(self, rate):
        # check rates
        self.assertEqual(self._type, rate.type)
        self.assertIn(self._cls.__name__, f"{rate}")
        value = self.eval(rate)
        self.assertIsFinite(value)
        self.assertNear(value, self.soln.forward_rate_constants[self._index])

    def test_from_parts(self):
        # check constructors (from argument list)
        self.check_rate(self.from_parts())

    def test_from_yaml(self):
        # check constructors (from yaml input)
        self.check_rate(self.from_yaml())

    def test_from_dict(self):
        # check constructors (from dictionary)
        self.check_rate(self.from_dict())

    def test_from_input(self):
        # check constructors (from 'input_data' argument)
        self.check_rate(self.from_input())

    def test_unconfigured(self):
        # check behavior of unconfigured rate object
        rate0 = self.from_input({})
        self.assertIsNaN(self.eval(rate0))
        input_data = rate0.input_data
        rate1 = self.from_dict(input_data)
        self.assertEqual(rate1.type, self._type)
        self.assertIsNaN(self.eval(rate1))

    def test_roundtrip(self):
        # check round-trip instantiation via input_data
        rate0 = self.from_yaml()
        input_data = rate0.input_data
        rate1 = self.from_dict(input_data)
        self.check_rate(rate1)

    def test_with_units(self):
        # test custom units. Sticking coefficients are dimensionless, so this is only
        # a concern for other rate types
        units = "units: {length: cm, quantity: mol}"
        yaml = f"{textwrap.dedent(self._yaml)}\n{units}"
        if "sticking" not in yaml:
            with self.assertRaisesRegex(Exception, "undefined units"):
                ct.ReactionRate.from_yaml(yaml)

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_derivative_ddT(self):
        # check temperature derivative against numerical derivative
        deltaT = self.soln.derivative_settings["rtol-delta"]
        deltaT *= self.soln.T
        rate = self.from_yaml()
        k0 = self.eval(rate)

        # derivative at constant pressure
        dcdt = - self.soln.density_mole / self.soln.T
        drate = self.soln.forward_rate_constants_ddT
        drate += self.soln.forward_rate_constants_ddC * dcdt
        self.soln.TP = self.soln.T + deltaT, self.soln.P
        k1 = self.eval(rate)
        self.assertNear((k1 - k0) / deltaT, drate[self._index], 1e-6)

    def test_derivative_ddP(self):
        # check pressure derivative against numerical derivative
        deltaP = self.soln.derivative_settings["rtol-delta"]
        deltaP *= self.soln.P
        rate = self.from_yaml()
        k0 = self.eval(rate)

        drate = self.soln.forward_rate_constants_ddP
        self.soln.TP = self.soln.T, self.soln.P + deltaP
        k1 = self.eval(rate)
        self.assertNear((k1 - k0) / deltaP, drate[self._index], 1e-6)


class TestArrheniusRate(ReactionRateTests, utilities.CanteraTest):
    # test Arrhenius rate expressions

    _cls = ct.ArrheniusRate
    _type = "Arrhenius"
    _index = 0
    _input = {"rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}}
    _parts = {"A": 38.7, "b": 2.7, "Ea": 26191840.0}
    _yaml = "rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}"

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        self.assertNear(self._parts["Ea"], rate.activation_energy)
        self.check_rate(rate)

    def test_negative_A(self):
        # test reaction rate property
        rate = self.from_parts()
        self.assertFalse(rate.allow_negative_pre_exponential_factor)
        rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(rate.allow_negative_pre_exponential_factor)

    def test_standalone(self):
        # test creation with unsupported alternative units
        yaml = "rate-constant: {A: 4.0e+21 cm^6/mol^2/s, b: 0.0, Ea: 1207.72688}"
        with self.assertRaisesRegex(Exception, "undefined units"):
            ct.ReactionRate.from_yaml(yaml)

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_derivative_ddT_exact(self):
        # check exact derivative against analytical and numerical derivatives
        rate = self.from_parts()
        T = 1000.
        self.soln.TP = T, self.soln.P

        R = ct.gas_constant
        Ea = rate.activation_energy
        b =  rate.temperature_exponent
        A = rate.pre_exponential_factor
        k0 = self.eval(rate)
        self.assertNear(k0, A * T**b * np.exp(-Ea/R/T))

        scaled_ddT = (Ea / R / T + b) / T

        dkdT = self.soln.forward_rate_constants_ddT[self._index]
        self.assertNear(dkdT, k0 * scaled_ddT) # exact

        dT = 1.e-6
        dkdT_numeric = (rate(T + dT) - rate(T)) / dT
        self.assertNear(dkdT, dkdT_numeric, 1.e-6)


class TestBlowersMaselRate(ReactionRateTests, utilities.CanteraTest):
    # test Blowers-Masel rate expressions

    _cls = ct.BlowersMaselRate
    _type = "Blowers-Masel"
    _index = 6
    _input = {"rate-constant": {"A": 38700, "b": 2.7, "Ea0": 1.0958665856e8, "w": 1.7505856e13}}
    _parts = {"A": 38700, "b": 2.7, "Ea0": 1.0958665856e8, "w": 1.7505856e13}
    _yaml = """
        type: Blowers-Masel
        rate-constant: {A: 38700, b: 2.7, Ea0: 1.0958665856e8, w: 1.7505856e13}
        """

    def eval(self, rate):
        rate.delta_enthalpy = self.soln.delta_enthalpy[self._index]
        with pytest.warns(UserWarning, match="BlowersMaselData::update"):
            return rate(self.soln.T)

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        self.assertNear(self._parts["Ea0"], rate.activation_energy)
        self.assertNear(self._parts["w"], rate.bond_energy)
        self.check_rate(rate)

    def test_negative_A(self):
        # test reaction rate property
        rate = self.from_parts()
        self.assertFalse(rate.allow_negative_pre_exponential_factor)
        rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(rate.allow_negative_pre_exponential_factor)

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    def test_derivative_ddT(self):
        super().test_derivative_ddT()


class TestTwoTempPlasmaRate(ReactionRateTests, utilities.CanteraTest):
    # test TwoTempPlasma rate expressions

    _cls = ct.TwoTempPlasmaRate
    _type = "two-temperature-plasma"
    _index = 11
    _input = {"rate-constant": {"A": 17283, "b": -3.1, "Ea-gas": -5820000, "Ea-electron": 1081000}}
    _parts = {"A": 17283, "b": -3.1, "Ea_gas": -5820000, "Ea_electron": 1081000}
    _yaml = """
        type: two-temperature-plasma
        rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
        """

    def eval(self, rate):
        # check evaluation as a function of temperature and electron temperature
        return rate(self.soln.T, self.soln.Te)

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        self.assertAlmostEqual(self._parts["Ea_gas"], rate.activation_energy)
        self.assertAlmostEqual(self._parts["Ea_electron"], rate.activation_electron_energy)
        self.check_rate(rate)

    def test_negative_A(self):
        # test reaction rate property
        rate = self.from_parts()
        self.assertFalse(rate.allow_negative_pre_exponential_factor)
        rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(rate.allow_negative_pre_exponential_factor)

    def test_derivative_ddT(self):
        # check temperature derivative against numerical derivative
        deltaT = self.soln.derivative_settings["rtol-delta"]
        deltaT *= self.soln.T
        rate = self.from_yaml()
        k0 = self.eval(rate)

        # derivative at constant pressure and constant electron temperature
        dcdt = - self.soln.density_mole / self.soln.T
        drate = self.soln.forward_rate_constants_ddT
        drate += self.soln.forward_rate_constants_ddC * dcdt
        self.soln.TP = self.soln.T + deltaT, self.soln.P
        # Due to Te changes automatically with T, the initial value is used instead.
        k1 = rate(self.soln.T, self.soln.Te - deltaT)
        self.assertNear((k1 - k0) / deltaT, drate[self._index], 1e-6)


class TestTwoTempPlasmaRateShort(TestTwoTempPlasmaRate, utilities.CanteraTest):
    # test TwoTempPlasma rate expressions

    _index = 12
    _input = {"rate-constant": {"A": 17283, "b": -3.1}}
    _parts = {"A": 17283, "b": -3.1}
    _yaml = """
        type: two-temperature-plasma
        rate-constant: {A: 17283, b: -3.1, Ea-gas: 0.0 J/mol, Ea-electron: 0.0 J/mol}
        """

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        self.assertAlmostEqual(rate.activation_energy, 0.)
        self.assertAlmostEqual(rate.activation_electron_energy, 0.)
        self.check_rate(rate)


class FalloffRateTests(ReactionRateTests):
    # test Falloff rate expressions
    _type = "falloff"
    _n_data = [0] # list of valid falloff coefficient array lengths

    @classmethod
    def setUpClass(cls):
        ReactionRateTests.setUpClass()
        param = cls._input["low-P-rate-constant"]
        cls._parts["low"] = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        param = cls._input["high-P-rate-constant"]
        cls._parts["high"] = ct.Arrhenius(param["A"], param["b"], param["Ea"])

    def eval(self, rate):
        concm = self.soln.third_body_concentrations[self._index]
        return rate(self.soln.T, concm)

    def test_data(self):
        rate = self.from_parts()
        for n in self._n_data:
            rate.falloff_coeffs = np.random.rand(n)

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_derivative_ddT(self):
        pert = self.soln.derivative_settings["rtol-delta"]
        deltaT = self.soln.T * pert
        TP = self.soln.TP
        rate = self.from_yaml()
        k0 = self.eval(rate)

        # derivative at constant volume
        drate = self.soln.forward_rate_constants_ddT
        self.soln.TP = self.soln.T * (1 + pert), self.soln.P * (1 + pert)
        k1 = self.eval(rate)
        self.assertNear((k1 - k0) / deltaT, drate[self._index], 1e-6)

        # derivative at constant pressure
        self.soln.TP = TP
        dcdt = - self.soln.density_mole / self.soln.T
        drate += self.soln.forward_rate_constants_ddC * dcdt
        self.soln.TP = self.soln.T * (1 + pert), self.soln.P
        k1 = self.eval(rate)
        self.assertNear((k1 - k0) / deltaT, drate[self._index], 1e-6)

    def test_derivative_ddP(self):
        pert = self.soln.derivative_settings["rtol-delta"]
        deltaP = self.soln.P * pert
        rate = self.from_yaml()
        k0 = self.eval(rate)

        # derivative at constant temperature
        drate = self.soln.forward_rate_constants_ddP
        dcdp = self.soln.density_mole / self.soln.P
        drate += self.soln.forward_rate_constants_ddC * dcdp
        self.soln.TP = self.soln.T, self.soln.P + deltaP
        k1 = self.eval(rate)
        self.assertNear((k1 - k0) / deltaP, drate[self._index], 1e-6)


class TestLindemannRate(FalloffRateTests, utilities.CanteraTest):
    # test Lindemann rate expressions

    _cls = ct.LindemannRate
    _index = 7
    _parts = {
        "falloff_coeffs": [],
        }
    _input = {
        "type": "falloff",
        "low-P-rate-constant": {"A": 2.3e+12, "b": -0.9, "Ea": -7112800.0},
        "high-P-rate-constant": {"A": 7.4e+10, "b": -0.37, "Ea": 0.0},
        }
    _yaml = """
        type: falloff
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -7112800.0}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0}
        """

    def test_falloff_function(self):
        rate = self.from_parts()
        # Falloff-function for Lindemann is unity by definition
        assert np.isclose(rate.falloff_function(300, 0.), 1.)


class TestTroeRate(FalloffRateTests, utilities.CanteraTest):
    # test Troe rate expressions

    _cls = ct.TroeRate
    _index = 2
    _parts = {"falloff_coeffs": [0.7346, 94.0, 1756.0, 5182.0]}
    _input = {
        "type": "falloff",
        "low-P-rate-constant": {"A": 2.3e+12, "b": -0.9, "Ea": -7112800.0},
        "high-P-rate-constant": {"A": 7.4e+10, "b": -0.37, "Ea": 0.0},
        "Troe": {"A": 0.7346, "T3": 94.0, "T1": 1756.0, "T2": 5182.0},
        }
    _yaml = """
        type: falloff
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -7112800.0}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0}
        Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
        """
    _n_data = [3, 4]

    def test_unexpected_parameter(self):
        yaml = """
            type: falloff
            low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -7112800.0}
            high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0}
            Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 0.}
            """

        with pytest.warns(UserWarning, match="Unexpected parameter value T2=0"):
            ct.ReactionRate.from_yaml(yaml)


class TestSriRate(FalloffRateTests, utilities.CanteraTest):
    # test SRI rate expressions

    _cls = ct.SriRate
    _index = 8
    _parts = {"falloff_coeffs": [1.1, 700.0, 1234.0, 56.0, 0.7]}
    _input = {
        "type": "falloff",
        "high-P-rate-constant": {"A": 4.0e+15, "b": -0.5, "Ea": 418400.0},
        "low-P-rate-constant": {"A": 7.0e+20, "b": -1.0, "Ea": 0.0},
        "SRI": {"A": 1.1, "B": 700.0, "C": 1234.0, "D": 56.0, "E": 0.7},
        }
    _yaml = """
        type: falloff
        high-P-rate-constant: {A: 4.0e+15, b: -0.5, Ea: 100.0 cal/mol}
        low-P-rate-constant: {A: 7.0e+20, b: -1.0, Ea: 0.0 cal/mol}
        SRI: {A: 1.1, B: 700.0, C: 1234.0, D: 56.0, E: 0.7}
        """
    _n_data = [3, 5]


class TestTsangRate(FalloffRateTests, utilities.CanteraTest):
    # test Tsang rate expressions

    _cls = ct.TsangRate
    _index = 9
    _parts = {"falloff_coeffs": [0.95, -1.0e-04]}
    _input = {
        "type": "falloff",
        "high-P-rate-constant": {"A": 4.0e+15, "b": -0.5, "Ea": 418400.0},
        "low-P-rate-constant": {"A": 7.0e+20, "b": -1.0, "Ea": 0.0},
        "Tsang": {"A": 0.95, "B": -1.0e-04}
        }
    _yaml = """
        type: falloff
        high-P-rate-constant: {A: 4.0e+15, b: -0.5, Ea: 100.0 cal/mol}
        low-P-rate-constant: {A: 7.0e+20, b: -1.0, Ea: 0.0 cal/mol}
        Tsang: {A: 0.95, B: -1.0e-04}
        """
    _n_data = [1, 2]


class TestPlogRate(ReactionRateTests, utilities.CanteraTest):
    # test Plog rate expressions

    _cls = ct.PlogRate
    _type = "pressure-dependent-Arrhenius"
    _index = 3
    _input = {"rate-constants": [
        {"P": 1013.25, "A": 1.2124e+16, "b": -0.5779, "Ea": 45491376.8},
        {"P": 101325., "A": 4.9108e+31, "b": -4.8507, "Ea": 103649395.2},
        {"P": 1013250., "A": 1.2866e+47, "b": -9.0246, "Ea": 166508556.0},
        {"P": 10132500., "A": 5.9632e+56, "b": -11.529, "Ea": 220076726.4}]}
    _yaml = """
        type: pressure-dependent-Arrhenius
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """

    @classmethod
    def setUpClass(cls):
        ReactionRateTests.setUpClass()
        cls._parts = {
            "rates": [(rc["P"], ct.Arrhenius(rc["A"], rc["b"], rc["Ea"]))
                      for rc in cls._input["rate-constants"]],
            }

    def eval(self, rate):
        # check evaluation as a function of temperature and pressure
        return rate(self.soln.T, self.soln.P)

    def test_get_rates(self):
        # test getter for property rates
        rate = self.from_parts()
        rates = rate.rates
        self.assertIsInstance(rates, list)

        other = self._input["rate-constants"]
        self.assertEqual(len(rates), len(other))
        for index, item in enumerate(rates):
            P, rate = item
            self.assertNear(P, other[index]["P"])
            self.assertNear(rate.pre_exponential_factor, other[index]["A"])
            self.assertNear(rate.temperature_exponent, other[index]["b"])
            self.assertNear(rate.activation_energy, other[index]["Ea"])

    def test_set_rates(self):
        # test setter for property rates
        other = [
            {"P": 100., "A": 1.2124e+16, "b": -1., "Ea": 45491376.8},
            {"P": 10000., "A": 4.9108e+31, "b": -2., "Ea": 103649395.2},
            {"P": 1000000., "A": 1.2866e+47, "b": -3., "Ea": 166508556.0}]
        rate = ct.PlogRate([(o["P"], ct.Arrhenius(o["A"], o["b"], o["Ea"]))
                            for o in other])
        rates = rate.rates
        self.assertEqual(len(rates), len(other))

        for index, item in enumerate(rates):
            P, rate = item
            self.assertNear(P, other[index]["P"])
            self.assertNear(rate.pre_exponential_factor, other[index]["A"])
            self.assertNear(rate.temperature_exponent, other[index]["b"])
            self.assertNear(rate.activation_energy, other[index]["Ea"])

    def test_no_rates(self):
        # test instantiation of empty rate
        rate = ct.PlogRate()
        self.assertIsInstance(rate.rates, list)

    def test_standalone(self):
        yaml = """
            type: pressure-dependent-Arrhenius
            rate-constants:
            - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
            - {P: 1.0 atm, A: 4.9108e+31 cm^6/mol^2/s, b: -4.8507, Ea: 2.47728e+04 cal/mol}
            - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
            - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
            """
        with self.assertRaisesRegex(Exception, "undefined units"):
            ct.ReactionRate.from_yaml(yaml)


class TestChebyshevRate(ReactionRateTests, utilities.CanteraTest):
    # test Chebyshev rate expressions

    _cls = ct.ChebyshevRate
    _type = "Chebyshev"
    _index = 4
    _input = {"data": [[8.2883, -1.1397, -0.12059, 0.016034],
                       [1.9764, 1.0037, 0.0072865, -0.030432],
                       [0.3177, 0.26889, 0.094806, -0.0076385]],
              "pressure-range": [1000.0, 10000000.0],
              "temperature-range": [290.0, 3000.0]}
    _yaml = """
        type: Chebyshev
        temperature-range: [290.0, 3000.0]
        pressure-range: [9.869232667160128e-03 atm, 98.69232667160128 atm]
        data:
        - [8.2883, -1.1397, -0.12059, 0.016034]
        - [1.9764, 1.0037, 7.2865e-03, -0.030432]
        - [0.3177, 0.26889, 0.094806, -7.6385e-03]
        """

    @classmethod
    def setUpClass(cls):
        ReactionRateTests.setUpClass()
        cls._parts = {
            "pressure_range": cls._input["pressure-range"],
            "temperature_range": cls._input["temperature-range"],
            "data": cls._input["data"],
        }

    def eval(self, rate):
        # check evaluation as a function of temperature and pressure
        return rate(self.soln.T, self.soln.P)

    def test_from_parts(self):
        rate = self.from_parts()
        temperature_range = self._parts["temperature_range"]
        self.assertEqual(temperature_range[0], rate.temperature_range[0])
        self.assertEqual(temperature_range[1], rate.temperature_range[1])
        pressure_range = self._parts["pressure_range"]
        self.assertEqual(pressure_range[0], rate.pressure_range[0])
        self.assertEqual(pressure_range[1], rate.pressure_range[1])
        self.assertTrue(np.all(self._parts["data"] == rate.data))
        self.assertEqual(rate.n_pressure, rate.data.shape[1])
        self.assertEqual(rate.n_temperature, rate.data.shape[0])


class TestPhotolysisRateParameterization(ReactionRateTests, utilities.CanteraTest):
    # test Photolysis rate expressions
    _cls = ct.PhotolysisRate
    _type = "photolysis"
    _index = 14
    _input = {"rate-constant": {"l": 4.775e-4, "m": 0.298, "n": 0.080}}
    _parts = {"l": 4.775e-4, "m": 0.298, "n": 0.080}
    _yaml = """
            type: photolysis
            rate-constant: {l: 4.775e-4, m: 0.298, n: 0.080}
            """
    zenith_angle = 0.785

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.soln = ct.Solution("kineticsfromscratch.yaml")
        cls.soln.zenith_angle = cls.zenith_angle

    def eval(self, rate):
        # evaluate rate expression
        return rate(self.zenith_angle)

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["l"], rate.l_param)
        self.assertEqual(self._parts["m"], rate.m_param)
        self.assertEqual(self._parts["n"], rate.n_param)
        self.check_rate(rate)

    def test_derivative_ddT(self):
        """ Photolysis rate is not temperature dependent. """
        drate = self.soln.forward_rate_constants_ddT
        self.assertEqual(0.0, drate[self._index])

    def test_derivative_ddP(self):
        """ Photolysis rate is not pressure dependent. """
        drate = self.soln.forward_rate_constants_ddP
        self.assertEqual(0.0, drate[self._index])


class TestPhotolysisRateConstant(ReactionRateTests, utilities.CanteraTest):
    # test Photolysis rate expressions
    _cls = ct.PhotolysisRate
    _type = "photolysis"
    _index = 15
    _input = {"rate-constant": {"J": 0.5}}
    _parts = {"J": 0.5}
    _yaml = """
            type: photolysis
            rate-constant: {J: 0.5}
            """

    def eval(self, rate):
        return rate(0)

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertNear(self._parts["J"], self.eval(rate))
        self.check_rate(rate)

    def test_derivative_ddT(self):
        """ Photolysis rate is not temperature dependent. """
        drate = self.soln.forward_rate_constants_ddT
        self.assertEqual(0.0, drate[self._index])

    def test_derivative_ddP(self):
        """ Photolysis rate is not pressure dependent. """
        drate = self.soln.forward_rate_constants_ddP
        self.assertEqual(0.0, drate[self._index])

class SurfaceReactionRateTests(ReactionRateTests):
    # test suite for surface reaction rate expressions

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.soln = ct.Interface("kineticsfromscratch.yaml", "Pt_surf", transport_model=None)
        cls.gas = cls.soln.adjacent["ohmech"]

    def setUp(self):
        self.soln.TP = 900, ct.one_atm
        self.gas.X = "H2:0.05, H2O:0.01, O:1e-4, OH: 1e5, H:2e-5, O2:0.21, AR:0.79"
        self.gas.TP = 900, ct.one_atm

    def eval(self, rate, species=True):
        # evaluate rate expression
        if species:
            rate.set_species(self.soln.species_names)
        rate.site_density = self.soln.site_density
        self.assertEqual(rate.site_density, self.soln.site_density)
        if "Blowers-Masel" in self._type:
            rate.delta_enthalpy = self.soln.delta_enthalpy[self._index]
        with pytest.warns(UserWarning, match="InterfaceData::update"):
            return rate(self.soln.T, self.soln.coverages)

    def from_parts(self):
        rate = super().from_parts()
        if "coverage-dependencies" in self._input:
            rate.coverage_dependencies = self._input["coverage-dependencies"]
        return rate

    def test_no_species(self):
        # evaluate rate expression without providing species
        rate = self.from_yaml()
        value = self.eval(rate, species=False)
        if "coverage-dependencies" in self._input:
            self.assertIsNaN(value)
        else:
            self.assertIsFinite(value)

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        if "Ea" in self._parts:
            self.assertNear(self._parts["Ea"], rate.activation_energy)
        else:
            self.assertNear(self._parts["w"], rate.bond_energy)
        self.check_rate(rate)

    @pytest.mark.skip("Derivative is not supported")
    def test_derivative_ddT(self):
        pass

    @pytest.mark.skip("Derivative is not supported")
    def test_derivative_ddP(self):
        pass


class StickingReactionRateTests(SurfaceReactionRateTests):
    # test suite for surface reaction rate expressions

    _sticking_species = None
    _sticking_order = None

    def finalize(self, rate):
        weight = self.gas.molecular_weights[self.gas.species_index(self._sticking_species)]
        rate.sticking_species = self._sticking_species
        rate.sticking_order = self._sticking_order
        rate.sticking_weight = weight
        return rate

    def from_parts(self):
        rate = super().from_parts()
        rate.motz_wise_correction = "Motz-Wise" in self._input
        return self.finalize(rate)

    def test_sticking_coeffs(self):
        rate = self.from_yaml()
        species = self._sticking_species
        weight = self.gas.molecular_weights[self.gas.species_index(self._sticking_species)]
        assert rate.sticking_species == self._sticking_species
        assert rate.sticking_order == self._sticking_order
        assert rate.sticking_weight == pytest.approx(weight)


class TestSurfaceArrheniusRate(SurfaceReactionRateTests, utilities.CanteraTest):
    # test interface-Arrhenius rate expressions without coverage dependency

    _cls = ct.InterfaceArrheniusRate
    _type = "interface-Arrhenius"
    _index = 0
    _input = {"rate-constant": {"A": 3.7e+20, "b": 0, "Ea": 1.15e7}}
    _parts = _input["rate-constant"]
    _yaml = """
        rate-constant: {A: 3.7e+20, b: 0, Ea: 11500 J/mol}
        type: interface-Arrhenius
        """


class TestInterfaceArrheniusRate(SurfaceReactionRateTests, utilities.CanteraTest):
    # test interface-Arrhenius rate expressions with coverage dependency

    _cls = ct.InterfaceArrheniusRate
    _type = "interface-Arrhenius"
    _index = 1
    _input = {
        "rate-constant": {"A": 3.7e+20, "b": 0, "Ea": 213200000.},
        "coverage-dependencies": {"O(S)": {"a": 0.0, "m": 0.0, "E": -60000000.}},
    }
    _parts = _input["rate-constant"]
    _yaml = """
        rate-constant: {A: 3.7e+20, b: 0, Ea: 213200 J/mol}
        coverage-dependencies:
          O(S): {a: 0.0, m: 0.0, E: -6.0e+04 J/mol}
        type: interface-Arrhenius
        """

class TestStickingRate(StickingReactionRateTests, utilities.CanteraTest):
    # test surface-sticking rate expressions without coverage dependency

    _cls = ct.StickingArrheniusRate
    _type = "sticking-Arrhenius"
    _index = 2
    _input = {"sticking-coefficient": {"A": 1., "b": 0, "Ea": 0}}
    _parts = _input["sticking-coefficient"]
    _sticking_species = "H"
    _sticking_order = 1.0
    _yaml = """
        sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
        type: sticking-Arrhenius
        """


class TestCoverageStickingRate(StickingReactionRateTests, utilities.CanteraTest):
    # test sticking rate expressions with coverage dependency

    _cls = ct.StickingArrheniusRate
    _type = "sticking-Arrhenius"
    _index = 3
    _input = {
        "sticking-coefficient": {"A": 0.046, "b": 0, "Ea": 0},
        "coverage-dependencies": {"PT(S)": {"a": 0.0, "m": -1.0, "E": 0.0}},
    }
    _parts = _input["sticking-coefficient"]
    _sticking_species = "H2"
    _sticking_order = 2.0
    _yaml = """
        sticking-coefficient: {A: 0.046, b: 0, Ea: 0}
        coverage-dependencies:
          PT(S): {a: 0.0, m: -1.0, E: 0.0}
        type: sticking-Arrhenius
        """


class TestMotzWiseStickingRate(StickingReactionRateTests, utilities.CanteraTest):
    # test interface reaction with coverages

    _cls = ct.StickingArrheniusRate
    _type = "sticking-Arrhenius"
    _index = 4
    _input = {
        "sticking-coefficient": {"A": 1., "b": 0, "Ea": 0.},
        "Motz-Wise": True,
    }
    _parts = _input["sticking-coefficient"]
    _sticking_species = "OH"
    _sticking_order = 1.0
    _yaml = """
        sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
        Motz-Wise: true
        type: sticking-Arrhenius
        """


class TestSurfaceBMRate(SurfaceReactionRateTests, utilities.CanteraTest):
    # test coverage-Blowers-Masel rate expressions with coverage dependency

    _cls = ct.InterfaceBlowersMaselRate
    _type = "interface-Blowers-Masel"
    _index = 5
    _input = {
        "rate-constant": {"A": 3.7e+20, "b": 0, "Ea0": 67400000.0, 'w': 1000000000.0},
    }
    _parts = _input["rate-constant"]
    _yaml = """
        rate-constant: {A: 3.7e+20, b: 0, Ea0: 67400000.0, w: 1000000000.0}
        type: interface-Blowers-Masel
        """


class TestSurfaceBMRate(SurfaceReactionRateTests, utilities.CanteraTest):
    # test coverage-Blowers-Masel rate expressions with coverage dependency

    _cls = ct.InterfaceBlowersMaselRate
    _type = "interface-Blowers-Masel"
    _index = 6
    _input = {
        "rate-constant": {"A": 3.7e+20, "b": 0, "Ea0": 67400000.0, 'w': 1000000000.0},
        "coverage-dependencies": {"H(S)": {"a": 0.0, "m": 0.0, "E": -6000000.0}},
    }
    _parts = _input["rate-constant"]
    _yaml = """
        rate-constant: {A: 3.7e+20, b: 0, Ea0: 67400000.0, w: 1000000000.0}
        coverage-dependencies:
          H(S): {a: 0.0, m: 0.0, E: -6000000.0}
        type: interface-Blowers-Masel
        """


class TestBMStickate(StickingReactionRateTests, utilities.CanteraTest):
    # test coverage-Blowers-Masel stick expressions with coverage dependency

    _cls = ct.StickingBlowersMaselRate
    _type = "sticking-Blowers-Masel"
    _index = 7
    _input = {
        "sticking-coefficient": {"A": 1., "b": 0., "Ea0": 0., 'w': 100000.},
        "Motz-Wise": True,
    }
    _parts = _input["sticking-coefficient"]
    _sticking_species = "OH"
    _sticking_order = 1.0
    _yaml = """
        sticking-coefficient: {A: 1.0, b: 0, Ea0: 0, w: 100000}
        Motz-Wise: true
        type: sticking-Blowers-Masel
        """


class ReactionTests:
    # test suite for reaction expressions

    _cls = ct.Reaction # reaction object to be tested
    _rate_cls = None # corresponding reaction rate type
    _equation = None # reaction equation string
    _rate = None # parameters for reaction rate object constructor
    _3rd_body = None # object representing third-body collider
    _rate_obj = None # reaction rate object
    _kwargs = {} # additional parameters required by constructor
    _index = None # index of reaction in "kineticsfromscratch.yaml"
    _rate_type = None # name of reaction rate type
    _yaml = None # YAML parameterization
    _rc_units = None # Units of the rate coefficient

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.soln = ct.Solution("kineticsfromscratch.yaml", transport_model=None)
        cls.species = cls.soln.species()

    def setUp(self):
        self.soln.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5"
        self.soln.TP = 900, 2*ct.one_atm
        self.adj = []

    def eval_rate(self, rate):
        # evaluate rate expression
        return rate(self.soln.T)

    def finalize(self, rxn):
        # perform additional setup after construction (whenever applicable)
        return rxn

    def from_yaml(self):
        # create reaction object from yaml
        rxn = self.finalize(ct.Reaction.from_yaml(self._yaml, kinetics=self.soln))
        return self.finalize(rxn)

    def from_dict(self):
        # create reaction rate object from input data
        input_data = self.from_yaml().input_data
        rxn = ct.Reaction.from_dict(input_data, kinetics=self.soln)
        return self.finalize(rxn)

    def from_empty(self):
        # create reaction object with an "empty" rate of the correct type
        rxn = ct.Reaction(equation=self._equation,
                          rate=self._rate_cls(), third_body=self._3rd_body)
        return self.finalize(rxn)

    def from_rate(self, rate):
        # create reaction object from dictionary
        rxn = ct.Reaction(equation=self._equation,
                          rate=rate, third_body=self._3rd_body)
        return self.finalize(rxn)

    def from_parts(self):
        # create reaction rate object from parts
        orig = self.soln.reaction(self._index)
        rxn = ct.Reaction(orig.reactants, orig.products,
                          rate=self._rate_obj, third_body=self._3rd_body)
        rxn.reversible = "<=>" in self._equation
        return self.finalize(rxn)

    def check_rate(self, rate_obj):
        rate = self.eval_rate(rate_obj)
        self.assertNear(rate, self.soln.forward_rate_constants[self._index])

    def check_rxn(self, rxn):
        # helper function that checks reaction configuration
        ix = self._index
        self.assertEqual(rxn.reactants, self.soln.reaction(ix).reactants)
        self.assertEqual(rxn.products, self.soln.reaction(ix).products)
        self.check_rate(rxn.rate)

        if self.soln.thermo_model.lower() == "ideal-surface":
            sol2 = ct.Interface(thermo="ideal-surface", kinetics="surface",
                                species=self.species, reactions=[rxn], adjacent=self.adj)
            sol2.site_density = self.soln.site_density
            sol2.coverages = self.soln.coverages
            sol2.TP = self.soln.TP
        else:
            sol2 = ct.Solution(thermo=self.soln.thermo_model, kinetics=self.soln.kinetics_model,
                               species=self.species, reactions=[rxn])
            sol2.TPX = self.soln.TPX
        self.check_solution(sol2)

    def check_solution(self, sol2):
        # helper function that checks evaluation of reaction rates
        ix = self._index
        self.assertNear(sol2.forward_rate_constants[0],
                        self.soln.forward_rate_constants[ix])
        self.assertNear(sol2.net_rates_of_progress[0],
                        self.soln.net_rates_of_progress[ix])

    def test_rate(self):
        # check consistency of reaction rate and forward rate constant
        self.check_rate(self._rate_obj)

    def test_from_rate(self):
        # check instantiation from keywords / rate defined by dictionary
        self.check_rxn(self.from_rate(self._rate))

    def test_from_rate_obj(self):
        # check instantiation from keywords / rate provided as object
        self.check_rxn(self.from_rate(self._rate_obj))

    def test_from_parts(self):
        # check instantiation from parts (reactants, products, rate expression)
        self.check_rxn(self.from_parts())

    def test_from_yaml(self):
        # check constructors (from yaml input)
        self.check_rxn(self.from_yaml())

    def test_from_dict(self):
        # check instantiation from a yaml dictionary (input_data)
        self.check_rxn(self.from_dict())

    def test_add_rxn(self):
        # check adding new reaction to solution
        if self._yaml is not None:
            rxn = self.from_yaml()
        else:
            rxn = self.from_rate(self._rate_obj)
        if self.soln.thermo_model.lower() == "ideal-surface":
            sol2 = ct.Interface(thermo="ideal-surface", kinetics="surface",
                                species=self.species, reactions=[], adjacent=self.adj)
            sol2.site_density = self.soln.site_density
            sol2.coverages = self.soln.coverages
            sol2.TP = self.soln.TP
        else:
            sol2 = ct.Solution(thermo=self.soln.thermo_model, kinetics=self.soln.kinetics_model,
                               species=self.species, reactions=[])
            sol2.TPX = self.soln.TPX
        sol2.add_reaction(rxn)
        self.check_solution(sol2)

    def test_raises_invalid_rate(self):
        # check exception for instantiation from keywords / invalid rate
        with self.assertRaises(TypeError):
            self.from_rate(tuple())
        with self.assertRaises(TypeError):
            self.from_rate("spam")
        with self.assertRaises(TypeError):
            self.from_rate(False)
        with self.assertRaises(TypeError):
            self.from_rate(1.)

    def test_no_rate(self):
        # check behavior for instantiation from keywords / no rate
        rxn = self.from_empty()
        self.assertIsNaN(self.eval_rate(rxn.rate))

        with pytest.raises(ct.CanteraError, match="validate"):
            ct.Solution(thermo=self.soln.thermo_model,
                        kinetics=self.soln.kinetics_model,
                        species=self.species, reactions=[rxn], adjacent=self.adj)

    def test_replace_rate(self):
        # check replacing reaction rate expression
        if self._yaml is not None:
            rxn = self.from_yaml()
        else:
            rxn = self.from_rate(self._rate_obj)
        rxn = self.from_empty()
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_roundtrip(self):
        # check round-trip instantiation via input_data
        rxn = self.from_rate(self._rate_obj)
        rate_input_data = rxn.rate.input_data
        rate_obj = rxn.rate.__class__(input_data=rate_input_data)
        rxn2 = self.from_rate(rate_obj)
        self.check_rxn(rxn2)

    def test_rate_coeff_units(self):
        rxn = self.from_yaml()
        assert str(rxn.rate_coeff_units) == str(self._rc_units)

    def test_rate_conversion_units(self):
        rxn = self.from_yaml()
        assert str(rxn.rate.conversion_units) == str(self._rc_units)

    def check_equal(self, one, two):
        # helper function for deprecation tests
        self.assertEqual(type(one), type(two))
        if isinstance(one, (list, tuple, np.ndarray)):
            self.assertArrayNear(one, two)
        elif isinstance(one, (dict, str)):
            assert one == two
        else:
            self.assertNear(one, two)


class TestElementary(ReactionTests, utilities.CanteraTest):
    # test elementary reaction

    _rate_cls = ct.ArrheniusRate
    _equation = "H2 + O <=> H + OH"
    _rate = {"A": 38.7, "b": 2.7, "Ea": 2.619184e+07}
    _index = 0
    _rate_type = "Arrhenius"
    _yaml = """
        equation: O + H2 <=> H + OH
        type: elementary
        rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        cls._rate_obj = ct.ArrheniusRate(**cls._rate)


class TestThreeBody(TestElementary):
    # test three-body reaction

    _equation = "2 O + M <=> O2 + M"
    _rate = {"A": 1.2e11, "b": -1.0, "Ea": 0.0}
    _3rd_body = ct.ThirdBody(efficiencies={"H2": 2.4, "H2O": 15.4, "AR": 0.83})
    _index = 1
    _yaml = """
        equation: 2 O + M <=> O2 + M
        # type: three-body # optional
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """
    _rc_units = ct.Units("m^6 / kmol^2 / s")

    def test_efficiencies(self):
        # check efficiencies
        rxn = self.from_rate(self._rate_obj)
        self.assertEqual(rxn.third_body.efficiencies, self._3rd_body.efficiencies)

    def test_serialization_type(self):
        # test serialization output
        rxn = self.from_yaml()
        assert "type" not in rxn.input_data
        orig = self.soln.reaction(self._index)
        assert "type" in orig.input_data
        assert orig.input_data["type"] == "three-body"


class TestImplicitThreeBody(TestThreeBody):
    # test three-body reactions with explicit collision partner

    _equation = "H + 2 O2 <=> HO2 + O2"
    _rate = {"A": 2.08e+19, "b": -1.24, "Ea": 0.0}
    _3rd_body = ct.ThirdBody("O2")
    _index = 5
    _yaml = """
        equation: H + 2 O2 <=> HO2 + O2
        rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
        """
    _rc_units = ct.Units("m^6 / kmol^2 / s")

    def test_efficiencies(self):
        # overload of default tester
        rxn = self.from_rate(self._rate_obj)
        self.assertEqual(rxn.third_body.efficiencies, {"O2": 1.})
        self.assertEqual(rxn.third_body.default_efficiency, 0.)

    def test_serialization_type(self):
        # test serialization output
        rxn = self.from_yaml()
        assert "type" not in rxn.input_data
        assert "efficiencies" not in rxn.input_data
        orig = self.soln.reaction(self._index)
        assert "type" not in orig.input_data
        assert "efficiencies" not in orig.input_data


class TestTwoTempPlasma(ReactionTests, utilities.CanteraTest):
    # test two-temperature plasma reaction

    _rate_cls = ct.TwoTempPlasmaRate
    _rate_type = "two-temperature-plasma"
    _equation = "O + H => O + H"
    _rate = {"A": 17283, "b": -3.1, "Ea_gas": -5820000, "Ea_electron": 1081000}
    _rate_obj = ct.TwoTempPlasmaRate(A=17283, b=-3.1, Ea_gas=-5820000, Ea_electron=1081000)
    _index = 11
    _yaml = """
        equation: O + H => O + H
        type: two-temperature-plasma
        rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    def eval_rate(self, rate):
        return rate(self.soln.T, self.soln.Te)

    def test_reversible(self):
        orig = self.soln.reaction(self._index)
        rxn = ct.Reaction(orig.reactants, orig.products, rate=self._rate_obj)
        with self.assertRaisesRegex(ct.CanteraError, "does not support reversible"):
            self.check_rxn(rxn)


class TestBlowersMasel(ReactionTests, utilities.CanteraTest):
    # test elementary version of Blowers-Masel reaction

    _rate_cls = ct.BlowersMaselRate
    _rate_type = "Blowers-Masel"
    _equation = "O + H2 <=> H + OH"
    _rate = {"A": 38700, "b": 2.7, "Ea0": 1.0958665856e8, "w": 1.7505856e13}
    _rate_obj = ct.BlowersMaselRate(A=38700, b=2.7, Ea0=1.0958665856e8, w=1.7505856e13)
    _index = 6
    _yaml = """
        equation: O + H2 <=> H + OH
        type: Blowers-Masel
        rate-constant: {A: 38700, b: 2.7, Ea0: 2.619184e4 cal/mol, w: 4.184e9 cal/mol}
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    def eval_rate(self, rate):
        rate.delta_enthalpy = self.soln.delta_enthalpy[self._index]
        with pytest.warns(UserWarning, match="BlowersMaselData::update"):
            return rate(self.soln.T)


class TestThreeBodyBlowersMasel(TestBlowersMasel):
    # test three-body version of Blowers-Masel reaction

    _equation = "O + H2 + M <=> H2O + M"
    _3rd_body = ct.ThirdBody()
    _index = 13
    _yaml = """
        equation: O + H2 + M <=> H2O + M
        type: Blowers-Masel
        rate-constant: {A: 38700, b: 2.7, Ea0: 2.619184e4 cal/mol, w: 4.184e9 cal/mol}
        """
    _rc_units = ct.Units("m^6 / kmol^2 / s")


class TestTroe(ReactionTests, utilities.CanteraTest):
    # test Troe falloff reaction

    _rate_cls = ct.TroeRate
    _equation = "2 OH (+ M) <=> H2O2 (+ M)"
    _rate = {
        "type": "falloff",
        "low_P_rate_constant": {"A": 2.3e+12, "b": -0.9, "Ea": -7112800.0},
        "high_P_rate_constant": {"A": 7.4e+10, "b": -0.37, "Ea": 0.0},
        "Troe": {"A": 0.7346, "T3": 94.0, "T1": 1756.0, "T2": 5182.0}
        }
    _3rd_body = ct.ThirdBody("(+M)", efficiencies={"AR": 0.7, "H2": 2.0, "H2O": 6.0})
    _index = 2
    _rate_type = "falloff"
    _yaml = """
        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
        type: falloff
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        param = cls._rate["low_P_rate_constant"]
        low = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        param = cls._rate["high_P_rate_constant"]
        high = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        param = cls._rate["Troe"]
        data = [param["A"], param["T3"], param["T1"], param["T2"]]
        cls._rate_obj = ct.TroeRate(low=low, high=high, falloff_coeffs=data)

    def eval_rate(self, rate):
        concm = self.soln.third_body_concentrations[self._index]
        return rate(self.soln.T, concm)


class TestLindemann(ReactionTests, utilities.CanteraTest):
    # test Lindemann falloff reaction

    _rate_cls = ct.LindemannRate
    _equation = "2 OH (+ M) <=> H2O2 (+ M)"
    _rate = {
        "type": "falloff",
        "low_P_rate_constant": {"A": 2.3e+12, "b": -0.9, "Ea": -7112800.0},
        "high_P_rate_constant": {"A": 7.4e+10, "b": -0.37, "Ea": 0.0}
        }
    _3rd_body = ct.ThirdBody("(+M)", efficiencies={"AR": 0.7, "H2": 2.0, "H2O": 6.0})
    _index = 7
    _rate_type = "falloff"
    _yaml = """
        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 8
        duplicate: true
        type: falloff
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        param = cls._rate["low_P_rate_constant"]
        low = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        param = cls._rate["high_P_rate_constant"]
        high = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        cls._rate_obj = ct.LindemannRate(low=low, high=high, falloff_coeffs=[])

    def eval_rate(self, rate):
        concm = self.soln.third_body_concentrations[self._index]
        return rate(self.soln.T, concm)


class TestChemicallyActivated(ReactionTests, utilities.CanteraTest):
    # test Chemically Activated falloff reaction

    _rate_cls = ct.LindemannRate
    _equation = "H2O + OH (+M) <=> HO2 + H2 (+M)"
    _rate = {
        "type": "chemically-activated",
        "low_P_rate_constant": {"A": 282.320078, "b": 1.46878, "Ea": -13684043.7508},
        "high_P_rate_constant": {"A": 5.88E-14, "b": 6.721, "Ea": -12644997.768}
        }
    _index = 10
    _rate_type = "chemically-activated"
    _yaml = """
        equation: H2O + OH (+M) <=> HO2 + H2 (+M)  # Reaction 11
        units: {length: cm, quantity: mol, activation-energy: cal/mol}
        type: chemically-activated
        low-P-rate-constant: [282320.078, 1.46878, -3270.56495]
        high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        param = cls._rate["low_P_rate_constant"]
        low = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        param = cls._rate["high_P_rate_constant"]
        high = ct.Arrhenius(param["A"], param["b"], param["Ea"])
        cls._rate_obj = ct.LindemannRate(low=low, high=high, falloff_coeffs=[])
        cls._rate_obj.chemically_activated = True

    def eval_rate(self, rate):
        concm = self.soln.third_body_concentrations[self._index]
        return rate(self.soln.T, concm)


class TestPlog(ReactionTests, utilities.CanteraTest):
    # test Plog reaction

    _rate_cls = ct.PlogRate
    _rate_type = "pressure-dependent-Arrhenius"
    _equation = "H2 + O2 <=> 2 OH"
    _rate = {
        "type": "pressure-dependent-Arrhenius",
        "rate-constants":
             [{"P": 1013.25, "A": 1.2124e+16, "b": -0.5779, "Ea": 45491376.8},
              {"P": 101325., "A": 4.9108e+31, "b": -4.8507, "Ea": 103649395.2},
              {"P": 1013250., "A": 1.2866e+47, "b": -9.0246, "Ea": 166508556.0},
              {"P": 10132500., "A": 5.9632e+56, "b": -11.529, "Ea": 220076726.4}]}
    _index = 3
    _yaml = """
        equation: H2 + O2 <=> 2 OH
        type: pressure-dependent-Arrhenius
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """
    _rc_units = ct.Units("m^3 / kmol / s")

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        cls._rate_obj = ct.ReactionRate.from_dict(cls._rate)

    def eval_rate(self, rate):
        return rate(self.soln.T, self.soln.P)

    def check_rates(self, rates, other):
        # helper function used by deprecation tests
        self.assertEqual(len(rates), len(other))
        for index, item in enumerate(rates):
            P, rate = item
            self.assertNear(P, other[index][0])
            self.assertNear(rate.pre_exponential_factor, other[index][1].pre_exponential_factor)
            self.assertNear(rate.temperature_exponent, other[index][1].temperature_exponent)
            self.assertNear(rate.activation_energy, other[index][1].activation_energy)


class TestChebyshev(ReactionTests, utilities.CanteraTest):
    # test Chebyshev reaction

    _rate_cls = ct.ChebyshevRate
    _rate_type = "Chebyshev"
    _equation = "HO2 <=> OH + O"
    _rate = {
        "type": "Chebyshev",
        "temperature-range": [290., 3000.],
        "pressure-range": [1000., 10000000.0],
        "data":
            [[8.2883, -1.1397, -0.12059, 0.016034],
             [1.9764, 1.0037, 7.2865e-03, -0.030432],
             [0.3177, 0.26889, 0.094806, -7.6385e-03]]}
    _rate_obj = ct.ChebyshevRate(
        temperature_range=(290., 3000.), pressure_range=(1000., 10000000.0),
        data=[[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
              [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
              [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]])
    _index = 4
    _yaml = """
        equation: HO2 <=> OH + O
        type: Chebyshev
        temperature-range: [290.0, 3000.0]
        pressure-range: [9.869232667160128e-03 atm, 98.69232667160128 atm]
        data:
        - [8.2883, -1.1397, -0.12059, 0.016034]
        - [1.9764, 1.0037, 7.2865e-03, -0.030432]
        - [0.3177, 0.26889, 0.094806, -7.6385e-03]
        """
    _rc_units = ct.Units("1 / s")

    def eval_rate(self, rate):
        return rate(self.soln.T, self.soln.P)


class TestCustom(ReactionTests, utilities.CanteraTest):
    # test Custom reaction

    # probe O + H2 <=> H + OH
    _rate_cls = ct.CustomRate
    _equation = "H2 + O <=> H + OH"
    _rate_obj = ct.CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
    _index = 0
    _rate_type = "custom-rate-function"
    _yaml = None
    _rc_units = ct.Units("m^3 / kmol / s")

    def setUp(self):
        # need to overwrite rate to ensure correct type ("method" is not compatible with Func1)
        super().setUp()
        self._rate = lambda T: 38.7 * T**2.7 * exp(-3150.15428/T)

    def from_yaml(self):
        pytest.skip("Reactions with CustomRate do not support YAML")

    def test_roundtrip(self):
        pytest.skip("Reactions with CustomRate do not support roundtrip conversion")

    def test_raises_invalid_rate(self):
        # check exception for instantiation from keywords / invalid rate
        with self.assertRaises(TypeError):
            self.from_rate(tuple())
        with self.assertRaises(TypeError):
            self.from_rate("spam")

    def test_from_func1(self):
        # check instantiation from keywords / rate provided as func1
        f = ct.Func1(self._rate)
        rxn = self.from_rate(f)
        self.check_rxn(rxn)

    def test_rate_func(self):
        # check result of rate expression
        f = ct.Func1(self._rate)
        rate = ct.CustomRate(f)
        self.assertNear(rate(self.soln.T), self.soln.forward_rate_constants[self._index])

    def test_custom_lambda(self):
        # check instantiation from keywords / rate provided as lambda function
        rxn = self.from_rate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
        self.check_rxn(rxn)

    def test_persistent(self):
        # ensure that temporary Python objects are buffered
        gas = ct.Solution(thermo=self.soln.thermo_model,
                         kinetics=self.soln.kinetics_model,
                         species=self.species, reactions=[])
        for i in range(5):
            gas.add_reaction(
                ct.Reaction(equation='H + OH <=> H2O', rate=ct.CustomRate(lambda T: T))
            )
            gas.add_reaction(ct.Reaction(equation='H + OH <=> H2O', rate=lambda T: T))
            def func(T):
                return T
            gas.add_reaction(ct.Reaction(equation='H + OH <=> H2O', rate=func))
        assert (gas.forward_rate_constants == gas.T).all()


class UserRate1Data(ct.ExtensibleRateData):
    def __init__(self):
        self.T = None

    def update(self, gas):
        T = gas.T
        if T != self.T:
            self.T = T
            return True
        else:
            return False


@ct.extension(name="user-rate-1", data=UserRate1Data)
class UserRate1(ct.ExtensibleRate):
    def __init__(self, *args, **kwargs):
        # Do default initialization before calling parent init since that init function
        # may call set_parameters and we don't want to overwrite those values
        self.A = np.nan
        self.eval_error = False
        super().__init__(*args, **kwargs)

    def set_parameters(self, params, units):
        self.A = params["A"]

    def get_parameters(self, params):
        params["A"] = self.A

    def validate(self, equation, soln):
        if np.isnan(self.A):
            raise ValueError("'A' is NaN")

    def eval(self, data):
        if self.eval_error:
            raise ValueError("Error evaluating rate")
        return self.A * data.T**2.7 * exp(-3150.15428/data.T)


class TestExtensible(ReactionTests, utilities.CanteraTest):
    # test general functionality of ExtensibleRate
    _phase_def = """
    phases:
    - name: gas
      thermo: ideal-gas
      species:
      - h2o2.yaml/species: [AR, O, H2, H, OH, O2, H2O, H2O2, HO2]
      - ozone-photolysis.yaml/species: [O3]
      kinetics: gas
      reactions: none
      state: {T: 300.0, P: 1 atm}
    """

    # probe O + H2 <=> H + OH
    _rate_cls = UserRate1
    _equation = "H2 + O <=> H + OH"
    _index = 0
    _rate_type = "user-rate-1"
    _rate = {
        "type": "user-rate-1",
        "A": 38.7
    }
    _yaml = """
        equation: H2 + O <=> H + OH
        type: user-rate-1
        A: 38.7
    """
    _rc_units = ct.Units("m^3 / kmol / s")

    def setUp(self):
        super().setUp()
        self._rate_obj = ct.ReactionRate.from_dict(self._rate)

    def eval_rate(self, rate):
        gas = ct.Solution(yaml=self._phase_def)
        gas.TDY = self.soln.TDY
        gas.add_reaction(ct.Reaction(equation=self._equation, rate=rate))
        return gas.forward_rate_constants[0]

    def test_no_rate(self):
        # Slightly different from the base case since we normally check evaluation via
        # a Kinetics object, which will fail validation
        rxn = self.from_empty()
        with pytest.raises(ct.CanteraError, match="validate"):
            self.eval_rate(rxn.rate)

    def test_parameter_access(self):
        gas = ct.Solution(yaml=self._phase_def)
        R = ct.Reaction.from_yaml(self._yaml, gas)
        assert R.rate.A == 38.7

    def test_set_parameters_error(self):
        with pytest.raises(KeyError):
            # Instantiate with no A factor
            ct.ReactionRate.from_dict({"type": "user-rate-1"})

    def test_standalone_rate(self):
        R = ct.ReactionRate.from_dict({"type": "user-rate-1", "A": 101})
        assert R.type == "user-rate-1"

    def test_eval_error(self):
        # Instantiate with non-numeric A factor to cause an exception during evaluation
        R = ct.ReactionRate.from_dict({"type": "user-rate-1", "A": 12})
        R.eval_error = True
        with pytest.raises(ValueError):
            self.eval_rate(R)


class TestExtensible2(utilities.CanteraTest):
    # Test handling of ExtensibleRate defined in a separate Python module

    _input_template = """
    extensions:
    - type: python
      name: {module}

    phases:
    - name: gas
      thermo: ideal-gas
      species: [{{h2o2.yaml/species: all}}]
      kinetics: gas
      reactions: none
      state: {{T: 300.0, P: 1 atm}}
    """

    @classmethod
    def setUpClass(cls):
        here = str(Path(__file__).parent)
        if here not in sys.path:
            sys.path.append(here)

    def test_load_module(self):
        gas = ct.Solution("extensible-reactions.yaml", transport_model=None)

        for T in np.linspace(300, 3000, 10):
            gas.TP = T, None
            assert gas.forward_rate_constants[0] == pytest.approx(3.14 * T**2)

    def test_missing_module(self):
        with pytest.raises(ct.CanteraError, match="No module named 'fake_ext'"):
            ct.Solution(yaml=self._input_template.format(module="fake_ext"))

    def test_invalid_module(self):
        with pytest.raises(ct.CanteraError, match="SyntaxError"):
            ct.Solution(yaml=self._input_template.format(module="user_ext_invalid"))

    def test_memory_management(self):
        # Make sure objects are being correctly cleaned up and not stuck in
        # mixed Python/C++ ownership cycles
        import user_ext

        for _ in range(3):
            gc.collect()

        initialRate = user_ext.SquareRate.use_count[0]
        initialData = user_ext.SquareRateData.use_count[0]

        def run():
            gas = ct.Solution("extensible-reactions.yaml", transport_model=None)
            assert gas.forward_rate_constants[0] > 0
            assert user_ext.SquareRate.use_count[0] == initialRate + 1
            assert user_ext.SquareRateData.use_count[0] == initialData + 1

            standalone = ct.ReactionRate.from_dict({"type": "square-rate", "A": 101})
            assert user_ext.SquareRate.use_count[0] == initialRate + 2
            assert user_ext.SquareRateData.use_count[0] == initialData + 1

        run()

        # The number of instances for both classes should go back to its previous value
        # after deleting the Solution (may not be zero due to other Solution instances)
        # held by other test classes. Takes a few GC runs to clean up reference cycles
        for _ in range(3):
            gc.collect()
        assert user_ext.SquareRate.use_count[0] == initialRate
        assert user_ext.SquareRateData.use_count[0] == initialData


@ct.extension(name="user-rate-2", data=UserRate1Data)
class UserRate2(ct.ExtensibleRate):
    def set_parameters(self, params, rc_units):
        self.A = params.convert_rate_coeff("A", rc_units)
        self.length = params.convert("L", "m")
        self.Ta = params.convert_activation_energy("Ea", "K")

    def get_parameters(self, params):
        params.set_quantity("A", self.A, self.conversion_units)
        params.set_quantity("L", self.length, "m")
        params.set_activation_energy("Ea", self.Ta, "K")

    def validate(self, equation, soln):
        if self.length < 0:
            raise ValueError(f"Negative length found in reaction {equation}")

    def eval(self, data):
        return self.A * (self.length / 2.0)**2 * exp(-self.Ta/data.T)

class TestExtensible3(utilities.CanteraTest):
    # Additional ExtensibleRate tests

    def setUp(self):
        self.gas = ct.Solution('h2o2.yaml', transport_model=None)

    def test_explicit_units(self):
        rxn = """
        equation: H2 + OH = H2O + H
        type: user-rate-2
        A: 1000 cm^3/kmol/s
        L: 200 cm
        Ea: 1000
        """
        rxn = ct.Reaction.from_yaml(rxn, kinetics=self.gas)
        assert rxn.rate.length == 2
        assert rxn.rate.Ta == pytest.approx(1000 / ct.gas_constant)

    def test_implicit_units(self):
        rxn = """
        equation: H2 + OH = H2O + H
        units: {length: cm}
        type: user-rate-2
        A: 1000
        L: 200
        Ea: 1000
        """
        rxn = ct.Reaction.from_yaml(rxn, kinetics=self.gas)
        assert rxn.rate.length == 2
        assert rxn.rate.Ta == pytest.approx(1000 / ct.gas_constant)

    def test_output_units(self):
        rxn = """
        equation: H2 + OH = H2O + H
        type: user-rate-2
        A: 1000
        L: 200
        Ea: 50
        """
        rxn = ct.Reaction.from_yaml(rxn, kinetics=self.gas)
        N = self.gas.n_reactions
        self.gas.add_reaction(rxn)

        self.gas.write_yaml(self.test_work_path / 'user-rate-units.yaml',
                            units={'length': 'mm', 'activation-energy': 'eV'})

        yml = utilities.load_yaml(self.test_work_path / 'user-rate-units.yaml')
        rxn = yml['reactions'][-1]
        assert rxn['type'] == 'user-rate-2'
        assert rxn['A'] == pytest.approx(1000 * 1000**3)
        assert rxn['L'] == pytest.approx(200 * 1000)
        assert rxn['Ea'] == pytest.approx(50 / ct.faraday)

    def test_validate_error(self):
        rxn = """
        equation: H2 + OH = H2O + H
        type: user-rate-2
        A: 1000
        L: -200
        Ea: 50
        """
        rxn = ct.Reaction.from_yaml(rxn, kinetics=self.gas)
        N = self.gas.n_reactions
        with pytest.raises(ct.CanteraError, match="Negative"):
            self.gas.add_reaction(rxn)


class InterfaceReactionTests(ReactionTests):
    # test suite for surface reaction expressions

    _value = np.NAN # reference value
    _coverage_deps = None

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.soln = ct.Interface("kineticsfromscratch.yaml", "Pt_surf", transport_model=None)
        cls.adj = [cls.soln.adjacent["ohmech"]]
        cls.gas = cls.adj[0]
        cls.species = cls.soln.species()
        cls.concentrations = cls.soln.concentrations
        if cls._rate_cls is not None:
            cls._rate_obj = cls._rate_cls(**cls._rate)
            if cls._coverage_deps:
                cls._rate_obj.coverage_dependencies = cls._coverage_deps

    def setUp(self):
        self.soln.TP = 900, ct.one_atm
        gas = self.adj[0]
        gas.X = "H2:0.05, H2O:0.01, O:1e-4, OH: 1e5, H:2e-5, O2:0.21, AR:0.79"
        gas.TP = 900, ct.one_atm

    def finalize(self, rxn):
        rxn = super().finalize(rxn)
        return rxn

    def eval_rate(self, rate):
        rate.set_species(self.soln.species_names)
        rate.site_density = self.soln.site_density
        self.assertEqual(rate.site_density, self.soln.site_density)
        if "Blowers-Masel" in self._rate_type:
            rate.delta_enthalpy = self.soln.delta_enthalpy[self._index]
        with pytest.warns(UserWarning, match="InterfaceData::update"):
            return rate(self.soln.T, self.soln.coverages)

    def check_rate(self, rate_obj):
        rate = self.eval_rate(rate_obj)
        self.assertNear(self._value, rate)
        self.assertNear(self._value, self.soln.forward_rate_constants[self._index])

    def from_rate(self, rate):
        if isinstance(rate, dict):
            pytest.skip("Detection of rate from dictionary is ambiguous")
        return self.finalize(super().from_rate(rate))

    def test_electrochemistry(self):
        rxn = self.from_yaml()
        sol2 = ct.Interface(thermo="ideal-surface", kinetics="surface",
                    species=self.species, reactions=[rxn], adjacent=self.adj)
        rate2 = sol2.reaction(0).rate
        assert not rate2.uses_electrochemistry
        assert np.isnan(rate2.beta)


class TestArrheniusInterfaceReaction(InterfaceReactionTests, utilities.CanteraTest):
    # test interface reaction without coverages

    _equation = "H(S) + O(S) <=> OH(S) + PT(S)"
    _rate = {"A": 3.7e+20, "b": 0, "Ea": 1.15e7}
    _index = 0
    _rate_type = "interface-Arrhenius"
    _rate_cls = ct.InterfaceArrheniusRate
    _yaml = """
        equation: H(S) + O(S) <=> OH(S) + PT(S)
        rate-constant: {A: 3.7e+20, b: 0, Ea: 11500 J/mol}
        type: interface-Arrhenius
        """
    _rc_units = ct.Units("m^2 / kmol / s")
    _value = 7.9574172975288e+19


class TestArrheniusCoverageReaction(InterfaceReactionTests, utilities.CanteraTest):
    # test interface reaction with coverages

    _equation = "2 O(S) => O2 + 2 PT(S)"
    _rate = {"A": 3.7e+20, "b": 0, "Ea": 213200000.}
    _index = 1
    _rate_type = "interface-Arrhenius"
    _rate_cls = ct.InterfaceArrheniusRate
    _coverage_deps = {"O(S)": {"a": 0.0, "m": 0.0, "E": -60000000.}}
    _yaml = """
        equation: 2 O(S) => O2 + 2 PT(S)
        rate-constant: {A: 3.7e+21, b: 0, Ea: 213200}
        coverage-dependencies:
          O(S): {a: 0.0, m: 0.0, E: -6.0e+04}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-Arrhenius
        """
    _rc_units = ct.Units("m^2 / kmol / s")
    _value = 349029090.19755


class TestBMInterfaceReaction(InterfaceReactionTests, utilities.CanteraTest):
    # test coverage-Blowers-Masel rate expressions with coverage dependency

    _equation = "2 H(S) => H2 + 2 PT(S)"
    _rate = {"A": 3.7e+20, "b": 0, "Ea0": 67400000.0, 'w': 1000000000.0}
    _index = 5
    _rate_type = "interface-Blowers-Masel"
    _rate_cls = ct.InterfaceBlowersMaselRate
    _yaml = """
        equation: 2 H(S) => H2 + 2 PT(S)
        rate-constant: {A: 3.7e+21, b: 0, Ea0: 67400, w: 1000000}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-Blowers-Masel
        """
    _rc_units = ct.Units("m^2 / kmol / s")
    _value = 1.2891970390741e+14


class TestBMCoverageReaction(InterfaceReactionTests, utilities.CanteraTest):
    # test coverage-Blowers-Masel rate expressions with coverage dependency

    _equation = "2 H(S) => H2 + 2 PT(S)"
    _rate = {"A": 3.7e+20, "b": 0, "Ea0": 67400000.0, 'w': 1000000000.0}
    _index = 6
    _coverage_deps = {"H(S)": {"a": 0.0, "m": 0.0, "E": -6000000.0}}
    _rate_type = "interface-Blowers-Masel"
    _rate_cls = ct.InterfaceBlowersMaselRate
    _yaml = """
        equation: 2 H(S) => H2 + 2 PT(S)
        rate-constant: {A: 3.7e+21, b: 0, Ea0: 67400, w: 1000000}
        coverage-dependencies:
          H(S): {a: 0.0, m: 0.0, E: -6000.0}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-Blowers-Masel
        """
    _rc_units = ct.Units("m^2 / kmol / s")
    _value = 1.7068593925679e+14


class StickReactionTests(InterfaceReactionTests):
    # test suite for surface reaction expressions

    _sticking_species = None
    _sticking_order = None

    def finalize(self, rxn):
        rxn = super().finalize(rxn)
        weight = self.gas.molecular_weights[self.gas.species_index(self._sticking_species)]
        rxn.rate.sticking_species = self._sticking_species
        rxn.rate.sticking_order = self._sticking_order
        rxn.rate.sticking_weight = weight
        return rxn

    def from_rate(self, rate):
        rxn = super().from_rate(rate)
        rxn.rate.motz_wise_correction = "Motz-Wise" in self._yaml
        return rxn

    def from_parts(self):
        rxn = super().from_parts()
        rxn.rate.motz_wise_correction = "Motz-Wise" in self._yaml
        return rxn

    def test_sticking_coeffs(self):
        rxn = self.from_yaml()
        if "Motz-Wise" in self._yaml:
            self.assertTrue(rxn.rate.motz_wise_correction)
        else:
            self.assertFalse(rxn.rate.motz_wise_correction)
        weight = self.gas.molecular_weights[self.gas.species_index(self._sticking_species)]
        assert rxn.rate.sticking_species == self._sticking_species
        assert rxn.rate.sticking_order == self._sticking_order
        assert rxn.rate.sticking_weight == pytest.approx(weight)

    def test_site_density(self):
        self.assertEqual(self.soln.site_density,
            self.soln.reaction(self._index).rate.site_density)

    def test_rate_conversion_units(self):
        rxn = self.from_yaml()
        assert str(rxn.rate.conversion_units) == str(ct.Units('1'))


class TestArrheniusStickReaction(StickReactionTests, utilities.CanteraTest):
    # test interface reaction without coverages

    _equation = "H + PT(S) => H(S)"
    _rate = {"A": 1., "b": 0, "Ea": 0.}
    _index = 2
    _rate_type = "sticking-Arrhenius"
    _rate_cls = ct.StickingArrheniusRate
    _sticking_species = "H"
    _sticking_order = 1.0
    _yaml = """
        equation: H + PT(S) => H(S)
        sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: sticking-Arrhenius
        """
    _rc_units = ct.Units("m^3 / kmol / s")
    _value = 401644856274.97


class TestArrheniusCovStickReaction(StickReactionTests, utilities.CanteraTest):
    # test interface reaction with coverages

    _equation = "H2 + 2 PT(S) => 2 H(S)"
    _rate = {"A": 0.046, "b": 0, "Ea": 0.}
    _index = 3
    _rate_type = "sticking-Arrhenius"
    _rate_cls = ct.StickingArrheniusRate
    _coverage_deps = {"PT(S)": {"a": 0.0, "m": -1.0, "E": 0.}}
    _sticking_species = "H2"
    _sticking_order = 2.0
    _yaml = """
        equation: H2 + 2 PT(S) => 2 H(S)
        sticking-coefficient: {A: 0.046, b: 0, Ea: 0}
        coverage-dependencies:
          PT(S): {a: 0.0, m: -1.0, E: 0.0}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: sticking-Arrhenius
        """
    _rc_units = ct.Units("m^5 / kmol^2 / s")
    _value = 1.3792438668539e+19


class TestArrheniusMotzStickReaction(StickReactionTests, utilities.CanteraTest):
    # test interface reaction with coverages

    _equation = "OH + PT(S) => OH(S)"
    _rate = {"A": 1., "b": 0, "Ea": 0.}
    _index = 4
    _rate_type = "sticking-Arrhenius"
    _rate_cls = ct.StickingArrheniusRate
    _sticking_species = "OH"
    _sticking_order = 1.0
    _yaml = """
        equation: OH + PT(S) => OH(S)
        sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
        Motz-Wise: true
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: sticking-Arrhenius
        """
    _rc_units = ct.Units("m^3 / kmol / s")
    _value = 195563866595.97


class TestBlowersMaselStickReaction(StickReactionTests, utilities.CanteraTest):
    # test interface reaction with coverages

    _equation = "OH + PT(S) => OH(S)"
    _rate = {"A": 1., "b": 0, "Ea0": 0., "w": 100000}
    _index = 7
    _rate_type = "sticking-Blowers-Masel"
    _rate_cls = ct.StickingBlowersMaselRate
    _sticking_species = "OH"
    _sticking_order = 1.0
    _yaml = """
        equation: OH + PT(S) => OH(S)
        sticking-coefficient: {A: 1.0, b: 0, Ea0: 0, w: 100000}
        Motz-Wise: true
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: sticking-Blowers-Masel
        """
    _rc_units = ct.Units("m^3 / kmol / s")
    _value = 195563866595.97
