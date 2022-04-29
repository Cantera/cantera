from math import exp
from pathlib import Path
import textwrap

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
        self.assertEqual(rxn1.reaction_type, "three-body")
        self.assertEqual(rxn1.default_efficiency, 0.)
        self.assertEqual(rxn1.efficiencies, {"O2": 1})

        yaml2 = """
            equation: H + O2 + M <=> HO2 + M
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            type: three-body
            default-efficiency: 0
            efficiencies: {O2: 1.0}
            """
        rxn2 = ct.Reaction.from_yaml(yaml2, self.gas)
        self.assertEqual(rxn1.efficiencies, rxn2.efficiencies)
        self.assertEqual(rxn1.default_efficiency, rxn2.default_efficiency)

    def test_duplicate(self):
        # @todo simplify this test
        #     duplicates are currently only checked for import from file
        gas1 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.gas.species(), reactions=[])

        yaml1 = """
            equation: H + O2 + H2O <=> HO2 + H2O
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            """
        rxn1 = ct.Reaction.from_yaml(yaml1, gas1)

        yaml2 = """
            equation: H + O2 + M <=> HO2 + M
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            type: three-body
            default-efficiency: 0
            efficiencies: {H2O: 1}
            """
        rxn2 = ct.Reaction.from_yaml(yaml2, gas1)

        self.assertEqual(rxn1.reaction_type, rxn2.reaction_type)
        self.assertEqual(rxn1.reactants, rxn2.reactants)
        self.assertEqual(rxn1.products, rxn2.products)
        self.assertEqual(rxn1.efficiencies, rxn2.efficiencies)
        self.assertEqual(rxn1.default_efficiency, rxn2.default_efficiency)

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
        self.assertEqual(rxn.reaction_type, "reaction")

    def test_not_three_body(self):
        # check that insufficient reactants prevent automatic conversion
        yaml = """
            equation: HCNO + H <=> H + HNCO  # Reaction 270
            rate-constant: {A: 2.1e+15, b: -0.69, Ea: 2850.0}
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "reaction")

    def test_user_override(self):
        # check that type specification prevents automatic conversion
        yaml = """
            equation: H + 2 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            type: elementary
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "reaction")


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
        ct.use_legacy_rate_constants(False)
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
        # test custom units
        units = "units: {length: cm, quantity: mol}"
        yaml = f"{textwrap.dedent(self._yaml)}\n{units}"
        with self.assertRaisesRegex(Exception, "not supported"):
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
        with self.assertRaisesRegex(Exception, "not supported"):
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
    _type = "Lindemann"
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
    _type = "Troe"
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
    _type = "SRI"
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
    _type = "Tsang"
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
        with self.assertRaisesRegex(Exception, "not supported"):
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


class SurfaceReactionRateTests(ReactionRateTests):
    # test suite for surface reaction rate expressions

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        ct.use_legacy_rate_constants(False)
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
    _rate_obj = None # reaction rate object
    _kwargs = {} # additional parameters required by constructor
    _index = None # index of reaction in "kineticsfromscratch.yaml"
    _rxn_type = "reaction" # name of reaction type
    _rate_type = None # name of reaction rate type
    _legacy = False # object uses legacy framework
    _legacy_uses_rate = True # legacy object implements rate property
    _yaml = None # YAML parameterization
    _deprecated_getters = {} # test legacy getters (old framework)
    _deprecated_setters = {} # test legacy setters (old framework)
    _deprecated_callers = {} # test legacy callers (old framework)

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

    def from_yaml(self, deprecated=False):
        # create reaction object from yaml
        if deprecated:
            with self.assertWarnsRegex(DeprecationWarning, "is renamed to 'from_yaml'"):
                rxn = ct.Reaction.fromYaml(self._yaml, kinetics=self.soln)
        else:
            rxn = self.finalize(ct.Reaction.from_yaml(self._yaml, kinetics=self.soln))
        return self.finalize(rxn)

    def from_dict(self):
        # create reaction rate object from input data
        input_data = self.from_yaml().input_data
        rxn = ct.Reaction.from_dict(input_data, kinetics=self.soln)
        return self.finalize(rxn)

    def from_empty(self):
        # create reaction object with uninitialized rate
        if self._rate_cls is None:
            rate = None
        else:
            # Create an "empty" rate of the correct type for merged reaction types where
            # the only way they are distinguished is by the rate type
            rate = self._rate_cls()
        rxn = self._cls(equation=self._equation, rate=rate, kinetics=self.soln,
                        legacy=self._legacy, **self._kwargs)
        return self.finalize(rxn)

    def from_rate(self, rate):
        if not self._legacy_uses_rate:
            pytest.skip("Legacy: rate object is not defined [1]")
        if rate is None and self._rate is None:
            # this does not work when no specialized reaction class exists
            pytest.skip("Creation from dictionary is not supported")
        rxn = self._cls(equation=self._equation, rate=rate, kinetics=self.soln,
                        legacy=self._legacy, **self._kwargs)
        return self.finalize(rxn)

    def from_parts(self):
        # create reaction rate object from parts
        if self._rate_obj is None:
            pytest.skip("Legacy: rate object is not defined [2]")
        orig = self.soln.reaction(self._index)
        if self._legacy:
            rxn = self._cls(orig.reactants, orig.products, legacy=self._legacy)
            rxn.rate = self._rate_obj
            rxn.reversible = "<=>" in self._equation
            return rxn

        rxn = self._cls(orig.reactants, orig.products, rate=self._rate_obj,
                        legacy=self._legacy)
        rxn.reversible = "<=>" in self._equation
        return self.finalize(rxn)

    def check_rate(self, rate_obj):
        if self._legacy:
            rate = rate_obj(self.soln.T)
        else:
            rate = self.eval_rate(rate_obj)
        self.assertNear(rate, self.soln.forward_rate_constants[self._index])

    def check_rxn(self, rxn, check_legacy=True):
        # helper function that checks reaction configuration
        ix = self._index
        self.assertEqual(rxn.reactants, self.soln.reaction(ix).reactants)
        self.assertEqual(rxn.products, self.soln.reaction(ix).products)
        if check_legacy:
            # self.assertEqual(rxn.uses_legacy, self._rxn_type.endswith("-legacy"))
            self.assertEqual(rxn.uses_legacy, self._legacy)
            self.assertEqual(rxn.reaction_type, self._rxn_type)
            if not rxn.uses_legacy:
                self.assertEqual(rxn.rate.type, self._rate_type)

        if not self._legacy:
            # legacy rate evaluation is not consistent
            self.check_rate(rxn.rate)
        if self.soln.thermo_model.lower() == "surf":
            sol2 = ct.Interface(thermo="Surface", kinetics="interface",
                                species=self.species, reactions=[rxn], adjacent=self.adj)
            sol2.site_density = self.soln.site_density
            sol2.coverages = self.soln.coverages
            sol2.TP = self.soln.TP
        else:
            sol2 = ct.Solution(thermo=self.soln.thermo_model, kinetics=self.soln.kinetics_model,
                               species=self.species, reactions=[rxn])
            sol2.TPX = self.soln.TPX
        self.check_solution(sol2, check_legacy)

    def check_solution(self, sol2, check_legacy=True):
        # helper function that checks evaluation of reaction rates
        ix = self._index
        if check_legacy:
            self.assertEqual(sol2.reaction(0).reaction_type, self._rxn_type)
        self.assertNear(sol2.forward_rate_constants[0],
                        self.soln.forward_rate_constants[ix])
        self.assertNear(sol2.net_rates_of_progress[0],
                        self.soln.net_rates_of_progress[ix])

    def test_rate(self):
        # check consistency of reaction rate and forward rate constant
        if self._rate_obj is None:
            pytest.skip("Legacy: rate object is not defined [3]")
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
        self.check_rxn(self.from_yaml(deprecated=True))

    def test_from_dict(self):
        # check instantiation from a yaml dictionary (input_data)
        # cannot compare types as input_data does not recreate legacy objects
        self.check_rxn(self.from_dict(), check_legacy=False)

    def test_add_rxn(self):
        # check adding new reaction to solution
        if self._yaml is not None:
            rxn = self.from_yaml()
        else:
            rxn = self.from_rate(self._rate_obj)
        if self.soln.thermo_model.lower() == "surf":
            sol2 = ct.Interface(thermo="Surface", kinetics="interface",
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
        if self._rate_obj is None:
            return
        rxn = self.from_empty()
        if self._legacy:
            self.assertNear(rxn.rate(self.soln.T), 0.)
        else:
            self.assertIsNaN(self.eval_rate(rxn.rate))

        if self._legacy:
            sol2 = ct.Solution(thermo=self.soln.thermo_model,
                               kinetics=self.soln.kinetics_model,
                               species=self.species,
                               reactions=[rxn], adjacent=self.adj)
            sol2.TPX = self.soln.TPX
            self.assertNear(sol2.forward_rate_constants[0], 0.)
            self.assertNear(sol2.net_rates_of_progress[0], 0.)
        else:
            with self.assertRaisesRegex(ct.CanteraError, "validate"):
                ct.Solution(thermo=self.soln.thermo_model,
                            kinetics=self.soln.kinetics_model,
                            species=self.species, reactions=[rxn], adjacent=self.adj)

    def test_replace_rate(self):
        # check replacing reaction rate expression
        if not self._legacy_uses_rate:
            pytest.skip("Legacy: rate property not implemented")
        elif self._yaml is not None:
            rxn = self.from_yaml()
        else:
            rxn = self.from_rate(self._rate_obj)
        rxn = self.from_empty()
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_roundtrip(self):
        # check round-trip instantiation via input_data
        rxn = self.from_yaml()
        if self._legacy:
            pytest.skip("Legacy: round-trip conversion is not supported")
        rxn = self.from_rate(self._rate_obj)
        rate_input_data = rxn.rate.input_data
        rate_obj = rxn.rate.__class__(input_data=rate_input_data)
        rxn2 = self.from_rate(rate_obj)
        self.check_rxn(rxn2)

    def check_equal(self, one, two):
        # helper function for deprecation tests
        self.assertEqual(type(one), type(two))
        if isinstance(one, (list, tuple, np.ndarray)):
            self.assertArrayNear(one, two)
        elif isinstance(one, (dict, str)):
            assert one == two
        else:
            self.assertNear(one, two)

    def test_deprecated_getters(self):
        # check property getters deprecated in new framework
        rxn = self.from_yaml()
        for attr, default in self._deprecated_getters.items():
            if self._legacy:
                self.check_equal(getattr(rxn, attr), default)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "This property is"):
                    try:
                        self.check_equal(getattr(rxn, attr), default)
                    except Exception as err:
                        print(f"Exception raised when testing getter for '{attr}'")
                        raise err

    def test_deprecated_setters(self):
        # check property setters deprecated in new framework
        rxn = self.from_yaml()
        for attr, new in self._deprecated_setters.items():
            if self._legacy:
                setattr(rxn, attr, new)
                self.check_equal(getattr(rxn, attr), new)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "This property is"):
                    setattr(rxn, attr, new)
                with self.assertWarnsRegex(DeprecationWarning, "This property is"):
                    self.check_equal(getattr(rxn, attr), new)

    def test_deprecated_callers(self):
        # check methods deprecated in new framework
        rxn = self.from_yaml()
        for state, value in self._deprecated_callers.items():
            T, P = state
            if self._legacy:
                self.check_equal(rxn(T, P), value)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "method is moved"):
                    self.check_equal(rxn(T, P), value)


class TestElementary2(ReactionTests, utilities.CanteraTest):
    # test legacy version of elementary reaction

    _cls = ct.ElementaryReaction
    _equation = "H2 + O <=> H + OH"
    _rate = {"A": 38.7, "b": 2.7, "Ea": 2.619184e+07}
    _index = 0
    _rxn_type = "elementary-legacy"
    _legacy = True
    _yaml = """
        equation: O + H2 <=> H + OH
        type: elementary-legacy
        rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
        """
    _deprecated_getters = {"allow_negative_pre_exponential_factor": False}
    _deprecated_setters = {"allow_negative_pre_exponential_factor": True}

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        if cls._legacy:
            args = list(cls._rate.values())
            cls._rate_obj = ct.Arrhenius(*args)
        else:
            cls._rate_obj = ct.ArrheniusRate(**cls._rate)

    def test_arrhenius(self):
        # test assigning Arrhenius rate
        rate = ct.Arrhenius(self._rate["A"], self._rate["b"], self._rate["Ea"])
        rxn = self.from_empty()
        if self._legacy:
            rxn.rate = rate
        else:
            with self.assertWarnsRegex(DeprecationWarning, "'Arrhenius' object is deprecated"):
                rxn.rate = rate
        self.check_rxn(rxn)


class TestElementary(TestElementary2):
    # test updated version of elementary reaction

    _cls = ct.Reaction
    _rxn_type = "reaction"
    _rate_type = "Arrhenius"
    _rate_cls = ct.ArrheniusRate
    _legacy = False
    _yaml = """
        equation: O + H2 <=> H + OH
        rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
        """


class TestThreeBody2(TestElementary2):
    # test legacy version of three-body reaction

    _cls = ct.ThreeBodyReaction
    _equation = "2 O + M <=> O2 + M"
    _rate = {"A": 1.2e11, "b": -1.0, "Ea": 0.0}
    _kwargs = {"efficiencies": {"H2": 2.4, "H2O": 15.4, "AR": 0.83}}
    _index = 1
    _rxn_type = "three-body-legacy"
    _legacy = True
    _yaml = """
        equation: 2 O + M <=> O2 + M
        type: three-body-legacy
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """

    def from_parts(self):
        rxn = ReactionTests.from_parts(self)
        rxn.efficiencies = self._kwargs["efficiencies"]
        return rxn

    def test_efficiencies(self):
        # check efficiencies
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.soln,
                        legacy=self._legacy, **self._kwargs)

        self.assertEqual(rxn.efficiencies, self._kwargs["efficiencies"])


class TestThreeBody(TestThreeBody2):
    # test updated version of three-body reaction

    _legacy = False
    _rxn_type = "three-body"
    _rate_type = "Arrhenius"
    _yaml = """
        equation: 2 O + M <=> O2 + M
        type: three-body
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """


class TestImplicitThreeBody(TestThreeBody):
    # test three-body reactions with explicit collision parther

    _equation = "H + 2 O2 <=> HO2 + O2"
    _rate = {"A": 2.08e+19, "b": -1.24, "Ea": 0.0}
    _kwargs = {}
    _index = 5
    _yaml = """
        equation: H + 2 O2 <=> HO2 + O2
        rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
        """

    def from_parts(self):
        rxn = ReactionTests.from_parts(self)
        rxn.efficiencies = {"O2": 1.}
        rxn.default_efficiency = 0
        return rxn

    def test_efficiencies(self):
        # overload of default tester
        rxn = self.from_rate(self._rate_obj)
        self.assertEqual(rxn.efficiencies, {"O2": 1.})
        self.assertEqual(rxn.default_efficiency, 0.)


class TestTwoTempPlasma(ReactionTests, utilities.CanteraTest):
    # test two-temperature plasma reaction

    _rate_cls = ct.TwoTempPlasmaRate
    _rate_type = "two-temperature-plasma"
    _equation = "O + H => O + H"
    _rate_obj = ct.TwoTempPlasmaRate(A=17283, b=-3.1, Ea_gas=-5820000, Ea_electron=1081000)
    _index = 11
    _yaml = """
        equation: O + H => O + H
        type: two-temperature-plasma
        rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
        """

    def eval_rate(self, rate):
        return rate(self.soln.T, self.soln.Te)

    def test_reversible(self):
        orig = self.soln.reaction(self._index)
        rxn = ct.Reaction(orig.reactants, orig.products, rate=self._rate_obj)
        with self.assertRaisesRegex(ct.CanteraError, "does not support reversible"):
            self.check_rxn(rxn)


class TestBlowersMasel(ReactionTests, utilities.CanteraTest):
    # test updated version of Blowers-Masel reaction

    _rate_cls = ct.BlowersMaselRate
    _rate_type = "Blowers-Masel"
    _equation = "O + H2 <=> H + OH"
    _rate_obj = ct.BlowersMaselRate(A=38700, b=2.7, Ea0=1.0958665856e8, w=1.7505856e13)
    _index = 6
    _yaml = """
        equation: O + H2 <=> H + OH
        type: Blowers-Masel
        rate-constant: {A: 38700, b: 2.7, Ea0: 2.619184e4 cal/mol, w: 4.184e9 cal/mol}
        """

    def eval_rate(self, rate):
        rate.delta_enthalpy = self.soln.delta_enthalpy[self._index]
        with pytest.warns(UserWarning, match="BlowersMaselData::update"):
            return rate(self.soln.T)


class TestTroe2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Troe falloff reaction

    _cls = ct.FalloffReaction
    _equation = "2 OH (+ M) <=> H2O2 (+ M)"
    _kwargs = {"efficiencies": {"AR": 0.7, "H2": 2.0, "H2O": 6.0}}
    _index = 2
    _rxn_type = "falloff-legacy"
    _legacy = True
    _legacy_uses_rate = False
    _yaml = """
        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
        type: falloff-legacy
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
        """


class TestTroe(ReactionTests, utilities.CanteraTest):

    _cls = ct.FalloffReaction
    _equation = "2 OH (+ M) <=> H2O2 (+ M)"
    _rate = {
        "type": "falloff",
        "low_P_rate_constant": {"A": 2.3e+12, "b": -0.9, "Ea": -7112800.0},
        "high_P_rate_constant": {"A": 7.4e+10, "b": -0.37, "Ea": 0.0},
        "Troe": {"A": 0.7346, "T3": 94.0, "T1": 1756.0, "T2": 5182.0}
        }
    _kwargs = {"efficiencies": {"AR": 0.7, "H2": 2.0, "H2O": 6.0}}
    _index = 2
    _rxn_type = "falloff"
    _rate_type = "Troe"
    _yaml = """
        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
        type: falloff
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
        """

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

    def from_parts(self):
        rxn = ReactionTests.from_parts(self)
        rxn.efficiencies = self._kwargs["efficiencies"]
        return rxn


class TestLindemann(ReactionTests, utilities.CanteraTest):
    # test Lindemann falloff reaction

    _cls = ct.FalloffReaction
    _equation = "2 OH (+ M) <=> H2O2 (+ M)"
    _rate = {
        "type": "falloff",
        "low_P_rate_constant": {"A": 2.3e+12, "b": -0.9, "Ea": -7112800.0},
        "high_P_rate_constant": {"A": 7.4e+10, "b": -0.37, "Ea": 0.0}
        }
    _kwargs = {"efficiencies": {"AR": 0.7, "H2": 2.0, "H2O": 6.0}}
    _index = 7
    _rxn_type = "falloff"
    _rate_type = "Lindemann"
    _legacy = False
    _yaml = """
        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 8
        duplicate: true
        type: falloff
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
        """

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

    def from_parts(self):
        rxn = ReactionTests.from_parts(self)
        rxn.efficiencies = self._kwargs["efficiencies"]
        return rxn


class TestChemicallyActivated2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Chemically Activated falloff reaction

    _cls = ct.ChemicallyActivatedReaction
    _equation = "H2O + OH (+M) <=> HO2 + H2 (+M)"
    _index = 10
    _rxn_type = "chemically-activated-legacy"
    _legacy = True
    _legacy_uses_rate = False
    _yaml = """
        equation: H2O + OH (+M) <=> HO2 + H2 (+M)  # Reaction 11
        units: {length: cm, quantity: mol, activation-energy: cal/mol}
        type: chemically-activated-legacy
        low-P-rate-constant: [282320.078, 1.46878, -3270.56495]
        high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
        """


class TestChemicallyActivated(ReactionTests, utilities.CanteraTest):
    # test Chemically Activated falloff reaction

    _cls = ct.FalloffReaction
    _equation = "H2O + OH (+M) <=> HO2 + H2 (+M)"
    _rate = {
        "type": "chemically-activated",
        "low_P_rate_constant": {"A": 282.320078, "b": 1.46878, "Ea": -13684043.7508},
        "high_P_rate_constant": {"A": 5.88E-14, "b": 6.721, "Ea": -12644997.768}
        }
    _index = 10
    _rxn_type = "chemically-activated"
    _rate_type = "Lindemann"
    _yaml = """
        equation: H2O + OH (+M) <=> HO2 + H2 (+M)  # Reaction 11
        units: {length: cm, quantity: mol, activation-energy: cal/mol}
        type: chemically-activated
        low-P-rate-constant: [282320.078, 1.46878, -3270.56495]
        high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
        """

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


class TestPlog2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Plog reaction

    _cls = ct.PlogReaction
    _equation = "H2 + O2 <=> 2 OH"
    _rate = [(1013.25, ct.Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
             (101325., ct.Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
             (1013250., ct.Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
             (10132500., ct.Arrhenius(5.9632e+56, -11.529, 220076726.4))]
    _index = 3
    _rxn_type = "pressure-dependent-Arrhenius-legacy"
    _legacy = True
    _legacy_uses_rate = False
    _yaml = """
        equation: H2 + O2 <=> 2 OH
        type: pressure-dependent-Arrhenius-legacy
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """
    _deprecated_callers = {(1000., ct.one_atm): 530968934612.9017}

    def check_rates(self, rates, other):
        # helper function used by deprecation tests
        self.assertEqual(len(rates), len(other))
        for index, item in enumerate(rates):
            P, rate = item
            self.assertNear(P, other[index][0])
            self.assertNear(rate.pre_exponential_factor, other[index][1].pre_exponential_factor)
            self.assertNear(rate.temperature_exponent, other[index][1].temperature_exponent)
            self.assertNear(rate.activation_energy, other[index][1].activation_energy)

    def test_deprecated_getters(self):
        # overload default tester for deprecated property getters
        rxn = self.from_yaml()
        if self._legacy:
            self.check_rates(rxn.rates, self._rate)
        else:
            with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                self.check_rates(rxn.rates, TestPlog2._rate)

    def test_deprecated_setters(self):
        # overload default tester for deprecated property setters
        rate = ct.PlogRate(TestPlog2._rate)
        rates = rate.rates

        rxn = self.from_yaml()
        if self._legacy:
            rxn.rates = rates
            self.check_rates(rxn.rates, self._rate)
        else:
            with self.assertWarnsRegex(DeprecationWarning, "Setter is replaceable"):
                rxn.rates = rates
            with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                self.check_rates(rxn.rates, TestPlog2._rate)


class TestPlog(TestPlog2):
    # test updated version of Plog reaction

    _cls = ct.Reaction
    _rate_cls = ct.PlogRate
    _rate = {
        "type": "pressure-dependent-Arrhenius",
        "rate-constants":
             [{"P": 1013.25, "A": 1.2124e+16, "b": -0.5779, "Ea": 45491376.8},
              {"P": 101325., "A": 4.9108e+31, "b": -4.8507, "Ea": 103649395.2},
              {"P": 1013250., "A": 1.2866e+47, "b": -9.0246, "Ea": 166508556.0},
              {"P": 10132500., "A": 5.9632e+56, "b": -11.529, "Ea": 220076726.4}]}
    _rxn_type = "reaction"
    _rate_type = "pressure-dependent-Arrhenius"
    _legacy = False
    _legacy_uses_rate = True
    _yaml = """
        equation: H2 + O2 <=> 2 OH
        type: pressure-dependent-Arrhenius
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        cls._rate_obj = ct.ReactionRate.from_dict(cls._rate)

    def eval_rate(self, rate):
        return rate(self.soln.T, self.soln.P)


class TestChebyshev2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Chebyshev reaction

    _cls = ct.ChebyshevReaction
    _equation = "HO2 <=> OH + O"
    _rate = {"temperature_range": (290., 3000.), "pressure_range": (1000., 10000000.0),
             "data": [[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
                      [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
                      [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]]}
    _index = 4
    _rxn_type = "Chebyshev-legacy"
    _legacy = True
    _legacy_uses_rate = False
    _yaml = """
        equation: HO2 <=> OH + O
        type: Chebyshev-legacy
        temperature-range: [290.0, 3000.0]
        pressure-range: [9.869232667160128e-03 atm, 98.69232667160128 atm]
        data:
        - [8.2883, -1.1397, -0.12059, 0.016034]
        - [1.9764, 1.0037, 7.2865e-03, -0.030432]
        - [0.3177, 0.26889, 0.094806, -7.6385e-03]
        """
    _deprecated_getters = {"nPressure": 4, "nTemperature": 3}
    _deprecated_callers = {(1000., ct.one_atm): 2858762454.1119065}

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        cls._deprecated_getters.update({"coeffs": np.array(TestChebyshev2._rate["data"])})
        cls._deprecated_getters.update(
            {k: v for k, v in TestChebyshev2._rate.items()
                if k not in ["data", "temperature_range", "pressure_range"]})


class TestChebyshev(TestChebyshev2):
    # test updated version of Chebyshev reaction

    _cls = ct.Reaction
    _rate_cls = ct.ChebyshevRate
    _rxn_type = "reaction"
    _rate_type = "Chebyshev"
    _rate = None
    _rate_obj = ct.ChebyshevRate(
        temperature_range=(290., 3000.), pressure_range=(1000., 10000000.0),
        data=[[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
              [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
              [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]])
    _legacy = False
    _legacy_uses_rate = True
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

    def eval_rate(self, rate):
        return rate(self.soln.T, self.soln.P)


class TestCustom(ReactionTests, utilities.CanteraTest):
    # test Custom reaction

    # probe O + H2 <=> H + OH
    _cls = ct.CustomReaction
    _equation = "H2 + O <=> H + OH"
    _rate_obj = ct.CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
    _index = 0
    _rxn_type = "custom-rate-function"
    _rate_type = "custom-rate-function"
    _legacy = False
    _yaml = None

    def setUp(self):
        # need to overwrite rate to ensure correct type ("method" is not compatible with Func1)
        super().setUp()
        self._rate = lambda T: 38.7 * T**2.7 * exp(-3150.15428/T)

    def from_yaml(self, deprecated=False):
        pytest.skip("CustomReaction does not support YAML")

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


class InterfaceReactionTests(ReactionTests):
    # test suite for surface reaction expressions

    _value = np.NAN # reference value (obtained from legacy framework)
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
        if cls._legacy:
            args = list(cls._rate.values())
            cls._rate_obj = ct.Arrhenius(*args)
            cls._cls = ct.InterfaceReaction
            cls._rxn_type = "interface-legacy"
            cls._rate_type = None
            cls._rate_cls = None
            cls._legacy_uses_rate = True

    def setUp(self):
        self.soln.TP = 900, ct.one_atm
        gas = self.adj[0]
        gas.X = "H2:0.05, H2O:0.01, O:1e-4, OH: 1e5, H:2e-5, O2:0.21, AR:0.79"
        gas.TP = 900, ct.one_atm

    def finalize(self, rxn):
        rxn = super().finalize(rxn)
        if self._legacy and self._coverage_deps:
            if not rxn.coverage_deps:
                # legacy coverage dependencies use tuples
                coverage_deps = {key: tuple(val.values()) for key, val in self._coverage_deps.items()}
                rxn.coverage_deps = coverage_deps
        return rxn

    def eval_rate(self, rate):
        if self._legacy:
            return super().eval_rate(rate)
        rate.set_species(self.soln.species_names)
        rate.site_density = self.soln.site_density
        self.assertEqual(rate.site_density, self.soln.site_density)
        if "Blowers-Masel" in self._rate_type:
            rate.delta_enthalpy = self.soln.delta_enthalpy[self._index]
        with pytest.warns(UserWarning, match="InterfaceData::update"):
            return rate(self.soln.T, self.soln.coverages)

    def check_rate(self, rate_obj):
        if self._legacy:
            rate = rate_obj(self.soln.T)
        else:
            rate = self.eval_rate(rate_obj)
        self.assertNear(self._value, rate)
        self.assertNear(self._value, self.soln.forward_rate_constants[self._index])

    def from_rate(self, rate):
        if isinstance(rate, dict) and not self._legacy:
            pytest.skip("Detection of rate from dictionary is ambiguous")
        return self.finalize(super().from_rate(rate))

    def test_replace_rate(self):
        if self._legacy:
            pytest.skip("Legacy: modifying reactions is not supported")
        super().test_replace_rate()

    def test_from_parts(self):
        if self._legacy and self._coverage_deps:
            pytest.skip("Legacy: construction from parts is not tested")
        super().test_from_parts()

    def test_rate(self):
        if self._legacy and self._coverage_deps:
            pytest.skip("Legacy: interface rate does not include coverage terms")
        super().test_rate()

    def test_electrochemistry(self):
        if self._legacy:
            pytest.skip("Legacy: property uses_electrochemisty not implemented")
        rxn = self.from_yaml()
        sol2 = ct.Interface(thermo="Surface", kinetics="interface",
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
    _value = 7.9574172975288e+19
    _deprecated_getters = {
        "coverage_deps": {},
        "is_sticking_coefficient": False}
    _deprecated_setters = {
        "coverage_deps": {"O(S)": {"a": 0.0, "m": 0.0, "E": -60000000.}}}


class TestArrheniusInterfaceReaction2(TestArrheniusInterfaceReaction):
    # test legacy version of interface reaction

    _legacy = True
    _yaml = """
        equation: H(S) + O(S) <=> OH(S) + PT(S)
        rate-constant: {A: 3.7e+20, b: 0, Ea: 11500 J/mol}
        type: interface-legacy
        """
    _deprecated_setters = {"coverage_deps": {"O(S)": (0.0, 0.0, -60000000.)}}


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
    _value = 349029090.19755
    _deprecated_getters = {
        "coverage_deps": {"O(S)": {"a": 0.0, "m": 0.0, "E": -60000000.}},
        "is_sticking_coefficient": False}
    _deprecated_setters = {"coverage_deps": {}}


class TestArrheniusCoverageReaction2(TestArrheniusCoverageReaction):
    # test interface reaction with coverages

    _legacy = True
    _yaml = """
        equation: 2 O(S) => O2 + 2 PT(S)
        rate-constant: {A: 3.7e+21, b: 0, Ea: 213200}
        coverage-dependencies:
          O(S): {a: 0.0, m: 0.0, E: -6.0e+04}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-legacy
        """
    _deprecated_getters = {"coverage_deps": {"O(S)": (0.0, 0.0, -60000000.)}}


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
    _value = 1.7068593925679e+14


class StickReactionTests(InterfaceReactionTests):
    # test suite for surface reaction expressions

    _sticking_species = None
    _sticking_order = None

    def finalize(self, rxn):
        rxn = super().finalize(rxn)
        if not self._legacy:
            weight = self.gas.molecular_weights[self.gas.species_index(self._sticking_species)]
            rxn.rate.sticking_species = self._sticking_species
            rxn.rate.sticking_order = self._sticking_order
            rxn.rate.sticking_weight = weight
        return rxn

    def from_rate(self, rate):
        rxn = super().from_rate(rate)
        if self._legacy:
            pytest.skip("Legacy: construction from rate is not tested")
        rxn.rate.motz_wise_correction = "Motz-Wise" in self._yaml
        return rxn

    def from_parts(self):
        rxn = super().from_parts()
        rxn.rate.motz_wise_correction = "Motz-Wise" in self._yaml
        return rxn

    def test_rate(self):
        if self._legacy:
            pytest.skip("Legacy: interface rate does not include sticking terms")
        super().test_rate()

    @pytest.mark.skip("Legacy: construction from parts is not tested")
    def test_from_parts(self):
        super().test_from_parts()

    def test_sticking_coeffs(self):
        if self._legacy:
            pytest.skip("Legacy: explicit sticking coefficients are not supported")
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
        if self._legacy:
            pytest.skip("Legacy: interface rate does not include site density")
        self.assertEqual(self.soln.site_density,
            self.soln.reaction(self._index).rate.site_density)


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
    _value = 401644856274.97
    _deprecated_getters = {
        "coverage_deps": {},
        "sticking_species": "H",
        "is_sticking_coefficient": True,
        "use_motz_wise_correction": False}
    _deprecated_setters = {
        "coverage_deps": {"O(S)": {"a": 0.0, "m": 0.0, "E": -60000000.}}}


class TestArrheniusStickReaction2(TestArrheniusStickReaction):
    # test legacy interface reaction without coverages

    _legacy = True
    _yaml = """
        equation: H + PT(S) => H(S)
        sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-legacy
        """
    _deprecated_getters = {
        "coverage_deps": {},
        "is_sticking_coefficient": True,
        "use_motz_wise_correction": False}
    _deprecated_setters = {"coverage_deps": {"O(S)": (0.0, 0.0, -60000000.)}}


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
    _value = 1.3792438668539e+19


class TestArrheniusCovStickReaction2(TestArrheniusCovStickReaction):
    # test interface reaction with coverages

    _legacy = True
    _yaml = """
        equation: H2 + 2 PT(S) => 2 H(S)
        sticking-coefficient: {A: 0.046, b: 0, Ea: 0}
        coverage-dependencies:
          PT(S): {a: 0.0, m: -1.0, E: 0.0}
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-legacy
        """


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
    _value = 195563866595.97
    _deprecated_getters = {
        "coverage_deps": {},
        "is_sticking_coefficient": True,
        "sticking_species": "OH",
        "use_motz_wise_correction": True}
    _deprecated_setters = {
        "coverage_deps": {"O(S)": {"a": 0.0, "m": 0.0, "E": -60000000.}}}


class TestArrheniusMotzStickReaction2(TestArrheniusMotzStickReaction):
    # test interface reaction with coverages

    _legacy = True
    _yaml = """
        equation: OH + PT(S) => OH(S)
        sticking-coefficient: {A: 1.0, b: 0, Ea: 0}
        Motz-Wise: true
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        type: interface-legacy
        """
    _deprecated_getters = {
        "coverage_deps": {},
        "is_sticking_coefficient": True,
        "use_motz_wise_correction": True}
    _deprecated_setters = {"coverage_deps": {"O(S)": (0.0, 0.0, -60000000.)}}


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
    _value = 195563866595.97
