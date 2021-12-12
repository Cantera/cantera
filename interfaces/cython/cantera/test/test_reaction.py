from math import exp
from pathlib import Path
import textwrap

import cantera as ct
import numpy as np
from . import utilities


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
            gas2 = ct.Solution(fname)

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
        self.assertEqual(rxn.reaction_type, "elementary")

    def test_not_three_body(self):
        # check that insufficient reactants prevent automatic conversion
        yaml = """
            equation: HCNO + H <=> H + HNCO  # Reaction 270
            rate-constant: {A: 2.1e+15, b: -0.69, Ea: 2850.0}
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "elementary")

    def test_user_override(self):
        # check that type specification prevents automatic conversion
        yaml = """
            equation: H + 2 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            type: elementary
            """
        rxn = ct.Reaction.from_yaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "elementary")


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
        cls.gas = ct.Solution("kineticsfromscratch.yaml")
        cls.gas.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5"
        cls.gas.TP = 900, 2*ct.one_atm

    def from_parts(self):
        # create reaction rate object from parts
        return self._cls(**self._parts)

    def from_input(self, input=None):
        # create reaction rate object from input_data
        if input is None:
            input = self._input
        else:
            self.assertIsInstance(input, dict)
        return self._cls(input_data=input)

    def from_yaml(self):
        # create reaction rate object from yaml
        return ct.ReactionRate.from_yaml(self._yaml)

    def from_dict(self, input=None):
        # create reaction rate object from dictionary
        if input is None:
            input = self.from_parts().input_data
        else:
            self.assertIsInstance(input, dict)
        return ct.ReactionRate.from_dict(input)

    def eval(self, rate):
        # evaluate rate expression
        return rate(self.gas.T)

    def check_rate(self, rate):
        # check rates
        self.assertEqual(self._type, rate.type)
        self.assertIn(self._cls.__name__, f"{rate}")
        value = self.eval(rate)
        self.assertIsFinite(value)
        self.assertNear(value, self.gas.forward_rate_constants[self._index])

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
        rate0 = self.from_parts()
        input_data = rate0.input_data
        rate1 = self.from_dict(input_data)
        self.assertNear(self.eval(rate0), self.eval(rate1))

    def test_with_units(self):
        units = "units: {length: cm, quantity: mol}"
        yaml = f"{textwrap.dedent(self._yaml)}\n{units}"
        with self.assertRaisesRegex(Exception, "not supported"):
            ct.ReactionRate.from_yaml(yaml)

    def test_third_body(self):
        concm = self.gas.third_body_concentrations
        self.assertIsNaN(concm[self._index])


class TestArrheniusRate(ReactionRateTests, utilities.CanteraTest):
    # test Arrhenius rate expressions

    _cls = ct.ArrheniusRate
    _type = "Arrhenius"
    _index = 0
    _input = {"rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}}
    _yaml = "rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}"

    @classmethod
    def setUpClass(cls):
        ReactionRateTests.setUpClass()
        cls._parts = cls._input["rate-constant"]

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        self.assertEqual(self._parts["Ea"], rate.activation_energy)
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


class TestBlowersMaselRate(ReactionRateTests, utilities.CanteraTest):
    # test Blowers-Masel rate expressions

    _cls = ct.BlowersMaselRate
    _type = "Blowers-Masel"
    _index = 6
    _input = {"rate-constant": {"A": 38700, "b": 2.7, "Ea0": 1.0958665856e8, "w": 1.7505856e13}}
    _yaml = """
        type: Blowers-Masel
        rate-constant: {A: 38700, b: 2.7, Ea0: 1.0958665856e8, w: 1.7505856e13}
        """

    @classmethod
    def setUpClass(cls):
        ReactionRateTests.setUpClass()
        cls._parts = cls._input["rate-constant"]

    def eval(self, rate):
        delta_enthalpy = self.gas.delta_enthalpy[self._index]
        return rate(self.gas.T, delta_enthalpy)

    def test_from_parts(self):
        rate = self.from_parts()
        self.assertEqual(self._parts["A"], rate.pre_exponential_factor)
        self.assertEqual(self._parts["b"], rate.temperature_exponent)
        self.assertEqual(self._parts["Ea0"], rate.intrinsic_activation_energy)
        self.assertEqual(self._parts["w"], rate.bond_energy)
        self.check_rate(rate)

    def test_negative_A(self):
        # test reaction rate property
        rate = self.from_parts()
        self.assertFalse(rate.allow_negative_pre_exponential_factor)
        rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(rate.allow_negative_pre_exponential_factor)


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
        concm = self.gas.third_body_concentrations[self._index]
        return rate(self.gas.T, concm)

    def test_data(self):
        rate = self.from_parts()
        for n in self._n_data:
            rate.falloff_coeffs = np.random.rand(n)

    def test_third_body(self):
        concm = self.gas.third_body_concentrations
        self.assertIsFinite(concm[self._index])


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
        return rate(self.gas.T, self.gas.P)

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
        return rate(self.gas.T, self.gas.P)

    def test_from_parts(self):
        rate = self.from_parts()
        temperature_range = self._parts["temperature_range"]
        self.assertEqual(temperature_range[0], rate.temperature_range[0])
        self.assertEqual(temperature_range[1], rate.temperature_range[1])
        pressure_range = self._parts["pressure_range"]
        self.assertEqual(pressure_range[0], rate.pressure_range[0])
        self.assertEqual(pressure_range[1], rate.pressure_range[1])
        self.assertTrue(np.all(self._parts["data"] == rate.data))


class ReactionTests:
    # test suite for reaction expressions

    _cls = None # reaction object to be tested
    _equation = None # reaction equation string
    _rate = None # parameters for reaction rate object constructor
    _rate_obj = None # reaction rate object
    _kwargs = {} # additional parameters required by constructor
    _index = None # index of reaction in "kineticsfromscratch.yaml"
    _type = None # name of reaction rate
    _legacy = False # object uses legacy framework
    _yaml = None # YAML parameterization
    _deprecated_getters = {} # test legacy getters (old framework)
    _deprecated_setters = {} # test legacy setters (old framework)
    _deprecated_callers = {} # test legacy callers (old framework)

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution("kineticsfromscratch.yaml", transport_model=None)
        cls.species = cls.gas.species()

    def setUp(self):
        self.gas.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5"
        self.gas.TP = 900, 2*ct.one_atm

    def eval_rate(self, rate):
        # evaluate rate expression
        return rate(self.gas.T)

    def from_yaml(self, deprecated=False):
        # create reaction object from yaml
        if deprecated:
            with self.assertWarnsRegex(DeprecationWarning, "is renamed to 'from_yaml'"):
                return ct.Reaction.fromYaml(self._yaml, kinetics=self.gas)
        return ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)

    def from_dict(self):
        # create reaction rate object from input data
        input_data = self.from_yaml().input_data
        return ct.Reaction.from_dict(input_data, kinetics=self.gas)

    def from_rate(self, rate):
        # create reaction object from keywords / rate
        return self._cls(equation=self._equation, rate=rate, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)

    def from_parts(self):
        # create reaction rate object from parts
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products, legacy=self._legacy)
        rxn.rate = self._rate_obj
        return rxn

    def check_rate(self, rate_obj):
        if self._legacy:
            rate = rate_obj(self.gas.T)
        else:
            rate = self.eval_rate(rate_obj)
        self.assertNear(rate, self.gas.forward_rate_constants[self._index])

    def check_rxn(self, rxn, check_legacy=True):
        # helper function that checks reaction configuration
        ix = self._index
        self.assertEqual(rxn.reactants, self.gas.reaction(ix).reactants)
        self.assertEqual(rxn.products, self.gas.reaction(ix).products)
        if check_legacy:
            self.assertEqual(rxn.reaction_type, self._type)
            self.assertEqual(rxn.uses_legacy, self._type.endswith("-legacy"))
            self.assertEqual(rxn.uses_legacy, self._legacy)

        if not self._legacy:
            # legacy rate evaluation is not consistent
            self.check_rate(rxn.rate)
        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX
        self.check_solution(gas2, check_legacy)

    def check_solution(self, gas2, check_legacy=True):
        # helper function that checks evaluation of reaction rates
        ix = self._index
        if check_legacy:
            self.assertEqual(gas2.reaction_type_str(0), self._type)
        self.assertNear(gas2.forward_rate_constants[0],
                        self.gas.forward_rate_constants[ix])
        self.assertNear(gas2.net_rates_of_progress[0],
                        self.gas.net_rates_of_progress[ix])

    def test_rate(self):
        # check consistency of reaction rate and forward rate constant
        if self._rate_obj is not None:
            self.check_rate(self._rate_obj)

    def test_from_rate(self):
        # check instantiation from keywords / rate defined by dictionary
        self.check_rxn(self.from_rate(self._rate))

    def test_from_rate_obj(self):
        # check instantiation from keywords / rate provided as object
        if self._rate_obj is not None:
            self.check_rxn(self.from_rate(self._rate_obj))

    def test_from_parts(self):
        # check instantiation from parts (reactants, products, rate expression)
        if self._rate_obj is not None:
            self.check_rxn(self.from_parts())

    def test_from_yaml(self):
        # check constructors (from yaml input)
        if self._yaml is not None:
            # check instantiation from yaml string
            self.check_rxn(self.from_yaml())
            self.check_rxn(self.from_yaml(deprecated=True))

    def test_from_dict(self):
        # check instantiation from a yaml dictionary (input_data)
        if self._yaml is not None:
            # cannot compare types as input_data does not recreate legacy objects
            self.check_rxn(self.from_dict(), check_legacy=False)

    def test_add_rxn(self):
        # check adding new reaction to solution
        if self._rate_obj is None:
            return
        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[])
        gas2.TPX = self.gas.TPX
        rxn = self.from_rate(self._rate_obj)
        gas2.add_reaction(rxn)
        self.check_solution(gas2)

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
        rxn = self.from_rate(None)
        if self._legacy:
            self.assertNear(rxn.rate(self.gas.T), 0.)
        else:
            self.assertIsNaN(self.eval_rate(rxn.rate))

        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX
        if self._legacy:
            self.assertNear(gas2.forward_rate_constants[0], 0.)
            self.assertNear(gas2.net_rates_of_progress[0], 0.)
        elif not ct.debug_mode_enabled():
            self.assertIsNaN(gas2.forward_rate_constants[0])
            self.assertIsNaN(gas2.net_rates_of_progress[0])
        else:
            with self.assertRaisesRegex(ct.CanteraError, "not finite"):
                gas2.net_rates_of_progress

    def test_replace_rate(self):
        # check replacing reaction rate expression
        if self._rate_obj is None:
            return
        rxn = self.from_rate(None)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_roundtrip(self):
        # check round-trip instantiation via input_data
        if self._legacy:
            return
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
        else:
            self.assertNear(one, two)

    def test_deprecated_getters(self):
        # check property getters deprecated in new framework
        if self._yaml is None:
            return

        rxn = self.from_yaml()
        for attr, default in self._deprecated_getters.items():
            if self._legacy:
                self.check_equal(getattr(rxn, attr), default)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                    try:
                        self.check_equal(getattr(rxn, attr), default)
                    except Exception as err:
                        print(f"Exception raised when testing getter for '{attr}'")
                        raise err

    def test_deprecated_setters(self):
        # check property setters deprecated in new framework
        if self._yaml is None:
            return

        rxn = self.from_yaml()
        for attr, new in self._deprecated_setters.items():
            if self._legacy:
                setattr(rxn, attr, new)
                self.check_equal(getattr(rxn, attr), new)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                    setattr(rxn, attr, new)
                with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                    self.check_equal(getattr(rxn, attr), new)

    def test_deprecated_callers(self):
        # check methods deprecated in new framework
        if self._yaml is None:
            return

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
    _type = "elementary-legacy"
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
        rxn = self.from_rate(None)
        if self._legacy:
            rxn.rate = rate
        else:
            with self.assertWarnsRegex(DeprecationWarning, "'Arrhenius' object is deprecated"):
                rxn.rate = rate
        self.check_rxn(rxn)


class TestElementary(TestElementary2):
    # test updated version of elementary reaction

    _type = "elementary"
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
    _type = "three-body-legacy"
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

    def test_rate(self):
        # rate constant contains third-body concentration
        pass

    def test_efficiencies(self):
        # check efficiencies
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)

        self.assertEqual(rxn.efficiencies, self._kwargs["efficiencies"])


class TestThreeBody(TestThreeBody2):
    # test updated version of three-body reaction

    _legacy = False
    _type = "three-body"
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


class TestBlowersMasel(ReactionTests, utilities.CanteraTest):
    # test updated version of Blowers-Masel reaction

    _cls = ct.BlowersMaselReaction
    _type = "Blowers-Masel"
    _equation = "O + H2 <=> H + OH"
    _rate = {"A": 38700, "b": 2.7, "Ea0": 1.0958665856e8, "w": 1.7505856e13}
    _index = 6
    _yaml = """
        equation: O + H2 <=> H + OH
        type: Blowers-Masel
        rate-constant: {A: 38700, b: 2.7, Ea0: 2.619184e4 cal/mol, w: 4.184e9 cal/mol}
        """

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        cls._rate_obj = ct.BlowersMaselRate(**cls._rate)

    def eval_rate(self, rate):
        delta_enthalpy = self.gas.delta_enthalpy[self._index]
        return rate(self.gas.T, delta_enthalpy)


class TestTroe2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Troe falloff reaction

    _cls = ct.FalloffReaction
    _equation = "2 OH (+ M) <=> H2O2 (+ M)"
    _kwargs = {"efficiencies": {"AR": 0.7, "H2": 2.0, "H2O": 6.0}}
    _index = 2
    _type = "falloff-legacy"
    _legacy = True
    _yaml = """
        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
        type: falloff-legacy
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
        """

    def from_parts(self):
        rxn = ReactionTests.from_parts(self)
        rxn.efficiencies = self._kwargs["efficiencies"]
        return rxn

    def test_from_rate(self):
        # do not port creation from legacy Fallout objects
        pass

    def test_from_yaml(self):
        # check constructors (from yaml input)
        self.check_rxn(self.from_yaml(), check_legacy=False)


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
    _type = "falloff"
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
        concm = self.gas.third_body_concentrations[self._index]
        return rate(self.gas.T, concm)

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
    _type = "falloff"
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
        concm = self.gas.third_body_concentrations[self._index]
        return rate(self.gas.T, concm)

    def from_parts(self):
        rxn = ReactionTests.from_parts(self)
        rxn.efficiencies = self._kwargs["efficiencies"]
        return rxn


class TestChemicallyActivated2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Chemically Activated falloff reaction

    _cls = ct.ChemicallyActivatedReaction
    _equation = "H2O + OH (+M) <=> HO2 + H2 (+M)"
    _index = 10
    _type = "chemically-activated-legacy"
    _legacy = True
    _yaml = """
        equation: H2O + OH (+M) <=> HO2 + H2 (+M)  # Reaction 11
        units: {length: cm, quantity: mol, activation-energy: cal/mol}
        type: chemically-activated-legacy
        low-P-rate-constant: [282320.078, 1.46878, -3270.56495]
        high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
        """

    def test_from_rate(self):
        # do not port creation from legacy Fallout objects
        pass

    def test_from_yaml(self):
        # check constructors (from yaml input)
        self.check_rxn(self.from_yaml(), check_legacy=False)


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
    _type = "chemically-activated"
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
        concm = self.gas.third_body_concentrations[self._index]
        return rate(self.gas.T, concm)


class TestPlog2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Plog reaction

    _cls = ct.PlogReaction
    _equation = "H2 + O2 <=> 2 OH"
    _rate = [(1013.25, ct.Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
             (101325., ct.Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
             (1013250., ct.Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
             (10132500., ct.Arrhenius(5.9632e+56, -11.529, 220076726.4))]
    _index = 3
    _type = "pressure-dependent-Arrhenius-legacy"
    _legacy = True
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

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        if not cls._legacy:
            cls._rate_obj = ct.PlogRate(cls._rate)

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
                self.check_rates(rxn.rates, self._rate)

    def test_deprecated_setters(self):
        # overload default tester for deprecated property setters
        rate = ct.PlogRate(self._rate)
        rates = rate.rates

        rxn = self.from_yaml()
        if self._legacy:
            rxn.rates = rates
            self.check_rates(rxn.rates, self._rate)
        else:
            with self.assertWarnsRegex(DeprecationWarning, "Setter is replaceable"):
                rxn.rates = rates
            with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                self.check_rates(rxn.rates, self._rate)


class TestPlog(TestPlog2):
    # test updated version of Plog reaction

    _type = "pressure-dependent-Arrhenius"
    _legacy = False
    _yaml = """
        equation: H2 + O2 <=> 2 OH
        type: pressure-dependent-Arrhenius
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """

    def eval_rate(self, rate):
        return rate(self.gas.T, self.gas.P)


class TestChebyshev2(ReactionTests, utilities.CanteraTest):
    # test legacy version of Chebyshev reaction

    _cls = ct.ChebyshevReaction
    _equation = "HO2 <=> OH + O"
    _rate = {"temperature_range": (290., 3000.), "pressure_range": (1000., 10000000.0),
             "data": [[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
                      [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
                      [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]]}
    _index = 4
    _type = "Chebyshev-legacy"
    _legacy = True
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
        if not cls._legacy:
            cls._rate_obj = ct.ChebyshevRate(**cls._rate)
        cls._deprecated_getters.update({"coeffs": np.array(cls._rate["data"])})
        cls._deprecated_getters.update(
            {k: v for k, v in cls._rate.items()
                if k not in ["data", "temperature_range", "pressure_range"]})


class TestChebyshev(TestChebyshev2):
    # test updated version of Chebyshev reaction

    _type = "Chebyshev"
    _legacy = False
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
        return rate(self.gas.T, self.gas.P)


class TestCustom(ReactionTests, utilities.CanteraTest):
    # test Custom reaction

    # probe O + H2 <=> H + OH
    _cls = ct.CustomReaction
    _equation = "H2 + O <=> H + OH"
    _rate_obj = ct.CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
    _index = 0
    _type = "custom-rate-function"
    _legacy = False
    _yaml = None

    def setUp(self):
        # need to overwrite rate to ensure correct type ("method" is not compatible with Func1)
        super().setUp()
        self._rate = lambda T: 38.7 * T**2.7 * exp(-3150.15428/T)

    def test_roundtrip(self):
        # overload default tester for round trip
        pass

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
        self.assertNear(rate(self.gas.T), self.gas.forward_rate_constants[self._index])

    def test_custom_lambda(self):
        # check instantiation from keywords / rate provided as lambda function
        rxn = self.from_rate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
        self.check_rxn(rxn)
