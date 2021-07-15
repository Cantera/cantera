from math import exp
from pathlib import Path

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
    _uses_pressure = False # rate evaluation requires pressure
    _index = None # index of reaction in "kineticsfromscratch.yaml"
    _input = None # input parameters (dict corresponding to YAML)

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution("kineticsfromscratch.yaml")
        cls.gas.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5"
        cls.gas.TP = 900, 2*ct.one_atm

    def test_type(self):
        # check reaction type
        self.assertIn(self._type, "{}".format(self.rate))

    def test_rate_T(self):
        # check evaluation as a function of temperature only
        if self._uses_pressure:
            with self.assertRaisesRegex(ct.CanteraError, "reaction type requires pressure"):
                self.assertNear(self.rate(self.gas.T),
                                self.gas.forward_rate_constants[self._index])
        else:
            self.assertNear(self.rate(self.gas.T),
                            self.gas.forward_rate_constants[self._index])

    def test_rate_TP(self):
        # check evaluation as a function of temperature and pressure
        self.assertNear(self.rate(self.gas.T, self.gas.P),
                        self.gas.forward_rate_constants[self._index])

    def test_input(self):
        # check instantiation based on input_data
        rate = self._cls(input_data=self._input)
        self.assertIn(self._type, "{}".format(rate))
        self.assertNear(rate(self.gas.T, self.gas.P),
                        self.rate(self.gas.T, self.gas.P))

    def test_unconfigured(self):
        # check behavior of unconfigured rate object
        rate = self._cls(input_data={})
        self.assertTrue(np.isnan(rate(self.gas.T, self.gas.P)))
        input_data = rate.input_data
        self.assertIsInstance(input_data, dict)
        if input_data:
            self.assertEqual(input_data, {"type": self._type})

    def test_roundtrip(self):
        # check round-trip instantiation via input_data
        if self._index is None:
            return
        input_data = self.rate.input_data
        rate = self._cls(input_data=input_data)
        self.assertNear(rate(self.gas.T, self.gas.P),
                        self.rate(self.gas.T, self.gas.P))

    def test_from_dict(self):
        # check round-trip instantiation via input_data
        if self._index is None:
            return
        input_data = self.rate.input_data
        rate = ct.ReactionRate.from_dict(input_data, self.gas)
        self.assertNear(rate(self.gas.T, self.gas.P),
                        self.rate(self.gas.T, self.gas.P))


class TestArrheniusRate(ReactionRateTests, utilities.CanteraTest):
    # test Arrhenius rate expressions

    _cls = ct.ArrheniusRate
    _type = "ArrheniusRate"
    _uses_pressure = False
    _index = 0
    _input = {"rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}}

    def setUp(self):
        self.A = self.gas.reaction(self._index).rate.pre_exponential_factor
        self.b = self.gas.reaction(self._index).rate.temperature_exponent
        self.Ea = self.gas.reaction(self._index).rate.activation_energy
        self.rate = ct.ArrheniusRate(self.A, self.b, self.Ea)

    def test_parameters(self):
        # test reaction rate parameters
        self.assertEqual(self.A, self.rate.pre_exponential_factor)
        self.assertEqual(self.b, self.rate.temperature_exponent)
        self.assertEqual(self.Ea, self.rate.activation_energy)

    def test_allow_negative_pre_exponential_factor(self):
        # test reaction rate property
        self.assertFalse(self.rate.allow_negative_pre_exponential_factor)
        self.rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(self.rate.allow_negative_pre_exponential_factor)


class TestPlogRate(ReactionRateTests, utilities.CanteraTest):
    # test Plog rate expressions

    _cls = ct.PlogRate
    _type = "PlogRate"
    _uses_pressure = True
    _index = 3
    _input = {"rate-constants": [
        {"P": 1013.25, "A": 1.2124e+16, "b": -0.5779, "Ea": 45491376.8},
        {"P": 101325., "A": 4.9108e+31, "b": -4.8507, "Ea": 103649395.2},
        {"P": 1013250., "A": 1.2866e+47, "b": -9.0246, "Ea": 166508556.0},
        {"P": 10132500., "A": 5.9632e+56, "b": -11.529, "Ea": 220076726.4}]}

    def setUp(self):
        self.rate = ct.PlogRate([(1013.25, ct.Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
                                (101325., ct.Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
                                (1013250., ct.Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
                                (10132500., ct.Arrhenius(5.9632e+56, -11.529, 220076726.4))])

    def test_get_rates(self):
        # test getter for property rates
        rates = self.rate.rates
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


class TestChebyshevRate(ReactionRateTests, utilities.CanteraTest):
    # test Chebyshev rate expressions

    _cls = ct.ChebyshevRate
    _type = "ChebyshevRate"
    _uses_pressure = True
    _index = 4
    _input = {"data": [[8.2883, -1.1397, -0.12059, 0.016034],
                       [1.9764, 1.0037, 0.0072865, -0.030432],
                       [0.3177, 0.26889, 0.094806, -0.0076385]],
              "pressure-range": [1000.0, 10000000.0],
              "temperature-range": [290.0, 3000.0]}

    def setUp(self):
        self.Tmin = self.gas.reaction(self._index).rate.Tmin
        self.Tmax = self.gas.reaction(self._index).rate.Tmax
        self.Pmin = self.gas.reaction(self._index).rate.Pmin
        self.Pmax = self.gas.reaction(self._index).rate.Pmax
        self.coeffs = self.gas.reaction(self._index).rate.coeffs
        self.rate = ct.ChebyshevRate(self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.coeffs)

    def test_parameters(self):
        # test getters for rate properties
        self.assertEqual(self.Tmin, self.rate.Tmin)
        self.assertEqual(self.Tmax, self.rate.Tmax)
        self.assertEqual(self.Pmin, self.rate.Pmin)
        self.assertEqual(self.Pmax, self.rate.Pmax)
        self.assertTrue(np.all(self.coeffs == self.rate.coeffs))


class ReactionTests:
    # test suite for reaction expressions

    _cls = None # reaction object to be tested
    _type = None # name of reaction rate
    _legacy = False # object uses legacy framework
    _equation = None # reaction equation string
    _rate = None # parameters for reaction rate object constructor
    _rate_obj = None # reaction rate object
    _kwargs = {} # additional parameters required by contructor
    _index = None # index of reaction in "kineticsfromscratch.yaml"
    _yaml = None # YAML parameterization
    _input = None # input parameters (dict corresponding to YAML)
    _deprecated_getters = {} # test legacy getters (old framework)
    _deprecated_setters = {} # test legacy setters (old framework)
    _deprecated_callers = {} # test legacy callers (old framework)

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution("kineticsfromscratch.yaml", transport_model=None)
        cls.gas.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5"
        cls.gas.TP = 900, 2*ct.one_atm
        cls.species = cls.gas.species()

    def check_rxn(self, rxn, check_legacy=True):
        # helper function that checks reaction configuration
        ix = self._index
        self.assertEqual(rxn.reactants, self.gas.reaction(ix).reactants)
        self.assertEqual(rxn.products, self.gas.reaction(ix).products)
        if check_legacy:
            self.assertEqual(rxn.reaction_type, self._type)
            self.assertEqual(rxn.uses_legacy, self._type.endswith("-legacy"))
            self.assertEqual(rxn.uses_legacy, self._legacy)

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
        if self._rate_obj is None:
            return
        if self._legacy:
            self.assertNear(self._rate_obj(self.gas.T), self.gas.forward_rate_constants[self._index])
        else:
            self.assertNear(self._rate_obj(self.gas.T, self.gas.P),
                            self.gas.forward_rate_constants[self._index])

    def test_from_parts(self):
        # check instantiation from parts (reactants, products, rate expression)
        if not self._rate_obj:
            return
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products, legacy=self._legacy)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_from_dict1(self):
        # check instantiation from keywords / rate defined by dictionary
        rxn = self._cls(equation=self._equation, rate=self._rate, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        self.check_rxn(rxn)

    def test_from_yaml(self):
        # check instantiation from yaml string
        if self._yaml is None:
            return
        rxn = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_deprecated_yaml(self):
        # check instantiation from yaml string
        if self._yaml is None:
            return
        with self.assertWarnsRegex(DeprecationWarning, "is renamed to 'from_yaml'"):
            rxn = ct.Reaction.fromYaml(self._yaml, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_from_dict2(self):
        # check instantiation from a yaml dictionary
        if self._yaml is None:
            return
        rxn1 = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
        input_data = rxn1.input_data
        rxn2 = ct.Reaction.from_dict(input_data, kinetics=self.gas)
        # cannot compare types as input_data does not recreate legacy objects
        self.check_rxn(rxn2, check_legacy=False)

    def test_from_rate(self):
        # check instantiation from keywords / rate provided as object
        if self._rate_obj is None:
            return
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        self.check_rxn(rxn)

    def test_add_rxn(self):
        # check adding new reaction to solution
        if self._rate_obj is None:
            return
        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[])
        gas2.TPX = self.gas.TPX

        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        gas2.add_reaction(rxn)
        self.check_solution(gas2)

    def test_raises_invalid_rate(self):
        # check exception for instantiation from keywords / invalid rate
        with self.assertRaises(TypeError):
            rxn = self._cls(equation=self._equation, rate=(), kinetics=self.gas,
                            legacy=self._legacy, **self._kwargs)

    def test_no_rate(self):
        # check behavior for instantiation from keywords / no rate
        if self._rate_obj is None:
            return
        rxn = self._cls(equation=self._equation, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        if self._legacy:
            self.assertNear(rxn.rate(self.gas.T), 0.)
        else:
            self.assertTrue(np.isnan(rxn.rate(self.gas.T, self.gas.P)))

        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX

        self.assertNear(gas2.forward_rate_constants[0], 0.)
        self.assertNear(gas2.net_rates_of_progress[0], 0.)

    def test_replace_rate(self):
        # check replacing reaction rate expression
        if self._rate_obj is None:
            return
        rxn = self._cls(equation=self._equation, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_roundtrip(self):
        # check round-trip instantiation via input_data
        if self._legacy:
            return
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        input_data = rxn.rate.input_data
        rate_obj = rxn.rate.__class__(input_data=input_data)
        rxn2 = self._cls(equation=self._equation, rate=rate_obj, kinetics=self.gas,
                         legacy=self._legacy, **self._kwargs)
        self.check_rxn(rxn2)

    def check_equal(self, one, two):
        # helper function for deprecation tests
        self.assertEqual(type(one), type(two))
        if isinstance(one, (list, tuple, np.ndarray)):
            self.assertArrayNear(one, two)
        else:
            self.assertEqual(one, two)

    def test_deprecated_getters(self):
        # check property getters deprecated in new framework
        if self._yaml is None:
            return

        rxn = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
        for attr, default in self._deprecated_getters.items():
            if self._legacy:
                self.check_equal(getattr(rxn, attr), default)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                    self.check_equal(getattr(rxn, attr), default)

    def test_deprecated_setters(self):
        # check property setters deprecated in new framework
        if self._yaml is None:
            return

        rxn = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
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

        rxn = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
        for state, value in self._deprecated_callers.items():
            T, P = state
            if self._legacy:
                self.check_equal(rxn(T, P), value)
            else:
                with self.assertWarnsRegex(DeprecationWarning, "method is moved"):
                    self.check_equal(rxn(T, P), value)


class TestElementary(ReactionTests, utilities.CanteraTest):
    # test updated version of elementary reaction

    _cls = ct.ElementaryReaction
    _type = "elementary"
    _equation = "H2 + O <=> H + OH"
    _rate = {"A": 38.7, "b": 2.7, "Ea": 2.619184e+07}
    _index = 0
    _yaml = """
        equation: O + H2 <=> H + OH
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
        rxn = self._cls(equation=self._equation, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)
        if self._legacy:
            rxn.rate = rate
        else:
            with self.assertWarnsRegex(DeprecationWarning, "'Arrhenius' object is deprecated"):
                rxn.rate = rate
        self.check_rxn(rxn)


class TestElementary2(TestElementary):
    # test legacy version of elementary reaction

    _cls = ct.ElementaryReaction
    _type = "elementary-legacy"
    _legacy = True
    _yaml = """
        equation: O + H2 <=> H + OH
        type: elementary-legacy
        rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
        """


class TestThreeBody(TestElementary):
    # test updated version of three-body reaction

    _cls = ct.ThreeBodyReaction
    _type = "three-body"
    _equation = "2 O + M <=> O2 + M"
    _rate = {"A": 1.2e11, "b": -1.0, "Ea": 0.0}
    _kwargs = {"efficiencies": {"H2": 2.4, "H2O": 15.4, "AR": 0.83}}
    _index = 1
    _yaml = """
        equation: 2 O + M <=> O2 + M
        type: three-body
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """

    def test_from_parts(self):
        # overload default reaction creation from parts
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products, legacy=self._legacy)
        rxn.rate = self._rate_obj
        rxn.efficiencies = self._kwargs["efficiencies"]
        self.check_rxn(rxn)

    def test_rate(self):
        # rate constant contains third-body concentration
        pass

    def test_efficiencies(self):
        # check efficiencies
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas,
                        legacy=self._legacy, **self._kwargs)

        self.assertEqual(rxn.efficiencies, self._kwargs["efficiencies"])


class TestThreeBody2(TestThreeBody):
    # test legacy version of three-body reaction

    _cls = ct.ThreeBodyReaction
    _type = "three-body-legacy"
    _legacy = True
    _yaml = """
        equation: 2 O + M <=> O2 + M
        type: three-body-legacy
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """


class TestImplicitThreeBody(TestThreeBody):
    # test three-body reactions with explicit collision parther

    _cls = ct.ThreeBodyReaction
    _type = "three-body"
    _equation = "H + 2 O2 <=> HO2 + O2"
    _rate = {"A": 2.08e+19, "b": -1.24, "Ea": 0.0}
    _index = 5
    _yaml = """
        equation: H + 2 O2 <=> HO2 + O2
        rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
        """

    def test_efficiencies(self):
        # overload of default tester
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas,
                        legacy=self._legacy)
        self.assertEqual(rxn.efficiencies, {"O2": 1.})
        self.assertEqual(rxn.default_efficiency, 0.)

    def test_from_parts(self):
        # overload of default tester
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products)
        rxn.rate = self._rate_obj
        rxn.efficiencies = {"O2": 1.}
        rxn.default_efficiency = 0
        self.check_rxn(rxn)


class TestPlog(ReactionTests, utilities.CanteraTest):
    # test updated version of Plog reaction

    _cls = ct.PlogReaction
    _type = "pressure-dependent-Arrhenius"
    _legacy = False
    _equation = "H2 + O2 <=> 2 OH"
    _rate = [(1013.25, ct.Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
             (101325., ct.Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
             (1013250., ct.Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
             (10132500., ct.Arrhenius(5.9632e+56, -11.529, 220076726.4))]
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
        rxn = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
        if self._legacy:
            self.check_rates(rxn.rates, self._rate)
        else:
            with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                self.check_rates(rxn.rates, self._rate)

    def test_deprecated_setters(self):
        # overload default tester for deprecated property setters
        rate = ct.PlogRate(self._rate)
        rates = rate.rates

        rxn = ct.Reaction.from_yaml(self._yaml, kinetics=self.gas)
        if self._legacy:
            rxn.rates = rates
            self.check_rates(rxn.rates, self._rate)
        else:
            with self.assertWarnsRegex(DeprecationWarning, "Setter is replaceable"):
                rxn.rates = rates
            with self.assertWarnsRegex(DeprecationWarning, "property is moved"):
                self.check_rates(rxn.rates, self._rate)


class TestPlog2(TestPlog):
    # test legacy version of Plog reaction

    _cls = ct.PlogReaction
    _type = "pressure-dependent-Arrhenius-legacy"
    _legacy = True
    _rate_obj = None
    _yaml = """
        equation: H2 + O2 <=> 2 OH
        type: pressure-dependent-Arrhenius-legacy
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """


class TestChebyshev(ReactionTests, utilities.CanteraTest):
    # test updated version of Chebyshev reaction

    _cls = ct.ChebyshevReaction
    _type = "Chebyshev"
    _equation = "HO2 <=> OH + O"
    _rate = {"Tmin": 290., "Tmax": 3000., "Pmin": 1000., "Pmax": 10000000.0,
             "data": [[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
                      [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
                      [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]]}
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
    _deprecated_getters = {"nPressure": 4, "nTemperature": 3}
    _deprecated_callers = {(1000., ct.one_atm): 2858762454.1119065}

    @classmethod
    def setUpClass(cls):
        ReactionTests.setUpClass()
        if not cls._legacy:
            cls._rate_obj = ct.ChebyshevRate(**cls._rate)
        cls._deprecated_getters.update({"coeffs": np.array(cls._rate["data"])})
        cls._deprecated_getters.update(
            {k: v for k, v in cls._rate.items() if k != "data"})


class TestChebyshev2(TestChebyshev):
    # test legacy version of Chebyshev reaction

    _cls = ct.ChebyshevReaction
    _type = "Chebyshev-legacy"
    _legacy = True
    _rate_obj = None
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


class TestCustom(ReactionTests, utilities.CanteraTest):
    # test Custom reaction

    # probe O + H2 <=> H + OH
    _cls = ct.CustomReaction
    _type = "custom-rate-function"
    _legacy = False
    _equation = "H2 + O <=> H + OH"
    _rate_obj = ct.CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
    _index = 0
    _yaml = None

    def setUp(self):
        # need to overwrite rate to ensure correct type ("method" is not compatible with Func1)
        self._rate = lambda T: 38.7 * T**2.7 * exp(-3150.15428/T)

    def test_roundtrip(self):
        # overload default tester for round trip
        pass

    def test_from_func1(self):
        # check instantiation from keywords / rate provided as func1
        f = ct.Func1(self._rate)
        rxn = ct.CustomReaction(equation=self._equation, rate=f, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_rate_func(self):
        # check result of rate expression
        f = ct.Func1(self._rate)
        rate = ct.CustomRate(f)
        self.assertNear(rate(self.gas.T), self.gas.forward_rate_constants[self._index])

    def test_custom_lambda(self):
        # check instantiation from keywords / rate provided as lambda function
        rxn = ct.CustomReaction(equation=self._equation,
                                rate=lambda T: 38.7 * T**2.7 * exp(-3150.15428/T),
                                kinetics=self.gas)
        self.check_rxn(rxn)
