from math import exp
from pathlib import Path

import cantera as ct
import numpy as np
from . import utilities


class TestImplicitThirdBody(utilities.CanteraTest):

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution("gri30.yaml")

    def test_implicit_three_body(self):
        yaml1 = """
            equation: H + 2 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            """
        rxn1 = ct.Reaction.fromYaml(yaml1, self.gas)
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
        rxn2 = ct.Reaction.fromYaml(yaml2, self.gas)
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
        rxn1 = ct.Reaction.fromYaml(yaml1, gas1)

        yaml2 = """
            equation: H + O2 + M <=> HO2 + M
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            type: three-body
            default-efficiency: 0
            efficiencies: {H2O: 1}
            """
        rxn2 = ct.Reaction.fromYaml(yaml2, gas1)

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
        yaml = """
            equation: H + O2 + H2O <=> HO2 + H2O
            rate-constant: {A: 1.126e+19, b: -0.76, Ea: 0.0}
            """
        rxn = ct.Reaction.fromYaml(yaml, self.gas)
        input_data = rxn.input_data

        self.assertNotIn("type", input_data)
        self.assertNotIn("default-efficiency", input_data)
        self.assertNotIn("efficiencies", input_data)

    def test_non_integer_stoich(self):
        yaml = """
            equation: H + 1.5 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            """
        rxn = ct.Reaction.fromYaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "elementary")

    def test_not_three_body(self):
        yaml = """
            equation: HCNO + H <=> H + HNCO  # Reaction 270
            rate-constant: {A: 2.1e+15, b: -0.69, Ea: 2850.0}
            """
        rxn = ct.Reaction.fromYaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "elementary")

    def test_user_override(self):
        yaml = """
            equation: H + 2 O2 <=> HO2 + O2
            rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
            type: elementary
            """
        rxn = ct.Reaction.fromYaml(yaml, self.gas)
        self.assertEqual(rxn.reaction_type, "elementary")


class TestReactionRate(utilities.CanteraTest):

    _cls = None
    _type = None
    _uses_pressure = False
    _index = None
    _input = None
    _cls = None

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution("kineticsfromscratch.yaml")
        cls.gas.X = "H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5"
        cls.gas.TP = 900, 2*ct.one_atm

    def test_type(self):
        if self._index is None:
            return
        self.assertIn(self._type, "{}".format(self.rate))

    def test_rate_T(self):
        if self._index is None:
            return
        if self._uses_pressure:
            with self.assertRaisesRegex(ct.CanteraError, "reaction type requires pressure"):
                self.assertNear(self.rate(self.gas.T),
                                self.gas.forward_rate_constants[self._index])
        else:
            self.assertNear(self.rate(self.gas.T),
                            self.gas.forward_rate_constants[self._index])

    def test_rate_TP(self):
        if self._index is None:
            return
        self.assertNear(self.rate(self.gas.T, self.gas.P),
                        self.gas.forward_rate_constants[self._index])

    def test_input(self):
        if self._input is None or self._cls is None:
            return
        rate = self._cls(input_data=self._input)
        self.assertIn(self._type, "{}".format(rate))
        self.assertNear(rate(self.gas.T, self.gas.P),
                        self.rate(self.gas.T, self.gas.P))

    def test_roundtrip(self):
        if self._index is None:
            return
        input_data = self.rate.input_data
        rate = self._cls(input_data=input_data)
        self.assertNear(rate(self.gas.T, self.gas.P),
                        self.rate(self.gas.T, self.gas.P))


class TestArrheniusRate(TestReactionRate):

    _type = "Arrhenius"
    _uses_pressure = False
    _index = 0
    _input = {"rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}}
    _cls = ct.ArrheniusRate

    def setUp(self):
        self.A = self.gas.reaction(self._index).rate.pre_exponential_factor
        self.b = self.gas.reaction(self._index).rate.temperature_exponent
        self.Ea = self.gas.reaction(self._index).rate.activation_energy
        self.rate = ct.ArrheniusRate(self.A, self.b, self.Ea)

    def test_parameters(self):
        self.assertEqual(self.A, self.rate.pre_exponential_factor)
        self.assertEqual(self.b, self.rate.temperature_exponent)
        self.assertEqual(self.Ea, self.rate.activation_energy)

    def test_allow_negative_pre_exponential_factor1(self):
        self.assertFalse(self.rate.allow_negative_pre_exponential_factor)
        self.rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(self.rate.allow_negative_pre_exponential_factor)

    def test_allow_negative_pre_exponential_factor2(self):
        # modify value in memory
        self.assertFalse(self.gas.reaction(self._index).rate.allow_negative_pre_exponential_factor)
        self.gas.reaction(self._index).rate.allow_negative_pre_exponential_factor = True
        self.assertTrue(self.gas.reaction(self._index).rate.allow_negative_pre_exponential_factor)


class TestPlogRate(TestReactionRate):

    _type = "Plog"
    _uses_pressure = True
    _index = 3
    _input = {"rate-constants": [
        {"P": 1013.25, "A": 1.2124e+16, "b": -0.5779, "Ea": 45491376.8},
        {"P": 101325., "A": 4.9108e+31, "b": -4.8507, "Ea": 103649395.2},
        {"P": 1013250., "A": 1.2866e+47, "b": -9.0246, "Ea": 166508556.0},
        {"P": 10132500., "A": 5.9632e+56, "b": -11.529, "Ea": 220076726.4}]}
    _cls = ct.PlogRate

    def setUp(self):
        self.rate = ct.PlogRate([(1013.25, ct.Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
                                (101325., ct.Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
                                (1013250., ct.Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
                                (10132500., ct.Arrhenius(5.9632e+56, -11.529, 220076726.4))])


class TestChebyshevRate(TestReactionRate):

    _type = "Chebyshev"
    _uses_pressure = True
    _index = 4
    _input = {"data": [[8.2883, -1.1397, -0.12059, 0.016034],
                       [1.9764, 1.0037, 0.0072865, -0.030432],
                       [0.3177, 0.26889, 0.094806, -0.0076385]],
              "pressure-range": [1000.0, 10000000.0],
              "temperature-range": [290.0, 3000.0]}
    _cls = ct.ChebyshevRate

    def setUp(self):
        self.Tmin = self.gas.reaction(self._index).rate.Tmin
        self.Tmax = self.gas.reaction(self._index).rate.Tmax
        self.Pmin = self.gas.reaction(self._index).rate.Pmin
        self.Pmax = self.gas.reaction(self._index).rate.Pmax
        self.coeffs = self.gas.reaction(self._index).rate.coeffs
        self.rate = ct.ChebyshevRate(self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.coeffs)

    def test_parameters(self):
        self.assertEqual(self.Tmin, self.rate.Tmin)
        self.assertEqual(self.Tmax, self.rate.Tmax)
        self.assertEqual(self.Pmin, self.rate.Pmin)
        self.assertEqual(self.Pmax, self.rate.Pmax)
        self.assertTrue(np.all(self.coeffs == self.rate.coeffs))


class TestReaction(utilities.CanteraTest):

    _cls = None
    _equation = None
    _rate = None
    _rate_obj = None
    _kwargs = {}
    _index = None
    _type = None
    _yaml = None
    _input = None

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('kineticsfromscratch.yaml', transport_model=None)
        cls.gas.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'
        cls.gas.TP = 900, 2*ct.one_atm
        cls.species = cls.gas.species()

    def check_rxn(self, rxn):
        ix = self._index
        self.assertEqual(rxn.reaction_type, self._type)
        self.assertEqual(rxn.reactants, self.gas.reaction(ix).reactants)
        self.assertEqual(rxn.products, self.gas.reaction(ix).products)

        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX
        self.check_solution(gas2)

    def check_solution(self, gas2):
        ix = self._index
        self.assertEqual(gas2.reaction_type_str(0), self._type)
        self.assertNear(gas2.forward_rate_constants[0],
                        self.gas.forward_rate_constants[ix])
        self.assertNear(gas2.net_rates_of_progress[0],
                        self.gas.net_rates_of_progress[ix])

    def test_rate(self):
        if self._rate_obj is None:
            return
        if "Rate" in type(self._rate_obj).__name__:
            self.assertNear(self._rate_obj(self.gas.T, self.gas.P),
                            self.gas.forward_rate_constants[self._index])
        else:
            self.assertNear(self._rate_obj(self.gas.T), self.gas.forward_rate_constants[self._index])

    def test_from_parts(self):
        if self._cls is None or not hasattr(self._cls, "rate"):
            return
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_from_dict(self):
        if self._cls is None:
            return
        rxn = self._cls(equation=self._equation, rate=self._rate, kinetics=self.gas, **self._kwargs)
        self.check_rxn(rxn)

    def test_from_yaml(self):
        if self._yaml is None:
            return
        rxn = ct.Reaction.fromYaml(self._yaml, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_from_rate(self):
        if self._cls is None or self._rate_obj is None:
            return
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas, **self._kwargs)
        self.check_rxn(rxn)

    def test_add_rxn(self):
        if self._cls is None or self._rate_obj is None:
            return
        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[])
        gas2.TPX = self.gas.TPX

        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas, **self._kwargs)
        gas2.add_reaction(rxn)
        self.check_solution(gas2)

    def test_wrong_rate(self):
        if self._cls is None:
            return
        with self.assertRaises(TypeError):
            rxn = self._cls(equation=self._equation, rate=(), kinetics=self.gas, **self._kwargs)

    def test_no_rate(self):
        if self._cls is None or not hasattr(self._cls, "rate"):
            return
        rxn = self._cls(equation=self._equation, kinetics=self.gas, **self._kwargs)
        if "Rate" in type(self._rate_obj).__name__:
            self.assertNear(rxn.rate(self.gas.T, self.gas.P), 0.)
        else:
            self.assertNear(rxn.rate(self.gas.T), 0.)

        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX

        self.assertNear(gas2.forward_rate_constants[0], 0.)
        self.assertNear(gas2.net_rates_of_progress[0], 0.)

    def test_replace_rate(self):
        if self._cls is None or self._rate_obj is None:
            return
        rxn = self._cls(equation=self._equation, kinetics=self.gas, **self._kwargs)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_roundtrip(self):
        if self._cls is None or self._type.endswith("-old"):
            return
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas, **self._kwargs)
        input_data = rxn.rate.input_data
        rate_obj = rxn.rate.__class__(input_data=input_data)
        rxn2 = self._cls(equation=self._equation, rate=rate_obj, kinetics=self.gas, **self._kwargs)
        self.check_rxn(rxn2)


class TestElementary(TestReaction):

    _cls = ct.ElementaryReaction
    _equation = "H2 + O <=> H + OH"
    _rate = {"A": 38.7, "b": 2.7, "Ea": 2.619184e+07}
    _rate_obj = ct.Arrhenius(38.7, 2.7, 2.619184e+07)
    _kwargs = {}
    _index = 0
    _type = "elementary-old"
    _yaml = """
        equation: O + H2 <=> H + OH
        type: elementary-old
        rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
        """


class TestElementary3(TestElementary):

    _cls = ct.ElementaryReaction3
    _rate_obj = ct.ArrheniusRate(38.7, 2.7, 2.619184e+07)
    _type = "elementary"
    _yaml = """
        equation: O + H2 <=> H + OH
        rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
        """


class TestThreeBody(TestReaction):

    _cls = ct.ThreeBodyReaction
    _equation = "2 O + M <=> O2 + M"
    _rate = {"A": 1.2e11, "b": -1.0, "Ea": 0.0}
    _rate_obj = ct.Arrhenius(1.2e11, -1., 0.)
    _kwargs = {"efficiencies": {"H2": 2.4, "H2O": 15.4, "AR": 0.83}}
    _index = 1
    _type = "three-body-old"
    _yaml = """
        equation: 2 O + M <=> O2 + M
        type: three-body-old
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """

    def test_from_parts(self):
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products)
        rxn.rate = self._rate_obj
        rxn.efficiencies = self._kwargs["efficiencies"]
        self.check_rxn(rxn)

    def test_rate(self):
        # rate constant contains third-body concentration
        pass

    def test_efficiencies(self):
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas, **self._kwargs)

        self.assertEqual(rxn.efficiencies, self._kwargs["efficiencies"])


class TestThreeBody3(TestThreeBody):

    _cls = ct.ThreeBodyReaction3
    _rate_obj = ct.ArrheniusRate(1.2e11, -1., 0.)
    _type = "three-body"
    _yaml = """
        equation: 2 O + M <=> O2 + M
        type: three-body
        rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        """


class TestImplicitThreeBody3(TestThreeBody):

    _cls = ct.ThreeBodyReaction3
    _equation = "H + 2 O2 <=> HO2 + O2"
    _rate = {"A": 2.08e+19, "b": -1.24, "Ea": 0.0}
    _rate_obj = ct.ArrheniusRate(2.08e+19, -1.24, 0.)
    _index = 5
    _type = "three-body"
    _yaml = """
        equation: H + 2 O2 <=> HO2 + O2
        rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
        """

    def test_efficiencies(self):
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas)
        self.assertEqual(rxn.efficiencies, {"O2": 1.})
        self.assertEqual(rxn.default_efficiency, 0.)

    def test_from_parts(self):
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products)
        rxn.rate = self._rate_obj
        rxn.efficiencies = {"O2": 1.}
        rxn.default_efficiency = 0
        self.check_rxn(rxn)


class TestPlog(TestReaction):

    _cls = ct.PlogReaction
    _equation = "H2 + O2 <=> 2 OH"
    _rate = [{"P": 1013.25, "rate-constant": {"A": 1.2124e+16, "b": -0.5779, "Ea": 45491376.8}},
             {"P": 101325., "rate-constant": {"A": 4.9108e+31, "b": -4.8507, "Ea": 103649395.2}},
             {"P": 1013250., "rate-constant": {"A": 1.2866e+47, "b": -9.0246, "Ea": 166508556.0}},
             {"P": 10132500., "rate-constant": {"A": 5.9632e+56, "b": -11.529, "Ea": 220076726.4}}]
    _type = "pressure-dependent-Arrhenius-old"
    _index = 3
    _yaml = """
        equation: H2 + O2 <=> 2 OH  # Reaction 4
        type: pressure-dependent-Arrhenius-old
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """


class TestPlog3(TestPlog):

    _cls = ct.PlogReaction3
    _rate_obj = ct.PlogRate([(1013.25, ct.Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
                             (101325., ct.Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
                             (1013250., ct.Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
                             (10132500., ct.Arrhenius(5.9632e+56, -11.529, 220076726.4))])
    _type = "pressure-dependent-Arrhenius"
    _yaml = """
        equation: H2 + O2 <=> 2 OH  # Reaction 4
        type: pressure-dependent-Arrhenius
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}
        """


class TestChebyshev(TestReaction):

    _cls = ct.ChebyshevReaction
    _equation = "HO2 <=> OH + O"
    _rate = {"Tmin": 290., "Tmax": 3000., "Pmin": 1000., "Pmax": 10000000.0,
             "data": [[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
                      [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
                      [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]]}
    _type = "Chebyshev-old"
    _index = 4
    _yaml = """
        equation: HO2 <=> OH + O  # Reaction 5
        type: Chebyshev-old
        temperature-range: [290.0, 3000.0]
        pressure-range: [9.869232667160128e-03 atm, 98.69232667160128 atm]
        data:
        - [8.2883, -1.1397, -0.12059, 0.016034]
        - [1.9764, 1.0037, 7.2865e-03, -0.030432]
        - [0.3177, 0.26889, 0.094806, -7.6385e-03]
        """


class TestChebyshev3(TestChebyshev):

    _cls = ct.ChebyshevReaction3
    _rate_obj = ct.ChebyshevRate(Tmin=290., Tmax=3000., Pmin=1000., Pmax=10000000.0,
                                 data=[[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
                                       [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
                                       [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]])
    _type = "Chebyshev"
    _yaml = """
        equation: HO2 <=> OH + O  # Reaction 5
        type: Chebyshev
        temperature-range: [290.0, 3000.0]
        pressure-range: [9.869232667160128e-03 atm, 98.69232667160128 atm]
        data:
        - [8.2883, -1.1397, -0.12059, 0.016034]
        - [1.9764, 1.0037, 7.2865e-03, -0.030432]
        - [0.3177, 0.26889, 0.094806, -7.6385e-03]
        """


class TestCustom(TestReaction):

    # probe O + H2 <=> H + OH
    _cls = ct.CustomReaction
    _equation = "H2 + O <=> H + OH"
    _rate_obj = ct.CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
    _index = 0
    _type = "custom-rate-function"
    _yaml = None

    def setUp(self):
        # need to overwrite rate to ensure correct type ("method" is not compatible with Func1)
        self._rate = lambda T: 38.7 * T**2.7 * exp(-3150.15428/T)

    def test_no_rate(self):
        rxn = self._cls(equation=self._equation, kinetics=self.gas)
        with self.assertRaisesRegex(ct.CanteraError, "Custom rate function is not initialized."):
            rxn.rate(self.gas.T)

        gas2 = ct.Solution(thermo="IdealGas", kinetics="GasKinetics",
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX

        with self.assertRaisesRegex(ct.CanteraError, "Custom rate function is not initialized."):
            gas2.forward_rate_constants

    def test_roundtrip(self):
        pass

    def test_from_func(self):
        f = ct.Func1(self._rate)
        rxn = ct.CustomReaction(equation=self._equation, rate=f, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_rate_func(self):
        f = ct.Func1(self._rate)
        rate = ct.CustomRate(f)
        self.assertNear(rate(self.gas.T), self.gas.forward_rate_constants[self._index])

    def test_custom(self):
        rxn = ct.CustomReaction(equation=self._equation,
                                rate=lambda T: 38.7 * T**2.7 * exp(-3150.15428/T),
                                kinetics=self.gas)
        self.check_rxn(rxn)

