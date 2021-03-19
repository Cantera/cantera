from math import exp

import cantera as ct
from . import utilities


class TestElementary(utilities.CanteraTest):

    _cls = ct.ElementaryReaction
    _equation = 'H2 + O <=> H + OH'
    _rate = {'A': 38.7, 'b': 2.7, 'Ea': 2.619184e+07}
    _rate_obj = ct.Arrhenius(38.7, 2.7, 2.619184e+07)
    _index = 2
    _type = "elementary"

    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('h2o2.xml')
        cls.gas.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'
        cls.gas.TP = 900, 2*ct.one_atm
        cls.species = ct.Species.listFromFile('h2o2.xml')

    def check_rxn(self, rxn):
        ix = self._index
        self.assertEqual(rxn.reaction_type, self._type)
        self.assertNear(rxn.rate(self.gas.T), self.gas.forward_rate_constants[ix])
        self.assertEqual(rxn.reactants, self.gas.reaction(ix).reactants)
        self.assertEqual(rxn.products, self.gas.reaction(ix).products)

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX
        self.check_sol(gas2)

    def check_sol(self, gas2):
        ix = self._index
        self.assertEqual(gas2.reaction_type_str(0), self._type)
        self.assertNear(gas2.forward_rate_constants[0],
                        self.gas.forward_rate_constants[ix])
        self.assertNear(gas2.net_rates_of_progress[0],
                        self.gas.net_rates_of_progress[ix])

    def test_rate(self):
        self.assertNear(self._rate_obj(self.gas.T), self.gas.forward_rate_constants[self._index])

    def test_from_parts(self):
        orig = self.gas.reaction(self._index)
        rxn = self._cls(orig.reactants, orig.products)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)

    def test_from_dict(self):
        rxn = self._cls(equation=self._equation, rate=self._rate, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_from_rate(self):
        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas)
        self.check_rxn(rxn)

    def test_add_rxn(self):
        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[])
        gas2.TPX = self.gas.TPX

        rxn = self._cls(equation=self._equation, rate=self._rate_obj, kinetics=self.gas)
        gas2.add_reaction(rxn)
        self.check_sol(gas2)

    def test_wrong_rate(self):
        with self.assertRaises(TypeError):
            rxn = self._cls(equation=self._equation, rate=[], kinetics=self.gas)

    def test_no_rate(self):
        rxn = self._cls(equation=self._equation, kinetics=self.gas)
        self.assertNear(rxn.rate(self.gas.T), 0.)

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX

        self.assertNear(gas2.forward_rate_constants[0], 0.)
        self.assertNear(gas2.net_rates_of_progress[0], 0.)

    def test_replace_rate(self):
        rxn = self._cls(equation=self._equation, kinetics=self.gas)
        rxn.rate = self._rate_obj
        self.check_rxn(rxn)


class TestCustom(TestElementary):

    # probe O + H2 <=> H + OH
    _cls = ct.CustomReaction
    _equation = 'H2 + O <=> H + OH'
    _rate_obj = ct.CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15428/T))
    _index = 2
    _type = "custom-rate-function"

    def setUp(self):
        # need to overwrite rate to ensure correct type ('method' is not compatible with Func1)
        self._rate = lambda T: 38.7 * T**2.7 * exp(-3150.15428/T)

    def test_no_rate(self):
        rxn = self._cls(equation=self._equation, kinetics=self.gas)
        with self.assertRaisesRegex(ct.CanteraError, "Custom rate function is not initialized."):
            rxn.rate(self.gas.T)

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[rxn])
        gas2.TPX = self.gas.TPX

        with self.assertRaisesRegex(ct.CanteraError, "Custom rate function is not initialized."):
            gas2.forward_rate_constants

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


class TestElementaryNew(TestElementary):

    _cls = ct.TestReaction
    _equation = 'H2 + O <=> H + OH'
    _rate = {'A': 38.7, 'b': 2.7, 'Ea': 2.619184e+07}
    _rate_obj = ct.ArrheniusRate(38.7, 2.7, 2.619184e+07)
    _index = 2
    _type = "elementary-new"
