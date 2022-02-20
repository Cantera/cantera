import numpy as np
import pytest

import cantera as ct
from . import utilities
from .utilities import has_temperature_derivative_warnings


class RateExpressionTests:
    # Generic test class to check derivatives evaluated for a single reaction within
    # a reaction mechanism

    rxn_idx = None # index of reaction to be tested
    phase = None
    rtol = 1e-5
    orders = None
    ix3b = [] # three-body indices
    equation = None
    rate_type = None

    @classmethod
    def setUpClass(cls):
        ct.use_legacy_rate_constants(False)
        cls.tpx = cls.gas.TPX

        cls.r_stoich = cls.gas.reactant_stoich_coeffs3
        cls.p_stoich = cls.gas.product_stoich_coeffs3

        cls.rxn = cls.gas.reactions()[cls.rxn_idx]
        cls.rix = [cls.gas.species_index(k) for k in cls.rxn.reactants.keys()]
        cls.pix = [cls.gas.species_index(k) for k in cls.rxn.products.keys()]

    def setUp(self):
        self.gas.TPX = self.tpx
        self.gas.set_multiplier(0.)
        self.gas.set_multiplier(1., self.rxn_idx)
        self.gas.derivative_settings = {} # reset defaults

        # check stoichiometric coefficient output
        for k, v in self.rxn.reactants.items():
            ix = self.gas.species_index(k)
            self.assertEqual(self.r_stoich[ix, self.rxn_idx], v)
        for k, v in self.rxn.products.items():
            ix = self.gas.species_index(k)
            self.assertEqual(self.p_stoich[ix, self.rxn_idx], v)

    def test_input(self):
        # ensure that correct equation is referenced
        self.assertEqual(self.equation, self.rxn.equation)
        self.assertEqual(self.rate_type, self.rxn.rate.type)

    def rop_ddX(self, spc_ix, mode, const_t=True, rtol_deltac=1e-5, atol_deltac=1e-20):
        # numerical derivative for rates-of-progress with respect to mole fractions
        def calc():
            if mode == "forward":
                return self.gas.forward_rates_of_progress
            if mode == "reverse":
                return self.gas.reverse_rates_of_progress
            if mode == "net":
                return self.gas.net_rates_of_progress

        self.gas.TPX = self.tpx
        rop0 = calc()
        conc = self.gas.concentrations
        ctot0 = conc.sum()

        # perturb concentration
        dconc = conc[spc_ix] * rtol_deltac + atol_deltac
        conc[spc_ix] += dconc
        ctot1 = conc.sum()
        if const_t:
            # adjust pressure to compensate for concentration change
            pnew = self.gas.P * ctot1 / ctot0
            self.gas.TPX = self.gas.T, pnew, conc / ctot1
        else:
            # adjust temperature to compensate for concentration change
            tnew = self.gas.T * ctot1 / ctot0
            self.gas.TPX = tnew, self.gas.P, conc / ctot1
        drop = (calc() - rop0) / dconc
        self.gas.TPX = self.tpx
        return drop * self.gas.density_mole

    def test_forward_rop_ddX(self):
        # check derivatives of forward rates of progress with respect to mole fractions
        # against analytic result
        dropm = self.gas.forward_rates_of_progress_ddX
        dropp = self.gas.forward_rates_of_progress_ddP

        self.gas.derivative_settings = {"skip-third-bodies": True}
        drop = self.gas.forward_rates_of_progress_ddX
        rop = self.gas.forward_rates_of_progress
        for spc_ix in self.rix:
            if self.orders is None:
                order = self.r_stoich[spc_ix, self.rxn_idx]
            else:
                order = self.orders[self.gas.species_names[spc_ix]]
            self.assertNear(rop[self.rxn_idx],
                drop[self.rxn_idx, spc_ix] * self.gas.X[spc_ix] / order)

            drop_num = self.rop_ddX(spc_ix, mode="forward")
            self.assertArrayNear(dropm[:, spc_ix] + dropp * self.gas.P, drop_num, self.rtol)

        if isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.rix + self.ix3b):
            self.assertTrue(dropm[self.rxn_idx, spc_ix]) # non-zero
            dropm[self.rxn_idx, spc_ix] = 0
        self.assertFalse(dropm.any())

    def test_reverse_rop_ddX(self):
        # check derivatives of reverse rates of progress with respect to mole fractions
        # against analytic result
        dropm = self.gas.reverse_rates_of_progress_ddX
        dropp = self.gas.reverse_rates_of_progress_ddP

        self.gas.derivative_settings = {"skip-third-bodies": True}
        drop = self.gas.reverse_rates_of_progress_ddX
        rop = self.gas.reverse_rates_of_progress
        for spc_ix in self.pix:
            order = self.p_stoich[spc_ix, self.rxn_idx]
            self.assertNear(rop[self.rxn_idx],
                drop[self.rxn_idx, spc_ix] * self.gas.X[spc_ix] / order)

            drop_num = self.rop_ddX(spc_ix, mode="reverse")
            self.assertArrayNear(dropm[:, spc_ix] + dropp * self.gas.P, drop_num, self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.pix + self.ix3b):
            self.assertTrue(dropm[self.rxn_idx, spc_ix]) # non-zero
            dropm[self.rxn_idx, spc_ix] = 0
        self.assertFalse(dropm.any())

    def test_net_rop_ddX(self):
        # check derivatives of net rates of progress with respect to mole fractions
        # against numeric result
        drop = self.gas.net_rates_of_progress_ddX
        dropp = self.gas.net_rates_of_progress_ddP

        for spc_ix in self.rix + self.pix:
            drop_num = self.rop_ddX(spc_ix, mode="net")
            ix = drop[:, spc_ix] != 0
            drop_ = drop[:, spc_ix] + dropp * self.gas.P
            self.assertArrayNear(drop_[ix], drop_num[ix], self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.rix + self.pix + self.ix3b):
            self.assertTrue(drop[self.rxn_idx, spc_ix]) # non-zero
            drop[self.rxn_idx, spc_ix] = 0
        self.assertFalse(drop.any())

    def rop_ddT(self, mode=None, const_p=False, rtol=1e-6):
        # numerical derivative for rates-of-progress at constant pressure
        def calc():
            if mode == "forward":
                return self.gas.forward_rates_of_progress
            if mode == "reverse":
                return self.gas.reverse_rates_of_progress
            if mode == "net":
                return self.gas.net_rates_of_progress

        dt = self.tpx[0] * rtol
        dp = 0 if const_p else self.tpx[1] * rtol
        self.gas.TP = self.tpx[0] + dt, self.tpx[1] + dp
        rop1 = calc()
        self.gas.TP = self.tpx[:2]
        rop0 = calc()
        return (rop1[self.rxn_idx] - rop0[self.rxn_idx]) / dt

    def test_forward_rop_ddT(self):
        # check derivatives of forward rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.forward_rates_of_progress_ddT
        drop += self.gas.forward_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="forward", const_p=True)
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = self.gas.P / self.gas.T
        drop = self.gas.forward_rates_of_progress_ddT
        drop += self.gas.forward_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(mode="forward")
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

        if isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        self.assertTrue(drop[self.rxn_idx]) # non-zero
        drop[self.rxn_idx] = 0
        self.assertFalse(drop.any())

    def test_reverse_rop_ddT(self):
        # check derivatives of reverse rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.reverse_rates_of_progress_ddT
        drop += self.gas.reverse_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="reverse", const_p=True)
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = self.gas.P / self.gas.T
        drop = self.gas.reverse_rates_of_progress_ddT
        drop += self.gas.reverse_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(mode="reverse")
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        self.assertTrue(drop[self.rxn_idx]) # non-zero
        drop[self.rxn_idx] = 0
        self.assertFalse(drop.any())

    def test_net_rop_ddT(self):
        # check derivatives of net rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.net_rates_of_progress_ddT
        drop += self.gas.net_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="net", const_p=True)
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = self.gas.P / self.gas.T
        drop = self.gas.net_rates_of_progress_ddT
        drop += self.gas.net_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(mode="forward") - self.rop_ddT(mode="reverse")
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        self.assertTrue(drop[self.rxn_idx]) # non-zero
        drop[self.rxn_idx] = 0
        self.assertFalse(drop.any())

    def rop_ddP(self, mode=None, rtol=1e-6):
        # numerical derivative for rates-of-progress at constant pressure
        def calc():
            if mode == "forward":
                return self.gas.forward_rates_of_progress
            if mode == "reverse":
                return self.gas.reverse_rates_of_progress
            if mode == "net":
                return self.gas.net_rates_of_progress

        dp = self.tpx[1] * rtol
        self.gas.TP = self.tpx[0], self.tpx[1] + dp
        rop1 = calc()
        self.gas.TP = self.tpx[:2]
        rop0 = calc()
        return (rop1[self.rxn_idx] - rop0[self.rxn_idx]) / dp

    def test_forward_rop_ddP(self):
        # check derivatives of forward rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = self.gas.density_mole / self.gas.P
        drop = self.gas.forward_rates_of_progress_ddP
        drop += self.gas.forward_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(mode="forward")
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

    def test_reverse_rop_ddP(self):
        # check derivatives of reverse rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = self.gas.density_mole / self.gas.P
        drop = self.gas.reverse_rates_of_progress_ddP
        drop += self.gas.reverse_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(mode="reverse")
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

    def test_net_rop_ddP(self):
        # check derivatives of net rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = self.gas.density_mole / self.gas.P
        drop = self.gas.net_rates_of_progress_ddP
        drop += self.gas.net_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(mode="net")
        self.assertNear(drop[self.rxn_idx], drop_num, self.rtol)

    def rate_ddX(self, spc_ix, mode=None, const_t=True, rtol_deltac=1e-6, atol_deltac=1e-20):
        # numerical derivative for production rates with respect to mole fractions
        def calc(mode):
            if mode == "creation":
                return self.gas.creation_rates
            if mode == "destruction":
                return self.gas.destruction_rates
            if mode == "net":
                return self.gas.net_production_rates

        self.gas.TPX = self.tpx
        rate0 = calc(mode)
        conc = self.gas.concentrations
        ctot0 = conc.sum()

        # perturb concentration
        dconc = conc[spc_ix] * rtol_deltac + atol_deltac
        conc[spc_ix] += dconc
        ctot1 = conc.sum()
        if const_t:
            # adjust pressure to compensate for concentration change
            pnew = self.gas.P * ctot1 / ctot0
            self.gas.TPX = self.gas.T, pnew, conc / ctot1
        else:
            # adjust temperature to compensate for concentration change
            tnew = self.gas.T * ctot1 / ctot0
            self.gas.TPX = tnew, self.gas.P, conc / ctot1
        drate = (calc(mode) - rate0) / dconc
        self.gas.TPX = self.tpx
        return drate * self.gas.density_mole

    def test_creation_ddX(self):
        # check derivatives of creation rates with respect to mole fractions
        drate = self.gas.creation_rates_ddX
        dratep = self.gas.creation_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "creation")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            self.assertArrayNear(drate[ix, spc_ix], drate_num[ix], self.rtol)

    def test_destruction_ddX(self):
        # check derivatives of destruction rates with respect to mole fractions
        drate = self.gas.destruction_rates_ddX
        dratep = self.gas.destruction_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "destruction")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            self.assertArrayNear(drate[ix, spc_ix], drate_num[ix], self.rtol)

    def test_net_production_ddX(self):
        # check derivatives of destruction rates with respect to mole fractions
        drate = self.gas.net_production_rates_ddX
        dratep = self.gas.net_production_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "net")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            self.assertArrayNear(drate[ix, spc_ix], drate_num[ix], self.rtol)


class HydrogenOxygen(RateExpressionTests):

    @classmethod
    def setUpClass(cls):
        cls.gas = ct.Solution("h2o2.yaml", transport_model=None)
        #   species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2]
        cls.gas.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.3, 0.1]
        cls.gas.TP = 800, 2 * ct.one_atm
        super().setUpClass()


class TestElementaryRev(HydrogenOxygen, utilities.CanteraTest):
    # Standard elementary reaction with two reactants
    rxn_idx = 2
    equation = "H2 + O <=> H + OH"
    rate_type = "Arrhenius"


class TestElementarySelf(HydrogenOxygen, utilities.CanteraTest):
    # Elementary reaction with reactant reacting with itself
    rxn_idx = 27
    equation = "2 HO2 <=> H2O2 + O2"
    rate_type = "Arrhenius"


class TestFalloff(HydrogenOxygen, utilities.CanteraTest):
    # Fall-off reaction
    rxn_idx = 21
    equation = "2 OH (+M) <=> H2O2 (+M)"
    rate_type = "Troe"
    rtol = 1e-4


class TestThreeBody(HydrogenOxygen, utilities.CanteraTest):
    # Three body reaction with default efficiency
    rxn_idx = 1
    equation = "H + O + M <=> OH + M"
    rate_type = "Arrhenius"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ix3b = list(range(cls.gas.n_species))

    def test_thirdbodies_forward(self):
        drop = self.gas.forward_rates_of_progress_ddX
        self.gas.derivative_settings = {"skip-third-bodies": True}
        drops = self.gas.forward_rates_of_progress_ddX
        dropm = drop - drops
        rop = self.gas.forward_rates_of_progress
        self.assertNear(rop[self.rxn_idx], (dropm[self.rxn_idx] * self.gas.X).sum())

    def test_thirdbodies_reverse(self):
        drop = self.gas.reverse_rates_of_progress_ddX
        self.gas.derivative_settings = {"skip-third-bodies": True}
        drops = self.gas.reverse_rates_of_progress_ddX
        dropm = drop - drops
        rop = self.gas.reverse_rates_of_progress
        self.assertNear(rop[self.rxn_idx], (dropm[self.rxn_idx] * self.gas.X).sum())


class EdgeCases(RateExpressionTests):

    @classmethod
    def setUpClass(cls):
        cls.gas = ct.Solution("jacobian-tests.yaml", transport_model=None)
        #   species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR]
        cls.gas.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        cls.gas.TP = 800, 2 * ct.one_atm
        super().setUpClass()


class TestElementaryIrr(EdgeCases, utilities.CanteraTest):
    # Irreversible elementary reaction with two reactants
    rxn_idx = 0
    equation = "HO2 + O => O2 + OH"
    rate_type = "Arrhenius"


class TestElementaryOne(EdgeCases, utilities.CanteraTest):
    # Three-body reaction with single reactant species
    rxn_idx = 1
    equation = "H2 <=> 2 H"
    rate_type = "Arrhenius"


class TestElementaryThree(EdgeCases, utilities.CanteraTest):
    # Elementary reaction with three reactants
    rxn_idx = 2
    equation = "2 H + O <=> H2O"
    rate_type = "Arrhenius"


class TestElementaryFrac(EdgeCases, utilities.CanteraTest):
    # Elementary reaction with specified reaction order
    rxn_idx = 3
    orders = {"H2": 0.8, "O2": 1.0, "OH": 2.0}
    equation = "0.7 H2 + 0.2 O2 + 0.6 OH => H2O"
    rate_type = "Arrhenius"


class TestThreeBodyNoDefault(EdgeCases, utilities.CanteraTest):
    # Three body reaction without default efficiency
    rxn_idx = 4
    equation = "H + O + M <=> OH + M"
    rate_type = "Arrhenius"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        efficiencies = {"H2": 2.0, "H2O": 6.0, "AR": 0.7}
        cls.ix3b = [cls.gas.species_index(k) for k in efficiencies.keys()]


class FromScratchCases(RateExpressionTests):

    @classmethod
    def setUpClass(cls):
        cls.gas = ct.Solution("kineticsfromscratch.yaml", transport_model=None)
        #   species: [AR, O, H2, H, OH, O2, H2O, H2O2, HO2]
        cls.gas.X = [0.1, 3e-4, 5e-5, 6e-6, 3e-3, 0.6, 0.25, 1e-6, 2e-5]
        cls.gas.TP = 2000, 5 * ct.one_atm
        super().setUpClass()

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_forward_rop_ddT(self):
        super().test_forward_rop_ddT()

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_reverse_rop_ddT(self):
        super().test_reverse_rop_ddT()

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_net_rop_ddT(self):
        super().test_net_rop_ddT()


class TestPlog(FromScratchCases, utilities.CanteraTest):
    # Plog reaction
    rxn_idx = 3
    equation = "H2 + O2 <=> 2 OH"
    rate_type = "pressure-dependent-Arrhenius"


class TestChebyshev(FromScratchCases, utilities.CanteraTest):
    # Chebyshev reaction
    rxn_idx = 4
    equation = "HO2 <=> O + OH"
    rate_type = "Chebyshev"


class TestBlowersMasel(FromScratchCases, utilities.CanteraTest):
    # Blowers-Masel
    rxn_idx = 6
    equation = "H2 + O <=> H + OH"
    rate_type = "Blowers-Masel"

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    def test_forward_rop_ddT(self):
        super().test_forward_rop_ddT()

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    def test_reverse_rop_ddT(self):
        super().test_reverse_rop_ddT()

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    def test_net_rop_ddT(self):
        super().test_net_rop_ddT()


class FullTests:
    # Generic test class to check derivatives evaluated for an entire reaction mechanisms
    rtol = 1e-4

    @classmethod
    def setUpClass(cls):
        ct.use_legacy_rate_constants(False)
        cls.tpx = cls.gas.TPX

    def setUp(self):
        self.gas.TPX = self.tpx
        self.gas.derivative_settings = {} # reset

    def rop_ddX(self, mode, rtol_deltac=1e-9, atol_deltac=1e-20):
        # numerical derivative for rates-of-progress with respect to mole fractions
        def calc():
            if mode == "forward":
                return self.gas.forward_rates_of_progress
            if mode == "reverse":
                return self.gas.reverse_rates_of_progress
            if mode == "net":
                return self.gas.net_rates_of_progress

        n_spc, n_rxn = self.gas.n_species, self.gas.n_reactions
        drop = np.zeros((n_rxn, n_spc))

        self.gas.TPX = self.tpx
        rop0 = calc()
        ctot0 = self.gas.density_mole
        ctot0 = self.gas.concentrations.sum()

        for spc_ix in range(n_spc):
            conc = self.gas.concentrations
            dconc = conc[spc_ix] * rtol_deltac + atol_deltac
            conc[spc_ix] += dconc
            ctot1 = conc.sum()
            self.gas.TPX = self.tpx[0], self.tpx[1] * ctot1 / ctot0, conc / ctot1
            drop[:, spc_ix] = (calc() - rop0) / dconc
            self.gas.TPX = self.tpx

        return drop * self.gas.density_mole

    def test_forward_rop_ddX(self):
        # check forward rop against numerical derivative with respect to mole fractions
        drop = self.gas.forward_rates_of_progress_ddX
        dropp = self.gas.forward_rates_of_progress_ddP
        drop_num = self.rop_ddX(mode="forward")
        stoich = self.gas.reactant_stoich_coeffs3
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                self.assertArrayNear(drop_, drop_num[i, ix], self.rtol)
            except AssertionError as err:
                print(i, self.gas.reaction(i).rate.type)
                print(self.gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_reverse_rop_ddX(self):
        # check reverse rop against numerical derivative with respect to mole fractions
        drop = self.gas.reverse_rates_of_progress_ddX
        dropp = self.gas.reverse_rates_of_progress_ddP
        drop_num = self.rop_ddX(mode="reverse")
        stoich = self.gas.product_stoich_coeffs3
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                self.assertArrayNear(drop_, drop_num[i, ix], self.rtol)
            except AssertionError as err:
                print(i, self.gas.reaction(i).rate.type)
                print(self.gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_net_rop_ddX(self):
        # check net rop against numerical derivative with respect to mole fractions
        drop = self.gas.net_rates_of_progress_ddX
        dropp = self.gas.net_rates_of_progress_ddP
        drop_num = self.rop_ddX(mode="net")
        stoich = self.gas.product_stoich_coeffs3 - self.gas.reactant_stoich_coeffs3
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                self.assertArrayNear(drop_, drop_num[i, ix], self.rtol)
            except AssertionError as err:
                if self.gas.reaction(i).reversible:
                    print(i, self.gas.reaction(i).rate.type)
                    print(self.gas.reaction(i))
                    print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                    raise err

    def rop_ddT(self, mode=None, dt=1e-6):
        # numerical derivative for rates-of-progress at constant pressure
        def calc():
            if mode == "forward":
                return self.gas.forward_rates_of_progress
            if mode == "reverse":
                return self.gas.reverse_rates_of_progress
            if mode == "net":
                return self.gas.net_rates_of_progress
            return None

        self.gas.TP = self.tpx[0] + dt, self.tpx[1]
        rop1 = calc()
        self.gas.TP = self.tpx[:2]
        rop0 = calc()
        return (rop1 - rop0) / dt

    def test_forward_rop_ddT(self):
        # check forward rop against numerical derivative with respect to temperature
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.forward_rates_of_progress_ddT
        drop += self.gas.forward_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="forward")
        self.assertArrayNear(drop, drop_num, self.rtol)

    def test_reverse_rop_ddT(self):
        # check reverse rop against numerical derivative with respect to temperature
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.reverse_rates_of_progress_ddT
        drop += self.gas.reverse_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="reverse")
        self.assertArrayNear(drop, drop_num, self.rtol)

    def test_net_rop_ddT(self):
        # check net rop against numerical derivative with respect to temperature
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.net_rates_of_progress_ddT
        drop += self.gas.net_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="net")
        try:
            self.assertArrayNear(drop, drop_num, self.rtol)
        except AssertionError as err:
            i = np.argmax(2 * (drop - drop_num) / (drop + drop_num + 2e-4))
            print(i, self.gas.reaction(i).rate.type)
            print(self.gas.reaction(i))
            print(drop[i])
            print(drop_num[i])
            raise err


class FullHydrogenOxygen(FullTests, utilities.CanteraTest):

    @classmethod
    def setUpClass(cls):
        cls.gas = ct.Solution("h2o2.yaml", transport_model=None)
        cls.gas.TPX = 300, 5 * ct.one_atm, "H2:1, O2:3"
        cls.gas.equilibrate("HP")
        super().setUpClass()


class FullGriMech(FullTests, utilities.CanteraTest):

    @classmethod
    def setUpClass(cls):
        cls.gas = ct.Solution("gri30.yaml", transport_model=None)
        cls.gas.TPX = 300, 5 * ct.one_atm, "CH4:1, C3H8:.1, O2:1, N2:3.76"
        cls.gas.equilibrate("HP")
        super().setUpClass()


class FullEdgeCases(FullTests, utilities.CanteraTest):

    @classmethod
    def setUpClass(cls):
        cls.gas = ct.Solution("jacobian-tests.yaml", transport_model=None)
        #   species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR]
        cls.gas.TPX = 300, 2 * ct.one_atm, "H2:1, O2:3, AR:0.4"
        cls.gas.equilibrate("HP")
        super().setUpClass()
