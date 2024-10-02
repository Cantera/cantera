import numpy as np
import pytest
from pytest import approx

import cantera as ct


@pytest.fixture(scope='class')
def setup_rate_expression_tests(request):
    # Sets the TPX state for the gas object and the stoichiometric coefficients,
    # once per class.
    request.cls.tpx = request.cls.gas.TPX

    request.cls.r_stoich = request.cls.gas.reactant_stoich_coeffs
    request.cls.p_stoich = request.cls.gas.product_stoich_coeffs

    request.cls.rxn = request.cls.gas.reactions()[request.cls.rxn_idx]
    request.cls.rix = [request.cls.gas.species_index(k) for k in request.cls.rxn.reactants.keys()]
    request.cls.pix = [request.cls.gas.species_index(k) for k in request.cls.rxn.products.keys()]

@pytest.fixture(scope='function')
def setup_rate_expression_data(request, setup_rate_expression_tests):
    request.cls.gas.TPX = request.cls.tpx
    request.cls.gas.set_multiplier(0.)
    request.cls.gas.set_multiplier(1., request.cls.rxn_idx)
    request.cls.gas.derivative_settings = {} # reset defaults

    # check stoichiometric coefficient output
    for k, v in request.cls.rxn.reactants.items():
        ix = request.cls.gas.species_index(k)
        assert request.cls.r_stoich[ix, request.cls.rxn_idx] == v
    for k, v in request.cls.rxn.products.items():
        ix = request.cls.gas.species_index(k)
        assert request.cls.p_stoich[ix, request.cls.rxn_idx] == v

class RateExpressionTests:
    """
    Generic test class to check derivatives evaluated for a single reaction within
    a reaction mechanism. It is the responsibility of the deriving class to provide,
    at the minimum, the following class attributes:

    gas
    tpx

    """
    rxn_idx = None # index of reaction to be tested
    phase = None
    rtol = 1e-5
    orders = None
    ix3b = [] # three-body indices
    equation = None
    rate_type = None
    gas = None
    tpx = None

    def test_input(self):
        # ensure that correct equation is referenced
        assert self.equation == self.rxn.equation
        assert self.rate_type == self.rxn.rate.type

    def rop_derivs(self, spc_ix, mode, const_t=True, rtol_deltac=1e-5, atol_deltac=1e-20, ddX=True):
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
        if ddX:
            return drop * self.gas.density_mole
        else:
            return drop

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
            assert rop[self.rxn_idx] == approx(
                drop[self.rxn_idx, spc_ix] * self.gas.X[spc_ix] / order)

            drop_num = self.rop_derivs(spc_ix, mode="forward")
            assert dropm[:, spc_ix] + dropp * self.gas.P == approx(drop_num, rel=self.rtol)

        if isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.rix + self.ix3b):
            assert dropm[self.rxn_idx, spc_ix] # non-zero
            dropm[self.rxn_idx, spc_ix] = 0
        assert not dropm.any()

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
            assert rop[self.rxn_idx] == approx(
                drop[self.rxn_idx, spc_ix] * self.gas.X[spc_ix] / order)

            drop_num = self.rop_derivs(spc_ix, mode="reverse")
            assert dropm[:, spc_ix] + dropp * self.gas.P == approx(drop_num, rel=self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.pix + self.ix3b):
            assert dropm[self.rxn_idx, spc_ix] # non-zero
            dropm[self.rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_net_rop_ddX(self):
        # check derivatives of net rates of progress with respect to mole fractions
        # against numeric result
        drop = self.gas.net_rates_of_progress_ddX
        dropp = self.gas.net_rates_of_progress_ddP

        for spc_ix in self.rix + self.pix:
            drop_num = self.rop_derivs(spc_ix, mode="net")
            ix = drop[:, spc_ix] != 0
            drop_ = drop[:, spc_ix] + dropp * self.gas.P
            assert drop_[ix] == approx(drop_num[ix], rel=self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.rix + self.pix + self.ix3b):
            assert drop[self.rxn_idx, spc_ix] # non-zero
            drop[self.rxn_idx, spc_ix] = 0
        assert not drop.any()

    def test_forward_rop_ddCi(self):
        # check derivatives of forward rates of progress with respect to species
        # concentrations against analytic result
        dropm = self.gas.forward_rates_of_progress_ddCi
        dropp = self.gas.forward_rates_of_progress_ddP

        self.gas.derivative_settings = {"skip-third-bodies": True}
        drop = self.gas.forward_rates_of_progress_ddCi
        rop = self.gas.forward_rates_of_progress
        for spc_ix in self.rix:
            if self.orders is None:
                order = self.r_stoich[spc_ix, self.rxn_idx]
            else:
                order = self.orders[self.gas.species_names[spc_ix]]
            assert rop[self.rxn_idx] == approx(
                drop[self.rxn_idx, spc_ix] * self.gas.concentrations[spc_ix] / order)

            drop_num = self.rop_derivs(spc_ix, mode="forward", ddX=False)
            assert dropm[:, spc_ix] + dropp * self.gas.P == approx(drop_num, rel=1e-3)

        if isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.rix + self.ix3b):
            assert dropm[self.rxn_idx, spc_ix] # non-zero
            dropm[self.rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_reverse_rop_ddCi(self):
        # check derivatives of reverse rates of progress with respect to species
        # concentrations against analytic result
        dropm = self.gas.reverse_rates_of_progress_ddCi
        dropp = self.gas.reverse_rates_of_progress_ddP

        self.gas.derivative_settings = {"skip-third-bodies": True}
        drop = self.gas.reverse_rates_of_progress_ddCi
        rop = self.gas.reverse_rates_of_progress
        for spc_ix in self.pix:
            order = self.p_stoich[spc_ix, self.rxn_idx]
            assert rop[self.rxn_idx] == approx(
                drop[self.rxn_idx, spc_ix] * self.gas.concentrations[spc_ix] / order)

            drop_num = self.rop_derivs(spc_ix, mode="reverse", ddX=False)
            assert dropm[:, spc_ix] + dropp * self.gas.P == approx(drop_num, rel=1e-3)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(self.pix + self.ix3b):
            assert dropm[self.rxn_idx, spc_ix] # non-zero
            dropm[self.rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_net_rop_ddCi(self):
        # check derivatives of net rates of progress with respect to species
        # concentrations against numeric result
        drop = self.gas.net_rates_of_progress_ddCi
        dropp = self.gas.net_rates_of_progress_ddP

        for spc_ix in self.rix + self.pix:
            drop_num = self.rop_derivs(spc_ix, mode="net", ddX=False)
            ix = drop[:, spc_ix] != 0
            drop_ = drop[:, spc_ix] + dropp * self.gas.P
            assert drop_[ix] == approx(drop_num[ix], rel=1e-4)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

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
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = self.gas.P / self.gas.T
        drop = self.gas.forward_rates_of_progress_ddT
        drop += self.gas.forward_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(mode="forward")
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

        if isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        assert drop[self.rxn_idx] # non-zero
        drop[self.rxn_idx] = 0
        assert not drop.any()

    def test_reverse_rop_ddT(self):
        # check derivatives of reverse rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.reverse_rates_of_progress_ddT
        drop += self.gas.reverse_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="reverse", const_p=True)
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = self.gas.P / self.gas.T
        drop = self.gas.reverse_rates_of_progress_ddT
        drop += self.gas.reverse_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(mode="reverse")
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        assert drop[self.rxn_idx] # non-zero
        drop[self.rxn_idx] = 0
        assert not drop.any()

    def test_net_rop_ddT(self):
        # check derivatives of net rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.net_rates_of_progress_ddT
        drop += self.gas.net_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="net", const_p=True)
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = self.gas.P / self.gas.T
        drop = self.gas.net_rates_of_progress_ddT
        drop += self.gas.net_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(mode="forward") - self.rop_ddT(mode="reverse")
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

        if not self.rxn.reversible or isinstance(self.rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        assert drop[self.rxn_idx] # non-zero
        drop[self.rxn_idx] = 0
        assert not drop.any()

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
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

    def test_reverse_rop_ddP(self):
        # check derivatives of reverse rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = self.gas.density_mole / self.gas.P
        drop = self.gas.reverse_rates_of_progress_ddP
        drop += self.gas.reverse_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(mode="reverse")
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

    def test_net_rop_ddP(self):
        # check derivatives of net rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = self.gas.density_mole / self.gas.P
        drop = self.gas.net_rates_of_progress_ddP
        drop += self.gas.net_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(mode="net")
        assert drop[self.rxn_idx] == approx(drop_num, rel=self.rtol)

    def rate_ddT(self, mode=None, const_p=False, rtol=1e-6):
        # numerical derivative for production rates with respect to temperature
        def calc():
            if mode == "creation":
                return self.gas.creation_rates
            if mode == "destruction":
                return self.gas.destruction_rates
            if mode == "net":
                return self.gas.net_production_rates

        dt = self.tpx[0] * rtol
        dp = 0 if const_p else self.tpx[1] * rtol
        self.gas.TP = self.tpx[0] + dt, self.tpx[1] + dp
        rate1 = calc()
        self.gas.TP = self.tpx[:2]
        rate0 = calc()
        return (rate1 - rate0) / dt

    def test_net_rate_ddT(self):
        # check equivalence of numerical and analytical derivatives of net creation
        # rates with respect to temperature

        # constant pressure - need to account for density change
        # numeric: d(omegadot)/dT =
        # analytic: d(omegadot)/dT + dC/dT d(omegadot)/dC
        dcdt = - self.gas.density_mole / self.gas.T
        drate = self.gas.net_production_rates_ddT
        drate += self.gas.net_production_rates_ddC * dcdt
        drate_num = self.rate_ddT(mode="net", const_p=True)
        for spc_ix in self.rix + self.pix:
            assert drate[spc_ix] == approx(drate_num[spc_ix], rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        # numeric: d(omegadot)/dT =
        # analytic: d(omegadot)/dT + dP/dT d(omegadot)/dP
        dpdt = self.gas.P / self.gas.T
        drate = self.gas.net_production_rates_ddT
        drate += self.gas.net_production_rates_ddP * dpdt
        drate_num = self.rate_ddT(mode="creation") - self.rate_ddT(mode="destruction")
        for spc_ix in self.rix + self.pix:
            assert drate[spc_ix] == approx(drate_num[spc_ix], rel=self.rtol)

    def rate_ddX(self, spc_ix, mode=None, const_t=True, rtol_deltac=1e-6,
                 atol_deltac=1e-20, ddX=True):
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
        # cantera calculates kinetics derivatives with respect to mole fractions
        # and concentrations, when ddX flag is true it will return the numerical
        # derivatives in the form of mole fractions but otherwise return concentrations
        if ddX:
            return drate * self.gas.density_mole
        else:
            return drate

    def test_creation_ddX(self):
        # check derivatives of creation rates with respect to mole fractions
        drate = self.gas.creation_rates_ddX
        dratep = self.gas.creation_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "creation")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=self.rtol)

    def test_destruction_ddX(self):
        # check derivatives of destruction rates with respect to mole fractions
        drate = self.gas.destruction_rates_ddX
        dratep = self.gas.destruction_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "destruction")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=self.rtol)

    def test_net_production_ddX(self):
        # check derivatives of destruction rates with respect to mole fractions
        drate = self.gas.net_production_rates_ddX
        dratep = self.gas.net_production_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "net")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=self.rtol)

    def test_creation_ddCi(self):
        # check derivatives of creation rates with respect to mole fractions
        drate = self.gas.creation_rates_ddCi
        dratep = self.gas.creation_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "creation", ddX=False)
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=1e-3)

    def test_destruction_ddCi(self):
        # check derivatives of destruction rates with respect to mole fractions
        drate = self.gas.destruction_rates_ddCi
        dratep = self.gas.destruction_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "destruction", ddX=False)
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * self.gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=1e-3)

    def test_net_production_ddCi(self):
        # check derivatives of destruction rates with respect to mole fractions
        drate = self.gas.net_production_rates_ddCi
        dratep = self.gas.net_production_rates_ddP
        for spc_ix in self.rix + self.pix:
            drate_num = self.rate_ddX(spc_ix, "net", ddX=False)
            ix = drate[:, spc_ix] != 0
            # drate[:, spc_ix] += dratep * self.gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=1e-3)


@pytest.mark.usefixtures("setup_rate_expression_data")
class HydrogenOxygen(RateExpressionTests):

    gas = ct.Solution("h2o2.yaml", transport_model=None)
    #   species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2]
    gas.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.3, 0.1]
    gas.TP = 800, 2 * ct.one_atm

class TestElementaryRev(HydrogenOxygen):
    # Standard elementary reaction with two reactants
    rxn_idx = 2
    equation = "H2 + O <=> H + OH"
    rate_type = "Arrhenius"

class TestElementarySelf(HydrogenOxygen):
    # Elementary reaction with reactant reacting with itself
    rxn_idx = 27
    equation = "2 HO2 <=> H2O2 + O2"
    rate_type = "Arrhenius"

class TestFalloff(HydrogenOxygen):
    # Fall-off reaction
    rxn_idx = 21
    equation = "2 OH (+M) <=> H2O2 (+M)"
    rate_type = "falloff"
    rtol = 1e-4

@pytest.fixture(scope='class')
def setup_three_body_tests(request, setup_rate_expression_tests):
    request.cls.ix3b = list(range(request.cls.gas.n_species))

@pytest.mark.usefixtures("setup_three_body_tests")
class TestThreeBody(HydrogenOxygen):
    # Three body reaction with default efficiency
    rxn_idx = 1
    equation = "H + O + M <=> OH + M"
    rate_type = "Arrhenius"

    def test_thirdbodies_forward(self):
        drop = self.gas.forward_rates_of_progress_ddX
        self.gas.derivative_settings = {"skip-third-bodies": True}
        drops = self.gas.forward_rates_of_progress_ddX
        dropm = drop - drops
        rop = self.gas.forward_rates_of_progress
        assert rop[self.rxn_idx] == approx((dropm[self.rxn_idx] * self.gas.X).sum())

    def test_thirdbodies_reverse(self):
        drop = self.gas.reverse_rates_of_progress_ddX
        self.gas.derivative_settings = {"skip-third-bodies": True}
        drops = self.gas.reverse_rates_of_progress_ddX
        dropm = drop - drops
        rop = self.gas.reverse_rates_of_progress
        assert rop[self.rxn_idx] == approx((dropm[self.rxn_idx] * self.gas.X).sum())

@pytest.mark.usefixtures("setup_rate_expression_data")
class EdgeCases(RateExpressionTests):
    gas = ct.Solution("jacobian-tests.yaml", transport_model=None)
    #   species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR]
    gas.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
    gas.TP = 800, 2 * ct.one_atm

class TestElementaryIrr(EdgeCases):
    # Irreversible elementary reaction with two reactants
    rxn_idx = 0
    equation = "HO2 + O => O2 + OH"
    rate_type = "Arrhenius"

class TestElementaryOne(EdgeCases):
    # Three-body reaction with single reactant species
    rxn_idx = 1
    equation = "H2 <=> 2 H"
    rate_type = "Arrhenius"

class TestElementaryThree(EdgeCases):
    # Elementary reaction with three reactants
    rxn_idx = 2
    equation = "2 H + O <=> H2O"
    rate_type = "Arrhenius"

class TestElementaryFrac(EdgeCases):
    # Elementary reaction with specified reaction order
    rxn_idx = 3
    orders = {"H2": 0.8, "O2": 1.0, "OH": 2.0}
    equation = "0.7 H2 + 0.2 O2 + 0.6 OH => H2O"
    rate_type = "Arrhenius"


@pytest.fixture(scope='class')
def setup_three_body_no_default_tests(request, setup_rate_expression_tests):
    efficiencies = {"H2": 2.0, "H2O": 6.0, "AR": 0.7}
    request.cls.ix3b = [request.cls.gas.species_index(k) for k in efficiencies.keys()]
@pytest.mark.usefixtures("setup_three_body_no_default_tests")
class TestThreeBodyNoDefault(EdgeCases):
    # Three body reaction without default efficiency
    rxn_idx = 4
    equation = "H + O + M <=> OH + M"
    rate_type = "Arrhenius"

class FromScratchCases(RateExpressionTests):

    gas = ct.Solution("kineticsfromscratch.yaml", transport_model=None)
    #   species: [AR, O, H2, H, OH, O2, H2O, H2O2, HO2]
    gas.X = [0.1, 3e-4, 5e-5, 6e-6, 3e-3, 0.6, 0.25, 1e-6, 2e-5]
    gas.TP = 2000, 5 * ct.one_atm

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_forward_rop_ddT(self):
        super().test_forward_rop_ddT()

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_reverse_rop_ddT(self):
        super().test_reverse_rop_ddT()

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_net_rop_ddT(self):
        super().test_net_rop_ddT()

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_net_rate_ddT(self):
        super().test_net_rate_ddT()

@pytest.mark.usefixtures("setup_rate_expression_data")
class TestPlog(FromScratchCases):
    # Plog reaction
    rxn_idx = 3
    equation = "H2 + O2 <=> 2 OH"
    rate_type = "pressure-dependent-Arrhenius"

@pytest.mark.usefixtures("setup_rate_expression_data")
class TestChebyshev(FromScratchCases):
    # Chebyshev reaction
    rxn_idx = 4
    equation = "HO2 <=> O + OH"
    rate_type = "Chebyshev"

@pytest.mark.usefixtures("setup_rate_expression_data")
class TestBlowersMasel(FromScratchCases):
    # Blowers-Masel
    rxn_idx = 6
    equation = "H2 + O <=> H + OH"
    rate_type = "Blowers-Masel"

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_forward_rop_ddT(self):
        super().test_forward_rop_ddT()

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_reverse_rop_ddT(self):
        super().test_reverse_rop_ddT()

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_net_rop_ddT(self):
        super().test_net_rop_ddT()

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_net_rate_ddT(self):
        super().test_net_rate_ddT()


@pytest.fixture(scope='class')
def setup_full_tests(request):
    request.cls.tpx = request.cls.gas.TPX
    request.cls.gas.equilibrate("HP")

@pytest.fixture(scope='function')
def full_tests_data(request, setup_full_tests):
    request.cls.gas.TPX = request.cls.tpx
    request.cls.gas.derivative_settings = {} # reset

class FullTests:
    # Generic test class to check derivatives evaluated for an entire reaction mechanisms
    rtol = 1e-4
    gas = None

    def rop_derivs(self, mode, rtol_deltac=1e-9, atol_deltac=1e-20, ddX=True):
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

        if ddX:
            return drop * self.gas.density_mole
        else:
            return drop

    def test_forward_rop_ddX(self):
        # check forward rop against numerical derivative with respect to mole fractions
        drop = self.gas.forward_rates_of_progress_ddX
        dropp = self.gas.forward_rates_of_progress_ddP
        drop_num = self.rop_derivs(mode="forward")
        stoich = self.gas.reactant_stoich_coeffs
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, self.gas.reaction(i).rate.type)
                print(self.gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_reverse_rop_ddX(self):
        # check reverse rop against numerical derivative with respect to mole fractions
        drop = self.gas.reverse_rates_of_progress_ddX
        dropp = self.gas.reverse_rates_of_progress_ddP
        drop_num = self.rop_derivs(mode="reverse")
        stoich = self.gas.product_stoich_coeffs
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, self.gas.reaction(i).rate.type)
                print(self.gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_net_rop_ddX(self):
        # check net rop against numerical derivative with respect to mole fractions
        drop = self.gas.net_rates_of_progress_ddX
        dropp = self.gas.net_rates_of_progress_ddP
        drop_num = self.rop_derivs(mode="net")
        stoich = self.gas.product_stoich_coeffs - self.gas.reactant_stoich_coeffs
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                if self.gas.reaction(i).reversible:
                    print(i, self.gas.reaction(i).rate.type)
                    print(self.gas.reaction(i))
                    print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                    raise err

    def test_forward_rop_ddCi(self):
        # check forward rop against numerical derivative with respect to species
        # concentrations
        drop = self.gas.forward_rates_of_progress_ddCi
        dropp = self.gas.forward_rates_of_progress_ddP
        drop_num = self.rop_derivs(mode="forward", ddX=False)
        stoich = self.gas.reactant_stoich_coeffs
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, self.gas.reaction(i).rate.type)
                print(self.gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_reverse_rop_ddCi(self):
        # check reverse rop against numerical derivative with respect to species
        # concentrations
        drop = self.gas.reverse_rates_of_progress_ddCi
        dropp = self.gas.reverse_rates_of_progress_ddP
        drop_num = self.rop_derivs(mode="reverse", ddX=False)
        stoich = self.gas.product_stoich_coeffs
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, self.gas.reaction(i).rate.type)
                print(self.gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_net_rop_ddCi(self):
        # check net rop against numerical derivative with respect to species
        # concentrations
        drop = self.gas.net_rates_of_progress_ddCi
        dropp = self.gas.net_rates_of_progress_ddP
        drop_num = self.rop_derivs(mode="net", ddX=False)
        stoich = self.gas.product_stoich_coeffs - self.gas.reactant_stoich_coeffs
        for i in range(self.gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * self.gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
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
        assert drop == approx(drop_num, rel=self.rtol)

    def test_reverse_rop_ddT(self):
        # check reverse rop against numerical derivative with respect to temperature
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.reverse_rates_of_progress_ddT
        drop += self.gas.reverse_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="reverse")
        assert drop == approx(drop_num, rel=self.rtol)

    def test_net_rop_ddT(self):
        # check net rop against numerical derivative with respect to temperature
        dcdt = - self.gas.density_mole / self.gas.T
        drop = self.gas.net_rates_of_progress_ddT
        drop += self.gas.net_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(mode="net")
        try:
            assert drop == approx(drop_num, rel=self.rtol)
        except AssertionError as err:
            i = np.argmax(2 * (drop - drop_num) / (drop + drop_num + 2e-4))
            print(i, self.gas.reaction(i).rate.type)
            print(self.gas.reaction(i))
            print(drop[i])
            print(drop_num[i])
            raise err

@pytest.mark.usefixtures("full_tests_data")
class TestFullHydrogenOxygen(FullTests):

    gas = ct.Solution("h2o2.yaml", transport_model=None)
    gas.TPX = 300, 5 * ct.one_atm, "H2:1, O2:3"


@pytest.mark.usefixtures("full_tests_data")
class TestFullGriMech(FullTests):

    gas = ct.Solution("gri30.yaml", transport_model=None)
    gas.TPX = 300, 5 * ct.one_atm, "CH4:1, C3H8:.1, O2:1, N2:3.76"

@pytest.mark.usefixtures("full_tests_data")
class TestFullEdgeCases(FullTests):

    gas = ct.Solution("jacobian-tests.yaml", transport_model=None)
    #   species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR]
    gas.TPX = 300, 2 * ct.one_atm, "H2:1, O2:3, AR:0.4"


@pytest.fixture(scope='class')
def setup_surface_rate_expression_tests(request):

    # gas state information
    request.cls.gas_tpx = request.cls.gas.TPX
    request.cls.surf_tpx = request.cls.surf.TPX

    # all species indices
    all_species = request.cls.surf.species() + request.cls.gas.species()
    request.cls.sidxs  = {spec.name:i for i, spec in enumerate(all_species)}

    # kinetics objects
    request.cls.r_stoich = request.cls.surf.reactant_stoich_coeffs
    request.cls.p_stoich = request.cls.surf.product_stoich_coeffs
    request.cls.rxn = request.cls.surf.reactions()[request.cls.rxn_idx]
    request.cls.rix = [request.cls.sidxs[k] for k in request.cls.rxn.reactants.keys()]
    request.cls.pix = [request.cls.sidxs[k] for k in request.cls.rxn.products.keys()]

@pytest.fixture(scope='function')
def surface_rate_expression_data(request, setup_surface_rate_expression_tests):
    # gas phase
    request.cls.gas.TPX = request.cls.gas_tpx
    request.cls.gas.set_multiplier(0)
    request.cls.gas.derivative_settings = {} # reset defaults

    # surface phase
    request.cls.surf.TPX = request.cls.surf_tpx
    request.cls.surf.set_multiplier(0.)
    request.cls.surf.set_multiplier(1., request.cls.rxn_idx)
    request.cls.surf.derivative_settings = {"skip-coverage-dependence": True, "skip-electrochemistry": True}

    # check stoichiometric coefficient output
    for k, v in request.cls.rxn.reactants.items():
        ix = request.cls.sidxs[k]
        assert request.cls.r_stoich[ix, request.cls.rxn_idx] == v
    for k, v in request.cls.rxn.products.items():
        ix = request.cls.sidxs[k]
        assert request.cls.p_stoich[ix, request.cls.rxn_idx] == v

class SurfaceRateExpressionTests:
    """
    Generic test class to check derivatives evaluated for a single reaction within
    a reaction mechanism for surfaces.
    """

    rxn_idx = None # index of reaction to be tested
    phase = None
    rtol = 1e-5
    orders = None
    ix3b = [] # three-body reaction indices
    equation = None
    rate_type = None

    def test_input(self):
        # ensure that correct equation is referenced
        assert self.equation == self.rxn.equation
        assert self.rate_type == self.rxn.rate.type

    def test_forward_rop_ddCi(self):
        # check derivatives of forward rates of progress with respect to species
        # concentrations against analytic result
        drop = self.surf.forward_rates_of_progress_ddCi
        rop = self.surf.forward_rates_of_progress
        concentrations = np.concatenate((self.surf.concentrations, self.gas.concentrations))
        specs = self.surf.species_names + self.gas.species_names
        for spc_ix in self.rix:
            if self.orders is None:
                order = self.r_stoich[spc_ix, self.rxn_idx]
            else:
                order = self.orders[specs[spc_ix]]
            assert rop[self.rxn_idx] == approx(
                drop[self.rxn_idx, spc_ix] * concentrations[spc_ix] / order)

    def test_reverse_rop_ddCi(self):
        # check derivatives of forward rates of progress with respect to species
        # concentrations against analytic result
        drop = self.surf.reverse_rates_of_progress_ddCi
        rop = self.surf.reverse_rates_of_progress
        concentrations = np.concatenate((self.surf.concentrations, self.gas.concentrations))
        specs = self.surf.species_names + self.gas.species_names
        for spc_ix in self.pix:
            if self.orders is None:
                order = self.p_stoich[spc_ix, self.rxn_idx]
            else:
                order = self.orders[specs[spc_ix]]
            assert rop[self.rxn_idx] == approx(
                       drop[self.rxn_idx, spc_ix] * concentrations[spc_ix] / order)

    def test_net_rop_ddCi(self):
        # check derivatives of net rates of progress with respect to species
        # concentrations against analytic
        rop = self.surf.net_rates_of_progress
        drop = self.surf.net_rates_of_progress_ddCi
        concentrations = np.concatenate((self.surf.concentrations, self.gas.concentrations))

        drop *= concentrations
        for spc_ix in self.rix + self.pix:
            ix = np.abs(drop[:, spc_ix]) > 1
            drop_ = drop[:, spc_ix]
            assert drop_[ix] == approx(rop[ix], rel=1e-3)

@pytest.mark.usefixtures("surface_rate_expression_data")
class PlatinumHydrogen(SurfaceRateExpressionTests):

    phase_defs = """
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        phases:
        - name: gas
          thermo: ideal-gas
          species:
          - gri30.yaml/species: [H2, H2O, H2O2, O2]
          kinetics: gas
          reactions:
          - gri30.yaml/reactions: declared-species
          skip-undeclared-third-bodies: true
        - name: Pt_surf
          thermo: ideal-surface
          species:
          - ptcombust.yaml/species: [PT(S), H(S), H2O(S), OH(S), O(S)]
          kinetics: surface
          reactions: [ptcombust.yaml/reactions: declared-species]
          site-density: 3e-09
    """
    # create phase objects
    gas = ct.Solution(yaml=phase_defs, name="gas")
    surf = ct.Interface(yaml=phase_defs, name="Pt_surf", adjacent=[gas])
    gas.TPX = 800, 2*ct.one_atm, "H2:1.5, O2:1.0, H2O2:0.75, H2O:0.3"
    surf.TPX = 800, 2*ct.one_atm , "PT(S):4.0, H(S):0.5, H2O(S):0.1, OH(S):0.2, O(S):0.8"


class TestSurfInterfaceArrhenius(PlatinumHydrogen):
    rxn_idx = 7
    equation = "H(S) + O(S) <=> OH(S) + PT(S)"
    rate_type = "interface-Arrhenius"

class TestSurfGasFwdStickingArrhenius(PlatinumHydrogen):
    rxn_idx = 5
    equation = "H2O + PT(S) => H2O(S)"
    rate_type = "sticking-Arrhenius"

class TestSurfGasInterfaceArrhenius(PlatinumHydrogen):
    rxn_idx = 0
    equation = "H2 + 2 PT(S) => 2 H(S)"
    rate_type = "interface-Arrhenius"
    orders = {"PT(S)": 1, "H2": 1, "H(S)": 2}

class TestGasSurfInterfaceArrhenius(PlatinumHydrogen):
    rxn_idx = 6
    equation = "H2O(S) => H2O + PT(S)"
    rate_type = "interface-Arrhenius"


@pytest.fixture(scope='class')
def setup_surface_full_tests(request):
    # all species indices
    all_species = request.cls.surf.species() + request.cls.gas.species()
    request.cls.sidxs  = {spec.name:i for i, spec in enumerate(all_species)}
    request.cls.gas_tpx = request.cls.gas.TPX
    request.cls.surf_tpx = request.cls.surf.TPX


@pytest.fixture(scope='function')
def setup_surface_full_data(request, setup_surface_full_tests):
    # gas phase
    request.cls.gas.TPX = request.cls.gas_tpx
    request.cls.gas.derivative_settings = {} # reset defaults

    # surface phase
    request.cls.surf.TPX = request.cls.surf_tpx
    request.cls.surf.derivative_settings = {"skip-coverage-dependence": True, "skip-electrochemistry": True}

class SurfaceFullTests:
    # Generic test class to check derivatives evaluated for an entire reaction mechanisms
    rtol = 1e-4

    gas = None
    surf = None

    # closure to get concentrations vector
    def get_concentrations(self):
        return np.concatenate((self.surf.concentrations, self.gas.concentrations))

    def test_forward_rop_ddCi(self):
        # matrix multiplication of the forward rates of progress derivatives w.r.t
        # concentration and the concentrations should provide the rate of progress
        # for each species and can be compared to the directly calculated rate
        drop = self.surf.forward_rates_of_progress_ddCi
        rop = self.surf.forward_rates_of_progress
        conc = self.get_concentrations()
        # multiply derivatives with concentrations
        drop = drop @ conc
        # get total reactant reaction orders
        total_orders = []
        for rxn in self.surf.reactions():
            orders = rxn.orders
            curr_order = 0
            for k, v in rxn.reactants.items():
                if k in orders:
                    curr_order += orders[k]
                else:
                    curr_order += v
            total_orders.append(curr_order)
        total_orders = np.array(total_orders)
        # rates of progress do not factor in reaction order it must be accounted for
        drop /= total_orders
        # compare the rate of progress vectors produced in different ways
        assert drop == approx(rop, rel=self.rtol)

    def test_reverse_rop_ddCi(self):
        # matrix multiplication of the reverse rate of progress derivatives  w.r.t
        # concentration and the concentrations should provide the rate of progress
        # for each species and can be compared to the directly calculated rate
        drop = self.surf.reverse_rates_of_progress_ddCi
        rop = self.surf.reverse_rates_of_progress
        conc = self.get_concentrations()
        # multiply derivatives with concentrations
        drop = drop @ conc
        # get total reactant reaction orders
        total_orders = []
        for rxn in self.surf.reactions():
            orders = rxn.orders
            curr_order = 0
            for k, v in rxn.products.items():
                if k in orders:
                    curr_order += orders[k]
                else:
                    curr_order += v
            total_orders.append(curr_order)
        total_orders = np.array(total_orders)
        # rates of progress do not factor in reaction order it must be accounted for
        drop /= total_orders
        # compare the rate of progress vectors produced in different ways
        assert drop == approx(rop, rel=self.rtol)

    def test_net_rop_ddCi(self):
        # check derivatives of net rates of progress with respect to species
        # concentrations against analytic
        ropf = self.surf.forward_rates_of_progress
        ropr = self.surf.reverse_rates_of_progress
        drop = self.surf.net_rates_of_progress_ddCi
        conc = self.get_concentrations()
        # multiply derivatives with concentrations
        drop = drop @ conc
        # reaction orders are not yet accounted for in rates of progress
        # so they must be included manually
        for i, rxn in enumerate(self.surf.reactions()):
            orders = rxn.orders
            curr_order = 0
            # adjust forward rates by reactant order
            for k, v in rxn.reactants.items():
                curr_order += orders[k] if k in orders else v
            ropf[i] *= curr_order
            curr_order = 0
            # adjust reverse rates by product order
            for k, v in rxn.products.items():
                curr_order += orders[k] if k in orders else v
            ropr[i] *= curr_order
        # compare the rate of progress vectors produced in different ways
        assert drop == approx(ropf - ropr, rel=self.rtol)

@pytest.mark.usefixtures("setup_surface_full_data")
class TestFullPlatinumHydrogen(SurfaceFullTests):

    phase_defs = """
        units: {length: cm, quantity: mol, activation-energy: J/mol}
        phases:
        - name: gas
          thermo: ideal-gas
          species:
          - gri30.yaml/species: [H2, H2O, H2O2, O2]
          kinetics: gas
          reactions:
          - gri30.yaml/reactions: declared-species
          skip-undeclared-third-bodies: true
        - name: Pt_surf
          thermo: ideal-surface
          species:
          - ptcombust.yaml/species: [PT(S), H(S), H2O(S), OH(S), O(S)]
          kinetics: surface
          reactions: [ptcombust.yaml/reactions: declared-species]
          site-density: 3e-09
    """
    # create phase objects
    gas = ct.Solution(yaml=phase_defs, name="gas")
    surf = ct.Interface(yaml=phase_defs, name="Pt_surf", adjacent=[gas])
    gas.TPX = 800, 2*ct.one_atm, "H2:1.5, O2:1.0, H2O2:0.75, H2O:0.3"
    surf.TPX = 800, 2*ct.one_atm , "PT(S):1.0, H(S):0.5, H2O(S):0.1, OH(S):0.2, O(S):0.8"
