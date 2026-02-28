import numpy as np
import pytest
from pytest import approx

import cantera as ct


class RateExpressionTests:
    """
    Generic test class to check derivatives evaluated for a single reaction within
    a reaction mechanism. Derived classes must provide the following fixtures:

    - gas: Cantera Solution object representing the gas phase.
    - rxn_idx: Index of the reaction to be tested.
    - equation: String containing the reaction equation.
    - rate_type: String indicating the type of rate expression.
    - orders: Dictionary specifying reaction orders (if applicable).
    - ix3b: List of species indices involved in three-body interactions (if any).
    """
    rtol = 1e-5 # relative tolerance for comparisons

    rxn_idx = None
    equation = None
    rate_type = None
    orders = None
    ix3b = []

    @pytest.fixture(scope='class')
    def rxn(self, gas, rxn_idx):
        """Fixture to retrieve the reaction object."""
        return gas.reactions()[rxn_idx]

    @pytest.fixture(scope='class')
    def r_stoich(self, gas):
        """Fixture to retrieve reactant stoichiometric coefficients."""
        return gas.reactant_stoich_coeffs

    @pytest.fixture(scope='class')
    def p_stoich(self, gas):
        """Fixture to retrieve product stoichiometric coefficients."""
        return gas.product_stoich_coeffs

    @pytest.fixture(scope='class')
    def rix(self, gas, rxn):
        """Fixture to retrieve reactant species indices."""
        return [gas.species_index(k) for k in rxn.reactants.keys()]

    @pytest.fixture(scope='class')
    def pix(self, gas, rxn):
        """Fixture to retrieve product species indices."""
        return [gas.species_index(k) for k in rxn.products.keys()]

    @pytest.fixture(scope='class')
    def rxn_idx(self, test_data):
        return test_data["rxn_idx"]

    @pytest.fixture(scope='class')
    def equation(self, test_data):
        return test_data["equation"]

    @pytest.fixture(scope='class')
    def rate_type(self, test_data):
        return test_data["rate_type"]

    @pytest.fixture(scope='class')
    def ix3b(self, test_data):
        """Fixture to retrieve species indices involved in three-body interactions."""
        return test_data.get("ix3b", [])

    @pytest.fixture(scope='class')
    def orders(self, test_data):
        """Fixture to retrieve reaction orders."""
        return test_data.get("orders", None)

    @pytest.fixture(scope='class')
    def setup_rate_expression_tests(self, gas, r_stoich, p_stoich, rxn):
        """
        Sets the TPX state for the gas object and provides stoichiometric coefficients
        and species indices as well as setting test-specific attributes. Runs once
        per class.
        """
        tpx = gas.TPX
        return {
            "tpx": tpx,
            "r_stoich": r_stoich,
            "p_stoich": p_stoich
        }

    @pytest.fixture(scope='function', autouse=True)
    def setup_rate_expression_data(self, gas, setup_rate_expression_tests, rxn_idx, rxn):
        """
        Resets the gas state before each test and configures reaction multipliers.
        Also verifies stoichiometric coefficients.
        """
        data = setup_rate_expression_tests
        tpx = data["tpx"]
        r_stoich = data["r_stoich"]
        p_stoich = data["p_stoich"]

        # Reset gas phase
        gas.TPX = tpx
        gas.set_multiplier(0.0)
        gas.set_multiplier(1.0, rxn_idx)
        gas.derivative_settings = {}

        # Check stoichiometric coefficients for reactants
        for k, v in rxn.reactants.items():
            ix = gas.species_index(k)
            actual = r_stoich[ix, rxn_idx]
            assert actual == v, f"Reactant stoich mismatch for species '{k}'"

        # Check stoichiometric coefficients for products
        for k, v in rxn.products.items():
            ix = gas.species_index(k)
            actual = p_stoich[ix, rxn_idx]
            assert actual == v, f"Product stoich mismatch for species '{k}'"

    def test_input(self, equation, rate_type, rxn):
        """Ensure that correct equation is referenced"""
        assert equation == rxn.equation
        assert rate_type == rxn.rate.type

    def rop_derivs(self, gas, spc_ix, mode, const_t=True, rtol_deltac=1e-5, atol_deltac=1e-20, ddX=True):
        """
        Numerical derivative for rates-of-progress with respect to mole fractions
        """
        def calc():
            if mode == "forward":
                return gas.forward_rates_of_progress
            if mode == "reverse":
                return gas.reverse_rates_of_progress
            if mode == "net":
                return gas.net_rates_of_progress

        tpx = gas.TPX
        rop0 = calc()
        conc = gas.concentrations
        ctot0 = conc.sum()

        # perturb concentration
        dconc = conc[spc_ix] * rtol_deltac + atol_deltac
        conc[spc_ix] += dconc
        ctot1 = conc.sum()
        if const_t:
            # adjust pressure to compensate for concentration change
            pnew = gas.P * ctot1 / ctot0
            gas.TPX = gas.T, pnew, conc / ctot1
        else:
            # adjust temperature to compensate for concentration change
            tnew = gas.T * ctot1 / ctot0
            gas.TPX = tnew, gas.P, conc / ctot1
        drop = (calc() - rop0) / dconc

        gas.TPX = tpx # reset state
        if ddX:
            return drop * gas.density_mole
        else:
            return drop

    def test_forward_rop_ddX(self, gas, r_stoich, rxn_idx, rix, rxn, ix3b, orders):
        # check derivatives of forward rates of progress with respect to mole fractions
        # against analytic result
        dropm = gas.forward_rates_of_progress_ddX
        dropp = gas.forward_rates_of_progress_ddP

        gas.derivative_settings = {"skip-third-bodies": True}
        drop = gas.forward_rates_of_progress_ddX
        rop = gas.forward_rates_of_progress
        for spc_ix in rix:
            if orders is None:
                order = r_stoich[spc_ix, rxn_idx]
            else:
                order = orders[gas.species_names[spc_ix]]
            assert rop[rxn_idx] == approx(
                drop[rxn_idx, spc_ix] * gas.X[spc_ix] / order)

            drop_num = self.rop_derivs(gas, spc_ix, mode="forward")
            assert dropm[:, spc_ix] + dropp * gas.P == approx(drop_num, rel=self.rtol)

        if isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(rix + ix3b):
            assert dropm[rxn_idx, spc_ix] # non-zero
            dropm[rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_reverse_rop_ddX(self, gas, rxn_idx, p_stoich, pix, rxn, ix3b):
        # check derivatives of reverse rates of progress with respect to mole fractions
        # against analytic result
        dropm = gas.reverse_rates_of_progress_ddX
        dropp = gas.reverse_rates_of_progress_ddP

        gas.derivative_settings = {"skip-third-bodies": True}
        drop = gas.reverse_rates_of_progress_ddX
        rop = gas.reverse_rates_of_progress
        for spc_ix in pix:
            order = p_stoich[spc_ix, rxn_idx]
            assert rop[rxn_idx] == approx(
                drop[rxn_idx, spc_ix] * gas.X[spc_ix] / order)

            drop_num = self.rop_derivs(gas, spc_ix, mode="reverse")
            assert dropm[:, spc_ix] + dropp * gas.P == approx(drop_num, rel=self.rtol)

        if not rxn.reversible or isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(pix + ix3b):
            assert dropm[rxn_idx, spc_ix] # non-zero
            dropm[rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_net_rop_ddX(self, gas, rxn_idx, rix, pix, rxn, ix3b):
        # check derivatives of net rates of progress with respect to mole fractions
        # against numeric result
        drop = gas.net_rates_of_progress_ddX
        dropp = gas.net_rates_of_progress_ddP

        for spc_ix in rix + pix:
            drop_num = self.rop_derivs(gas, spc_ix, mode="net")
            ix = drop[:, spc_ix] != 0
            drop_ = drop[:, spc_ix] + dropp * gas.P
            assert drop_[ix] == approx(drop_num[ix], rel=self.rtol)

        if not rxn.reversible or isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(rix + pix + ix3b):
            assert drop[rxn_idx, spc_ix] # non-zero
            drop[rxn_idx, spc_ix] = 0
        assert not drop.any()

    def test_forward_rop_ddCi(self, gas, rxn_idx, r_stoich, rix, rxn, ix3b, orders):
        # check derivatives of forward rates of progress with respect to species
        # concentrations against analytic result
        dropm = gas.forward_rates_of_progress_ddCi
        dropp = gas.forward_rates_of_progress_ddP

        gas.derivative_settings = {"skip-third-bodies": True}
        drop = gas.forward_rates_of_progress_ddCi
        rop = gas.forward_rates_of_progress
        for spc_ix in rix:
            if orders is None:
                order = r_stoich[spc_ix, rxn_idx]
            else:
                order = orders[gas.species_names[spc_ix]]
            assert rop[rxn_idx] == approx(
                drop[rxn_idx, spc_ix] * gas.concentrations[spc_ix] / order)

            drop_num = self.rop_derivs(gas, spc_ix, mode="forward", ddX=False)
            assert dropm[:, spc_ix] + dropp * gas.P == approx(drop_num, rel=1e-3)

        if isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(rix + ix3b):
            assert dropm[rxn_idx, spc_ix] # non-zero
            dropm[rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_reverse_rop_ddCi(self, gas, rxn_idx, r_stoich, pix, p_stoich, rix, rxn, ix3b, orders):
        # check derivatives of reverse rates of progress with respect to species
        # concentrations against analytic result
        dropm = gas.reverse_rates_of_progress_ddCi
        dropp = gas.reverse_rates_of_progress_ddP

        gas.derivative_settings = {"skip-third-bodies": True}
        drop = gas.reverse_rates_of_progress_ddCi
        rop = gas.reverse_rates_of_progress
        for spc_ix in pix:
            order = p_stoich[spc_ix, rxn_idx]
            assert rop[rxn_idx] == approx(
                drop[rxn_idx, spc_ix] * gas.concentrations[spc_ix] / order)

            drop_num = self.rop_derivs(gas, spc_ix, mode="reverse", ddX=False)
            assert dropm[:, spc_ix] + dropp * gas.P == approx(drop_num, rel=1e-3)

        if not rxn.reversible or isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        for spc_ix in set(pix + ix3b):
            assert dropm[rxn_idx, spc_ix] # non-zero
            dropm[rxn_idx, spc_ix] = 0
        assert not dropm.any()

    def test_net_rop_ddCi(self, gas, rxn_idx, r_stoich, rix, pix, rxn, ix3b, orders):
        # check derivatives of net rates of progress with respect to species
        # concentrations against numeric result
        drop = gas.net_rates_of_progress_ddCi
        dropp = gas.net_rates_of_progress_ddP

        for spc_ix in rix + pix:
            drop_num = self.rop_derivs(gas, spc_ix, mode="net", ddX=False)
            ix = drop[:, spc_ix] != 0
            drop_ = drop[:, spc_ix] + dropp * gas.P
            assert drop_[ix] == approx(drop_num[ix], rel=1e-4)

        if not rxn.reversible or isinstance(rxn.rate, ct.FalloffRate):
            return

    def rop_ddT(self, gas, rxn_idx, mode=None, const_p=False, rtol=1e-6):
        # numerical derivative for rates-of-progress at constant pressure
        def calc():
            if mode == "forward":
                return gas.forward_rates_of_progress
            if mode == "reverse":
                return gas.reverse_rates_of_progress
            if mode == "net":
                return gas.net_rates_of_progress

        tpx = gas.TPX
        dt = tpx[0] * rtol
        dp = 0 if const_p else tpx[1] * rtol
        gas.TP = tpx[0] + dt, tpx[1] + dp
        rop1 = calc()
        gas.TP = tpx[:2]
        rop0 = calc()
        gas.TPX = tpx
        return (rop1[rxn_idx] - rop0[rxn_idx]) / dt

    def test_forward_rop_ddT(self, gas, rxn_idx, rxn):
        # check derivatives of forward rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - gas.density_mole / gas.T
        drop = gas.forward_rates_of_progress_ddT
        drop += gas.forward_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(gas, rxn_idx, mode="forward", const_p=True)
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = gas.P / gas.T
        drop = gas.forward_rates_of_progress_ddT
        drop += gas.forward_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(gas, rxn_idx, mode="forward")
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

        if isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        assert drop[rxn_idx] # non-zero
        drop[rxn_idx] = 0
        assert not drop.any()

    def test_reverse_rop_ddT(self, gas, rxn_idx, rxn):
        # check derivatives of reverse rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - gas.density_mole / gas.T
        drop = gas.reverse_rates_of_progress_ddT
        drop += gas.reverse_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(gas, rxn_idx, mode="reverse", const_p=True)
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = gas.P / gas.T
        drop = gas.reverse_rates_of_progress_ddT
        drop += gas.reverse_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(gas, rxn_idx, mode="reverse")
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

        if not rxn.reversible or isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        assert drop[rxn_idx] # non-zero
        drop[rxn_idx] = 0
        assert not drop.any()

    def test_net_rop_ddT(self, gas, rxn_idx, rxn):
        # check derivatives of net rop with respect to temperature

        # constant pressure - need to account for density change
        dcdt = - gas.density_mole / gas.T
        drop = gas.net_rates_of_progress_ddT
        drop += gas.net_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(gas, rxn_idx, mode="net", const_p=True)
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        dpdt = gas.P / gas.T
        drop = gas.net_rates_of_progress_ddT
        drop += gas.net_rates_of_progress_ddP * dpdt
        drop_num = self.rop_ddT(gas, rxn_idx, mode="forward") - self.rop_ddT(gas, rxn_idx, mode="reverse")
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

        if not rxn.reversible or isinstance(rxn.rate, ct.FalloffRate):
            return

        # ensure all zeros are in the correct spots
        assert drop[rxn_idx] # non-zero
        drop[rxn_idx] = 0
        assert not drop.any()

    def rop_ddP(self, gas, rxn_idx, mode=None, rtol=1e-6):
        # numerical derivative for rates-of-progress at constant pressure
        def calc():
            if mode == "forward":
                return gas.forward_rates_of_progress
            if mode == "reverse":
                return gas.reverse_rates_of_progress
            if mode == "net":
                return gas.net_rates_of_progress

        tpx = gas.TPX
        dp = tpx[1] * rtol
        gas.TP = tpx[0], tpx[1] + dp
        rop1 = calc()
        gas.TP = tpx[:2]
        rop0 = calc()
        return (rop1[rxn_idx] - rop0[rxn_idx]) / dp

    def test_forward_rop_ddP(self, gas, rxn_idx):
        # check derivatives of forward rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = gas.density_mole / gas.P
        drop = gas.forward_rates_of_progress_ddP
        drop += gas.forward_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(gas, rxn_idx, mode="forward")
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

    def test_reverse_rop_ddP(self, gas, rxn_idx):
        # check derivatives of reverse rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = gas.density_mole / gas.P
        drop = gas.reverse_rates_of_progress_ddP
        drop += gas.reverse_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(gas, rxn_idx, mode="reverse")
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

    def test_net_rop_ddP(self, gas, rxn_idx):
        # check derivatives of net rop with respect to pressure

        # constant temperature - need to account for density change
        dcdp = gas.density_mole / gas.P
        drop = gas.net_rates_of_progress_ddP
        drop += gas.net_rates_of_progress_ddC * dcdp
        drop_num = self.rop_ddP(gas, rxn_idx, mode="net")
        assert drop[rxn_idx] == approx(drop_num, rel=self.rtol)

    def rate_ddT(self, gas, mode=None, const_p=False, rtol=1e-6):
        # numerical derivative for production rates with respect to temperature
        def calc():
            if mode == "creation":
                return gas.creation_rates
            if mode == "destruction":
                return gas.destruction_rates
            if mode == "net":
                return gas.net_production_rates

        tpx = gas.TPX
        dt = tpx[0] * rtol
        dp = 0 if const_p else tpx[1] * rtol
        gas.TP = tpx[0] + dt, tpx[1] + dp
        rate1 = calc()
        gas.TP = tpx[:2]
        rate0 = calc()
        gas.TPX = tpx
        return (rate1 - rate0) / dt

    def test_net_rate_ddT(self, gas, rix, pix):
        # check equivalence of numerical and analytical derivatives of net creation
        # rates with respect to temperature

        # constant pressure - need to account for density change
        # numeric: d(omegadot)/dT =
        # analytic: d(omegadot)/dT + dC/dT d(omegadot)/dC
        dcdt = - gas.density_mole / gas.T
        drate = gas.net_production_rates_ddT
        drate += gas.net_production_rates_ddC * dcdt
        drate_num = self.rate_ddT(gas, mode="net", const_p=True)
        for spc_ix in rix + pix:
            assert drate[spc_ix] == approx(drate_num[spc_ix], rel=self.rtol)

        # constant density (volume) - need to account for pressure change
        # numeric: d(omegadot)/dT =
        # analytic: d(omegadot)/dT + dP/dT d(omegadot)/dP
        dpdt = gas.P / gas.T
        drate = gas.net_production_rates_ddT
        drate += gas.net_production_rates_ddP * dpdt
        drate_num = self.rate_ddT(gas, mode="creation") - self.rate_ddT(gas, mode="destruction")
        for spc_ix in rix + pix:
            assert drate[spc_ix] == approx(drate_num[spc_ix], rel=self.rtol)

    def rate_ddX(self, gas, spc_ix, mode=None, const_t=True, rtol_deltac=1e-6,
                 atol_deltac=1e-20, ddX=True):
        # numerical derivative for production rates with respect to mole fractions
        def calc(mode):
            if mode == "creation":
                return gas.creation_rates
            if mode == "destruction":
                return gas.destruction_rates
            if mode == "net":
                return gas.net_production_rates

        tpx = gas.TPX

        rate0 = calc(mode)
        conc = gas.concentrations
        ctot0 = conc.sum()

        # perturb concentration
        dconc = conc[spc_ix] * rtol_deltac + atol_deltac
        conc[spc_ix] += dconc
        ctot1 = conc.sum()
        if const_t:
            # adjust pressure to compensate for concentration change
            pnew = gas.P * ctot1 / ctot0
            gas.TPX = gas.T, pnew, conc / ctot1
        else:
            # adjust temperature to compensate for concentration change
            tnew = gas.T * ctot1 / ctot0
            gas.TPX = tnew, gas.P, conc / ctot1
        drate = (calc(mode) - rate0) / dconc

        gas.TPX = tpx # reset state
        # cantera calculates kinetics derivatives with respect to mole fractions
        # and concentrations, when ddX flag is true it will return the numerical
        # derivatives in the form of mole fractions but otherwise return concentrations
        if ddX:
            return drate * gas.density_mole
        else:
            return drate

    def test_creation_ddX(self, gas, rix, pix):
        # check derivatives of creation rates with respect to mole fractions
        drate = gas.creation_rates_ddX
        dratep = gas.creation_rates_ddP
        for spc_ix in rix + pix:
            drate_num = self.rate_ddX(gas, spc_ix, "creation")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=self.rtol)

    def test_destruction_ddX(self, gas, rix, pix):
        # check derivatives of destruction rates with respect to mole fractions
        drate = gas.destruction_rates_ddX
        dratep = gas.destruction_rates_ddP
        for spc_ix in rix + pix:
            drate_num = self.rate_ddX(gas, spc_ix, "destruction")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=self.rtol)

    def test_net_production_ddX(self, gas, rix, pix):
        # check derivatives of destruction rates with respect to mole fractions
        drate = gas.net_production_rates_ddX
        dratep = gas.net_production_rates_ddP
        for spc_ix in rix + pix:
            drate_num = self.rate_ddX(gas, spc_ix, "net")
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=self.rtol)

    def test_creation_ddCi(self, gas, rix, pix):
        # check derivatives of creation rates with respect to mole fractions
        drate = gas.creation_rates_ddCi
        dratep = gas.creation_rates_ddP
        for spc_ix in rix + pix:
            drate_num = self.rate_ddX(gas, spc_ix, "creation", ddX=False)
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=1e-3)

    def test_destruction_ddCi(self, gas, rix, pix):
        # check derivatives of destruction rates with respect to mole fractions
        drate = gas.destruction_rates_ddCi
        dratep = gas.destruction_rates_ddP
        for spc_ix in rix + pix:
            drate_num = self.rate_ddX(gas, spc_ix, "destruction", ddX=False)
            ix = drate[:, spc_ix] != 0
            drate[:, spc_ix] += dratep * gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=1e-3)

    def test_net_production_ddCi(self, gas, rix, pix):
        # check derivatives of destruction rates with respect to mole fractions
        drate = gas.net_production_rates_ddCi
        dratep = gas.net_production_rates_ddP
        for spc_ix in rix + pix:
            drate_num = self.rate_ddX(gas, spc_ix, "net", ddX=False)
            ix = drate[:, spc_ix] != 0
            # drate[:, spc_ix] += dratep * gas.P
            assert drate[ix, spc_ix] == approx(drate_num[ix], rel=1e-3)


class HydrogenOxygen(RateExpressionTests):

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("h2o2.yaml", transport_model=None)
        # species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2]
        gas.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.3, 0.1]
        gas.TP = 800, 2 * ct.one_atm
        return gas

class TestElementaryRev(HydrogenOxygen):
    """Standard elementary reaction with two reactants"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 2,
            "equation": "H2 + O <=> H + OH",
            "rate_type": "Arrhenius"
        }

class TestElementarySelf(HydrogenOxygen):
    """Elementary reaction with reactant reacting with itself"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 27,
            "equation": "2 HO2 <=> H2O2 + O2",
            "rate_type": "Arrhenius"
        }


class TestFalloff(HydrogenOxygen):
    """ Fall-off reaction"""
    rtol = 2e-4

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 21,
            "equation": "2 OH (+M) <=> H2O2 (+M)",
            "rate_type": "falloff"
        }

class TestThreeBody(HydrogenOxygen):
    """ Three-body reaction with default efficiency"""

    @pytest.fixture(scope='class')
    def test_data(self, gas):
        return {
            "rxn_idx": 1,
            "equation": "H + O + M <=> OH + M",
            "rate_type": "Arrhenius",
            "ix3b": list(range(gas.n_species))
        }

    def test_thirdbodies_forward(self, gas, rxn_idx):
        drop = gas.forward_rates_of_progress_ddX
        gas.derivative_settings = {"skip-third-bodies": True}
        drops = gas.forward_rates_of_progress_ddX
        dropm = drop - drops
        rop = gas.forward_rates_of_progress
        assert rop[rxn_idx] == approx((dropm[rxn_idx] * gas.X).sum())

    def test_thirdbodies_reverse(self, gas, rxn_idx):
        drop = gas.reverse_rates_of_progress_ddX
        gas.derivative_settings = {"skip-third-bodies": True}
        drops = gas.reverse_rates_of_progress_ddX
        dropm = drop - drops
        rop = gas.reverse_rates_of_progress
        assert rop[rxn_idx] == approx((dropm[rxn_idx] * gas.X).sum())

class EdgeCases(RateExpressionTests):

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("jacobian-tests.yaml", transport_model=None)
        # species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR]
        gas.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        gas.TP = 800, 2 * ct.one_atm
        return gas

class TestElementaryIrr(EdgeCases):
    """Irreversible elementary reaction with two reactants"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 0,
            "equation": "HO2 + O => O2 + OH",
            "rate_type": "Arrhenius"
        }

class TestElementaryOne(EdgeCases):
    """Three-body reaction with single reactant species"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 1,
            "equation": "H2 <=> 2 H",
            "rate_type": "Arrhenius"
        }

class TestElementaryThree(EdgeCases):
    """Elementary reaction with three reactants"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 2,
            "equation": "2 H + O <=> H2O",
            "rate_type": "Arrhenius"
        }


class TestElementaryFrac(EdgeCases):
    """Elementary reaction with specified reaction order"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 3,
            "equation": "0.7 H2 + 0.2 O2 + 0.6 OH => H2O",
            "rate_type": "Arrhenius",
            "orders": {"H2": 0.8, "O2": 1.0, "OH": 2.0}
        }


class TestThreeBodyNoDefault(EdgeCases):
    """Three body reaction without default efficiency"""

    @pytest.fixture(scope='class')
    def test_data(self, gas):
        efficiencies = {"H2": 2.0, "H2O": 6.0, "AR": 0.7}
        return {
            "rxn_idx": 4,
            "equation": "H + O + M <=> OH + M",
            "rate_type": "Arrhenius",
            "ix3b": [gas.species_index(k) for k in efficiencies.keys()]
        }


class FromScratchCases(RateExpressionTests):

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("kineticsfromscratch.yaml", transport_model=None)
        # species: [AR, O, H2, H, OH, O2, H2O, H2O2, HO2]
        gas.X = [0.1, 3e-4, 5e-5, 6e-6, 3e-3, 0.6, 0.25, 1e-6, 2e-5]
        gas.TP = 2000, 5 * ct.one_atm
        return gas

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_forward_rop_ddT(self, gas, rxn_idx, rxn):
        """Override to handle temperature derivative warnings."""
        super().test_forward_rop_ddT(gas, rxn_idx, rxn)

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_reverse_rop_ddT(self, gas, rxn_idx, rxn):
        """Override to handle temperature derivative warnings."""
        super().test_reverse_rop_ddT(gas, rxn_idx, rxn)

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_net_rop_ddT(self, gas, rxn_idx, rxn):
        """Override to handle temperature derivative warnings."""
        super().test_net_rop_ddT(gas, rxn_idx, rxn)

    @pytest.mark.usefixtures("has_temperature_derivative_warnings")
    def test_net_rate_ddT(self, gas, rix, pix):
        """Override to handle temperature derivative warnings."""
        super().test_net_rate_ddT(gas, rix, pix)


class TestPlog(FromScratchCases):
    """ Plog reaction"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 3,
            "equation": "H2 + O2 <=> 2 OH",
            "rate_type": "pressure-dependent-Arrhenius"
        }


class TestChebyshev(FromScratchCases):
    """Chebyshev reaction"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 4,
            "equation": "HO2 <=> O + OH",
            "rate_type": "Chebyshev"
        }


class TestBlowersMasel(FromScratchCases):
    """Blowers-Masel"""

    @pytest.fixture(scope='class')
    def test_data(self):
        return {
            "rxn_idx": 6,
            "equation": "H2 + O <=> H + OH",
            "rate_type": "Blowers-Masel"
        }

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_forward_rop_ddT(self, gas, rxn_idx, rxn):
        """Override to handle issues with Blowers-Masel derivatives."""
        super().test_forward_rop_ddT(gas, rxn_idx, rxn)

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_reverse_rop_ddT(self, gas, rxn_idx, rxn):
        """Override to handle issues with Blowers-Masel derivatives."""
        super().test_reverse_rop_ddT(gas, rxn_idx, rxn)

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_net_rop_ddT(self, gas, rxn_idx, rxn):
        """Override to handle issues with Blowers-Masel derivatives."""
        super().test_net_rop_ddT(gas, rxn_idx, rxn)

    @pytest.mark.xfail(reason="Change of reaction enthalpy is not considered")
    @pytest.mark.filterwarnings("ignore:.*does not consider.*(electron|enthalpy).*:UserWarning")
    def test_net_rate_ddT(self, gas, rix, pix):
        """Override to handle issues with Blowers-Masel derivatives."""
        super().test_net_rate_ddT(gas, rix, pix)



class FullTests:
    """
    Generic test class to check derivatives evaluated for an entire reaction mechanisms
    """
    rtol = 1e-4

    @pytest.fixture(scope='class')
    def initial_conditions(self, gas):
        """Store initial TPX conditions for gas."""
        return gas.TPX

    @pytest.fixture(scope='class', autouse=True)
    def setup_full_tests(self, gas, initial_conditions):
        """
        Equilibrate the gas phase at constant enthalpy and pressure.
        """
        gas.equilibrate("HP")

    @pytest.fixture(scope='function', autouse=True)
    def full_tests_data(self, gas, initial_conditions):
        """
        Reset TPX before each test.
        """
        gas.TPX = initial_conditions

    def rop_derivs(self, gas, mode, rtol_deltac=1e-9, atol_deltac=1e-20, ddX=True):
        # numerical derivative for rates-of-progress with respect to mole fractions
        def calc():
            if mode == "forward":
                return gas.forward_rates_of_progress
            if mode == "reverse":
                return gas.reverse_rates_of_progress
            if mode == "net":
                return gas.net_rates_of_progress

        n_spc, n_rxn = gas.n_species, gas.n_reactions
        drop = np.zeros((n_rxn, n_spc))

        tpx = gas.TPX
        rop0 = calc()
        ctot0 = gas.density_mole
        ctot0 = gas.concentrations.sum()

        for spc_ix in range(n_spc):
            conc = gas.concentrations
            dconc = conc[spc_ix] * rtol_deltac + atol_deltac
            conc[spc_ix] += dconc
            ctot1 = conc.sum()
            gas.TPX = tpx[0], tpx[1] * ctot1 / ctot0, conc / ctot1
            drop[:, spc_ix] = (calc() - rop0) / dconc
            gas.TPX = tpx

        if ddX:
            return drop * gas.density_mole
        else:
            return drop

    def test_forward_rop_ddX(self, gas):
        """
        Check forward rop against numerical derivative with respect to mole fractions
        """
        drop = gas.forward_rates_of_progress_ddX
        dropp = gas.forward_rates_of_progress_ddP
        drop_num = self.rop_derivs(gas, mode="forward")
        stoich = gas.reactant_stoich_coeffs
        for i in range(gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, gas.reaction(i).rate.type)
                print(gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_reverse_rop_ddX(self, gas):
        """
        Check reverse rop against numerical derivative with respect to mole fractions
        """
        drop = gas.reverse_rates_of_progress_ddX
        dropp = gas.reverse_rates_of_progress_ddP
        drop_num = self.rop_derivs(gas, mode="reverse")
        stoich = gas.product_stoich_coeffs
        for i in range(gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, gas.reaction(i).rate.type)
                print(gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_net_rop_ddX(self, gas):
        """
        Check net rop against numerical derivative with respect to mole fractions
        """
        drop = gas.net_rates_of_progress_ddX
        dropp = gas.net_rates_of_progress_ddP
        drop_num = self.rop_derivs(gas, mode="net")
        stoich = gas.product_stoich_coeffs - gas.reactant_stoich_coeffs
        for i in range(gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                if gas.reaction(i).reversible:
                    print(i, gas.reaction(i).rate.type)
                    print(gas.reaction(i))
                    print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                    raise err

    def test_forward_rop_ddCi(self, gas):
        """
        Check forward rop against numerical derivative with respect to species
        concentrations
        """
        drop = gas.forward_rates_of_progress_ddCi
        dropp = gas.forward_rates_of_progress_ddP
        drop_num = self.rop_derivs(gas, mode="forward", ddX=False)
        stoich = gas.reactant_stoich_coeffs
        for i in range(gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, gas.reaction(i).rate.type)
                print(gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_reverse_rop_ddCi(self, gas):
        """
        Check reverse rop against numerical derivative with respect to species
        concentrations
        """
        drop = gas.reverse_rates_of_progress_ddCi
        dropp = gas.reverse_rates_of_progress_ddP
        drop_num = self.rop_derivs(gas, mode="reverse", ddX=False)
        stoich = gas.product_stoich_coeffs
        for i in range(gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                print(i, gas.reaction(i).rate.type)
                print(gas.reaction(i))
                print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                raise err

    def test_net_rop_ddCi(self, gas):
        """
        Check net rop against numerical derivative with respect to species
        concentrations
        """
        drop = gas.net_rates_of_progress_ddCi
        dropp = gas.net_rates_of_progress_ddP
        drop_num = self.rop_derivs(gas, mode="net", ddX=False)
        stoich = gas.product_stoich_coeffs - gas.reactant_stoich_coeffs
        for i in range(gas.n_reactions):
            try:
                # test entries that are not spurious
                ix = np.abs((stoich[:, i] != 0) * drop[i, :]) > 1e-6
                drop_ = drop[i, ix] + dropp[i] * gas.P
                assert drop_ == approx(drop_num[i, ix], rel=self.rtol)
            except AssertionError as err:
                if gas.reaction(i).reversible:
                    print(i, gas.reaction(i).rate.type)
                    print(gas.reaction(i))
                    print(np.vstack([drop[i, ix], drop_num[i, ix]]).T)
                    raise err

    def rop_ddT(self, gas, mode=None, dt=1e-6):
        """Numerical derivative for rates-of-progress at constant pressure"""
        def calc():
            if mode == "forward":
                return gas.forward_rates_of_progress
            if mode == "reverse":
                return gas.reverse_rates_of_progress
            if mode == "net":
                return gas.net_rates_of_progress
            return None

        tpx = gas.TPX
        gas.TP = tpx[0] + dt, tpx[1]
        rop1 = calc()
        gas.TP = tpx[:2]
        rop0 = calc()
        gas.TPX = tpx
        return (rop1 - rop0) / dt

    def test_forward_rop_ddT(self, gas):
        """
        Check forward rop against numerical derivative with respect to temperature
        """
        dcdt = - gas.density_mole / gas.T
        drop = gas.forward_rates_of_progress_ddT
        drop += gas.forward_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(gas, mode="forward")
        assert drop == approx(drop_num, rel=self.rtol)

    def test_reverse_rop_ddT(self, gas):
        """
        Check reverse rop against numerical derivative with respect to temperature
        """
        dcdt = - gas.density_mole / gas.T
        drop = gas.reverse_rates_of_progress_ddT
        drop += gas.reverse_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(gas, mode="reverse")
        assert drop == approx(drop_num, rel=self.rtol)

    def test_net_rop_ddT(self, gas):
        """
        Check net rop against numerical derivative with respect to temperature
        """
        dcdt = - gas.density_mole / gas.T
        drop = gas.net_rates_of_progress_ddT
        drop += gas.net_rates_of_progress_ddC * dcdt
        drop_num = self.rop_ddT(gas, mode="net")
        try:
            assert drop == approx(drop_num, rel=self.rtol)
        except AssertionError as err:
            i = np.argmax(2 * (drop - drop_num) / (drop + drop_num + 2e-4))
            print(i, gas.reaction(i).rate.type)
            print(gas.reaction(i))
            print(drop[i])
            print(drop_num[i])
            raise err

class TestFullHydrogenOxygen(FullTests):

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("h2o2.yaml", transport_model=None)
        gas.TPX = 300, 5 * ct.one_atm, "H2:1, O2:3"
        return gas


class TestFullHydrogenOxygenRedlichKwong(FullTests):
    # Use a loose tolerance because derivatives currently do not account for any
    # non-ideal effects.
    rtol = 2e-3

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("h2o2.yaml", "ohmech-RK", transport_model=None)
        gas.derivative_settings = {"skip-nonideal": True}
        gas.TPX = 500, 2 * ct.one_atm, "H2:1, O2:3"
        return gas


class TestFullGriMech(FullTests):

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("gri30.yaml", transport_model=None)
        gas.TPX = 300, 5 * ct.one_atm, "CH4:1, C3H8:.1, O2:1, N2:3.76"
        return gas

class TestFullEdgeCases(FullTests):

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution("jacobian-tests.yaml", transport_model=None)
        # species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR]
        gas.TPX = 300, 2 * ct.one_atm, "H2:1, O2:3, AR:0.4"
        return gas


class SurfaceRateExpressionTests:
    """
    Generic test class to check derivatives evaluated for a single reaction within
    a reaction mechanism for surfaces.
    """
    rtol = 1e-5  # Relative tolerance for approximate comparisons

    @pytest.fixture(scope='class')
    def all_species(self, gas, surf):
        """Combine species from both gas and surface phases."""
        return surf.species() + gas.species()

    @pytest.fixture(scope='class')
    def sidxs(self, all_species):
        """Create a dictionary mapping species names to indices."""
        return {spec.name: i for i, spec in enumerate(all_species)}

    @pytest.fixture(scope='class')
    def rxn_data(self, surf, sidxs, rxn_idx):
        """
        Retrieve reaction-specific data based on the reaction index.
        Returns a dictionary containing reaction, reactant indices, and product indices.
        """
        rxn = surf.reactions()[rxn_idx]
        r_stoich = surf.reactant_stoich_coeffs
        p_stoich = surf.product_stoich_coeffs
        rix = [sidxs[k] for k in rxn.reactants.keys()]
        pix = [sidxs[k] for k in rxn.products.keys()]
        return {
            "rxn": rxn,
            "r_stoich": r_stoich,
            "p_stoich": p_stoich,
            "rix": rix,
            "pix": pix
        }

    @pytest.fixture(scope='class')
    def initial_conditions(self, gas, surf):
        """Store initial TPX conditions for gas and surface."""
        return {
            "gas_tpx": gas.TPX,
            "surf_tpx": surf.TPX
        }

    @pytest.fixture(scope='function', autouse=True)
    def setup_surface_rate_expression_data(self, gas, surf, initial_conditions, rxn_data, rxn_idx, sidxs):
        """
        Reset TPX and derivative settings before each test and configure reaction multipliers.
        Also verify stoichiometric coefficients.
        """
        # Reset gas phase
        gas.TPX = initial_conditions["gas_tpx"]
        gas.set_multiplier(0)
        gas.derivative_settings = {}  # Reset to defaults

        # Reset surface phase
        surf.TPX = initial_conditions["surf_tpx"]
        surf.set_multiplier(0.0)
        surf.set_multiplier(1.0, rxn_idx)  # Enable multiplier for specific reaction
        surf.derivative_settings = {
            "skip-coverage-dependence": True,
            "skip-electrochemistry": True
        }

        # Retrieve reaction data
        rxn = rxn_data["rxn"]
        r_stoich = rxn_data["r_stoich"]
        p_stoich = rxn_data["p_stoich"]

        # Check stoichiometric coefficients for reactants
        for k, v in rxn.reactants.items():
            ix = sidxs[k]
            actual = r_stoich[ix, rxn_idx]
            assert actual == v, f"Reactant stoich mismatch for species '{k}'"

        # Check stoichiometric coefficients for products
        for k, v in rxn.products.items():
            ix = sidxs[k]
            actual = p_stoich[ix, rxn_idx]
            assert actual == v, f"Product stoich mismatch for species '{k}'"

    def get_concentrations(self, surf, gas):
        """Concatenate concentrations from surface and gas phases."""
        return np.concatenate((surf.concentrations, gas.concentrations))

    # Test methods
    def test_input(self, rxn_data, equation, rate_type):
        """Ensure that the correct equation and rate type are referenced."""
        rxn = rxn_data["rxn"]
        assert equation == rxn.equation
        assert rate_type == rxn.rate.type

    def test_forward_rop_ddCi(self, surf, gas, rxn_data, orders, rxn_idx):
        """Test forward rates of progress derivatives."""
        drop = surf.forward_rates_of_progress_ddCi
        rop = surf.forward_rates_of_progress
        concentrations = self.get_concentrations(surf, gas)
        specs = surf.species_names + gas.species_names

        rix = rxn_data["rix"]
        r_stoich = rxn_data["r_stoich"]

        for spc_ix in rix:
            if orders is None:
                order = r_stoich[spc_ix, rxn_idx]
            else:
                order = orders.get(specs[spc_ix], 1)
            expected_rop = drop[rxn_idx, spc_ix] * concentrations[spc_ix] / order
            assert rop[rxn_idx] == pytest.approx(expected_rop, rel=self.rtol), (
                f"Forward ROP derivative mismatch for species index {spc_ix}"
            )

    def test_reverse_rop_ddCi(self, surf, gas, rxn_data, orders, rxn_idx):
        """Test reverse rates of progress derivatives."""
        drop = surf.reverse_rates_of_progress_ddCi
        rop = surf.reverse_rates_of_progress
        concentrations = self.get_concentrations(surf, gas)
        specs = surf.species_names + gas.species_names

        pix = rxn_data["pix"]
        p_stoich = rxn_data["p_stoich"]

        for spc_ix in pix:
            if orders is None:
                order = p_stoich[spc_ix, rxn_idx]
            else:
                order = orders.get(specs[spc_ix], 1)
            expected_rop = drop[rxn_idx, spc_ix] * concentrations[spc_ix] / order
            assert rop[rxn_idx] == pytest.approx(expected_rop, rel=self.rtol), (
                f"Reverse ROP derivative mismatch for species index {spc_ix}"
            )

    def test_net_rop_ddCi(self, surf, gas, rxn_data, rxn_idx, orders):
        """Test net rates of progress derivatives."""
        rop = surf.net_rates_of_progress
        drop = surf.net_rates_of_progress_ddCi
        concentrations = self.get_concentrations(surf, gas)

        drop_result = drop @ concentrations

        # Adjust forward and reverse rates by reaction orders
        ropf = surf.forward_rates_of_progress.copy()
        ropr = surf.reverse_rates_of_progress.copy()

        rxn = rxn_data["rxn"]
        orders_rxn = rxn.orders

        # Calculate total reactant orders
        total_orders_fwd = sum(orders_rxn.get(k, v) for k, v in rxn.reactants.items())
        ropf[rxn_idx] *= total_orders_fwd

        # Calculate total product orders
        total_orders_rev = sum(orders_rxn.get(k, v) for k, v in rxn.products.items())
        ropr[rxn_idx] *= total_orders_rev

        expected_net_rop = ropf - ropr

        assert drop_result[rxn_idx] == pytest.approx(expected_net_rop[rxn_idx], rel=self.rtol), (
            f"Net ROP derivative mismatch for reaction index {rxn_idx}"
        )


class PlatinumHydrogen(SurfaceRateExpressionTests):
    """
    Derived test class for Platinum-Hydrogen system.
    Provides specific fixtures for gas and surface phases.
    """

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

    @pytest.fixture(scope='class')
    def gas(self):
        """Create and configure the gas phase."""
        gas = ct.Solution(yaml=self.phase_defs, name="gas")
        gas.TPX = 800, 2 * ct.one_atm, "H2:1.5, O2:1.0, H2O2:0.75, H2O:0.3"
        return gas

    @pytest.fixture(scope='class')
    def surf(self, gas):
        """Create and configure the surface phase."""
        surf = ct.Interface(yaml=self.phase_defs, name="Pt_surf", adjacent=[gas])
        surf.TPX = 800, 2 * ct.one_atm, "PT(S):4.0, H(S):0.5, H2O(S):0.1, OH(S):0.2, O(S):0.8"
        return surf

class TestSurfInterfaceArrhenius(PlatinumHydrogen):
    """Test case for interface-Arrhenius reaction H(S) + O(S) <=> OH(S) + PT(S)."""

    @pytest.fixture(scope='class')
    def rxn_idx(self):
        return 7

    @pytest.fixture(scope='class')
    def equation(self):
        return "H(S) + O(S) <=> OH(S) + PT(S)"

    @pytest.fixture(scope='class')
    def rate_type(self):
        return "interface-Arrhenius"

    @pytest.fixture(scope='class')
    def orders(self):
        return None

class TestSurfGasFwdStickingArrhenius(PlatinumHydrogen):
    """Test case for sticking-Arrhenius reaction H2O + PT(S) => H2O(S)."""

    @pytest.fixture(scope='class')
    def rxn_idx(self):
        return 5

    @pytest.fixture(scope='class')
    def equation(self):
        return "H2O + PT(S) => H2O(S)"

    @pytest.fixture(scope='class')
    def rate_type(self):
        return "sticking-Arrhenius"

    @pytest.fixture(scope='class')
    def orders(self):
        return None

class TestSurfGasInterfaceArrhenius(PlatinumHydrogen):
    """Test case for interface-Arrhenius reaction H2 + 2 PT(S) => 2 H(S)."""

    @pytest.fixture(scope='class')
    def rxn_idx(self):
        return 0

    @pytest.fixture(scope='class')
    def equation(self):
        return "H2 + 2 PT(S) => 2 H(S)"

    @pytest.fixture(scope='class')
    def rate_type(self):
        return "interface-Arrhenius"

    @pytest.fixture(scope='class')
    def orders(self):
        return {"PT(S)": 1, "H2": 1, "H(S)": 2}

class TestGasSurfInterfaceArrhenius(PlatinumHydrogen):
    """Test case for interface-Arrhenius reaction H2O(S) => H2O + PT(S)."""

    @pytest.fixture(scope='class')
    def rxn_idx(self):
        return 6

    @pytest.fixture(scope='class')
    def equation(self):
        return "H2O(S) => H2O + PT(S)"

    @pytest.fixture(scope='class')
    def rate_type(self):
        return "interface-Arrhenius"

    @pytest.fixture(scope='class')
    def orders(self):
        return None


class SurfaceFullTests:
    # Generic test class to check derivatives evaluated for an entire reaction
    # mechanisms
    rtol = 1e-4

    @pytest.fixture(scope='class')
    def all_species(self, gas, surf):
        """Combine species from both gas and surface phases."""
        return surf.species() + gas.species()

    @pytest.fixture(scope='class')
    def sidxs(self, all_species):
        """Create a dictionary mapping species names to indices."""
        return {spec.name: i for i, spec in enumerate(all_species)}

    @pytest.fixture(scope='class')
    def initial_conditions(self, gas, surf):
        """Store initial TPX conditions for gas and surface."""
        return {
            "gas_tpx": gas.TPX,
            "surf_tpx": surf.TPX
        }

    @pytest.fixture(scope='function', autouse=True)
    def setup_surface_full_data(self, gas, surf, initial_conditions):
        """Reset TPX and derivative settings before each test."""
        # Initialize gas phase
        gas.TPX = initial_conditions["gas_tpx"]
        gas.derivative_settings = {}

        # Initialize surface phase
        surf.TPX = initial_conditions["surf_tpx"]
        surf.derivative_settings = {
            "skip-coverage-dependence": True,
            "skip-electrochemistry": True
        }

    def get_concentrations(self, surf, gas):
        """Concatenate concentrations from surface and gas phases."""
        return np.concatenate((surf.concentrations, gas.concentrations))

    def test_forward_rop_ddCi(self, gas, surf):
        # matrix multiplication of the forward rates of progress derivatives w.r.t
        # concentration and the concentrations should provide the rate of progress
        # for each species and can be compared to the directly calculated rate
        drop = surf.forward_rates_of_progress_ddCi
        rop = surf.forward_rates_of_progress
        conc = self.get_concentrations(surf, gas)
        # multiply derivatives with concentrations
        drop = drop @ conc
        # get total reactant reaction orders
        total_orders = []
        for rxn in surf.reactions():
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

    def test_reverse_rop_ddCi(self, gas, surf):
        # matrix multiplication of the reverse rate of progress derivatives  w.r.t
        # concentration and the concentrations should provide the rate of progress
        # for each species and can be compared to the directly calculated rate
        drop = surf.reverse_rates_of_progress_ddCi
        rop = surf.reverse_rates_of_progress
        conc = self.get_concentrations(surf, gas)
        # multiply derivatives with concentrations
        drop = drop @ conc
        # get total reactant reaction orders
        total_orders = []
        for rxn in surf.reactions():
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

    def test_net_rop_ddCi(self, gas, surf):
        # check derivatives of net rates of progress with respect to species
        # concentrations against analytic
        ropf = surf.forward_rates_of_progress
        ropr = surf.reverse_rates_of_progress
        drop = surf.net_rates_of_progress_ddCi
        conc = self.get_concentrations(surf, gas)
        # multiply derivatives with concentrations
        drop = drop @ conc
        # reaction orders are not yet accounted for in rates of progress
        # so they must be included manually
        for i, rxn in enumerate(surf.reactions()):
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

    @pytest.fixture(scope='class')
    def gas(self):
        """Fixture to create and configure the gas phase."""
        gas = ct.Solution(yaml=self.phase_defs, name="gas")
        gas.TPX = 800, 2 * ct.one_atm, "H2:1.5, O2:1.0, H2O2:0.75, H2O:0.3"
        return gas

    @pytest.fixture(scope='class')
    def surf(self, gas):
        """Fixture to create and configure the surface phase."""
        surf = ct.Interface(yaml=self.phase_defs, name="Pt_surf", adjacent=[gas])
        surf.TPX = 800, 2 * ct.one_atm, "PT(S):1.0, H(S):0.5, H2O(S):0.1, OH(S):0.2, O(S):0.8"
        return surf
