import importlib.metadata
import itertools
import numpy as np
import pytest
from pytest import approx
import re

import cantera as ct

from .utilities import (
    compare
)

# avoid explicit dependence of cantera on scipy
try:
    importlib.metadata.version('scipy')
except importlib.metadata.PackageNotFoundError:
    _scipy_sparse = ImportError('Method requires a working scipy installation.')
else:
    from scipy import sparse as _scipy_sparse


class TestKinetics:

    @pytest.fixture
    def phase(self):
        phase = ct.Solution('h2o2.yaml', transport_model=None)
        phase.TPX = 800, 2 * ct.one_atm, [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4, 0]
        return phase

    def test_counts(self, phase):
        assert phase.n_reactions == 29
        assert phase.n_total_species == 10
        assert phase.n_phases == 1

    def test_is_reversible(self, phase):
        for i in range(phase.n_reactions):
            assert phase.reaction(i).reversible

    def test_multiplier(self, phase):
        fwd_rates0 = phase.forward_rates_of_progress
        rev_rates0 = phase.reverse_rates_of_progress

        phase.set_multiplier(2.0, 0)
        phase.set_multiplier(0.1, 6)

        fwd_rates1 = phase.forward_rates_of_progress
        rev_rates1 = phase.reverse_rates_of_progress

        assert 2 * fwd_rates0[0] == approx(fwd_rates1[0])
        assert 0.1 * fwd_rates0[6] == approx(fwd_rates1[6])
        assert 2 * rev_rates0[0] == approx(rev_rates1[0])
        assert 0.1 * rev_rates0[6] == approx(rev_rates1[6])
        for i in range(phase.n_reactions):
            if i not in (0, 6):
                assert fwd_rates0[i] == approx(fwd_rates1[i])
                assert rev_rates0[i] == approx(rev_rates1[i])

        phase.set_multiplier(0.5)
        fwd_rates2 = phase.forward_rates_of_progress
        rev_rates2 = phase.reverse_rates_of_progress
        assert 0.5 * fwd_rates0 == approx(fwd_rates2)
        assert 0.5 * rev_rates0 == approx(rev_rates2)

    def test_legacy_reaction_rate(self, phase):
        ct.use_legacy_rate_constants(True)
        fwd_rates_legacy = phase.forward_rate_constants
        ct.use_legacy_rate_constants(False)
        fwd_rates = phase.forward_rate_constants
        ix_3b = np.array([r.reaction_type == "three-body-Arrhenius" for r in phase.reactions()])
        ix_other = ix_3b == False

        assert fwd_rates_legacy[ix_other] == approx(fwd_rates[ix_other])
        assert not (fwd_rates_legacy[ix_3b] == fwd_rates[ix_3b]).any()

    def test_reaction_type(self, phase):
        assert phase.reaction(0).reaction_type == "three-body-Arrhenius"
        assert phase.reaction(2).reaction_type == "Arrhenius"
        assert phase.reaction(2).rate.type == "Arrhenius"
        assert phase.reaction(21).reaction_type == "falloff-Troe"

        with pytest.raises(ct.CanteraError, match="outside valid range"):
            phase.reaction(33).reaction_type
        with pytest.raises(ct.CanteraError, match="outside valid range"):
            phase.reaction(-2).reaction_type

    def test_reaction_equations(self, phase):
        assert phase.n_reactions == len(phase.reaction_equations())
        r, p = [x.split() for x in phase.reaction(18).equation.split('<=>')]
        assert 'H' in r
        assert 'H2O2' in r
        assert 'HO2' in p
        assert 'H2' in p

    def test_reactants_products(self, phase):
        for i in range(phase.n_reactions):
            R = phase.reaction(i).reactant_string
            P = phase.reaction(i).product_string
            assert phase.reaction(i).equation.startswith(R)
            assert phase.reaction(i).equation.endswith(P)
            for k in range(phase.n_species):
                if phase.reactant_stoich_coeff(k, i) != 0:
                    assert phase.species_name(k) in R
                if phase.product_stoich_coeff(k, i) != 0:
                    assert phase.species_name(k) in P

    def test_stoich_coeffs(self, phase):
        nu_r = phase.reactant_stoich_coeffs
        nu_p = phase.product_stoich_coeffs

        def check_reactant(s, i, value):
            k = phase.kinetics_species_index(s)
            assert phase.reactant_stoich_coeff(s, i) == value
            assert phase.reactant_stoich_coeff(k, i) == value
            assert nu_r[k, i] == value

        def check_product(s, i, value):
            k = phase.kinetics_species_index(s)
            assert phase.product_stoich_coeff(k, i) == value
            assert phase.product_stoich_coeff(s, i) == value
            assert nu_p[k, i] == value

        # H + H2O2 <=> HO2 + H2
        check_reactant('H', 18, 1)
        check_reactant('H2O2', 18, 1)
        check_reactant('HO2', 18, 0)
        check_reactant('H2', 18, 0)

        check_product('H', 18, 0)
        check_product('H2O2', 18, 0)
        check_product('HO2', 18, 1)
        check_product('H2', 18, 1)

        # 2 O + M <=> O2 + M
        check_reactant('O', 0, 2)
        check_reactant('O2', 0, 0)
        check_product('O', 0, 0)
        check_product('O2', 0, 1)

    @pytest.mark.skipif(isinstance(_scipy_sparse, ImportError), reason="scipy is not installed")
    def test_stoich_coeffs_sparse(self, phase):
        nu_r_dense = phase.reactant_stoich_coeffs
        nu_p_dense = phase.product_stoich_coeffs

        ct.use_sparse(True)
        nu_r_sparse = phase.reactant_stoich_coeffs
        nu_p_sparse = phase.product_stoich_coeffs

        assert (nu_r_sparse.toarray() == nu_r_dense).all()
        assert (nu_p_sparse.toarray() == nu_p_dense).all()

        ct.use_sparse(False)

    def test_rates_of_progress(self, phase):
        assert len(phase.net_rates_of_progress) == phase.n_reactions
        assert (phase.forward_rates_of_progress - phase.reverse_rates_of_progress
                == approx(phase.net_rates_of_progress))

    def test_heat_release(self, phase):
        hrr = - phase.partial_molar_enthalpies.dot(phase.net_production_rates)
        assert hrr == approx(phase.heat_release_rate)
        assert hrr == approx(sum(phase.heat_production_rates))

    def test_rate_constants(self, phase):
        assert len(phase.forward_rate_constants) == phase.n_reactions
        ix = phase.reverse_rate_constants != 0.
        assert (phase.forward_rate_constants[ix] / phase.reverse_rate_constants[ix]
                == approx(phase.equilibrium_constants[ix]))

    def test_species_rates(self, phase):
        nu_p = phase.product_stoich_coeffs
        nu_r = phase.reactant_stoich_coeffs
        creation = (np.dot(nu_p, phase.forward_rates_of_progress)
                    + np.dot(nu_r, phase.reverse_rates_of_progress))
        destruction = (np.dot(nu_r, phase.forward_rates_of_progress)
                       + np.dot(nu_p, phase.reverse_rates_of_progress))

        assert phase.creation_rates == approx(creation)
        assert phase.destruction_rates == approx(destruction)
        assert phase.net_production_rates == approx(creation - destruction)

    def test_reaction_deltas(self, phase):
        assert (phase.delta_enthalpy - phase.delta_entropy * phase.T
                == approx(phase.delta_gibbs))

        assert (phase.delta_standard_enthalpy - phase.delta_standard_entropy * phase.T
                == approx(phase.delta_standard_gibbs))


class TestKineticsFromReactions:
    """
    Test for Kinetics objects which are constructed directly from Reaction
    objects instead of from input files.
    """
    def test_idealgas(self):
        gas1 = ct.Solution("h2o2.yaml")

        S = ct.Species.list_from_file("h2o2.yaml")
        R = ct.Reaction.list_from_file("h2o2.yaml", gas1)
        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=S, reactions=R)

        assert gas1.n_reactions == gas2.n_reactions
        gas1.TPY = 800, 2*ct.one_atm, 'H2:0.3, O2:0.7, OH:2e-4, O:1e-3, H:5e-5'
        gas2.TPY = gas1.TPY

        assert (gas1.reactant_stoich_coeffs == gas2.reactant_stoich_coeffs).all()
        assert (gas1.product_stoich_coeffs == gas2.product_stoich_coeffs).all()

        assert gas1.delta_gibbs == approx(gas2.delta_gibbs)
        assert gas1.reverse_rate_constants == approx(gas2.reverse_rate_constants)
        assert gas1.net_production_rates == approx(gas2.net_production_rates)

    def test_surface(self):
        surf1 = ct.Interface("ptcombust.yaml", "Pt_surf")
        gas = surf1.adjacent["gas"]
        surf_species = ct.Species.list_from_file("ptcombust.yaml")
        reactions = ct.Reaction.list_from_file("ptcombust.yaml", surf1)

        surf2 = ct.Interface(thermo='ideal-surface', kinetics='surface',
                             species=surf_species, reactions=reactions,
                             adjacent=[gas])
        surf1.site_density = surf2.site_density = 5e-9
        gas.TP = surf2.TP = surf1.TP = 900, 2*ct.one_atm
        surf2.coverages = surf1.coverages

        assert surf1.n_reactions == surf2.n_reactions

        for k, i in itertools.product(['PT(S)', 'H2', 'OH', 'OH(S)'], range(surf1.n_species)):
            assert surf1.reactant_stoich_coeff(k, i) == surf2.reactant_stoich_coeff(k, i)
            assert surf1.product_stoich_coeff(k, i) == surf2.product_stoich_coeff(k, i)

        for i in range(surf1.n_reactions):
            r1 = surf1.reaction(i)
            r2 = surf2.reaction(i)
            assert r1.reactants == r2.reactants
            assert r1.products == r2.products
            assert r1.rate.pre_exponential_factor == r2.rate.pre_exponential_factor
            assert r1.rate.temperature_exponent == r2.rate.temperature_exponent
            assert r1.rate.activation_energy == approx(r2.rate.activation_energy)

        assert surf1.delta_enthalpy == approx(surf2.delta_enthalpy)
        assert surf1.forward_rate_constants == approx(surf2.forward_rate_constants)
        assert surf1.reverse_rate_constants == approx(surf2.reverse_rate_constants)

        rop1 = surf1.net_production_rates
        rop2 = surf2.net_production_rates
        for k in gas.species_names + surf1.species_names:
            k1 = surf1.kinetics_species_index(k)
            k2 = surf2.kinetics_species_index(k)
            assert rop1[k1] == approx(rop2[k2])

    def test_add_reaction(self):
        gas1 = ct.Solution('h2o2.yaml', transport_model=None)

        S = ct.Species.list_from_file("h2o2.yaml")
        R = ct.Reaction.list_from_file("h2o2.yaml", gas1)
        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=S, reactions=R[:5])

        gas1.TPY = 800, 2*ct.one_atm, 'H2:0.3, O2:0.7, OH:2e-4, O:1e-3, H:5e-5'
        gas2.TPY = gas1.TPY

        for r in R[5:]:
            gas2.add_reaction(r)

        assert gas1.n_reactions == gas2.n_reactions

        assert (gas1.reactant_stoich_coeffs == gas2.reactant_stoich_coeffs).all()
        assert (gas1.product_stoich_coeffs == gas2.product_stoich_coeffs).all()

        assert gas1.delta_gibbs == approx(gas2.delta_gibbs)
        assert gas1.reverse_rate_constants == approx(gas2.reverse_rate_constants)
        assert gas1.net_production_rates == approx(gas2.net_production_rates)

    def test_coverage_dependence_flags(self):
        surf = ct.Interface("ptcombust.yaml", "Pt_surf")
        surf.TP = 900, ct.one_atm
        surf.coverages = {"PT(S)": 1}
        with pytest.raises(NotImplementedError):
            surf.net_rates_of_progress_ddCi
        # set skip and try to get jacobian again
        surf.derivative_settings = {"skip-coverage-dependence": True}
        surf.net_rates_of_progress_ddCi

    def test_electrochemistry_flags(self):
        anode_int = ct.Interface("lithium_ion_battery.yaml", "edge_anode_electrolyte")
        with pytest.raises(NotImplementedError):
            anode_int.net_rates_of_progress_ddCi
        # set skip and try to get jacobian again
        anode_int.derivative_settings = {"skip-electrochemistry": True}
        anode_int.net_rates_of_progress_ddCi

    def test_submechanism(self):
        # Simplified samples/python/kinetics/extract_submechanism.py
        gri30 = ct.Solution('gri30.yaml', transport_model=None)
        h2o2 = ct.Solution('h2o2.yaml', transport_model=None)

        reactions_plus = []
        reactions = []
        dest_species = set(h2o2.species_names)
        colliders = dest_species.union([None, "M"])
        for R in gri30.reactions():
            if not all(S in dest_species for S in R.reactants):
                continue
            if not all(S in dest_species for S in R.products):
                continue
            reactions_plus.append(R)
            if R.third_body_name not in colliders:
                continue
            reactions.append(R)

        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=h2o2.species(), reactions=reactions_plus)
        # there is one third-body reaction with an undeclared third body species
        assert gas.n_reactions < len(reactions_plus)
        assert gas.n_species == len(h2o2.species_names)

        gas = ct.Solution(thermo='ideal-gas', kinetics='gas', species=h2o2.species(),
                          reactions=reactions)
        assert gas.n_reactions == len(reactions)
        assert gas.n_species == len(h2o2.species_names)

        yaml_file = self.test_work_path / "reduced.yaml"
        gas.write_yaml(yaml_file)
        restored = ct.Solution(yaml_file)
        assert gas.species_names == restored.species_names
        assert gas.reaction_equations() == restored.reaction_equations()


class TestKineticsRepeatability:
    """
    Tests to make sure that lazily evaluated terms in the rate expression
    are always updated correctly.
    """
    T0 = 1200
    T1 = 1300
    rho0 = 2.4
    rho1 = 3.1
    P0 = 1.4e5
    P1 = 3.7e6

    def setup_gas(self, mech):
        gas = ct.Solution(mech)
        self.X0 = 1 + np.sin(range(1, gas.n_species + 1))
        self.X1 = 1 + np.sin(range(2, gas.n_species + 2))
        return gas

    def check_rates_composition(self, mech):
        gas = self.setup_gas(mech)
        gas.TDX = self.T0, self.rho0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TDX = self.T1, self.rho1, self.X1
        w2 = gas.net_production_rates

        gas.TDX = self.T0, self.rho0, self.X1
        w3 = gas.net_production_rates

        # change only composition, and make sure the rates match
        gas.TDX = self.T0, self.rho0, self.X0
        w4 = gas.net_production_rates

        assert w1 == approx(w4)

    def check_rates_temperature1(self, mech):
        gas = self.setup_gas(mech)
        gas.TDX = self.T0, self.rho0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TDX = self.T1, self.rho1, self.X1
        w2 = gas.net_production_rates

        gas.TDX = self.T1, self.rho0, self.X0
        w3 = gas.net_production_rates

        # change only temperature, and make sure the rates match
        gas.TDX = self.T0, self.rho0, self.X0
        w4 = gas.net_production_rates

        assert w1 == approx(w4)

    def check_rates_temperature2(self, mech):
        gas = self.setup_gas(mech)
        gas.TPX = self.T0, self.P0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TPX = self.T1, self.P1, self.X1
        w2 = gas.net_production_rates

        gas.TPX = self.T1, self.P0, self.X0
        w3 = gas.net_production_rates

        # change only temperature, and make sure the rates match
        gas.TPX = self.T0, self.P0, self.X0
        w4 = gas.net_production_rates

        assert w1 == approx(w4)

    def check_rates_pressure(self, mech):
        gas = self.setup_gas(mech)
        gas.TPX = self.T0, self.P0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TPX = self.T1, self.P1, self.X1
        w2 = gas.net_production_rates

        gas.TPX = self.T0, self.P1, self.X0
        w3 = gas.net_production_rates

        # change only pressure, and make sure the rates match
        gas.TPX = self.T0, self.P0, self.X0
        w4 = gas.net_production_rates

        assert w1 == approx(w4)

    def test_gri30_composition(self):
        self.check_rates_composition("gri30.yaml")

    def test_gri30_temperature(self):
        self.check_rates_temperature1("gri30.yaml")
        self.check_rates_temperature2("gri30.yaml")

    def test_gri30_pressure(self):
        self.check_rates_pressure("gri30.yaml")

    def test_h2o2_composition(self):
        self.check_rates_composition("h2o2.yaml")

    def test_h2o2_temperature(self):
        self.check_rates_temperature1("h2o2.yaml")
        self.check_rates_temperature2("h2o2.yaml")

    def test_h2o2_pressure(self):
        self.check_rates_pressure("h2o2.yaml")

    def test_pdep_composition(self):
        self.check_rates_composition('pdep-test.yaml')

    def test_pdep_temperature(self):
        self.check_rates_temperature1('pdep-test.yaml')
        self.check_rates_temperature2('pdep-test.yaml')

    def test_pdep_pressure(self):
        self.check_rates_pressure('pdep-test.yaml')

    def test_modify_thermo(self):
        # Make sure that thermo modifications propagate through to Kinetics
        gas = self.setup_gas("gri30.yaml")
        gas.TPX = self.T0, self.P0, self.X0
        gas.equilibrate('TP')
        gas.TP = gas.T + 20, None

        S = {sp.name: sp for sp in ct.Species.list_from_file("gri30.yaml")}
        w1 = gas.net_rates_of_progress

        OH = gas.species('OH')
        OH.thermo = S['CO2'].thermo
        gas.modify_species(gas.species_index('OH'), OH)
        w2 = gas.net_rates_of_progress

        for i, R in enumerate(gas.reactions()):
            if ('OH' in R.reactants or 'OH' in R.products) and R.reversible:
                # Rate should be different if reaction involves OH
                assert not 1 == approx(w2[i] / w1[i], rel=1e-5)
            else:
                # Rate should be the same if reaction does not involve OH
                assert 1 == approx(w2[i] / w1[i], rel=1e-5)

    def test_pdep_err(self):
        err_msg = ("InputFileError thrown by PlogRate::validate:",
                   "Invalid rate coefficient for reaction 'CH2CHOO <=> CH3O + CO'",
                   "at P = 32019, T = 500.0",
                   "at P = 32019, T = 1000.0",
                   "at P = 1.0132e+05, T = 500.0",
                   "at P = 1.0132e+05, T = 1000.0")

        if not ct.debug_mode_enabled():
            err_msg += (
                "Invalid rate coefficient for reaction 'CH2CHOO <=> CH3 + CO2'",
                "at P = 1.0132e+05, T = 500.0",
                "InputFileError thrown by Reaction::checkBalance:",
                "The following reaction is unbalanced: H2O2 + OH <=> 2 H2O + HO2"
            )
        with pytest.raises(ct.CanteraError) as e:
            ct.Solution('addReactions_err_test.yaml')
        err_msg_list = str(e.value).splitlines()
        for msg in err_msg:
            assert msg in err_msg_list

    def test_sticking_coeff_err(self):
        err_msg = (r"Sticking coefficient is greater than 1 for reaction 'O2 \+ 2 PT\(S\) => 2 O\(S\)'",
                   "at T = 200.0",
                   "at T = 500.0",
                   "at T = 1000.0",
                   "at T = 2000.0",
                   "at T = 5000.0",
                   "at T = 10000.0",
                   "Sticking coefficient is greater than 1 for reaction",
                   "StickingRate::validate:")

        for err in err_msg:
            with pytest.warns(UserWarning, match=err):
                ct.Interface("sticking_coeff_check.yaml", "Pt_surf")


def check_raises(yaml, err_msg, line):
    """Helper function for InputFileErrors"""
    err_msg = [err_msg]
    err_msg.append("InputFileError thrown by ")
    err_msg.append(f"Error on line {line} of ")
    with pytest.raises(ct.CanteraError) as ex:
        ct.Solution(yaml=yaml)
    msg = str(ex.value)
    for err in err_msg:
        assert err in msg


class TestUndeclared:

    _gas_def = """
            phases:
            - name: gas
              species:
              - h2o2.yaml/species: [H, O2, H2O, HO2]
              thermo: ideal-gas
              kinetics: gas
            """

    def test_raise_undeclared_species(self):

        gas_def = self._gas_def + """
            reactions:
            - equation: O + H2 <=> H + OH  # Reaction 3
              rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
            """
        check_raises(gas_def, "contains undeclared species", line=10)

    def test_raise_undeclared_third_bodies(self):

        gas_def = self._gas_def + """
            reactions:
            - equation: H + O2 + AR <=> HO2 + AR  # Reaction 10
              rate-constant: {A: 7.0e+17, b: -0.8, Ea: 0.0}
            """
        check_raises(gas_def, "three-body reaction with undeclared", line=10)

    def test_skip_undeclared_third_bodies1(self):

        gas_def = self._gas_def + """
              reactions:
              - h2o2.yaml/reactions: declared-species
              skip-undeclared-third-bodies: true
            """

        gas = ct.Solution(yaml=gas_def)
        assert gas.n_reactions == 3

    def test_skip_undeclared_third_bodies2(self):

        gas_def = """
            phases:
            - name: gas
              species:
              - h2o2.yaml/species: [H, O2, HO2]
              thermo: ideal-gas
              skip-undeclared-third-bodies: true
              kinetics: gas
              reactions:
              - h2o2.yaml/reactions: declared-species
            """

        gas = ct.Solution(yaml=gas_def)
        found = False
        for rxn in gas.reactions():
            if rxn.equation == "H + O2 + M <=> HO2 + M":
                found = True
                break
        assert found

    def test_skip_undeclared_orders(self):

        gas_def = self._gas_def + """
            reactions:
            - equation: H + O2 => HO2
              rate-constant: {A: 1.255943e+13, b: -2.0, Ea: 5000.0 cal/mol}
              orders:
                H2O: 0.2
              nonreactant-orders: true
            """

        gas = ct.Solution(yaml=gas_def)
        assert gas.n_reactions == 1

    def test_raise_nonreactant_orders(self):

        gas_def = self._gas_def + """
            reactions:
            - equation: H + O2 => HO2
              rate-constant: {A: 1.255943e+13, b: -2.0, Ea: 5000.0 cal/mol}
              orders:
                H2O: 0.2
            """
        check_raises(gas_def, "Reaction order specified", line=10)

    def test_raise_undeclared_orders(self):

        gas_def = self._gas_def + """
            reactions:
            - equation: H + O2 => HO2
              rate-constant: {A: 1.255943e+13, b: -2.0, Ea: 5000.0 cal/mol}
              orders:
                N2: 0.2
              nonreactant-orders: true
            """
        check_raises(gas_def, "reaction orders for undeclared", line=13)

    def test_skip_undeclared_surf_species(self):
        phase_defs = """
            phases:
            - name: gas
              thermo: ideal-gas
              species:
              - gri30.yaml/species: [H2, H, O, OH, H2O, CO2]
            - name: Pt_surf
              thermo: ideal-surface
              species:
              - ptcombust.yaml/species: [PT(S), H(S), H2O(S), OH(S), CO2(S), CH2(S)s,
                  CH(S), C(S), O(S)]
              kinetics: surface
              reactions: [ptcombust.yaml/reactions: declared-species]
              site-density: 2.7063e-09
            """
        gas = ct.Solution(yaml=phase_defs, name="gas")
        surf = ct.Interface(yaml=phase_defs, name="Pt_surf", adjacent=[gas])
        assert surf.n_reactions == 14


class TestInvalidInput:

    _gas_def = """
            phases:
            - name: gas
              thermo: ideal-gas
              species: [{h2o2.yaml/species: all}]
              kinetics: gas
              reactions: all
              state: {T: 300.0, P: 1 atm, X: {O2: 0.21, N2: 0.79}}
        """

    def test_failing_convert1(self):
        # invalid preexponential factor units
        gas_def = self._gas_def + """
            reactions:
            - equation: O + H2 <=> H + OH  # Reaction 3
              rate-constant: {A: 3.87e+04 cm^6/mol^2/s, b: 2.7, Ea: 6260.0}
            """
        check_raises(gas_def, "UnitSystem::convert:", line=12)

    def test_failing_convert2(self):
        # invalid activation energy units
        gas_def = self._gas_def + """
            reactions:
            - equation: O + H2 <=> H + OH  # Reaction 3
              rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0 m}
            """
        check_raises(gas_def, "UnitSystem::convertActivationEnergy:", line=12)

    def test_failing_unconfigured1(self):
        # missing reaction rate definition
        gas_def = self._gas_def + """
            reactions:
            - equation: O + H2 <=> H + OH
            """
        check_raises(gas_def, "is not configured", line=11)

    def test_failing_unconfigured2(self):
        # missing reaction rate definition
        gas_def = self._gas_def + """
            reactions:
            - equation: 2 OH (+M) <=> H2O2 (+M)
              type: falloff
            """
        check_raises(gas_def, "is not configured", line=11)

    def test_failing_unconfigured3(self):
        # missing reaction rate definition
        gas_def = self._gas_def + """
            reactions:
            - equation: O + H2 <=> H + OH
              type: pressure-dependent-Arrhenius
            """
        check_raises(gas_def, "is not configured", line=11)

    def test_failing_invalid_duplicate(self):
        # Make sure that the error we get is for the invalid reaction, not for the
        # missing duplicate
        gas_def = self._gas_def + """
            reactions:
            - equation: O + H2 <=> H + OH
              rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
              duplicate: true
            - equation: O + H2 <=> H + OH
              rate-constant: {A: -3.87e+04, b: 2.7, Ea: 6260.0}
              duplicate: true
            """
        check_raises(gas_def, "negative pre-exponential factor", line=14)


class TestEmptyKinetics:

    def test_empty(self):
        gas = ct.Solution("air-no-reactions.yaml")

        assert gas.n_reactions == 0
        assert gas.creation_rates == approx(np.zeros(gas.n_species))
        assert gas.destruction_rates == approx(np.zeros(gas.n_species))
        assert gas.net_production_rates == approx(np.zeros(gas.n_species))


class TestReactionPath:

    @pytest.fixture(scope='class')
    def gas(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 1300.0, ct.one_atm, 'CH4:0.4, O2:1, N2:3.76'

        # Advance the reactor
        r = ct.IdealGasReactor(gas)
        net = ct.ReactorNet([r])
        T = r.T
        while T < 1900:
            net.step()
            T = r.T
        return gas

    def check_dot(self, gas, diagram, element):
        diagram.label_threshold = 0
        diagram.threshold = 0
        dot = diagram.get_dot()
        dot = dot.replace('\n', '')
        nodes1 = set()
        nodes2 = set()
        species = set()
        for line in dot.split(';'):
            m = re.match(r'(.*)\[(.*)\]', line)
            if not m:
                continue
            A, B = m.groups()
            if '->' in A:
                # edges
                nodes1.update(s.strip() for s in A.split('->'))
            else:
                # nodes
                nodes2.add(A.strip())
                spec = re.search('label="(.*?)"', B).group(1)
                assert spec not in species
                species.add(spec)

        # Make sure that the output was actually parsable and that we
        # found some nodes
        assert nodes1
        assert species

        # All nodes should be connected to some edge (this does not
        # require the graph to be connected)
        assert nodes1 == nodes2

        # All of the species in the graph should contain the element whose
        # flux we're looking at
        for spec in species:
            assert gas.n_atoms(spec, element) > 0

        # return fluxes from the dot file for further tests
        return [float(re.search('label *= *"(.*?)"', line).group(1))
                for line in dot.split(';')
                if line.startswith('s') and 'arrowsize' in line]

    def get_fluxes(self, diagram):
        directional = {}
        net = {}
        for line in diagram.get_data().strip().split('\n')[1:]:
            s = line.split()
            fwd = float(s[2])
            rev = float(s[3])
            directional[s[0], s[1]] = fwd
            directional[s[1], s[0]] = -rev
            if fwd + rev > 0:
                net[s[0], s[1]] = fwd + rev
            else:
                net[s[1], s[0]] = - fwd - rev
        return directional, net

    def test_dot_net_autoscaled(self, gas):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(gas, element)
            dot_fluxes = self.check_dot(gas, diagram, element)
            assert max(dot_fluxes) == 1.0

    def test_dot_net_unscaled(self, gas):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(gas, element)
            diagram.scale = 1.0
            dot_fluxes = sorted(self.check_dot(gas, diagram, element))
            _, fluxes = self.get_fluxes(diagram)
            fluxes = sorted(fluxes.values())

            for i in range(1, 20):
                assert dot_fluxes[-i] == approx(fluxes[-i], rel=1e-2)

    def test_dot_oneway_autoscaled(self, gas):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(gas, element)
            diagram.flow_type = 'OneWayFlow'
            dot_fluxes = self.check_dot(gas, diagram, element)
            assert max(dot_fluxes) == 1.0

    def test_dot_oneway_unscaled(self, gas):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(gas, element)
            diagram.scale = 1.0
            diagram.flow_type = 'OneWayFlow'
            dot_fluxes = sorted(self.check_dot(gas, diagram, element))
            fluxes, _ = self.get_fluxes(diagram)
            fluxes = sorted(fluxes.values())

            for i in range(1, 20):
                assert dot_fluxes[-i] == approx(fluxes[-i], rel=1e-2)

    def test_fluxes(self):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        gas.TPX = 1100, 10*ct.one_atm, 'H:0.1, HO2:0.1, AR:6'
        diagram = ct.ReactionPathDiagram(gas, 'H')
        ropf = gas.forward_rates_of_progress
        ropr = gas.reverse_rates_of_progress
        fluxes, _ = self.get_fluxes(diagram)

        assert fluxes['HO2','H'] == approx(ropr[5] + ropr[9], rel=1e-5)
        assert fluxes['H', 'H2'] == approx(2*ropf[11] + ropf[16], rel=1e-5)
        assert fluxes['H', 'H2O'] == approx(ropf[15], rel=1e-5)
        assert fluxes['H', 'OH'] == approx(ropf[17], rel=1e-5)
        assert fluxes['HO2','H2'] == approx(ropf[16], rel=1e-5)
        assert fluxes['HO2', 'H2O'] == approx(ropf[15], rel=1e-5)
        assert fluxes['HO2', 'OH'] == approx(ropf[17], rel=1e-5)
        assert fluxes['HO2', 'H2O2'] == approx(2*ropf[26] + 2*ropf[27], rel=1e-5)


class TestChemicallyActivated:
    def test_rate_evaluation(self):
        gas = ct.Solution("chemically-activated-reaction.yaml")
        P = [2026.5, 202650.0, 10132500.0] # pressure

        # forward rate of progress, computed using Chemkin
        Rf = [2.851022e-04, 2.775924e+00, 2.481792e+03]

        for i in range(len(P)):
            gas.TPX = 900.0, P[i], [0.01, 0.01, 0.04, 0.10, 0.84]
            assert gas.forward_rates_of_progress[0] == approx(Rf[i], rel=2e-5)


class TestExplicitForwardOrder:

    @pytest.fixture
    def gas(self):
        gas = ct.Solution("explicit-forward-order.yaml")
        gas.TPX = 800, 101325, [0.01, 0.90, 0.02, 0.03, 0.04]
        return gas

    def test_irreversibility(self, gas):
        # Reactions are irreversible
        Rr = gas.reverse_rate_constants
        for i in range(3):
            assert Rr[i] == 0.0

    def test_rateConstants(self, gas):
        # species order: [H, AR, R1A, R1B, P1]
        C = gas.concentrations
        Rf = gas.forward_rates_of_progress
        kf = gas.forward_rate_constants
        assert Rf[0] == approx(kf[0] * C[2]**1.5 * C[3]**0.5)
        assert Rf[1] == approx(kf[1] * C[0]**1.0 * C[4]**0.2)
        assert Rf[2] == approx(kf[2] * C[2]**3.0)

    def test_ratio1(self, gas):
        rop1 = gas.forward_rates_of_progress
        # Double concentration of H and R1A
        gas.TPX = None, None, [0.02, 0.87, 0.04, 0.03, 0.04]
        rop2 = gas.forward_rates_of_progress
        ratio = rop2/rop1
        assert ratio[0] == approx(2**1.5) # order of R1A is 1.5
        assert ratio[1] == approx(2**1.0) # order of H is 1.0
        assert ratio[2] == approx(2**3) # order of R1A is 3

    def test_ratio2(self, gas):
        rop1 = gas.forward_rates_of_progress
        # Double concentration of P1 and R1B
        gas.TPX = None, None, [0.01, 0.83, 0.02, 0.06, 0.08]
        rop2 = gas.forward_rates_of_progress
        ratio = rop2/rop1
        assert ratio[0] == approx(2**0.5) # order of R1B is 0.5
        assert ratio[1] == approx(2**0.2) # order of P1 is 1.0
        assert ratio[2] == approx(2**0.0) # order of R1B is 0


class TestSofcKinetics:
    """ Test based on sofc.py """
    _mech = "sofc.yaml"

    def test_sofc(self, test_data_path):
        mech = self._mech
        T = 1073.15  # T in K
        P = ct.one_atm
        TPB_length_per_area = 1.0e7  # TPB length per unit area [1/m]

        def newton_solve(f, xstart, C=0.0):
            """ Solve f(x) = C by Newton iteration. """
            x0 = xstart
            dx = 1.0e-6

            n = 0
            while True:
                n += 1
                f0 = f(x0) - C
                x0 -= f0/(f(x0 + dx) - C - f0)*dx
                if n > 1000:
                    raise Exception('No convergence in Newton solve')
                if abs(f0) < 0.00001:
                    return x0

        # Anode-side phases
        tpb_a = ct.Interface(mech, "tpb")
        anode_surf = tpb_a.adjacent["metal_surface"]
        gas_a = anode_surf.adjacent["gas"]
        oxide_surf_a = tpb_a.adjacent["oxide_surface"]
        oxide_a = oxide_surf_a.adjacent["oxide_bulk"]
        anode_bulk = tpb_a.adjacent["metal"]

        # Cathode-side phases
        tpb_c = ct.Interface(mech, "tpb")
        cathode_surf = tpb_c.adjacent["metal_surface"]
        gas_c = cathode_surf.adjacent["gas"]
        oxide_surf_c = tpb_c.adjacent["oxide_surface"]
        oxide_c = oxide_surf_c.adjacent["oxide_bulk"]
        cathode_bulk = tpb_c.adjacent["metal"]

        kElectron_a = tpb_a.kinetics_species_index("electron")
        def anode_curr(E):
            anode_bulk.electric_potential = E
            w = tpb_a.net_production_rates
            return ct.faraday * w[kElectron_a] * TPB_length_per_area

        kElectron_c = tpb_c.kinetics_species_index("electron")
        def cathode_curr(E):
            cathode_bulk.electric_potential = E + oxide_c.electric_potential
            w = tpb_c.net_production_rates
            return -ct.faraday * w[kElectron_c] * TPB_length_per_area

        # initialization
        gas_a.TPX = T, P, 'H2:0.97, H2O:0.03'
        gas_a.equilibrate('TP')
        gas_c.TPX = T, P, 'O2:1.0, H2O:0.001'
        gas_c.equilibrate('TP')

        for p in [anode_bulk, anode_surf, oxide_surf_a, oxide_a, cathode_bulk,
                  cathode_surf, oxide_surf_c, oxide_c, tpb_a, tpb_c]:
            p.TP = T, P

        for s in [anode_surf, oxide_surf_a, cathode_surf, oxide_surf_c]:
            s.advance_coverages(50.0)

        # These values are just a regression test with no theoretical basis
        assert anode_surf.coverages == approx(
               [6.18736878e-01, 3.81123655e-01, 8.6303646e-05,
               2.59274203e-06, 5.05700981e-05], rel=1e-7)

        assert oxide_surf_a.coverages == approx(
               [4.99435780e-02, 9.48927983e-01, 1.12840418e-03,
               3.35936530e-08], rel=1e-7)

        assert cathode_surf.coverages == approx(
               [1.48180380e-07, 7.57234727e-14, 9.99999827e-01,
               2.49235513e-08, 4.03296469e-13], rel=1e-7)

        assert oxide_surf_c.coverages == approx(
               [4.99896947e-02, 9.49804199e-01, 2.06104679e-04,
               1.11970271e-09], rel=1e-7)

        Ea0 = newton_solve(anode_curr, xstart=-0.51)
        Ec0 = newton_solve(cathode_curr, xstart=0.51)

        data = []

        # vary the anode overpotential, from cathodic to anodic polarization
        for Ea in np.linspace(Ea0 - 0.25, Ea0 + 0.25, 20):
            anode_bulk.electric_potential = Ea
            curr = anode_curr(Ea)
            delta_V = curr * 5.0e-5 / 2.0
            phi_oxide_c = -delta_V
            oxide_c.electric_potential = phi_oxide_c
            oxide_surf_c.electric_potential = phi_oxide_c
            Ec = newton_solve(cathode_curr, xstart=Ec0+0.1, C=curr)
            cathode_bulk.electric_potential = phi_oxide_c + Ec
            data.append([Ea - Ea0, 0.1*curr, Ec - Ec0, delta_V,
                         cathode_bulk.electric_potential
                         - anode_bulk.electric_potential])

        compare(data, test_data_path / "sofc-test.csv", rtol=1e-7)


class TestLithiumIonBatteryKinetics:
    """ Test based on lithium_ion_battery.py """
    _mech = "lithium_ion_battery.yaml"

    def test_lithium_ion_battery(self, test_data_path):
        mech = self._mech
        samples = 11
        soc = np.linspace(0., 1., samples)  # [-] Input state of charge (0...1)
        current = -1  # [A] Externally-applied current, negative for discharge
        T = 293  # T in K
        P = ct.one_atm
        R_electrolyte = 0.0384  # [Ohm] Electrolyte resistance
        area_cathode = 1.1167  # [m^2] Cathode total active material surface area
        area_anode = 0.7824  # [m^2] Anode total active material surface area

        # Calculate mole fractions from SOC
        X_Li_anode_0 = 0.01  # [-] anode Li mole fraction at SOC = 0
        X_Li_anode_1 = 0.75  # [-] anode Li mole fraction at SOC = 100
        X_Li_cathode_0 = 0.99  # [-] cathode Li mole fraction at SOC = 0
        X_Li_cathode_1 = 0.49  # [-] cathode Li mole fraction at SOC = 100
        X_Li_anode = (X_Li_anode_1 - X_Li_anode_0) * soc + X_Li_anode_0
        X_Li_cathode = (X_Li_cathode_0 - X_Li_cathode_1) * (1 - soc) + X_Li_cathode_1

        def newton_solve(f, xstart, C=0.0):
            """ Solve f(x) = C by Newton iteration. """
            x0 = xstart
            dx = 1.0e-6
            n = 0
            while True:
                n += 1
                f0 = f(x0) - C
                x0 -= f0 / (f(x0 + dx) - C - f0) * dx
                if n > 1000:
                    raise Exception('No convergence in Newton solve')
                if abs(f0) < 0.00001:
                    return x0

        # Phases
        anode, cathode, metal, electrolyte = ct.import_phases(
            mech, ["anode", "cathode", "electron", "electrolyte"])
        anode_int = ct.Interface(
            mech, "edge_anode_electrolyte", adjacent=[anode, metal, electrolyte])
        cathode_int = ct.Interface(
            mech, "edge_cathode_electrolyte", adjacent=[cathode, metal, electrolyte])

        # initialization
        for phase in [anode, cathode, metal, electrolyte, anode_int, cathode_int]:
            phase.TP = T, P

        # This function returns the Cantera calculated anode current
        def anode_current(phi_s, phi_l, X_Li_anode):
            # Set mole fraction and electrode and electrolyte potential
            anode.X = {"Li[anode]": X_Li_anode, "V[anode]": 1 - X_Li_anode}
            metal.electric_potential = phi_s
            electrolyte.electric_potential = phi_l

            # Calculate the current.
            return ct.faraday * anode_int.net_rates_of_progress[0] * area_anode

        # This function returns the Cantera calculated cathode current
        def cathode_current(phi_s, phi_l, X_Li_cathode):
            # Set mole fraction and electrode and electrolyte potential
            cathode.X = {"Li[cathode]": X_Li_cathode, "V[cathode]": 1 - X_Li_cathode}
            metal.electric_potential = phi_s
            electrolyte.electric_potential = phi_l

            # Calculate the current. Should be negative for cell discharge.
            return - ct.faraday * cathode_int.net_rates_of_progress[0] * area_cathode

        # Calculate cell voltage, separately for each entry of the input vectors
        data = []
        phi_l_anode = 0
        phi_s_cathode = 0
        for i in range(samples):
            # Calculate anode electrolyte potential
            phi_s_anode = 0
            phi_l_anode = newton_solve(
                lambda E: anode_current(phi_s_anode, E, X_Li_anode[i]),
                phi_l_anode, C=current)

            # Calculate cathode electrode potential
            phi_l_cathode = phi_l_anode + current * R_electrolyte
            phi_s_cathode = newton_solve(
                lambda E: cathode_current(E, phi_l_cathode, X_Li_cathode[i]),
                phi_s_cathode, C=current)

            # Calculate cell voltage
            data.append(phi_s_cathode - phi_s_anode)

        data = np.array(data).ravel()
        ref = np.genfromtxt(test_data_path / "lithium-ion-battery-test.csv")
        assert data == approx(ref, rel=1e-7)

    def test_interface_current(self):
        anode_int = ct.Interface("lithium_ion_battery.yaml", "edge_anode_electrolyte")
        anode = anode_int.adjacent["anode"]
        elect = anode_int.adjacent["electron"]
        elyte = anode_int.adjacent["electrolyte"]

        anode.X = [0.9, 0.1]
        elyte.X = [0.4, 0.3, 0.15, 0.15]

        anode.electric_potential = 0.
        elyte.electric_potential = 3.

        phases = [anode_int, anode, elect, elyte]

        for p in phases:
            net_prod_rates = anode_int.get_net_production_rates(p)
            charges = p.charges

            method = anode_int.interface_current(p)
            manual = sum(net_prod_rates * charges) * ct.faraday

            assert method == manual


class TestDuplicateReactions:
    infile = 'duplicate-reactions.yaml'

    def check(self, name):
        with pytest.raises(ct.CanteraError, match='duplicate reaction'):
            ct.Solution(self.infile, name)

    def test_forward_multiple(self):
        self.check('A')

    def test_opposite_direction1(self):
        self.check('B')

    def test_opposite_direction2(self):
        self.check('C')

    def test_opposite_direction3(self):
        self.check('D')

    def test_opposite_direction4(self):
        gas = ct.Solution(self.infile, 'E')
        assert gas.n_reactions == 2

    def test_common_efficiencies(self):
        self.check('F')

    def test_disjoint_efficiencies(self):
        gas = ct.Solution(self.infile, 'G')
        assert gas.n_reactions == 2

    def test_different_type(self):
        gas = ct.Solution(self.infile, 'H')
        assert gas.n_reactions == 2

    def test_declared_duplicate(self):
        gas = ct.Solution(self.infile, 'I')
        assert gas.n_reactions == 2

    def test_unmatched_duplicate(self):
        self.check('J')

    def test_nonreacting_species(self):
        gas = ct.Solution(self.infile, 'K')
        assert gas.n_reactions == 3


@pytest.fixture(scope='class')
def gas():
    gas = ct.Solution('h2o2.yaml', transport_model=None)
    gas.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'
    gas.TP = 900, 2*ct.one_atm
    return gas

@pytest.fixture(scope='class')
def species():
    return ct.Species.list_from_file("h2o2.yaml")

class TestReaction:

    def test_from_yaml(self, gas):
        r = ct.Reaction.from_yaml(
                "{equation: 2 O + M <=> O2 + M,"
                " type: three-body,"
                " rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0},"
                " efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}}",
                gas)

        assert r.third_body is not None
        assert r.reactants['O'] == 2
        assert r.products['O2'] == 1
        assert r.third_body.efficiencies['H2O'] == 15.4
        assert r.rate.temperature_exponent == -1.0
        assert 'O' in r
        assert 'O2' in r
        assert 'H2O' not in r

    def test_list_from_file(self, gas):
        R = ct.Reaction.list_from_file("h2o2.yaml", gas)
        eq1 = [r.equation for r in R]
        eq2 = [r.equation for r in gas.reactions()]
        assert eq1 == eq2

    def test_list_from_yaml(self, gas):
        yaml = """
            - equation: O + H2 <=> H + OH  # Reaction 3
              rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
            - equation: O + HO2 <=> OH + O2  # Reaction 4
              rate-constant: {A: 2.0e+13, b: 0.0, Ea: 0.0}
            - equation: O + H2O2 <=> OH + HO2  # Reaction 5
              rate-constant: {A: 9.63e+06, b: 2.0, Ea: 4000.0}
        """
        R = ct.Reaction.list_from_yaml(yaml, gas)
        assert len(R) == 3
        assert 'HO2' in R[2].products
        assert R[0].rate.temperature_exponent == 2.7

    def test_input_data_from_file(self, gas):
        R = gas.reaction(0)
        data = R.input_data
        assert data['type'] == 'three-body'
        assert data['efficiencies'] == {'H2': 2.4, 'H2O': 15.4, 'AR': 0.83}
        assert data['equation'] == R.equation

    def test_input_data_from_scratch(self):
        r = ct.Reaction({"O": 1, "H2": 1}, {"H": 1, "OH": 1},
                        ct.ArrheniusRate(3.87e1, 2.7, 2.6e7))
        data = r.input_data
        assert data['rate-constant']['A'] == approx(3.87e1)
        assert data['rate-constant']['b'] == approx(2.7)
        assert data['rate-constant']['Ea'] == approx(2.6e7)
        terms = data['equation'].split()
        assert 'O' in terms
        assert 'OH' in terms

    def test_custom_from_scratch(self, gas):
        species = gas.species()
        custom_reactions = gas.reactions()

        L = lambda T: 38.7 * T**2.7 * np.exp(-3150.15/T)
        rate1 = ct.CustomRate(L)
        custom_reactions[2] = ct.Reaction(
            equation='H2 + O <=> H + OH',
            rate=rate1)

        gas1 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=custom_reactions)
        gas1.TPX = gas.TPX

        # remove references to Python objects
        del custom_reactions
        del rate1

        assert gas1.reaction(2).rate.type == 'custom-rate-function'
        assert gas1.net_production_rates[2] == approx(gas.net_production_rates[2],
                                                      rel=1e-5)

    def test_modify_invalid(self, gas):
        # different reaction type
        tbr = gas.reaction(0)
        R2 = ct.Reaction(tbr.reactants, tbr.products, tbr.rate)
        with pytest.raises(ct.CanteraError, match='types are different'):
            gas.modify_reaction(0, R2)

        # different reactants
        R = gas.reaction(4)
        with pytest.raises(ct.CanteraError, match='Reactants are different'):
            gas.modify_reaction(24, R)

        # different products
        R = gas.reaction(15)
        with pytest.raises(ct.CanteraError, match='Products are different'):
            gas.modify_reaction(16, R)

class TestElementaryReaction:

    def test_elementary(self, gas, species):
        r = ct.Reaction({"O":1, "H2":1}, {"H":1, "OH":1},
                        ct.ArrheniusRate(3.87e1, 2.7, 6260*1000*4.184))

        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=[r])
        gas2.TPX = gas.TPX

        assert (gas2.forward_rate_constants[0]
                == approx(gas.forward_rate_constants[2]))
        assert (gas2.net_rates_of_progress[0]
                == approx(gas.net_rates_of_progress[2]))

    def test_arrhenius_rate(self, gas):
        R = gas.reaction(2)
        assert R.rate(gas.T) == approx(gas.forward_rate_constants[2])

    def test_negative_A(self):
        species = ct.Species.list_from_file("gri30.yaml")
        rate = ct.ArrheniusRate(-2.16e13, -0.23, 0)
        assert rate.allow_negative_pre_exponential_factor is False

        with pytest.raises(ct.CanteraError, match='negative pre-exponential'):
            r = ct.Reaction("NH:1, NO:1", "N2O:1, H:1", rate)

        rate.allow_negative_pre_exponential_factor = True
        r = ct.Reaction("NH:1, NO:1", "N2O:1, H:1", rate)
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=species, reactions=[r])

    def test_modify_elementary(self, gas):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        gas.TPX = gas.TPX
        R = gas.reaction(2)
        A1 = R.rate.pre_exponential_factor
        b1 = R.rate.temperature_exponent
        Ta1 = R.rate.activation_energy / ct.gas_constant
        T = gas.T
        assert A1*T**b1*np.exp(-Ta1/T) == approx(gas.forward_rate_constants[2])

        A2 = 1.5 * A1
        b2 = b1 + 0.1
        Ta2 = Ta1 * 1.2
        R.rate = ct.ArrheniusRate(A2, b2, Ta2 * ct.gas_constant)
        gas.modify_reaction(2, R)
        assert A2*T**b2*np.exp(-Ta2/T) == approx(gas.forward_rate_constants[2])

class TestFalloffReaction:

    def test_negative_A_falloff(self):
        species = ct.Species.list_from_file("gri30.yaml")
        low_rate = ct.Arrhenius(2.16e13, -0.23, 0)
        high_rate = ct.Arrhenius(-8.16e12, -0.5, 0)

        with pytest.raises(ct.CanteraError, match='pre-exponential'):
            ct.LindemannRate(low_rate, high_rate, ())

        rate = ct.LindemannRate()
        assert not rate.allow_negative_pre_exponential_factor
        rate.allow_negative_pre_exponential_factor = True
        rate.high_rate = high_rate
        # Should still fail because of mixed positive and negative A factors
        with pytest.raises(ct.CanteraError, match='pre-exponential'):
            rate.low_rate = low_rate

        rate.low_rate = ct.Arrhenius(-2.16e13, -0.23, 0)
        rxn = ct.Reaction("NH:1, NO:1", "N2O:1, H:1", rate)
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=species, reactions=[rxn])
        assert gas.forward_rate_constants < 0

    def test_falloff(self, gas, species):
        high_rate = ct.Arrhenius(7.4e10, -0.37, 0.0)
        low_rate = ct.Arrhenius(2.3e12, -0.9, -1700 * 1000 * 4.184)
        tb = ct.ThirdBody(efficiencies={"AR":0.7, "H2":2.0, "H2O":6.0})
        r = ct.Reaction("OH:2", "H2O2:1",
                        ct.TroeRate(low_rate, high_rate, [0.7346, 94, 1756, 5182]),
                        third_body=tb)
        assert r.rate.type == "falloff"

        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=[r])
        gas2.TPX = gas.TPX

        assert (gas2.forward_rate_constants[0]
                == approx(gas.forward_rate_constants[21]))
        assert (gas2.net_rates_of_progress[0]
                == approx(gas.net_rates_of_progress[21]))

    def test_modify_falloff(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 1100, 3 * ct.one_atm, 'CH4:1.0, O2:0.4, CO2:0.1, H2O:0.05'
        r0 = gas.reaction(11)
        assert r0.rate.type == "falloff"
        # these two reactions happen to have the same third-body efficiencies
        r1 = gas.reaction(49)
        r2 = gas.reaction(53)
        assert r2.rate.type == "falloff"
        assert r1.third_body.efficiencies == r2.third_body.efficiencies
        r2.rate = r1.rate

        gas.modify_reaction(53, r2)
        kf = gas.forward_rate_constants
        assert kf[49] == approx(kf[53])

class TestThreebodyReaction:

    def test_threebody(self, gas, species):
        tb = ct.ThirdBody(efficiencies={"AR":0.7, "H2":2.0, "H2O":6.0})
        r = ct.Reaction({"O":1, "H":1}, {"OH":1},
                        ct.ArrheniusRate(5e11, -1.0, 0.0), third_body=tb)

        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=[r])
        gas2.TPX = gas.TPX

        assert gas2.forward_rate_constants[0] == approx(gas.forward_rate_constants[1])
        assert gas2.net_rates_of_progress[0] == approx(gas.net_rates_of_progress[1])

    def test_modify_third_body(self, gas):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        gas.TPX = gas.TPX
        R = gas.reaction(5)
        A1 = R.rate.pre_exponential_factor
        b1 = R.rate.temperature_exponent
        T = gas.T
        kf1 = gas.forward_rate_constants[5]

        A2 = 1.7 * A1
        b2 = b1 - 0.1
        R.rate = ct.ArrheniusRate(A2, b2, 0.0)
        gas.modify_reaction(5, R)
        kf2 = gas.forward_rate_constants[5]
        assert (A2*T**b2) / (A1*T**b1) == approx(kf2/kf1)

class TestPlogReaction:

    def test_plog(self):
        gas1 = ct.Solution('pdep-test.yaml')
        species = ct.Species.list_from_file("pdep-test.yaml")

        rate = ct.PlogRate([
            (0.01*ct.one_atm, ct.Arrhenius(1.2124e13, -0.5779, 10872.7*4184)),
            (1.0*ct.one_atm, ct.Arrhenius(4.9108e28, -4.8507, 24772.8*4184)),
            (10.0*ct.one_atm, ct.Arrhenius(1.2866e44, -9.0246, 39796.5*4184)),
            (100.0*ct.one_atm, ct.Arrhenius(5.9632e53, -11.529, 52599.6*4184))
        ])
        r = ct.Reaction({"R1A":1, "R1B":1}, {"P1":1, "H":1}, rate)

        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=[r])

        gas2.X = gas1.X = 'R1A:0.3, R1B:0.6, P1:0.1'

        for P in [0.001, 0.01, 0.2, 1.0, 1.1, 9.0, 10.0, 99.0, 103.0]:
            gas1.TP = gas2.TP = 900, P * ct.one_atm
            assert (gas2.forward_rate_constants[0]
                    == approx(gas1.forward_rate_constants[0]))
            assert (gas2.net_rates_of_progress[0]
                    == approx(gas1.net_rates_of_progress[0]))

    def test_plog_rate(self):
        gas1 = ct.Solution('pdep-test.yaml')
        gas1.TP = 800, 2*ct.one_atm
        for i in range(4):
            assert (gas1.reaction(i).rate(gas1.T, gas1.P)
                    == approx(gas1.forward_rate_constants[i]))

    def test_plog_invalid_third_body(self):
        with pytest.raises(ct.CanteraError, match="Found superfluous"):
            gas = ct.Solution("pdep-test.yaml", "plog-invalid")

    def test_modify_plog(self):
        gas = ct.Solution('pdep-test.yaml')
        gas.TPX = 1010, 0.12 * ct.one_atm, 'R1A:0.3, R1B:0.2, H:0.1, R2:0.4'

        r0 = gas.reaction(0)
        r1 = gas.reaction(1)
        r0.rate = ct.PlogRate(r1.rate.rates)
        gas.modify_reaction(0, r0)
        kf = gas.forward_rate_constants
        assert kf[0] == approx(kf[1])

        # Removing the high-pressure rates should have no effect at low P...
        r1.rate = ct.PlogRate(rates=r1.rate.rates[:-4])
        gas.modify_reaction(1, r1)
        assert kf[1] == approx(gas.forward_rate_constants[1])

        # ... but should change the rate at higher pressures
        gas.TP = 1010, 12.0 * ct.one_atm
        kf = gas.forward_rates_of_progress
        assert kf[0] != approx(kf[1])

class TestLinearBurkeReaction:
    @pytest.fixture(scope='class')
    def gas_baseline(self):
        return ct.Solution('linearBurke-test.yaml', name='baseline_mechanism')

    @pytest.fixture(scope='class')
    def gas_linearBurke(self):
        return ct.Solution('linearBurke-test.yaml', name='linear-Burke_mechanism')

    @pytest.mark.parametrize("reaction", [
        pytest.param("H + OH <=> H2O", id="PLOG"),
        pytest.param("H + O2 (+M) <=> HO2 (+M)", id="Troe"),
        pytest.param("H2O2 <=> 2 OH", id="Chebyshev"),
    ])
    def test_efficiency(self, gas_baseline, gas_linearBurke, reaction):
        T = 1000 # [K]
        P_ls = [0.1,1,10,100] # [atm]
        for P in P_ls:
            def getK(gas, T, P, X):
                gas.TPX = T,P,X
                eqn = gas.reaction_equations().index(reaction)
                return gas.forward_rate_constants[eqn]
            # collider 'O2' treated as M in this test reaction
            k_baseline = getK(gas_baseline, T, P,'O2:1')
            k_linearBurke = getK(gas_linearBurke, T, P,'O2:1')
            assert k_baseline == approx(k_linearBurke)
            # collider 'H2O' must behave as 'M' if 'M' were eval. at 10x the pressure
            k_baseline = getK(gas_baseline, T, P*10,'H2O:1')
            k_linearBurke = getK(gas_linearBurke, T, P,'H2O:1')
            assert k_baseline == approx(k_linearBurke)


class TestChebyshevReaction:
    def test_chebyshev(self):
        gas1 = ct.Solution('pdep-test.yaml')
        species = ct.Species.list_from_file("pdep-test.yaml")

        rate = ct.ChebyshevRate(
            temperature_range=(300.0, 2000.0),
            pressure_range=(1000, 10000000),
            data=[[ 5.28830e+00, -1.13970e+00, -1.20590e-01,  1.60340e-02],
                  [ 1.97640e+00,  1.00370e+00,  7.28650e-03, -3.04320e-02],
                  [ 3.17700e-01,  2.68890e-01,  9.48060e-02, -7.63850e-03],
                  [-3.12850e-02, -3.94120e-02,  4.43750e-02,  1.44580e-02]])
        r = ct.Reaction("R5:1, H:1", "P5A:1, P5B:1", rate)

        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=[r])

        gas2.X = gas1.X = 'R5:0.3, P5A:0.6, H:0.1'

        for T,P in itertools.product([300, 500, 1500], [1e4, 4e5, 3e6]):
            gas1.TP = gas2.TP = T, P
            assert (gas2.forward_rate_constants[0]
                    == approx(gas1.forward_rate_constants[4]))
            assert (gas2.net_rates_of_progress[0]
                    == approx(gas1.net_rates_of_progress[4]))

    def test_chebyshev_single_P(self):
        species = ct.Species.list_from_file("pdep-test.yaml")
        rate = ct.ChebyshevRate(
            temperature_range=(300.0, 2000.0),
            pressure_range=(1000, 10000000),
            data=[[ 5.28830e+00],
                  [ 1.97640e+00],
                  [ 3.17700e-01],
                  [-3.12850e-02]])
        r = ct.Reaction("R5:1, H:1", "P5A:1, P5B:1", rate)

        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=species, reactions=[r])

        # rate constant should be pressure independent
        for T in [300, 500, 1500]:
            gas.TP = T, 1e4
            k1 = gas.forward_rate_constants[0]
            gas.TP = T, 1e6
            k2 = gas.forward_rate_constants[0]
            assert k1 == approx(k2)

    def test_chebyshev_single_T(self):
        species = ct.Species.list_from_file("pdep-test.yaml")
        rate = ct.ChebyshevRate(
            temperature_range=(300.0, 2000.0),
            pressure_range=(1000, 10000000),
            data=[[ 5.28830e+00, -1.13970e+00, -1.20590e-01,  1.60340e-02]])
        r = ct.Reaction("R5:1, H:1", "P5A:1, P5B:1", rate)

        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=species, reactions=[r])

        # rate constant should be temperature independent
        for P in [1e4, 2e5, 8e6]:
            gas.TP = 400, P
            k1 = gas.forward_rate_constants[0]
            gas.TP = 1700, P
            k2 = gas.forward_rate_constants[0]
            assert k1 == approx(k2)

    def test_chebyshev_rate(self):
        gas1 = ct.Solution('pdep-test.yaml')
        gas1.TP = 800, 2*ct.one_atm
        for i in range(4,6):
            assert (gas1.reaction(i).rate(gas1.T, gas1.P)
                    == approx(gas1.forward_rate_constants[i]))

    def test_chebyshev_bad_shape_yaml(self):
        species = ct.Species.list_from_file("pdep-test.yaml")
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=species, reactions=[])

        with pytest.raises(ct.CanteraError, match="Inconsistent"):
            r = ct.Reaction.from_yaml('''
                equation: R5 + H <=> P5A + P5B
                type: Chebyshev
                temperature-range: [300.0, 2000.0]
                pressure-range: [9.86e-03 atm, 98.6 atm]
                data:
                - [8.2883, -1.1397, -0.12059, 0.016034]
                - [1.9764, 1.0037, 7.2865e-03]
                - [0.3177, 0.26889, 0.094806, -7.6385e-03]
                - [-0.031285, -0.039412, 0.044375, 0.014458]''', gas)

    def test_chebyshev_deprecated_third_body(self):
        with pytest.raises(ct.CanteraError, match="in the reaction equation"):
            gas = ct.Solution("pdep-test.yaml", "chebyshev-deprecated")

    def test_modify_chebyshev(self):
        gas = ct.Solution('pdep-test.yaml')
        gas.TPX = 1010, 0.34 * ct.one_atm, 'R1A:0.3, R1B:0.2, H:0.1, R2:0.4'

        r1 = gas.reaction(4)
        r2 = gas.reaction(5)
        r1.rate = ct.ChebyshevRate(r2.rate.temperature_range, r2.rate.pressure_range,
                                   r2.rate.data)

        # rates should be different before calling 'modify_reaction'
        kf = gas.forward_rate_constants
        assert kf[4] != approx(kf[5])

        gas.modify_reaction(4, r1)
        kf = gas.forward_rate_constants
        assert kf[4] == approx(kf[5])

class TestBlowersMaselReaction:

    def test_BlowersMasel(self, gas):
        r = ct.Reaction({"O":1, "H2":1}, {"H":1, "OH":1},
            ct.BlowersMaselRate(3.87e1, 2.7, 6260*1000*4.184, 1e9*1000*4.184))

        gas1 = ct.Solution("blowers-masel.yaml", "gas")
        assert isinstance(gas1.reaction(0).rate, ct.BlowersMaselRate)

        gas2 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=gas1.species(), reactions=[r])

        gas1.TP = gas.TP
        gas2.TP = gas.TP
        gas1.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'
        gas2.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'

        assert (gas2.forward_rate_constants[0]
                == approx(gas1.forward_rate_constants[0], rel=1e-7))
        assert (gas2.net_rates_of_progress[0]
                == approx(gas1.net_rates_of_progress[0], rel=1e-7))

    def test_negative_A_blowersmasel(self):
        species = ct.Solution("blowers-masel.yaml").species()
        rate = ct.BlowersMaselRate(-3.87e1, 2.7, 6260*1000*4.184, 1e9)
        assert rate.allow_negative_pre_exponential_factor is False
        with pytest.raises(ct.CanteraError, match='negative pre-exponential'):
            r = ct.Reaction({'O':1, 'H2':1}, {'H':1, 'OH':1}, rate)

        rate.allow_negative_pre_exponential_factor = True
        r = ct.Reaction({'O':1, 'H2':1}, {'H':1, 'OH':1}, rate)
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=species, reactions=[r])

    def test_Blowers_Masel_change_enthalpy(self):
        gas = ct.Solution("blowers-masel.yaml")
        r = gas.reaction(0)
        E0 = r.rate.activation_energy
        w = r.rate.bond_energy
        A = r.rate.pre_exponential_factor
        b = r.rate.temperature_exponent
        vp = 2 * w * (w+E0) / (w - E0)
        deltaH = gas.delta_enthalpy[0]
        E = ((w + deltaH / 2) * (vp - 2 * w + deltaH) ** 2
             / (vp ** 2 - 4 * w ** 2 + deltaH ** 2))

        deltaH_high = 10 * gas.reaction(0).rate.activation_energy
        deltaH_low = -20 * gas.reaction(0).rate.activation_energy
        index = gas.species_index('OH')
        species = gas.species('OH')

        gas.reaction(0).rate.delta_enthalpy = deltaH
        assert gas.reaction(0).rate.delta_enthalpy == approx(deltaH)
        assert E == approx(gas.reaction(0).rate.activation_energy)

        perturbed_coeffs = species.thermo.coeffs.copy()
        perturbed_coeffs[6] += deltaH_high / ct.gas_constant
        perturbed_coeffs[13] += deltaH_high / ct.gas_constant
        species.thermo = ct.NasaPoly2(species.thermo.min_temp, species.thermo.max_temp,
                            species.thermo.reference_pressure, perturbed_coeffs)
        gas.modify_species(index, species)
        gas.reaction(0).rate.delta_enthalpy = deltaH_high
        assert gas.reaction(0).rate.delta_enthalpy == approx(deltaH_high)
        assert deltaH_high == approx(gas.reaction(0).rate.activation_energy)
        assert (A*gas.T**b*np.exp(-deltaH_high/ct.gas_constant/gas.T)
                == approx(gas.forward_rate_constants[0]))

        perturbed_coeffs = species.thermo.coeffs.copy()
        perturbed_coeffs[6] += deltaH_low / ct.gas_constant
        perturbed_coeffs[13] += deltaH_low / ct.gas_constant
        species.thermo = ct.NasaPoly2(species.thermo.min_temp, species.thermo.max_temp,
                            species.thermo.reference_pressure, perturbed_coeffs)
        gas.modify_species(index, species)
        gas.reaction(0).rate.delta_enthalpy = deltaH_low
        assert gas.reaction(0).rate.delta_enthalpy == approx(deltaH_low)
        assert gas.reaction(0).rate.activation_energy == 0
        assert (A*gas.T**b*np.exp(0/ct.gas_constant/gas.T)
                == approx(gas.forward_rate_constants[0]))

    def test_modify_BlowersMasel(self, gas):
        gas = ct.Solution("blowers-masel.yaml")
        gas.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'
        gas.TP = gas.TP
        R = gas.reaction(0)
        delta_enthalpy = gas.delta_enthalpy[0]
        A1 = R.rate.pre_exponential_factor
        b1 = R.rate.temperature_exponent
        R.rate.delta_enthalpy = delta_enthalpy
        Ta1 = R.rate.activation_energy / ct.gas_constant
        T = gas.T
        assert A1 * T**b1 * np.exp(-Ta1 / T) == approx(gas.forward_rate_constants[0])

        # randomly modify the rate parameters of a Blowers-Masel reaction
        A2 = 1.5 * A1
        b2 = b1 + 0.1
        Ta_intrinsic = R.rate.activation_energy * 1.2
        w = R.rate.bond_energy * 0.8
        R.rate = ct.BlowersMaselRate(A2, b2, Ta_intrinsic, w)
        delta_enthalpy = gas.delta_enthalpy[0]
        R.rate.delta_enthalpy = delta_enthalpy
        Ta2 = R.rate.activation_energy / ct.gas_constant
        gas.modify_reaction(0, R)
        assert A2 * T**b2 * np.exp(-Ta2 / T) == approx(gas.forward_rate_constants[0])

class TestInterfaceReaction:

    def test_interface(self):
        surf_species = ct.Species.list_from_file("ptcombust.yaml")
        surf1 = ct.Interface("ptcombust.yaml", "Pt_surf")
        gas = surf1.adjacent["gas"]

        rate = ct.InterfaceArrheniusRate(3.7e20, 0, 67.4e6)
        rate.coverage_dependencies = {'H(S)': (0, 0, -6e6)}
        assert rate.coverage_dependencies["H(S)"]["E"] == approx(-6e6)
        r1 = ct.Reaction(equation="2 H(S) <=> H2 + 2 PT(S)", rate=rate)

        surf2 = ct.Interface(thermo='ideal-surface', species=surf_species,
                             kinetics='surface', reactions=[r1], adjacent=[gas])

        surf2.site_density = surf1.site_density
        surf1.coverages = surf2.coverages = 'PT(S):0.7, H(S):0.3'
        gas.TP = surf2.TP = surf1.TP

        for T in [300, 500, 1500]:
            gas.TP = surf1.TP = surf2.TP = T, 5*ct.one_atm
            assert (surf1.forward_rate_constants[1]
                    == approx(surf2.forward_rate_constants[0]))
            assert (surf1.net_rates_of_progress[1]
                    == approx(surf2.net_rates_of_progress[0]))

    def test_BlowersMaselinterface(self):
        gas = ct.Solution("gri30.yaml", transport_model=None)
        gas.TPX = 300, ct.one_atm, {"CH4": 0.095, "O2": 0.21, "AR": 0.79}
        surf1 = ct.Interface("blowers-masel.yaml", "Pt_surf", [gas])
        rate = ct.InterfaceBlowersMaselRate(3.7e20, 0, 67.4e6, 1e9)
        rate.coverage_dependencies = {"H(S)": (0, 0, -6e6)}
        assert rate.coverage_dependencies["H(S)"]["E"] == approx(-6e6)

        r1 = ct.Reaction("H(S):2", "H2:1, PT(S):2", rate)

        surf_species = []
        for species in surf1.species():
            surf_species.append(species)
        surf2 = ct.Interface(thermo='ideal-surface', species=surf_species,
                             kinetics='surface', reactions=[r1], adjacent=[gas])

        surf2.site_density = surf1.site_density
        surf1.coverages = surf2.coverages = 'PT(S):0.7, H(S):0.3'
        gas.TP = surf2.TP = surf1.TP

        for T in [300, 500, 1500]:
            gas.TP = surf1.TP = surf2.TP = T, 5*ct.one_atm
            assert (surf1.forward_rate_constants[0]
                    == approx(surf2.forward_rate_constants[0]))
            assert (surf1.net_rates_of_progress[0]
                    == approx(surf2.net_rates_of_progress[0]))

    def test_modify_interface(self):
        surf = ct.Interface("ptcombust.yaml", "Pt_surf")
        surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'

        R = surf.reaction(1)
        R.rate.coverage_dependencies = {'O(S)': (0.0, 0.0, -3e6)}
        surf.modify_reaction(1, R)

        # Rate constant should now be independent of H(S) coverage, but
        # dependent on O(S) coverage
        k1 = surf.forward_rate_constants[1]
        surf.coverages = 'O(S):0.2, PT(S):0.4, H(S):0.4'
        k2 = surf.forward_rate_constants[1]
        surf.coverages = 'O(S):0.2, PT(S):0.6, H(S):0.2'
        k3 = surf.forward_rate_constants[1]
        assert k1 != approx(k2)
        assert k2 == approx(k3)

    def test_modify_BMinterface(self):
        gas = ct.Solution("gri30.yaml", transport_model=None)
        gas.TPX = 300, ct.one_atm, {"CH4": 0.095, "O2": 0.21, "AR": 0.79}
        surf = ct.Interface("blowers-masel.yaml", "Pt_surf", [gas])
        surf.coverages = "O(S):0.1, PT(S):0.5, H(S):0.4"
        gas.TP = surf.TP

        R = surf.reaction(0)
        R.rate.coverage_dependencies = {'O(S)': (0.0, 0.0, -3e6)}
        surf.modify_reaction(0, R)

        # Rate constant should now be independent of H(S) coverage, but
        # dependent on O(S) coverage
        Ek = surf.reaction(0).rate.coverage_dependencies["O(S)"]["E"]
        k1 = surf.forward_rate_constants[0]
        O2_theta_k1 = surf.coverages[-1]
        surf.coverages = "O(S):0.2, PT(S):0.4, H(S):0.4"
        k2 = surf.forward_rate_constants[0]
        O2_theta_k2 = surf.coverages[-1]
        O2_delta_theta_k = O2_theta_k1 - O2_theta_k2
        surf.coverages = "O(S):0.2, PT(S):0.6, H(S):0.2"
        k3 = surf.forward_rate_constants[0]

        assert (np.exp(-O2_delta_theta_k * Ek / ct.gas_constant / surf.T)
               == approx(k1 / k2))
        assert k2 == approx(k3)


class TestStickingCoefficient:

    def test_invalid_sticking(self):
        yaml = """
        equation: OH + Csoot-H + CB-CB3 + CO => Csoot-* + 2 CO + H2
        sticking-coefficient: {A: 0.13, b: 0.0, Ea: 0.0}"""
        surf = ct.Interface("haca2.yaml", "soot_interface")
        rxn = ct.Reaction.from_yaml(yaml, surf)
        with pytest.raises(ct.CanteraError, match="non-interface species"):
            surf.add_reaction(rxn)

    def test_modify_sticking(self):
        surf = ct.Interface("ptcombust.yaml", "Pt_surf")
        surf.coverages = "O(S):0.1, PT(S):0.5, H(S):0.4"

        R = surf.reaction(2)
        R.rate = ct.StickingArrheniusRate(0.25, 0, 0) # original sticking coefficient = 1.0

        k1 = surf.forward_rate_constants[2]
        surf.modify_reaction(2, R)
        k2 = surf.forward_rate_constants[2]
        assert k1 == approx(4*k2)

    def test_motz_wise(self):
        # Motz & Wise off for all reactions
        surf1 = ct.Interface("ptcombust.yaml", "Pt_surf")
        surf1.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'

        # Motz & Wise correction on for some reactions
        surf2 = ct.Interface("ptcombust-motzwise.yaml", "Pt_surf")
        surf2.TPY = surf1.TPY

        k1 = surf1.forward_rate_constants
        k2 = surf2.forward_rate_constants

        # M&W toggled on (globally) for reactions 2 and 7
        assert 2.0 * k1[2] == approx(k2[2]) # sticking coefficient = 1.0
        assert 1.6 * k1[7] == approx(k2[7]) # sticking coefficient = 0.75

        # M&W toggled off (locally) for reaction 4
        assert k1[4] == approx(k2[4])

        # M&W toggled on (locally) for reaction 9
        assert 2.0 * k1[9] == approx(k2[9]) # sticking coefficient = 1.0

    def test_modify_BMsticking(self):
        gas = ct.Solution("gri30.yaml", transport_model=None)
        gas.TPX = 300, ct.one_atm, {"CH4": 0.095, "O2": 0.21, "AR": 0.79}
        surf = ct.Interface("blowers-masel.yaml", "Pt_surf", [gas])
        surf.coverages = "O(S):0.1, PT(S):0.5, H(S):0.4"
        gas.TP = surf.TP

        R = surf.reaction(1)
        R.rate = ct.StickingBlowersMaselRate(0.25, 0, 0, 1000000) # original sticking coefficient = 1.0

        k1 = surf.forward_rate_constants[1]
        surf.modify_reaction(1, R)
        k2 = surf.forward_rate_constants[1]
        assert k1 == approx(4*k2)

    def test_BMmotz_wise(self):
        # Motz & Wise off for all reactions
        surf1 = ct.Interface("blowers-masel.yaml", "Pt_surf")
        surf1.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'

        # Motz & Wise correction on for some reactions
        surf2 = ct.Interface("blowers-masel.yaml", "Pt_motz_wise")
        surf2.TPY = surf1.TPY

        k1 = surf1.forward_rate_constants
        k2 = surf2.forward_rate_constants

        # M&W toggled on (globally) for reactions 1 and 2
        assert 2.0 * k1[1] == approx(k2[1]) # sticking coefficient = 1.0
        assert 1.6 * k1[2] == approx(k2[2]) # sticking coefficient = 0.75

        # M&W toggled off (locally) for reaction 3
        assert k1[3] == approx(k2[3])

        # M&W toggled on (locally) for reaction 4
        assert k1[4] == approx(k2[4]) # sticking coefficient = 1.0

class TestElectronCollisionPlasmaReaction:

    def setup_plasma(self):
        gas = ct.Solution('oxygen-plasma.yaml',
                          'isotropic-electron-energy-plasma',
                          transport_model=None)
        gas.TPX = 300, ct.one_atm, {"O2": 1.0, "E": 1.0}
        return gas

    def test_ECPR_electron_energy_level(self):
        # Consistent test of using different electron energy method
        # This test compares the rate constants calculated from
        # the method discretized-electron-energy-plasma to
        # the method isotropic-electron-energy-plasma
        gas1 = self.setup_plasma()
        electron_energy_levels = gas1.electron_energy_levels
        electron_energy_dist = gas1.electron_energy_distribution
        k1 = gas1.forward_rate_constants[1]

        gas2 = ct.Solution('oxygen-plasma.yaml',
                           'discretized-electron-energy-plasma',
                           transport_model=None)
        gas2.TPX = 300, ct.one_atm, {"O2": 1.0, "E": 1.0}
        gas2.set_discretized_electron_energy_distribution(electron_energy_levels,
                                                          electron_energy_dist)
        k2 = gas2.forward_rate_constants[1]
        assert k1 == approx(k2)

    def test_ECPR_rate_of_progress(self):
        # Test of forward and reverse rate of progress
        # The reverse rate of progress calculates the super-elastic rate coefficient
        # which equals to the forward rate coefficient when the threshold energy is
        # zero. Therefore, the reverse rate of progress equals the forward rate of
        # progress.
        gas = self.setup_plasma()

        frop = approx(gas.forward_rates_of_progress[1])
        rrop = approx(gas.reverse_rates_of_progress[1])

        assert frop == rrop