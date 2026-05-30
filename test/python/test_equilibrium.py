import itertools

import numpy as np
import pytest
from pytest import approx

import cantera as ct
from .utilities import (
    compare
)


def assert_at_equilibrium(mix, abs_tol=1e-6):
    """
    Verify that ``mix`` is at chemical equilibrium by checking that the chemical
    potential of each species satisfies ``mu_k = sum_e n_{e,k} * lambda_e`` where
    ``lambda_e`` is the element potential for element e and ``n_{e,k}`` is the
    number of atoms of element e in species k. We solve for ``lambda`` by weighted
    least squares (weighting by moles so trace species do not dominate the fit),
    then require ``|mu_k - sum_e n_{e,k} lambda_e| / (R T) < abs_tol`` for every
    species whose moles are more than 1e-10 of the total. Returns the worst relative
    residual.
    """
    R = ct.gas_constant
    T = mix.T
    n_species = mix.n_species
    n_elements = mix.n_elements
    mu = np.asarray(mix.chemical_potentials)
    moles = np.asarray(mix.species_moles)
    composition = np.array([[mix.n_atoms(k, e) for e in range(n_elements)]
                            for k in range(n_species)])

    w = np.sqrt(np.maximum(moles, 0.0))
    mask = w > 0.0
    Aw = composition[mask] * w[mask, None]
    muw = mu[mask] * w[mask]
    lam, *_ = np.linalg.lstsq(Aw, muw, rcond=None)

    total = moles.sum()
    worst = 0.0
    for k in range(n_species):
        if moles[k] > 1e-10 * total:
            residual = (mu[k] - composition[k] @ lam) / (R * T)
            assert abs(residual) < abs_tol, (
                f"species {mix.species_name(k)} (k={k}): "
                f"moles={moles[k]:.3g} of total {total:.3g}, "
                f"|mu_k - sum_e n_{{e,k}} lambda_e| / RT = {abs(residual):.3g}"
            )
            worst = max(worst, abs(residual))
    return worst


class EquilTestCases:
    """
    Base class for equilibrium test cases, parameterized by the solver to use.
    """
    solver = None # must be set by subclass

    def check(self, gas, **moles):
        """
        Check that the mole fractions in `gas` match the expected values in `moles`.
        """
        nTotal = sum(moles.values())
        for name, X in moles.items():
            assert gas[name].X[0] == approx(X/nTotal)

    def test_equil_complete_stoichiometric(self):
        """
        Equilibrium should correspond to complete combustion
        """
        gas = ct.Solution("equilibrium.yaml", "complete")
        gas.TPX = 298, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_complete_lean(self):
        """
        Equilibrium should correspond to complete combustion (with excess O2)
        CH4 + 3 O2 -> CO2 + 2 H2O + O2
        """
        gas = ct.Solution("equilibrium.yaml", "complete")
        gas.TPX = 298, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_incomplete_stoichiometric(self):
        gas = ct.Solution("equilibrium.yaml", "incomplete")
        gas.TPX = 301, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_incomplete_lean(self):
        gas = ct.Solution("equilibrium.yaml", "incomplete")
        gas.TPX = 301, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_gri_stoichiometric(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 301, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_gri_lean(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 301, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_overconstrained1(self):
        gas = ct.Solution("equilibrium.yaml", "overconstrained-1")
        gas.TPX = 301, 100000, 'CH4:1.0, O2:1.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=1, O2=1)

    def test_equil_overconstrained2(self):
        gas = ct.Solution("equilibrium.yaml", "overconstrained-2")
        gas.TPX = 301, 100000, 'CH4:1.0, O2:1.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=1, O2=1)


class TestChemEquil(EquilTestCases):
    """
    Tests using the 'element_potential' solver.
    """

    solver = 'element_potential'

    def test_relaxed_temperature_limit(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 299.999, 100002.0, {
            'C2H4': 0.19274682199314855,
            'N2': 0.32751737685952764,
            'O2': 0.47973580114732384,
        }
        h0 = gas.enthalpy_mass

        assert not gas.enforce_temperature_limits
        gas.enforce_temperature_limits = True
        with pytest.raises(ct.CanteraError, match="enthalpy.*temperature bounds"):
            gas.equilibrate('HP', solver='element_potential')

        gas.TPX = 299.999, 100002.0, {
            'C2H4': 0.19274682199314855,
            'N2': 0.32751737685952764,
            'O2': 0.47973580114732384,
        }
        gas.enforce_temperature_limits = False
        with pytest.warns(UserWarning, match="outside valid range"):
            gas.equilibrate('HP', solver='element_potential')
        assert gas.T > gas.max_temp
        assert gas.enthalpy_mass == approx(h0, abs=1e-3)


class TestMultiphaseEquil(EquilTestCases):
    """
    Tests using the 'gibbs' solver.
    """

    solver = 'gibbs'

    def test_trace_solution_species(self):
        s = """
        phases:
        - name: gas
          thermo: ideal-gas
          species: [{nasa_gas.yaml/species: [C, CO, CO2, C2]}]
          skip-undeclared-elements: true
          elements: [C, O]
        """
        gas = ct.Solution(yaml=s)
        gas.TPX = 300, 1e5, 'CO2:1'
        gas.equilibrate('TP', self.solver, max_steps=10000)
        assert gas['CO2'].X[0] == approx(1.0)

    def test_equil_gri_stoichiometric(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 301, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_gri_lean(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 301, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)


class TestMultiphase_H2O2_TwoPhase:
    """
    Stress test for the 'gibbs' (``MultiPhaseEquil``) solver on a two-phase setup where
    the second phase is a subset of the h2o2 mechanism species and the mixture
    composition is stoichiometric with respect to H2 and O2. Verifies element
    conservation and that the converged state satisfies the element-potential
    consistency condition, sampling a range of phase2 compositions, temperatures,
    pressures, and inert moles. Regression coverage for the cluster of MultiPhaseEquil
    fixes addressing issue #1023.
    """

    # Spot-check conditions chosen to exercise the regimes that previously had failures
    # (issue #1023): low T where the answer is fully-combusted, mid T where minor
    # species start to matter, high T where dissociation becomes significant, and a
    # high-T/high-P combination.
    _H2O2_TPN = [
        (400, 1.0, 0.1),
        (500, 2.0, 0.3),
        (1000, 5.0, 0.7),
        (1300, 5.0, 0.7),
        (1700, 1.0, 0.0),
        (2000, 20.0, 0.7),
    ]

    @staticmethod
    def _h2o2_phase2_samples(n=20):
        """Deterministic sample of ~n phase2 species combinations from h2o2.yaml,
        spanning all combination sizes from 2 to the full set."""
        species = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'AR', 'N2']
        combos = [c for r in range(2, len(species) + 1)
                    for c in itertools.combinations(species, r)]
        stride = max(1, len(combos) // n)
        return [list(c) for c in combos[::stride][:n]]

    @pytest.mark.parametrize("phase2_names", _h2o2_phase2_samples())
    @pytest.mark.parametrize("T,P_atm,N_inert", _H2O2_TPN)
    def test_equilibrate(self, T, P_atm, N_inert, phase2_names):
        S = ct.Species.list_from_file('h2o2.yaml')
        phase1 = ct.Solution('h2o2.yaml')
        phase2 = ct.Solution(thermo='ideal-gas',
                             species=[s for s in S if s.name in phase2_names])

        P = P_atm * ct.one_atm
        initial = {"O2": 1.0, "H2": 2.0, "N2": N_inert}
        phase1.TPX = T, P, initial
        mix = ct.Mixture([(phase1, sum(initial.values())), (phase2, 0.0)])
        elements_before = np.array([mix.element_moles(e)
                                    for e in range(mix.n_elements)])

        mix.T = T
        mix.P = P
        # The default max_steps=1000 is not enough for a handful of corner cases (e.g.
        # T=1700, N_inert=0 with phase2 = {H2, O, O2, OH, H2O, N2}) where atom transfer
        # to phase2 is constrained by a single very-small species.
        mix.equilibrate('TP', solver='gibbs', max_steps=5000)

        elements_after = np.array([mix.element_moles(e)
                                   for e in range(mix.n_elements)])
        np.testing.assert_allclose(elements_after, elements_before, rtol=1e-7,
                                   atol=1e-12,
                                   err_msg="element moles changed during equilibration")
        assert_at_equilibrium(mix, abs_tol=1e-6)


class TestMultiphase_CHO:
    """
    Regression tests for the 'gibbs' (``MultiPhaseEquil``) solver on C/H/O compositions
    which previously failed to converge (Issue #265). Each case uses the GRI 3.0 gas
    phase at 923 K and 1 atm, with the element amounts set by integer atom counts (C, H,
    O). The cases exercise the two failure mechanisms that were fixed:

    - A trace single-species condensed phase (graphite) pinning the step size to a
      negligible value (oxygen-rich cases such as C=23, H=79, O=98).
    - A depleted species left in the component basis, which throttles every reaction
      that must consume it and sends the solver into a non-converging limit cycle
      (carbon-rich, low-oxygen cases such as C=77, H=120, O=3).
    """

    T = 923.0
    P = ct.one_atm

    @staticmethod
    def _ids(cases):
        return [f"C{C}-H{H}-O{O}" for C, H, O in cases]

    def _check(self, phases, C, H, O, max_steps):
        gas = phases[0][0]
        gas.TPX = self.T, self.P, {'C': C, 'H': H, 'O': O}
        mix = ct.Mixture(phases)
        elements_before = np.array([mix.element_moles(e)
                                    for e in range(mix.n_elements)])
        mix.T = self.T
        mix.P = self.P
        mix.equilibrate('TP', solver='gibbs', max_steps=max_steps)
        elements_after = np.array([mix.element_moles(e)
                                   for e in range(mix.n_elements)])
        np.testing.assert_allclose(
            elements_after, elements_before, rtol=1e-7, atol=1e-12,
            err_msg="element moles changed during equilibration")
        assert_at_equilibrium(mix, abs_tol=1e-6)

    # (C, H, O) atom counts. Single-phase (gas only) cases are carbon-rich, low-oxygen
    # compositions that previously stalled in a component-basis limit cycle.
    _single_phase = [
        (40, 147, 13), (77, 120, 3), (78, 119, 4),
        (80, 117, 3), (77, 119, 4), (80, 116, 4),
    ]
    @pytest.mark.parametrize("C,H,O", _single_phase, ids=_ids(_single_phase))
    def test_single_phase(self, C, H, O):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        self._check([(gas, 1.0)], C, H, O, max_steps=2000)

    # Two-phase (gas + graphite) cases: oxygen-rich (trace condensed-phase pinning)
    # through carbon-rich (slow convergence with the condensed phase present).
    _two_phase = [
        (23, 79, 98), (124, 67, 9), (160, 36, 4),
        (150, 46, 4), (100, 90, 10), (40, 147, 13),
    ]
    @pytest.mark.parametrize("C,H,O", _two_phase, ids=_ids(_two_phase))
    def test_two_phase(self, C, H, O):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        carbon = ct.Solution('graphite.yaml')
        self._check([(gas, 1.0), (carbon, 0.0)], C, H, O, max_steps=5000)


class TestMultiphase_VCS_CHO:
    """
    Regression tests for the 'vcs' (``vcs_MultiPhaseEquil``) solver on C/H/O
    compositions which previously failed to converge (Issue #266). Each case uses
    the GRI 3.0 gas phase plus a graphite condensed phase at 923 K and 1 atm, with
    the element amounts set by integer atom counts (C, H, O).
    """

    T = 923.0
    P = ct.one_atm

    @staticmethod
    def _ids(cases):
        return [f"C{C}-H{H}-O{O}" for C, H, O in cases]

    def _check(self, C, H, O, max_steps):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        carbon = ct.Solution('graphite.yaml')
        gas.TPX = self.T, self.P, {'C': C, 'H': H, 'O': O}
        mix = ct.Mixture([(gas, 1.0), (carbon, 0.0)])
        elements_before = np.array([mix.element_moles(e)
                                    for e in range(mix.n_elements)])
        mix.T = self.T
        mix.P = self.P
        mix.equilibrate('TP', solver='vcs', max_steps=max_steps)
        elements_after = np.array([mix.element_moles(e)
                                   for e in range(mix.n_elements)])
        np.testing.assert_allclose(
            elements_after, elements_before, rtol=1e-7, atol=1e-12,
            err_msg="element moles changed during equilibration")
        assert_at_equilibrium(mix, abs_tol=1e-6)

    # A representative spread across the band of compositions (C=13..30) that
    # previously triggered the basis-swap oscillation, including its edges.
    _cases = [
        (13, 78, 9), (15, 70, 15), (16, 65, 19), (20, 60, 20),
        (22, 56, 22), (25, 50, 25), (27, 44, 29), (30, 36, 34),
    ]
    @pytest.mark.parametrize("C,H,O", _cases, ids=_ids(_cases))
    def test_two_phase(self, C, H, O):
        self._check(C, H, O, max_steps=10000)


@pytest.fixture(scope='function')
def extra_elements(request):
    s = """
    phases:
    - name: gas
      thermo: ideal-gas
      elements: [H, Ar, C, O, Cl, N]
      species: [{gri30.yaml/species: [AR, N2, CH4, O2, CO2, H2O, CO, H2, OH]}]
    """
    request.cls.gas = ct.Solution(yaml=s)
    request.cls.gas.TP = 300, 101325
    request.cls.gas.set_equivalence_ratio(0.8, 'CH4', 'O2:1.0, N2:3.76')

@pytest.mark.usefixtures('extra_elements')
class TestEquilExtraElements:
    """
    Tests equilibrium with extra elements that are not involved in the reactions.
    """

    def test_auto(self):
        # Succeeds after falling back to VCS
        self.gas.equilibrate('TP')
        assert self.gas['CH4'].X[0] == approx(0.0)

    @pytest.mark.xfail
    def test_element_potential(self):
        self.gas.equilibrate('TP', solver='element_potential')
        assert self.gas['CH4'].X[0] == approx(0.0)

    def test_gibbs(self):
        self.gas.equilibrate('TP', solver='gibbs')
        assert self.gas['CH4'].X[0] == approx(0.0)

    def test_vcs(self):
        self.gas.equilibrate('TP', solver='vcs')
        assert self.gas['CH4'].X[0] == approx(0.0)


class TestVCS_EquilTest(EquilTestCases):
    """
    Tests using the 'vcs' solver.
    """

    solver = 'vcs'


@pytest.fixture(scope='function')
def koh_equil(request):
    phase_names = ['K_solid', 'K_liquid', 'KOH_a', 'KOH_b', 'KOH_liquid', 'K2O2_solid',
                'K2O_solid', 'KO2_solid', 'ice', 'liquid_water', 'KOH_plasma']
    phases = [ct.Solution('KOH.yaml', name) for name in phase_names]
    request.cls.mix = ct.Mixture(phases)

@pytest.mark.usefixtures('koh_equil')
class TestKOH_Equil:
    """
    Test roughly based on examples/multiphase/plasma_equilibrium.py
    """

    def test_equil_TP(self, test_data_path):
        temperatures = range(350, 5000, 300)
        data = np.zeros((len(temperatures), self.mix.n_species+1))
        data[:,0] = temperatures

        for i,T in enumerate(temperatures):
            self.mix.T = T
            self.mix.P = ct.one_atm
            self.mix.species_moles = 'K:1.03, H2:2.12, O2:0.9'
            self.mix.equilibrate('TP', solver='vcs')

            data[i,1:] = self.mix.species_moles

        # The reference values for this test are all completely non-physical, due to the
        # VCS solver extrapolating thermo polynomials outside of their valid range. See
        # https://github.com/Cantera/cantera/issues/270. The results show ice at
        # temperatures of over 1000 K, and liquid water for temperatures of 2000-5000 K.
        compare(data, test_data_path / "koh-equil-TP.csv")

    @pytest.mark.slow_test
    def test_equil_HP(self, test_data_path):
        temperatures = range(350, 5000, 300)
        data = np.zeros((len(temperatures), self.mix.n_species+2))
        data[:,0] = temperatures

        # The outer iteration for the temperature *should* be able to
        # converge from further away, but in practice, it can't. (Of course,
        # changing this value requires replacing the reference output)
        dT = 1
        self.mix.P = ct.one_atm

        for i,T in enumerate(temperatures):
            self.mix.species_moles = 'K:1.03, H2:2.12, O2:0.9'
            self.mix.T = T - dT
            self.mix.equilibrate('TP', solver='vcs')
            self.mix.T = T
            self.mix.equilibrate('HP', solver='vcs')

            data[i,1] = self.mix.T # equilibrated temperature
            data[i,2:] = self.mix.species_moles

        compare(data, test_data_path / "koh-equil-HP.csv")


@pytest.fixture(scope='function')
def carbon_equil(request):
    request.cls.gas = ct.Solution('gri30.yaml', transport_model=None)
    request.cls.carbon = ct.Solution("graphite.yaml")
    request.cls.fuel = 'CH4'
    request.cls.mix_phases = [(request.cls.gas, 1.0), (request.cls.carbon, 0.0)]
    request.cls.n_species = request.cls.gas.n_species + request.cls.carbon.n_species

@pytest.mark.usefixtures('carbon_equil')
class TestEquil_GasCarbon:
    "Test roughly based on examples/multiphase/adiabatic.py"

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, test_data_path):
        self.test_data_path = test_data_path

    @pytest.mark.parametrize("solver,P,ref_file,estimate_equil", [
        ('gibbs', ct.one_atm, "gas-carbon-equil.csv", 0),
        ('vcs', ct.one_atm, "gas-carbon-equil.csv", 0),
        ('vcs', ct.one_atm, "gas-carbon-equil.csv", -1),
        ('gibbs', 5 * ct.one_atm, "gas-carbon-equil-5atm.csv", 0),
        ('vcs', 5 * ct.one_atm, "gas-carbon-equil-5atm.csv", 0),
    ])
    def test_equilibrate(self, solver, P, ref_file, estimate_equil):
        n_points = 12
        T = 300
        data = np.zeros((n_points, 2+self.n_species))
        phi = np.linspace(0.3, 3.5, n_points)
        for i in range(n_points):
            self.gas.set_equivalence_ratio(phi[i], self.fuel,
                                           {'O2': 1.0, 'N2': 3.76})
            mix = ct.Mixture(self.mix_phases)
            mix.T = T
            mix.P = P

            # equilibrate the mixture adiabatically at constant P
            mix.equilibrate('HP', solver=solver, max_steps=1000,
                            estimate_equil=estimate_equil)
            data[i,:2] = (phi[i], mix.T)
            data[i,2:] = mix.species_moles

        compare(data, self.test_data_path / ref_file, rtol=1e-7)


class Test_IdealSolidSolnPhase_Equil:
    def test_equil(self):
        gas = ct.ThermoPhase("IdealSolidSolnPhaseExample.yaml")
        gas.TPX = 500, ct.one_atm, 'C2H2-graph: 1.0'

        gas.equilibrate('TP', solver='element_potential')
        assert gas['C-graph'].X[0] == approx(2.0 / 3.0)
        assert gas['H2-solute'].X[0] == approx(1.0 / 3.0)


class TestMultiphase_VaporLiquidHydrocarbon:
    """
    Regression test for element conservation in a two-phase vapor/liquid hydrocarbon
    mixture (C8H18, C6H6, C7H8). Previously, elements were not conserved after
    equilibration due to errors in MultiPhaseEquil with "trace" species.
    Regression coverage for issue #425.
    """

    _gas_yaml = """
    phases:
    - name: vapor
      thermo: ideal-gas
      elements: [O, H, C]
      species: [{nasa_gas.yaml/species: ["C8H18,n-octane", C6H6, C7H8]}]
      state: {T: 300.0, P: 1 atm}
    """

    _liquid_yaml = """
    phases:
    - name: liquid
      elements: [O, H, C]
      thermo: ideal-condensed
      species: ["C8H18(L),n-octa", C6H6(L), C7H8(L)]
      state: {T: 300.0, P: 1 atm}

    # Species data from nasa_condensed.yaml, augmented with (fake) density data
    species:
    - name: C8H18(L),n-octa
      composition: {C: 8, H: 18}
      thermo:
        model: NASA7
        temperature-ranges: [220.0, 300.0]
        data:
        - [71.413393, -0.5020795, 1.834199e-03, -2.0450165e-06, 0.0, -4.1243725e+04,
          -277.2224]
      equation-of-state:
        model: constant-volume
        density: 0.8765 g/cm^3
    - name: C6H6(L)
      composition: {C: 6, H: 6}
      thermo:
        model: NASA7
        temperature-ranges: [278.68, 500.0]
        data:
        - [63.6690229, -0.600534398, 2.6679281e-03, -5.06308828e-06, 3.63955562e-09,
          -1670.85472, -243.891797]
      equation-of-state:
        model: constant-volume
        density: 0.8765 g/cm^3
    - name: C7H8(L)
      composition: {C: 7, H: 8}
      thermo:
        model: NASA7
        temperature-ranges: [178.15, 500.0]
        data:
        - [29.3676022, -0.194722686, 9.74773096e-04, -1.91472689e-06, 1.48097019e-09,
          -4163.18442, -112.019966]
      equation-of-state:
        model: constant-volume
        density: 0.8765 g/cm^3
    """

    @pytest.mark.parametrize("n_C7H8_vapor", [1, 2, 3, 5, 9, 11, 20])
    def test_element_conservation(self, n_C7H8_vapor):
        liquid = ct.Solution(yaml=self._liquid_yaml)
        vapor = ct.Solution(yaml=self._gas_yaml)

        liquid.Y = {'C7H8(L)': 1.0, 'C6H6(L)': 1.0}
        vapor.Y = {'C7H8': n_C7H8_vapor, 'C6H6': 1.0}
        mix = ct.Mixture([(vapor, 1.0), (liquid, 1.0)])

        mix.T = 373.15
        mix.P = 101325.0
        elements_before = [mix.element_moles(e) for e in range(mix.n_elements)]
        mix.equilibrate('TP', solver='gibbs', max_steps=5000)
        elements_after = [mix.element_moles(e) for e in range(mix.n_elements)]

        assert elements_after == approx(elements_before, rel=1e-7)


@pytest.mark.parametrize("solver,include_h2o",
                         [("gibbs", False), ("vcs", False),
                          ("gibbs", True), ("vcs", True)])
def test_excluded_species_extra_element(solver, include_h2o):
    """
    Test equilibration of a Pb/O/C mixture where one of the condensed phases (PbO(yw))
    has thermo data only valid above 762 K. When equilibrating at 400 K, the solver
    should exclude PbO(yw) and redistribute its elemental contribution to valid species,
    regardless of whether PbO(yw) appears in the initial composition.

    Also test a case when a mixture includes a phase (H2O_s) that introduces an
    element (H) with zero abundance and only one phase containing that element.

    Adapted from issue report https://github.com/Cantera/cantera/issues/160.
    """
    _pb_phases_yaml_with_H2O = """
    phases:
    - name: gas
      thermo: ideal-gas
      species: [{nasa_gas.yaml/species: [Pb, O, C]}]
    - name: H2O_s
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [H2O(s)]}]
    - name: C_gr
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [C(gr)]}]
    - name: Pb_s
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [Pb(cr)]}]
    - name: PbO_rd
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [PbO(rd)]}]
    - name: PbO_yw
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [PbO(yw)]}]
    - name: PbO2_s
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [PbO2(s)]}]
    - name: Pb3O4_s
      thermo: fixed-stoichiometry
      species: [{nasa_condensed.yaml/species: [Pb3O4(s)]}]
    """
    if include_h2o:
        phase_names = ['gas', 'H2O_s', 'C_gr', 'Pb_s', 'PbO_rd', 'PbO_yw', 'PbO2_s',
                       'Pb3O4_s']
    else:
        phase_names = ['gas', 'C_gr', 'Pb_s', 'PbO_rd', 'PbO_yw', 'PbO2_s', 'Pb3O4_s']
    phases = [ct.Solution(yaml=_pb_phases_yaml_with_H2O, name=n) for n in phase_names]
    mix = ct.Mixture(phases)
    mix.T = 400
    mix.P = ct.one_atm
    mix.species_moles = 'O:3, Pb:1, C:1'
    mix.equilibrate('TP', solver='gibbs', rtol=1e-10)
    phase_moles_ref = mix.phase_moles()

    # Initial composition includes PbO(yw), which has no valid thermo at 400 K and is
    # excluded from the calculation. The solver should redistribute its Pb and O atoms
    # to valid species and converge to the same equilibrium as when initialized with
    # elemental species (mix1).
    mix.T = 400
    mix.P = ct.one_atm
    mix.species_moles = 'O:2, PbO(yw):1, C:1'  # same element totals
    mix.equilibrate('TP', solver=solver, rtol=1e-10, max_steps=10000, max_iter=10000)

    # PbO2(s) and C(gr) should be present; PbO(yw) excluded and absent; H2O(s) absent
    if include_h2o:
        assert mix.phase_moles("H2O_s") == approx(0.0)
    assert mix.phase_moles("PbO_yw") == approx(0.0)
    assert mix.phase_moles("PbO2_s") == approx(1.0, rel=1e-6)
    assert mix.phase_moles() == approx(phase_moles_ref, rel=1e-6)
