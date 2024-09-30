import numpy as np
import pytest
from pytest import approx

import cantera as ct
from .utilities import (
    compare
)

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


class TestMultiphaseEquil(EquilTestCases):
    """
    Tests using the 'gibbs' solver.
    """

    solver = 'gibbs'

    @pytest.mark.xfail
    def test_equil_gri_stoichiometric(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 301, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    @pytest.mark.xfail
    def test_equil_gri_lean(self):
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 301, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)


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
    phases = ct.import_phases("KOH.yaml",
            ['K_solid', 'K_liquid', 'KOH_a', 'KOH_b', 'KOH_liquid',
             'K2O2_solid', 'K2O_solid', 'KO2_solid', 'ice', 'liquid_water',
             'KOH_plasma'])
    request.cls.mix = ct.Mixture(phases)

@pytest.mark.usefixtures('koh_equil')
class TestKOH_Equil:
    """
    Test roughly based on examples/multiphase/plasma_equilibrium.py
    """

    def test_equil_TP(self):
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
        compare(data, self.test_data_path / "koh-equil-TP.csv")

    @pytest.mark.slow_test
    def test_equil_HP(self):
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

        compare(data, self.test_data_path / "koh-equil-HP.csv")


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

    def solve(self, solver, **kwargs):
        n_points = 12
        T = 300
        P = 101325
        data = np.zeros((n_points, 2+self.n_species))
        phi = np.linspace(0.3, 3.5, n_points)
        for i in range(n_points):
            self.gas.set_equivalence_ratio(phi[i], self.fuel,
                                           {'O2': 1.0, 'N2': 3.76})
            mix = ct.Mixture(self.mix_phases)
            mix.T = T
            mix.P = P

            # equilibrate the mixture adiabatically at constant P
            mix.equilibrate('HP', solver=solver, max_steps=1000, **kwargs)
            data[i,:2] = (phi[i], mix.T)
            data[i,2:] = mix.species_moles

        compare(data, self.test_data_path / "gas-carbon-equil.csv")

    @pytest.mark.slow_test
    def test_gibbs(self):
        self.solve('gibbs')

    @pytest.mark.slow_test
    def test_vcs(self):
        self.solve('vcs')

    def test_vcs_est(self):
        self.solve('vcs', estimate_equil=-1)


class Test_IdealSolidSolnPhase_Equil:
    def test_equil(self):
        gas = ct.ThermoPhase("IdealSolidSolnPhaseExample.yaml")
        gas.TPX = 500, ct.one_atm, 'C2H2-graph: 1.0'

        gas.equilibrate('TP', solver='element_potential')
        assert gas['C-graph'].X[0] == approx(2.0 / 3.0)
        assert gas['H2-solute'].X[0] == approx(1.0 / 3.0)
