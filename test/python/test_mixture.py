import pytest
import cantera as ct
from . import utilities
from .utilities import (
    assertNear,
    assertArrayNear
)

@pytest.fixture(scope='class')
def phases(request):
    # Create class-level data for each test
    request.cls.phase1 = ct.Solution('h2o2.yaml', transport_model=None)
    request.cls.phase2 = ct.Solution('air.yaml')

@pytest.fixture(scope='function')
def mixture(request, phases):
    # Create instance-level data for each test
    request.cls.mix = ct.Mixture([(request.cls.phase1, 1.0), (request.cls.phase2, 2.0)])

@pytest.mark.usefixtures("mixture")
class TestMixture:

    def test_sizes(self):
        assert self.mix.n_phases == 2

        assert self.mix.n_species == self.phase1.n_species + self.phase2.n_species

        E = set(self.phase1.element_names) | set(self.phase2.element_names)
        assert len(E) == self.mix.n_elements

    def test_element_index(self):
        m_H = self.mix.element_index('H')
        assert m_H == self.mix.element_index(m_H)

        with pytest.raises(ValueError, match='No such element'):
            self.mix.element_index('W')

        with pytest.raises(ValueError, match='No such element'):
            self.mix.element_index(41)

        with pytest.raises(TypeError, match='must be a string or a number'):
            self.mix.element_index(None)

    def test_speciesIndex(self):
        names = self.mix.species_names
        kOH = names.index('OH')
        kN2 = names.index('N2O')
        assert self.mix.species_name(kOH) == 'OH'
        assert self.mix.species_name(kN2) == 'N2O'

        assert self.mix.species_index(0, 'OH') == kOH
        assert self.mix.species_index(self.phase1, 'OH') == kOH
        assert self.mix.species_index(self.phase1.name, 'OH') == kOH
        assert self.mix.species_index(0, self.phase1.species_index('OH')) == kOH
        assert self.mix.species_index(1, self.phase2.species_index('N2O')) == kN2
        assert self.mix.species_index(1, 'N2O') == kN2

        with pytest.raises(IndexError, match='out of range'):
            self.mix.species_index(3, 'OH')

        with pytest.raises(ValueError, match='No such species'):
            self.mix.species_index(1, 'OH')

        with pytest.raises(ValueError, match='out of range'):
            self.mix.species_index(0, -2)

        with pytest.raises(ValueError, match='No such species'):
            self.mix.species_index(1, 'CO2')

    def test_n_atoms(self):
        names = self.mix.species_names
        kOH = names.index('OH')
        kN2 = names.index('N2')
        mH = self.mix.element_index('H')
        mN = self.mix.element_index('N')

        assert self.mix.n_atoms(kOH, 'H') == 1
        assert self.mix.n_atoms(kOH, 'O') == 1
        assert self.mix.n_atoms(kOH, mH) == 1
        assert self.mix.n_atoms(kOH, mN) == 0

        assert self.mix.n_atoms(kN2, mN) == 2
        assert self.mix.n_atoms(kN2, mH) == 0

    def test_phase(self):
        assert self.phase1 == self.mix.phase(0)
        assert self.phase2 == self.mix.phase(1)

        phaseNames = self.mix.phase_names
        assert len(phaseNames) == self.mix.n_phases
        assert phaseNames[0] == self.phase1.name
        assert phaseNames[1] == self.phase2.name

    def test_phase_index(self):
        assert self.mix.phase_index(self.phase1) == 0
        assert self.mix.phase_index(self.phase2) == 1
        assert self.mix.phase_index(self.phase2.name) == 1
        assert self.mix.phase_index(1) == 1

        with pytest.raises(KeyError):
            self.mix.phase_index('foobar')

        with pytest.raises(IndexError):
            self.mix.phase_index(2)

    def test_properties(self):
        self.mix.T = 350
        assert self.mix.T == 350

        self.mix.P = 2e5
        assert self.mix.P == 2e5
        assert self.mix.T == 350

        assert self.mix.max_temp > self.mix.min_temp

    def test_charge(self):
        C = sum(self.mix.phase_charge(i) for i in range(self.mix.n_phases))
        assert self.mix.charge == C

    def test_phase_moles(self):
        M = self.mix.phase_moles()
        assert M[0] == self.mix.phase_moles(0)
        assert M[1] == self.mix.phase_moles('air')

        self.mix.set_phase_moles('air', 4)
        assert self.mix.phase_moles(1) == 4

    def test_species_moles(self):
        self.mix.species_moles = 'H2:1.0, NO2:4.0'
        P = self.mix.phase_moles()
        S = self.mix.species_moles

        assert P[0] == 1
        assert P[1] == 4

        assert S[self.mix.species_index(0, 'H2')] == 1
        assert S[self.mix.species_index(1, 'NO2')] == 4

        S[2] = 7
        self.mix.species_moles = S
        assertNear(self.mix.species_moles[2], S[2])
        assertNear(self.mix.phase_moles(0), sum(S[:self.phase1.n_species]))

        with pytest.raises(ValueError):
            self.mix.species_moles = (1, 2, 3)

        with pytest.raises(TypeError):
            self.mix.species_moles = 9

    def test_element_moles(self):
        self.mix.species_moles = 'H2:1.0, OH:4.0'

        assertNear(self.mix.element_moles('H'), 6)
        assertNear(self.mix.element_moles('O'), 4)
        assertNear(self.mix.element_moles('N'), 0)

    def test_chemical_potentials(self):
        C = self.mix.chemical_potentials
        C1 = self.phase1.chemical_potentials
        C2 = self.phase2.chemical_potentials

        assertArrayNear(C[:self.phase1.n_species], C1)
        assertArrayNear(C[self.phase1.n_species:], C2)

    def test_equilibrate1(self):
        self.mix.species_moles = 'H2:1.0, O2:0.5, N2:1.0'
        self.mix.T = 400
        self.mix.P = 2 * ct.one_atm

        E1 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        self.mix.equilibrate('TP', solver='vcs', estimate_equil=-1)

        E2 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        assertArrayNear(E1, E2)
        assertNear(self.mix.T, 400)
        assertNear(self.mix.P, 2 * ct.one_atm)

    @pytest.mark.xfail(reason="See https://github.com/Cantera/cantera/issues/1023")
    def test_equilibrate2(self):
        self.mix.species_moles = 'H2:1.0, O2:0.5, N2:1.0'
        self.mix.T = 400
        self.mix.P = 2 * ct.one_atm

        E1 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        self.mix.equilibrate('TP', solver='gibbs')

        E2 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        assertArrayNear(E1, E2)
        assertNear(self.mix.T, 400)
        assertNear(self.mix.P, 2 * ct.one_atm)

    def test_invalid_property(self):
        x = self.mix
        with pytest.raises(AttributeError):
            x.foobar = 300
        with pytest.raises(AttributeError):
            x.foobar

    def test_invalid_phase_type(self):
        water = ct.Water()
        with pytest.raises(ct.CanteraError, match='not compatible'):
            self.mix = ct.Mixture([(self.phase1, 1.0), (water, 2.0)])
