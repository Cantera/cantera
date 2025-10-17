import pytest
from pytest import approx

import cantera as ct


class TestMixture:

    @pytest.fixture(scope='class')
    def phase1(self):
        return ct.Solution('h2o2.yaml', transport_model=None)

    @pytest.fixture(scope='class')
    def phase2(self):
        return ct.Solution('air.yaml')

    @pytest.fixture(scope='function')
    def mix(request, phase1, phase2):
        # Create instance-level data for each test
        return ct.Mixture([(phase1, 1.0), (phase2, 2.0)])

    def test_sizes(self, mix, phase1, phase2):
        assert mix.n_phases == 2

        assert mix.n_species == phase1.n_species + phase2.n_species

        E = set(phase1.element_names) | set(phase2.element_names)
        assert len(E) == mix.n_elements

    def test_element_index(self, mix):
        with pytest.raises(ct.CanteraError, match="Element 'W' not found"):
            mix.element_index('W')

        with pytest.raises(ct.CanteraError, match="outside valid range"):
            mix.element_index(41)

        with pytest.raises(TypeError, match='must be a string or a number'):
            mix.element_index(None)

    def test_speciesIndex(self, mix, phase1, phase2):
        names = mix.species_names
        kOH = names.index('OH')
        kN2 = names.index('N2O')
        assert mix.species_name(kOH) == 'OH'
        assert mix.species_name(kN2) == 'N2O'

        assert mix.species_index(0, 'OH') == kOH
        assert mix.species_index(phase1, 'OH') == kOH
        assert mix.species_index(phase1.name, 'OH') == kOH
        assert mix.species_index(0, phase1.species_index('OH')) == kOH
        assert mix.species_index(1, phase2.species_index('N2O')) == kN2
        assert mix.species_index(1, 'N2O') == kN2

        with pytest.raises(IndexError, match='out of range'):
            mix.species_index(3, 'OH')

        with pytest.raises(ct.CanteraError, match="Species 'OH' not found"):
            mix.species_index(1, 'OH')

        with pytest.raises(ct.CanteraError, match="outside valid range"):
            mix.species_index(0, -2)

        with pytest.raises(ct.CanteraError, match="Species 'CO2' not found"):
            mix.species_index(1, 'CO2')

    def test_n_atoms(self, mix):
        names = mix.species_names
        kOH = names.index('OH')
        kN2 = names.index('N2')
        mH = mix.element_index('H')
        mN = mix.element_index('N')

        assert mix.n_atoms(kOH, 'H') == 1
        assert mix.n_atoms(kOH, 'O') == 1
        assert mix.n_atoms(kOH, mH) == 1
        assert mix.n_atoms(kOH, mN) == 0

        assert mix.n_atoms(kN2, mN) == 2
        assert mix.n_atoms(kN2, mH) == 0

    def test_phase(self, mix, phase1, phase2):
        assert phase1 == mix.phase(0)
        assert phase2 == mix.phase(1)

        phaseNames = mix.phase_names
        assert len(phaseNames) == mix.n_phases
        assert phaseNames[0] == phase1.name
        assert phaseNames[1] == phase2.name

    def test_phase_index(self, mix, phase1, phase2):
        assert mix.phase_index(phase1) == 0
        assert mix.phase_index(phase2) == 1
        assert mix.phase_index(phase2.name) == 1
        assert mix.phase_index(1) == 1

        with pytest.raises(KeyError):
            mix.phase_index('foobar')

        with pytest.raises(IndexError):
            mix.phase_index(2)

    def test_properties(self, mix):
        mix.T = 350
        assert mix.T == 350

        mix.P = 2e5
        assert mix.P == 2e5
        assert mix.T == 350

        assert mix.max_temp > mix.min_temp

    def test_charge(self, mix):
        C = sum(mix.phase_charge(i) for i in range(mix.n_phases))
        assert mix.charge == C

    def test_phase_moles(self, mix):
        M = mix.phase_moles()
        assert M[0] == mix.phase_moles(0)
        assert M[1] == mix.phase_moles('air')

        mix.set_phase_moles('air', 4)
        assert mix.phase_moles(1) == 4

    def test_species_moles(self, mix, phase1):
        mix.species_moles = 'H2:1.0, NO2:4.0'
        P = mix.phase_moles()
        S = mix.species_moles

        assert P[0] == 1
        assert P[1] == 4

        assert S[mix.species_index(0, 'H2')] == 1
        assert S[mix.species_index(1, 'NO2')] == 4

        S[2] = 7
        mix.species_moles = S
        assert mix.species_moles[2] == approx(S[2])
        assert mix.phase_moles(0) == approx(sum(S[:phase1.n_species]))

        with pytest.raises(ValueError):
            mix.species_moles = (1, 2, 3)

        with pytest.raises(TypeError):
            mix.species_moles = 9

    def test_element_moles(self, mix):
        mix.species_moles = 'H2:1.0, OH:4.0'

        assert mix.element_moles('H') == approx(6)
        assert mix.element_moles('O') == approx(4)
        assert mix.element_moles('N') == approx(0)

    def test_chemical_potentials(self, mix, phase1, phase2):
        C = mix.chemical_potentials
        C1 = phase1.chemical_potentials
        C2 = phase2.chemical_potentials

        assert C[:phase1.n_species] == approx(C1)
        assert C[phase1.n_species:] == approx(C2)

    def test_equilibrate1(self, mix):
        mix.species_moles = 'H2:1.0, O2:0.5, N2:1.0'
        mix.T = 400
        mix.P = 2 * ct.one_atm

        E1 = [mix.element_moles(m) for m in range(mix.n_elements)]
        mix.equilibrate('TP', solver='vcs', estimate_equil=-1)

        E2 = [mix.element_moles(m) for m in range(mix.n_elements)]
        assert E1 == approx(E2)
        assert mix.T == approx(400)
        assert mix.P == approx(2 * ct.one_atm)

    @pytest.mark.xfail(reason="See https://github.com/Cantera/cantera/issues/1023")
    def test_equilibrate2(self, mix):
        mix.species_moles = 'H2:1.0, O2:0.5, N2:1.0'
        mix.T = 400
        mix.P = 2 * ct.one_atm

        E1 = [mix.element_moles(m) for m in range(mix.n_elements)]
        mix.equilibrate('TP', solver='gibbs')

        E2 = [mix.element_moles(m) for m in range(mix.n_elements)]
        assert E1 == approx(E2)
        assert mix.T == approx(400)
        assert mix.P == approx(2 * ct.one_atm)

    def test_invalid_property(self, mix):
        x = mix
        with pytest.raises(AttributeError):
            x.foobar = 300
        with pytest.raises(AttributeError):
            x.foobar

    def test_invalid_phase_type(self, mix, phase1):
        water = ct.Water()
        with pytest.raises(ct.CanteraError, match='not compatible'):
            mix = ct.Mixture([(phase1, 1.0), (water, 2.0)])
