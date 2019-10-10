import cantera as ct
from . import utilities


class TestMixture(utilities.CanteraTest):
    @classmethod
    def setUpClass(self):
        utilities.CanteraTest.setUpClass()
        self.phase1 = ct.Solution('h2o2.xml')
        self.phase2 = ct.Solution('air.xml')

    def setUp(self):
        self.mix = ct.Mixture([(self.phase1, 1.0), (self.phase2, 2.0)])

    def test_sizes(self):
        self.assertEqual(self.mix.n_phases, 2)

        self.assertEqual(self.mix.n_species,
                         self.phase1.n_species + self.phase2.n_species)

        E = set(self.phase1.element_names) | set(self.phase2.element_names)
        self.assertEqual(len(E), self.mix.n_elements)

    def test_element_index(self):
        m_H = self.mix.element_index('H')
        self.assertEqual(m_H, self.mix.element_index(m_H))

        with self.assertRaisesRegex(ValueError, 'No such element'):
            self.mix.element_index('W')

        with self.assertRaisesRegex(ValueError, 'No such element'):
            self.mix.element_index(41)

        with self.assertRaisesRegex(TypeError, 'must be a string or a number'):
            self.mix.element_index(None)

    def test_speciesIndex(self):
        names = self.mix.species_names
        kOH = names.index('OH')
        kN2 = names.index('N2')
        self.assertEqual(self.mix.species_name(kOH), 'OH')
        self.assertEqual(self.mix.species_name(kN2), 'N2')

        self.assertEqual(self.mix.species_index(0, 'OH'), kOH)
        self.assertEqual(self.mix.species_index(self.phase1, 'OH'), kOH)
        self.assertEqual(self.mix.species_index(self.phase1.name, 'OH'), kOH)
        self.assertEqual(self.mix.species_index(0, self.phase1.species_index('OH')), kOH)
        self.assertEqual(self.mix.species_index(1, self.phase2.species_index('N2')), kN2)
        self.assertEqual(self.mix.species_index(1, 'N2'), kN2)

        with self.assertRaisesRegex(IndexError, 'out of range'):
            self.mix.species_index(3, 'OH')

        with self.assertRaisesRegex(ValueError, 'No such species'):
            self.mix.species_index(1, 'OH')

        with self.assertRaisesRegex(ValueError, 'out of range'):
            self.mix.species_index(0, -2)

        with self.assertRaisesRegex(ValueError, 'No such species'):
            self.mix.species_index(1, 'CO2')

    def test_n_atoms(self):
        names = self.mix.species_names
        kOH = names.index('OH')
        kN2 = names.index('N2')
        mH = self.mix.element_index('H')
        mN = self.mix.element_index('N')

        self.assertEqual(self.mix.n_atoms(kOH, 'H'), 1)
        self.assertEqual(self.mix.n_atoms(kOH, 'O'), 1)
        self.assertEqual(self.mix.n_atoms(kOH, mH), 1)
        self.assertEqual(self.mix.n_atoms(kOH, mN), 0)

        self.assertEqual(self.mix.n_atoms(kN2, mN), 2)
        self.assertEqual(self.mix.n_atoms(kN2, mH), 0)

    def test_phase(self):
        self.assertEqual(self.phase1, self.mix.phase(0))
        self.assertEqual(self.phase2, self.mix.phase(1))

        phaseNames = self.mix.phase_names
        self.assertEqual(len(phaseNames), self.mix.n_phases)
        self.assertEqual(phaseNames[0], self.phase1.name)
        self.assertEqual(phaseNames[1], self.phase2.name)

    def test_phase_index(self):
        self.assertEqual(self.mix.phase_index(self.phase1), 0)
        self.assertEqual(self.mix.phase_index(self.phase2), 1)
        self.assertEqual(self.mix.phase_index(self.phase2.name), 1)
        self.assertEqual(self.mix.phase_index(1), 1)

        with self.assertRaises(KeyError):
            self.mix.phase_index('foobar')

        with self.assertRaises(IndexError):
            self.mix.phase_index(2)

    def test_properties(self):
        self.mix.T = 350
        self.assertEqual(self.mix.T, 350)

        self.mix.P = 2e5
        self.assertEqual(self.mix.P, 2e5)
        self.assertEqual(self.mix.T, 350)

        self.assertGreater(self.mix.max_temp, self.mix.min_temp)

    def test_charge(self):
        C = sum(self.mix.phase_charge(i) for i in range(self.mix.n_phases))
        self.assertEqual(self.mix.charge, C)

    def test_phase_moles(self):
        M = self.mix.phase_moles()
        self.assertEqual(M[0], self.mix.phase_moles(0))
        self.assertEqual(M[1], self.mix.phase_moles('air'))

        self.mix.set_phase_moles('air', 4)
        self.assertEqual(self.mix.phase_moles(1), 4)

    def test_species_moles(self):
        self.mix.species_moles = 'H2:1.0, N2:4.0'
        P = self.mix.phase_moles()
        S = self.mix.species_moles

        self.assertEqual(P[0], 1)
        self.assertEqual(P[1], 4)

        self.assertEqual(S[self.mix.species_index(0, 'H2')], 1)
        self.assertEqual(S[self.mix.species_index(1, 'N2')], 4)

        S[2] = 7
        self.mix.species_moles = S
        self.assertNear(self.mix.species_moles[2], S[2])
        self.assertNear(self.mix.phase_moles(0), sum(S[:self.phase1.n_species]))

        with self.assertRaises(ValueError):
            self.mix.species_moles = (1,2,3)

        with self.assertRaises(TypeError):
            self.mix.species_moles = 9

    def test_element_moles(self):
        self.mix.species_moles = 'H2:1.0, OH:4.0'

        self.assertNear(self.mix.element_moles('H'), 6)
        self.assertNear(self.mix.element_moles('O'), 4)
        self.assertNear(self.mix.element_moles('N'), 0)

    def test_chemical_potentials(self):
        C = self.mix.chemical_potentials
        C1 = self.phase1.chemical_potentials
        C2 = self.phase2.chemical_potentials

        self.assertArrayNear(C[:self.phase1.n_species], C1)
        self.assertArrayNear(C[self.phase1.n_species:], C2)

    def test_equilibrate1(self):
        self.mix.species_moles = 'H2:1.0, O2:0.5, N2:1.0'
        self.mix.T = 400
        self.mix.P = 2 * ct.one_atm

        E1 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        self.mix.equilibrate('TP', solver='vcs', estimate_equil=-1)

        E2 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        self.assertArrayNear(E1, E2)
        self.assertNear(self.mix.T, 400)
        self.assertNear(self.mix.P, 2 * ct.one_atm)

    def test_equilibrate2(self):
        self.mix.species_moles = 'H2:1.0, O2:0.5, N2:1.0'
        self.mix.T = 400
        self.mix.P = 2 * ct.one_atm

        E1 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        self.mix.equilibrate('TP', solver='gibbs')

        E2 = [self.mix.element_moles(m) for m in range(self.mix.n_elements)]
        self.assertArrayNear(E1, E2)
        self.assertNear(self.mix.T, 400)
        self.assertNear(self.mix.P, 2 * ct.one_atm)

    def test_invalid_property(self):
        x = self.mix
        with self.assertRaises(AttributeError):
            x.foobar = 300
        with self.assertRaises(AttributeError):
            x.foobar

    def test_invalid_phase_type(self):
        water = ct.Water()
        with self.assertRaisesRegex(ct.CanteraError, 'not compatible'):
            self.mix = ct.Mixture([(self.phase1, 1.0), (water, 2.0)])
