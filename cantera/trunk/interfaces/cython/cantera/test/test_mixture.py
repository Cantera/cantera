import unittest
import cantera as ct
from . import utilities


class TestMixture(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        cls.phase1 = ct.Solution('h2o2.xml')
        cls.phase2 = ct.Solution('air.xml')

    def setUp(self):
        self.mix = ct.Mixture([(self.phase1, 1.0), (self.phase2, 2.0)])

    def test_sizes(self):
        self.assertEqual(self.mix.nPhases, 2)

        self.assertEqual(self.mix.nSpecies,
                         self.phase1.nSpecies + self.phase2.nSpecies)

        E = set(self.phase1.elementNames) | set(self.phase2.elementNames)
        self.assertEqual(len(E), self.mix.nElements)

    def test_elementIndex(self):
        m_H = self.mix.elementIndex('H')
        self.assertEqual(m_H, self.mix.elementIndex(m_H))

        with self.assertRaises(ValueError):
            self.mix.elementIndex('W')

        with self.assertRaises(ValueError):
            self.mix.elementIndex(41)

        with self.assertRaises(TypeError):
            self.mix.elementIndex(None)

    def test_speciesIndex(self):
        names = self.mix.speciesNames
        kOH = names.index('OH')
        kN2 = names.index('N2')
        self.assertEqual(self.mix.speciesName(kOH), 'OH')
        self.assertEqual(self.mix.speciesName(kN2), 'N2')

        self.assertEqual(self.mix.speciesIndex(0, 'OH'), kOH)
        self.assertEqual(self.mix.speciesIndex(self.phase1, 'OH'), kOH)
        self.assertEqual(self.mix.speciesIndex(self.phase1.name, 'OH'), kOH)
        self.assertEqual(self.mix.speciesIndex(0, self.phase1.speciesIndex('OH')), kOH)
        self.assertEqual(self.mix.speciesIndex(1, self.phase2.speciesIndex('N2')), kN2)
        self.assertEqual(self.mix.speciesIndex(1, 'N2'), kN2)

        with self.assertRaises(IndexError):
            self.mix.speciesIndex(3, 'OH')

        with self.assertRaises(ValueError):
            self.mix.speciesIndex(1, 'OH')

        with self.assertRaises(ValueError):
            self.mix.speciesIndex(0, -2)

        with self.assertRaises(ValueError):
            self.mix.speciesIndex(1, 'CO2')

    def test_nAtoms(self):
        names = self.mix.speciesNames
        kOH = names.index('OH')
        kN2 = names.index('N2')
        mH = self.mix.elementIndex('H')
        mN = self.mix.elementIndex('N')

        self.assertEqual(self.mix.nAtoms(kOH, 'H'), 1)
        self.assertEqual(self.mix.nAtoms(kOH, 'O'), 1)
        self.assertEqual(self.mix.nAtoms(kOH, mH), 1)
        self.assertEqual(self.mix.nAtoms(kOH, mN), 0)

        self.assertEqual(self.mix.nAtoms(kN2, mN), 2)
        self.assertEqual(self.mix.nAtoms(kN2, mH), 0)

    def test_phase(self):
        self.assertEqual(self.phase1, self.mix.phase(0))
        self.assertEqual(self.phase2, self.mix.phase(1))

        phaseNames = self.mix.phaseNames
        self.assertEqual(len(phaseNames), self.mix.nPhases)
        self.assertEqual(phaseNames[0], self.phase1.name)
        self.assertEqual(phaseNames[1], self.phase2.name)

    def test_phaseIndex(self):
        self.assertEqual(self.mix.phaseIndex(self.phase1), 0)
        self.assertEqual(self.mix.phaseIndex(self.phase2), 1)
        self.assertEqual(self.mix.phaseIndex(self.phase2.name), 1)
        self.assertEqual(self.mix.phaseIndex(1), 1)

        with self.assertRaises(KeyError):
            self.mix.phaseIndex('foobar')

        with self.assertRaises(IndexError):
            self.mix.phaseIndex(2)

    def test_properties(self):
        self.mix.T = 350
        self.assertEqual(self.mix.T, 350)

        self.mix.P = 2e5
        self.assertEqual(self.mix.P, 2e5)
        self.assertEqual(self.mix.T, 350)

        self.assertGreater(self.mix.maxTemp, self.mix.minTemp)

    def test_charge(self):
        C = sum(self.mix.phaseCharge(i) for i in range(self.mix.nPhases))
        self.assertEqual(self.mix.charge, C)

    def test_phaseMoles(self):
        M = self.mix.phaseMoles()
        self.assertEqual(M[0], self.mix.phaseMoles(0))
        self.assertEqual(M[1], self.mix.phaseMoles('air'))

        self.mix.setPhaseMoles('air', 4)
        self.assertEqual(self.mix.phaseMoles(1), 4)

    def test_speciesMoles(self):
        self.mix.setSpeciesMoles('H2:1.0, N2:4.0')
        P = self.mix.phaseMoles()
        S = self.mix.speciesMoles()

        self.assertEqual(P[0], 1)
        self.assertEqual(P[1], 4)

        self.assertEqual(S[self.mix.speciesIndex(0, 'H2')], 1)
        self.assertEqual(S[self.mix.speciesIndex(1, 'N2')], 4)

        S[2] = 7
        self.mix.setSpeciesMoles(S)
        self.assertNear(self.mix.speciesMoles(2), S[2])
        self.assertNear(self.mix.phaseMoles(0), sum(S[:self.phase1.nSpecies]))

        with self.assertRaises(ValueError):
            self.mix.setSpeciesMoles((1,2,3))

        with self.assertRaises(TypeError):
            self.mix.setSpeciesMoles(9)

    def test_elementMoles(self):
        self.mix.setSpeciesMoles('H2:1.0, OH:4.0')

        self.assertNear(self.mix.elementMoles('H'), 6)
        self.assertNear(self.mix.elementMoles('O'), 4)
        self.assertNear(self.mix.elementMoles('N'), 0)

    def test_chem_potentials(self):
        C = self.mix.chem_potentials
        C1 = self.phase1.chem_potentials
        C2 = self.phase2.chem_potentials

        self.assertArrayNear(C[:self.phase1.nSpecies], C1)
        self.assertArrayNear(C[self.phase1.nSpecies:], C2)

    def test_equilibrate1(self):
        self.mix.setSpeciesMoles('H2:1.0, O2:0.5, N2:1.0')
        self.mix.T = 400
        self.mix.P = 2 * ct.OneAtm

        E1 = [self.mix.elementMoles(m) for m in range(self.mix.nElements)]
        self.mix.equilibrate('TP')

        E2 = [self.mix.elementMoles(m) for m in range(self.mix.nElements)]
        self.assertArrayNear(E1, E2)
        self.assertNear(self.mix.T, 400)
        self.assertNear(self.mix.P, 2 * ct.OneAtm)

    def test_equilibrate2(self):
        self.mix.setSpeciesMoles('H2:1.0, O2:0.5, N2:1.0')
        self.mix.T = 400
        self.mix.P = 2 * ct.OneAtm

        E1 = [self.mix.elementMoles(m) for m in range(self.mix.nElements)]
        self.mix.equilibrate('TP', solver='gibbs')

        E2 = [self.mix.elementMoles(m) for m in range(self.mix.nElements)]
        self.assertArrayNear(E1, E2)
        self.assertNear(self.mix.T, 400)
        self.assertNear(self.mix.P, 2 * ct.OneAtm)
