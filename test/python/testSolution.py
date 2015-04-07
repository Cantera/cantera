import unittest

import numpy as np

import Cantera as ct

class ImportTest(unittest.TestCase):
    """
    Test the various ways of creating a Solution object
    """
    def check(self, gas, name, T, P, nSpec, nElem):
        self.assertEqual(gas.name(), name)
        self.assertAlmostEqual(gas.temperature(), T)
        self.assertAlmostEqual(gas.pressure(), P)
        self.assertEqual(gas.nSpecies(), nSpec)
        self.assertEqual(gas.nElements(), nElem)

    def test_importPhase_cti(self):
        gas1 = ct.importPhase('../data/air-no-reactions.cti', 'air')
        self.check(gas1, 'air', 300, 101325, 8, 3)

        gas2 = ct.importPhase('../data/air-no-reactions.cti', 'notair')
        self.check(gas2, 'notair', 900, 5*101325, 7, 2)

    def test_importPhase_cti2(self):
        # This should import the first phase, i.e. 'air'
        gas = ct.importPhase('../data/air-no-reactions.cti')
        self.check(gas, 'air', 300, 101325, 8, 3)

    def test_importPhase_xml(self):
        gas1 = ct.importPhase('../data/air-no-reactions.xml', 'air')
        self.check(gas1, 'air', 300, 101325, 8, 3)

        gas2 = ct.importPhase('../data/air-no-reactions.xml', 'notair')
        self.check(gas2, 'notair', 900, 5*101325, 7, 2)

    def test_import_GRI30(self):
        gas = ct.GRI30()
        self.check(gas, 'gri30', 300, 101325, 53, 5)

    def test_checkReactionBalance(self):
        self.assertRaises(Exception,
                          lambda: ct.IdealGasMix('../data/h2o2_unbalancedReaction.xml'))


class BasicTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gas = ct.importPhase('../data/air-no-reactions.xml', 'air')

    def setUp(self):
        # Workaround for Python 2.6 unittest, which does not call setUpClass
        try:
            self.gas
        except AttributeError:
            self.setUpClass()


    def test_counts(self):
        self.assertEqual(self.gas.nElements(), 3)
        self.assertEqual(self.gas.nSpecies(), 8)
        self.assertEqual(self.gas.nPhases(), 1)

    def test_elements(self):
        self.assertEqual(self.gas.nElements(), len(self.gas.elementNames()))
        for i,name in enumerate(self.gas.elementNames()):
            self.assertEqual(self.gas.elementName(i), name)
            self.assertEqual(self.gas.elementIndex(name), i)

    def test_species(self):
        self.assertEqual(self.gas.nSpecies(), len(self.gas.speciesNames()))
        for i,name in enumerate(self.gas.speciesNames()):
            self.assertEqual(self.gas.speciesName(i), name)
            self.assertEqual(self.gas.speciesIndex(name), i)

    def test_nAtoms(self):
        data = [(1, 'O', 'O'), (2, 'O', 'O2'), (1, 'N', 'NO'),
                (1, 'O', 'NO'), (1, 'N', 'NO2'), (2, 'O', 'NO2'),
                (0, 'O', 'N2'), (0, 'Ar', 'N2O')]
        for (n, elem, species) in data:
            self.assertEqual(self.gas.nAtoms(species, elem), n)
            mElem = self.gas.elementIndex(elem)
            kSpec = self.gas.speciesIndex(species)
            self.assertEqual(self.gas.nAtoms(kSpec, mElem), n)

    def test_weights(self):
        atomic_weights = self.gas.atomicWeights()
        molecular_weights = self.gas.molecularWeights()
        self.assertEqual(self.gas.nElements(), len(atomic_weights))
        self.assertEqual(self.gas.nSpecies(), len(molecular_weights))

        for i,mw in enumerate(molecular_weights):
            test_weight = 0.0
            for j,aw in enumerate(atomic_weights):
                test_weight += aw * self.gas.nAtoms(i,j)
            self.assertAlmostEqual(test_weight, mw)


class ThermoTest(unittest.TestCase):
    """
    Test the thermodynamic property accessor functions of a Solution
    """
    @classmethod
    def setUpClass(cls):
        cls.gas = ct.importPhase('../data/air-no-reactions.xml', 'air')
        cls.T0 = cls.gas.temperature()
        cls.P0 = cls.gas.pressure()
        cls.X0 = cls.gas.moleFractions()

    def setUp(self):
        # Workaround for Python 2.6 unittest, which does not call setUpClass
        try:
            self.gas
        except AttributeError:
            self.setUpClass()

        self.gas.set(T=self.T0, P=self.P0, X=self.X0)

    def test_volume(self):
        # This class should follow the ideal gas law
        g = self.gas
        self.assertAlmostEqual(
            g.pressure(),
            g.molarDensity() * ct.GasConstant * g.temperature())

        self.assertAlmostEqual(
            g.pressure() / g.density(),
            ct.GasConstant / g.meanMolecularWeight() * g.temperature())

        self.assertAlmostEqual(g.density(), 1.0 / g.volume_mass())

    def test_energy(self):
        g = self.gas
        mmw = g.meanMolecularWeight()
        self.assertAlmostEqual(g.enthalpy_mass(), g.enthalpy_mole() / mmw)
        self.assertAlmostEqual(g.intEnergy_mass(), g.intEnergy_mole() / mmw)
        self.assertAlmostEqual(g.gibbs_mass(), g.gibbs_mole() / mmw)
        self.assertAlmostEqual(g.entropy_mass(), g.entropy_mole() / mmw)

        self.assertAlmostEqual(g.cv_mass(), g.cv_mole() / mmw)
        self.assertAlmostEqual(g.cp_mass(), g.cp_mole() / mmw)
        self.assertAlmostEqual(g.cv_mole() + ct.GasConstant, g.cp_mole())

    def test_nondimensional(self):
        g = self.gas
        T = g.temperature()
        R = ct.GasConstant
        X = g.moleFractions()

        self.assertAlmostEqual(np.dot(g.cp_R(), X),
                               g.cp_mole() / R)
        self.assertAlmostEqual(np.dot(g.enthalpies_RT(), X),
                               g.enthalpy_mole() / (R*T))

        Smix_R = - np.dot(X, np.log(X+1e-20))
        self.assertAlmostEqual(np.dot(g.entropies_R(), X) + Smix_R,
                               g.entropy_mole() / R)
        self.assertAlmostEqual(np.dot(g.gibbs_RT(), X) - Smix_R,
                               g.gibbs_mole() / (R*T))
