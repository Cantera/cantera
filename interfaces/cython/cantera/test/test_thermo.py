import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestThermoPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_nSpecies(self):
        self.assertEqual(self.phase.nSpecies, 9)

    def test_setComposition(self):
        X = np.zeros(self.phase.nSpecies)
        X[2] = 1.0
        self.phase.X = X
        Y = self.phase.Y

        self.assertEqual(list(X), list(Y))

    def test_badLength(self):
        X = np.zeros(5)
        with self.assertRaises(ValueError):
            self.phase.X = X


    def test_getState(self):
        self.assertNear(self.phase.P, ct.OneAtm)
        self.assertNear(self.phase.T, 300)


class TestInterfacePhase(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('diamond.xml', 'gas')
        self.solid = ct.Solution('diamond.xml', 'diamond')
        self.interface = ct.Solution('diamond.xml', 'diamond_100',
                                      (self.gas, self.solid))

    def test_properties(self):
        self.interface.siteDensity = 100
        self.assertEqual(self.interface.siteDensity, 100)
