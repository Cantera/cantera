import unittest

import cantera as ct
from . import utilities

class TestMixture(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        cls.phase1 = ct.Solution('h2o2.xml')
        cls.phase2 = ct.Solution('air.xml')

    def test_properties(self):
        mix = ct.Mixture([(self.phase1, 1.0), (self.phase2, 2.0)])
        self.assertEqual(mix.nSpecies, self.phase1.nSpecies + self.phase2.nSpecies)

        mix.temperature = 350
        self.assertEqual(mix.temperature, 350)

        mix.pressure = 2e5
        self.assertEqual(mix.pressure, 2e5)
