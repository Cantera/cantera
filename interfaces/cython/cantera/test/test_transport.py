import unittest

import cantera as ct
from . import utilities

class TestTransport(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_viscosity(self):
        self.assertTrue(self.phase.viscosity > 0.0)
