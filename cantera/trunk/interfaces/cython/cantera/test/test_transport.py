import unittest

import cantera as ct
from . import utilities

class TestTransport(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_viscosity(self):
        self.assertTrue(self.phase.viscosity > 0.0)

class TestDustyGas(utilities.CanteraTest):
    def test_methods(self):
        self.phase = ct.DustyGas('h2o2.xml')

        self.phase.setPorosity(0.2)
