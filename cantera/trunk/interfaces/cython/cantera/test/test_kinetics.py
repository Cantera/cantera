import unittest

import cantera as ct
from . import utilities

class TestKinetics(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_nReactions(self):
        self.assertEqual(self.phase.nReactions, 27)
