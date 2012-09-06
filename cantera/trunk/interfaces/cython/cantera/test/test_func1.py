import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestFunc1(utilities.CanteraTest):
    def test_sin(self):
        f = ct.Sin1(3)
        self.assertNear(f(0), np.sin(0))
        self.assertNear(f(0.1), np.sin(3*0.1))
        self.assertNear(f(0.7), np.sin(3*0.7))
