import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestReactor(utilities.CanteraTest):
    def setUp(self):
        self.gas1 = ct.Solution('h2o2.xml')
        self.gas2 = ct.Solution('h2o2.xml')
        X = np.zeros(self.gas1.nSpecies)
        X[3] = 1.0
        self.gas2.setMoleFractions(X)

        self.r1 = ct.Reactor(self.gas1)
        self.r2 = ct.Reactor(self.gas2)

    def test_disjoint(self):
        T1 = self.gas1.temperature
        T2 = self.gas2.temperature
        P1 = self.gas1.pressure
        P2 = self.gas2.pressure

        net = ct.ReactorNet()
        net.addReactor(self.r1)
        net.addReactor(self.r2)
        net.advance(1.0)

        # Nothing should change from the initial condition
        self.assertNear(T1, self.gas1.temperature)
        self.assertNear(T2, self.gas2.temperature)
        self.assertNear(P1, self.gas1.pressure)
        self.assertNear(P2, self.gas2.pressure)

    def test_equalizePressure(self):
        w = ct.Wall()
        w.install(self.r1, self.r2)
        w.expansionRateCoeff = 0.1
        w.area = 1.0

        net = ct.ReactorNet()
        net.addReactor(self.r1)
        net.addReactor(self.r2)

        net.advance(1.0)
        self.assertNear(self.gas1.pressure, self.gas2.pressure)
