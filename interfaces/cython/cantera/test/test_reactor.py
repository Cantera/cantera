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
        self.gas2.X = X

        self.r1 = ct.Reactor(self.gas1)
        self.r2 = ct.Reactor(self.gas2)

    def test_disjoint(self):
        T1,P1 = self.gas1.TP
        T2,P2 = self.gas2.TP

        net = ct.ReactorNet()
        net.addReactor(self.r1)
        net.addReactor(self.r2)
        net.advance(1.0)

        # Nothing should change from the initial condition
        self.assertNear(T1, self.gas1.T)
        self.assertNear(T2, self.gas2.T)
        self.assertNear(P1, self.gas1.P)
        self.assertNear(P2, self.gas2.P)

    def test_equalizePressure(self):
        w = ct.Wall()
        w.install(self.r1, self.r2)
        w.expansionRateCoeff = 0.1
        w.area = 1.0

        net = ct.ReactorNet()
        net.addReactor(self.r1)
        net.addReactor(self.r2)

        net.advance(1.0)
        self.assertNear(self.gas1.P, self.gas2.P)
