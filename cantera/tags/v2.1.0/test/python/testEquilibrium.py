from __future__ import division

import unittest

import numpy as np

import Cantera as ct

class EquilTestCases(object):
    def __init__(self, solver):
        self.solver = solver

    """
    Test the ChemEquil equilibrium solver
    """
    def check(self, gas, **moles):
        nTotal = sum(moles.values())
        for name,X in moles.iteritems():
            self.assertAlmostEqual(gas.moleFraction(name), X/nTotal)

    def test_equil_complete_stoichiometric(self):
        """
        Equilibrium should correspond to complete combustion
        """
        gas = ct.importPhase('equilibrium.cti', 'complete')
        gas.set(X='CH4:1.0, O2:2.0', T=298, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_complete_lean(self):
        """
        Equilibrium should correspond to complete combustion (with excess O2)
        CH4 + 3 O2 -> CO2 + 2 H2O + O2
        """
        gas = ct.importPhase('equilibrium.cti', 'complete')
        gas.set(X='CH4:1.0, O2:3.0', T=298, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_incomplete_stoichiometric(self):
        gas = ct.importPhase('equilibrium.cti', 'incomplete')
        gas.set(X='CH4:1.0, O2:2.0', T=301, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_incomplete_lean(self):
        gas = ct.importPhase('equilibrium.cti', 'incomplete')
        gas.set(X='CH4:1.0, O2:3.0', T=301, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_gri_stoichiometric(self):
        gas = ct.importPhase('gri30.xml')
        gas.set(X='CH4:1.0, O2:2.0', T=301, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_gri_lean(self):
        gas = ct.importPhase('gri30.xml')
        gas.set(X='CH4:1.0, O2:3.0', T=301, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_overconstrained1(self):
        gas = ct.importPhase('equilibrium.cti', 'overconstrained-1')
        gas.set(X='CH4:1.0, O2:1.0', T=301, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=1, O2=1)

    def test_equil_overconstrained2(self):
        gas = ct.importPhase('equilibrium.cti', 'overconstrained-2')
        gas.set(X='CH4:1.0, O2:1.0', T=301, P=100000)
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=1, O2=1)


class ChemEquilTest(EquilTestCases, unittest.TestCase):
    def __init__(self, *args, **kwargs):
        EquilTestCases.__init__(self, 0)
        unittest.TestCase.__init__(self, *args, **kwargs)


class MultiphaseEquilTest(EquilTestCases, unittest.TestCase):
    def __init__(self, *args, **kwargs):
        EquilTestCases.__init__(self, 1)
        unittest.TestCase.__init__(self, *args, **kwargs)


class VCS_EquilTest(EquilTestCases, unittest.TestCase):
    def __init__(self, *args, **kwargs):
        EquilTestCases.__init__(self, 2)
        unittest.TestCase.__init__(self, *args, **kwargs)
