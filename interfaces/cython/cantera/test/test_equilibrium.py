from __future__ import division

import unittest

import numpy as np

import cantera as ct
from . import utilities

class EquilTestCases(object):
    def __init__(self, solver):
        self.solver = solver

    def check(self, gas, **moles):
        nTotal = sum(moles.values())
        for name,X in moles.items():
            self.assertAlmostEqual(gas[name].X[0], X/nTotal)

    def test_equil_complete_stoichiometric(self):
        """
        Equilibrium should correspond to complete combustion
        """
        gas = ct.Solution('equilibrium.cti', 'complete')
        gas.TPX = 298, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_complete_lean(self):
        """
        Equilibrium should correspond to complete combustion (with excess O2)
        CH4 + 3 O2 -> CO2 + 2 H2O + O2
        """
        gas = ct.Solution('equilibrium.cti', 'complete')
        gas.TPX = 298, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_incomplete_stoichiometric(self):
        gas = ct.Solution('equilibrium.cti', 'incomplete')
        gas.TPX = 301, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_incomplete_lean(self):
        gas = ct.Solution('equilibrium.cti', 'incomplete')
        gas.TPX = 301, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_gri_stoichiometric(self):
        gas = ct.Solution('gri30.xml')
        gas.TPX = 301, 100000, 'CH4:1.0, O2:2.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=0, H2O=2, CO2=1)

    def test_equil_gri_lean(self):
        gas = ct.Solution('gri30.xml')
        gas.TPX = 301, 100000, 'CH4:1.0, O2:3.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=0, O2=1, H2O=2, CO2=1)

    def test_equil_overconstrained1(self):
        gas = ct.Solution('equilibrium.cti', 'overconstrained-1')
        gas.TPX = 301, 100000, 'CH4:1.0, O2:1.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=1, O2=1)

    def test_equil_overconstrained2(self):
        gas = ct.Solution('equilibrium.cti', 'overconstrained-2')
        gas.TPX = 301, 100000, 'CH4:1.0, O2:1.0'
        gas.equilibrate('TP', self.solver)
        self.check(gas, CH4=1, O2=1)


class ChemEquilTest(EquilTestCases, utilities.CanteraTest):
    def __init__(self, *args, **kwargs):
        EquilTestCases.__init__(self, 'element_potential')
        unittest.TestCase.__init__(self, *args, **kwargs)


class MultiphaseEquilTest(EquilTestCases, utilities.CanteraTest):
    def __init__(self, *args, **kwargs):
        EquilTestCases.__init__(self, 'gibbs')
        unittest.TestCase.__init__(self, *args, **kwargs)


class VCS_EquilTest(EquilTestCases, utilities.CanteraTest):
    def __init__(self, *args, **kwargs):
        EquilTestCases.__init__(self, 'vcs')
        unittest.TestCase.__init__(self, *args, **kwargs)
