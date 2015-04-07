from __future__ import division

import unittest
import os
import warnings

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


class TestKOH_Equil(utilities.CanteraTest):
    "Test roughly based on examples/multiphase/plasma_equilibrium.py"
    @classmethod
    def setUpClass(cls):
        cls.phases = ct.import_phases('KOH.xml',
                ['K_solid', 'K_liquid', 'KOH_a', 'KOH_b', 'KOH_liquid',
                 'K2O2_solid', 'K2O_solid', 'KO2_solid', 'ice', 'liquid_water',
                 'KOH_plasma'])

    def setUp(self):
        self.mix = ct.Mixture(self.phases)

    def compare(self, data, reference_file):
        if os.path.exists(reference_file):
            # Compare with existing output file
            ref = np.genfromtxt(reference_file)
            self.assertEqual(data.shape, ref.shape)
            for i in range(ref.shape[0]):
                self.assertArrayNear(ref[i], data[i])
        else:
            # Generate the output file for the first time
            warnings.warn('Generating test data file:' +
                          os.path.abspath(reference_file))
            np.savetxt(reference_file, data)

    def test_equil_TP(self):
        temperatures = range(350, 5000, 300)
        data = np.zeros((len(temperatures), self.mix.n_species+1))
        data[:,0] = temperatures

        for i,T in enumerate(temperatures):
            self.mix.T = T
            self.mix.P = ct.one_atm
            self.mix.species_moles = 'K:1.03, H2:2.12, O2:0.9'
            self.mix.equilibrate('TP')

            data[i,1:] = self.mix.species_moles

        self.compare(data, '../data/koh-equil-TP.csv')

    def test_equil_HP(self):
        temperatures = range(350, 5000, 300)
        data = np.zeros((len(temperatures), self.mix.n_species+2))
        data[:,0] = temperatures

        # The outer iteration for the temperature *should* be able to
        # converge from further away, but in practice, it can't. (Of course,
        # changing this value requires replacing the reference output)
        dT = 1
        self.mix.P = ct.one_atm

        for i,T in enumerate(temperatures):
            self.mix.species_moles = 'K:1.03, H2:2.12, O2:0.9'
            self.mix.T = T - dT
            self.mix.equilibrate('TP')
            self.mix.T = T
            self.mix.equilibrate('HP')

            data[i,1] = self.mix.T # equilibrated temperature
            data[i,2:] = self.mix.species_moles

        self.compare(data, '../data/koh-equil-HP.csv')
