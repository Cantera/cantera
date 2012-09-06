import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestThermoPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_species(self):
        self.assertEqual(self.phase.nSpecies, 9)
        for i,name in enumerate(self.phase.speciesNames):
            self.assertEqual(name, self.phase.speciesName(i))
            self.assertEqual(i, self.phase.speciesIndex(name))
            self.assertEqual(i, self.phase.speciesIndex(i))

    def test_elements(self):
        self.assertEqual(self.phase.nElements, 3)
        for i,symbol in enumerate(self.phase.elementNames):
            self.assertEqual(symbol, self.phase.elementName(i))
            self.assertEqual(i, self.phase.elementIndex(symbol))
            self.assertEqual(i, self.phase.elementIndex(i))

    def test_nAtoms(self):
        self.assertEqual(self.phase.nAtoms('H2', 'H'), 2)
        self.assertEqual(self.phase.nAtoms('H2O2', 'H'), 2)
        self.assertEqual(self.phase.nAtoms('HO2', 'H'), 1)
        self.assertEqual(self.phase.nAtoms('AR', 'H'), 0)

        self.assertRaises(ValueError, lambda: self.phase.nAtoms('C', 'H2'))
        self.assertRaises(ValueError, lambda: self.phase.nAtoms('H', 'CH4'))

    def test_setComposition(self):
        X = np.zeros(self.phase.nSpecies)
        X[2] = 1.0
        self.phase.X = X
        Y = self.phase.Y

        self.assertEqual(list(X), list(Y))

    def test_setCompositionString(self):
        self.phase.X = 'H2:1.0, O2:1.0'
        X = self.phase.X
        self.assertEqual(X[0], 0.5)
        self.assertEqual(X[3], 0.5)

        def set_bad():
            self.phase.X = 'H2:1.0, CO2:1.5'

        self.assertRaises(Exception, set_bad)

    def test_report(self):
        report = self.phase.report()
        self.assertTrue(self.phase.name in report)
        self.assertTrue('temperature' in report)
        for name in self.phase.speciesNames:
            self.assertTrue(name in report)

    def test_name(self):
        self.assertEqual(self.phase.name, 'ohmech')

        self.phase.name = 'something'
        self.assertEqual(self.phase.name, 'something')
        self.assertTrue('something' in self.phase.report())

    def test_badLength(self):
        X = np.zeros(5)
        def set_X():
            self.phase.X = X
        def set_Y():
            self.phase.Y = X

        self.assertRaises(ValueError, set_X)
        self.assertRaises(ValueError, set_Y)

    def test_mass_basis(self):
        self.assertEqual(self.phase.basis, 'mass')
        self.assertEqual(self.phase.density_mass, self.phase.density)
        self.assertEqual(self.phase.enthalpy_mass, self.phase.h)
        self.assertEqual(self.phase.entropy_mass, self.phase.s)
        self.assertEqual(self.phase.intEnergy_mass, self.phase.u)
        self.assertEqual(self.phase.volume_mass, self.phase.v)
        self.assertEqual(self.phase.cv_mass, self.phase.cv)
        self.assertEqual(self.phase.cp_mass, self.phase.cp)

    def test_molar_basis(self):
        self.phase.basis = 'molar'
        self.assertEqual(self.phase.basis, 'molar')
        self.assertEqual(self.phase.density_mole, self.phase.density)
        self.assertEqual(self.phase.enthalpy_mole, self.phase.h)
        self.assertEqual(self.phase.entropy_mole, self.phase.s)
        self.assertEqual(self.phase.intEnergy_mole, self.phase.u)
        self.assertEqual(self.phase.volume_mole, self.phase.v)
        self.assertEqual(self.phase.cv_mole, self.phase.cv)
        self.assertEqual(self.phase.cp_mole, self.phase.cp)


    def check_setters(self, T1, rho1, Y1):
        T0, rho0, Y0 = self.phase.TDY
        self.phase.TDY = T1, rho1, Y1
        X1 = self.phase.X
        P1 = self.phase.P
        h1 = self.phase.h
        s1 = self.phase.s
        u1 = self.phase.u
        v1 = self.phase.v

        def check_state(T, rho, Y):
            self.assertNear(self.phase.T, T)
            self.assertNear(self.phase.density, rho)
            self.assertArrayNear(self.phase.Y, Y)

        self.phase.TDY = T0, rho0, Y0
        self.phase.TPY = T1, P1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.UVY = u1, v1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.HPY = h1, P1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SPY = s1, P1, Y1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.TPX = T1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.UVX = u1, v1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.HPX = h1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SPX = s1, P1, X1
        check_state(T1, rho1, Y1)


    def test_setState_mass(self):
        self.check_setters(T1 = 500.0, rho1 = 1.5,
                           Y1 = [0.1, 0.0, 0.0, 0.1, 0.4, 0.2, 0.0, 0.0, 0.2])

    def test_setState_mole(self):
        self.phase.basis = 'molar'
        self.check_setters(T1 = 750.0, rho1 = 0.02,
                           Y1 = [0.2, 0.1, 0.0, 0.3, 0.1, 0.0, 0.0, 0.2, 0.1])

    def check_getters(self):
        T,P,X = self.phase.TPX
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        S,P,Y = self.phase.SPY
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        U,V = self.phase.UV
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)

    def test_getState_mass(self):
        self.phase.TDY = 350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.check_getters()

    def test_getState_mole(self):
        self.phase.basis = 'molar'
        self.phase.TDX = 350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'
        self.check_getters()

    def test_getState(self):
        self.assertNear(self.phase.P, ct.OneAtm)
        self.assertNear(self.phase.T, 300)

    def test_partial_molar(self):
        self.phase.TDY = 350.0, 0.6, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.assertNear(sum(self.phase.partial_molar_enthalpies * self.phase.X),
                        self.phase.enthalpy_mole)

        self.assertNear(sum(self.phase.partial_molar_entropies * self.phase.X),
                        self.phase.entropy_mole)

        self.assertNear(sum(self.phase.partial_molar_int_energies * self.phase.X),
                        self.phase.intEnergy_mole)

        self.assertNear(sum(self.phase.chem_potentials * self.phase.X),
                        self.phase.gibbs_mole)

        self.assertNear(sum(self.phase.partial_molar_cp * self.phase.X),
                        self.phase.cp_mole)

    def test_nondimensional(self):
        self.phase.TDY = 850.0, 0.2, 'H2:0.1, H2O:0.6, AR:0.3'
        H = (sum(self.phase.standard_enthalpies_RT * self.phase.X) *
             ct.GasConstant * self.phase.T)
        self.assertNear(H, self.phase.enthalpy_mole)

        U = (sum(self.phase.standard_intEnergies_RT * self.phase.X) *
             ct.GasConstant * self.phase.T)
        self.assertNear(U, self.phase.intEnergy_mole)

        cp = sum(self.phase.standard_cp_R * self.phase.X) * ct.GasConstant
        self.assertNear(cp, self.phase.cp_mole)

    @unittest.skip("needs merge from 2.0 maintenance")
    def test_isothermalCompressibility(self):
        self.assertNear(self.phase.isothermalCompressibility, 1.0/self.phase.P)

    def test_thermalExpansionCoeff(self):
        self.assertNear(self.phase.thermalExpansionCoeff, 1.0/self.phase.T)

    def test_ref_info(self):
        self.assertEqual(self.phase.refPressure, ct.OneAtm)
        self.assertEqual(self.phase.minTemp, 300.0)
        self.assertEqual(self.phase.maxTemp, 3500.0)


class TestInterfacePhase(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('diamond.xml', 'gas')
        self.solid = ct.Solution('diamond.xml', 'diamond')
        self.interface = ct.Solution('diamond.xml', 'diamond_100',
                                      (self.gas, self.solid))

    def test_properties(self):
        self.interface.siteDensity = 100
        self.assertEqual(self.interface.siteDensity, 100)
