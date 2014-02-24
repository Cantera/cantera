from .utilities import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestThermoPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_species(self):
        self.assertEqual(self.phase.n_species, 9)
        for i,name in enumerate(self.phase.species_names):
            self.assertEqual(name, self.phase.species_name(i))
            self.assertEqual(i, self.phase.species_index(name))
            self.assertEqual(i, self.phase.species_index(i))

    def test_elements(self):
        self.assertEqual(self.phase.n_elements, 3)
        for i,symbol in enumerate(self.phase.element_names):
            self.assertEqual(symbol, self.phase.element_name(i))
            self.assertEqual(i, self.phase.element_index(symbol))
            self.assertEqual(i, self.phase.element_index(i))

    def test_n_atoms(self):
        self.assertEqual(self.phase.n_atoms('H2', 'H'), 2)
        self.assertEqual(self.phase.n_atoms('H2O2', 'H'), 2)
        self.assertEqual(self.phase.n_atoms('HO2', 'H'), 1)
        self.assertEqual(self.phase.n_atoms('AR', 'H'), 0)

        self.assertRaises(ValueError, lambda: self.phase.n_atoms('C', 'H2'))
        self.assertRaises(ValueError, lambda: self.phase.n_atoms('H', 'CH4'))

    def test_setComposition(self):
        X = np.zeros(self.phase.n_species)
        X[2] = 1.0
        self.phase.X = X
        Y = self.phase.Y

        self.assertEqual(list(X), list(Y))

    def test_setCompositionString(self):
        self.phase.X = 'H2:1.0, O2:1.0'
        X = self.phase.X
        self.assertNear(X[0], 0.5)
        self.assertNear(X[3], 0.5)

        def set_bad():
            self.phase.X = 'H2:1.0, CO2:1.5'

        self.assertRaises(Exception, set_bad)

    def test_report(self):
        report = self.phase.report()
        self.assertTrue(self.phase.name in report)
        self.assertTrue('temperature' in report)
        for name in self.phase.species_names:
            self.assertTrue(name in report)

    def test_name(self):
        self.assertEqual(self.phase.name, 'ohmech')

        self.phase.name = 'something'
        self.assertEqual(self.phase.name, 'something')
        self.assertTrue('something' in self.phase.report())

    def test_ID(self):
        self.assertEqual(self.phase.ID, 'ohmech')

        self.phase.ID = 'something'
        self.assertEqual(self.phase.ID, 'something')
        self.assertEqual(self.phase.name, 'ohmech')

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
        self.assertEqual(self.phase.int_energy_mass, self.phase.u)
        self.assertEqual(self.phase.volume_mass, self.phase.v)
        self.assertEqual(self.phase.cv_mass, self.phase.cv)
        self.assertEqual(self.phase.cp_mass, self.phase.cp)

    def test_molar_basis(self):
        self.phase.basis = 'molar'
        self.assertEqual(self.phase.basis, 'molar')
        self.assertEqual(self.phase.density_mole, self.phase.density)
        self.assertEqual(self.phase.enthalpy_mole, self.phase.h)
        self.assertEqual(self.phase.entropy_mole, self.phase.s)
        self.assertEqual(self.phase.int_energy_mole, self.phase.u)
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

        self.phase.TDY = T0, rho0, Y0
        self.phase.SVX = s1, v1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.SVY = s1, v1, Y1
        check_state(T1, rho1, Y1)


    def test_setState_mass(self):
        self.check_setters(T1 = 500.0, rho1 = 1.5,
                           Y1 = [0.1, 0.0, 0.0, 0.1, 0.4, 0.2, 0.0, 0.0, 0.2])

    def test_setState_mole(self):
        self.phase.basis = 'molar'
        self.check_setters(T1 = 750.0, rho1 = 0.02,
                           Y1 = [0.2, 0.1, 0.0, 0.3, 0.1, 0.0, 0.0, 0.2, 0.1])

    def test_setters_hold_constant(self):
        props = ('T','P','s','h','u','v','X','Y')
        pairs = [('TP', 'T', 'P'), ('SP', 's', 'P'),
                 ('UV', 'u', 'v')]

        self.phase.TDX = 1000, 1.5, 'H2O:0.1, O2:0.95, AR:3.0'
        values = {}
        for p in props:
            values[p] = getattr(self.phase, p)

        for pair, first, second in pairs:
            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            first_val = getattr(self.phase, first)
            second_val = getattr(self.phase, second)

            setattr(self.phase, pair, (values[first], None))
            self.assertNear(getattr(self.phase, first), values[first])
            self.assertNear(getattr(self.phase, second), second_val)

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair, (None, values[second]))
            self.assertNear(getattr(self.phase, first), first_val)
            self.assertNear(getattr(self.phase, second), values[second])

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'X', (None, None, values['X']))
            self.assertNear(getattr(self.phase, first), first_val)
            self.assertNear(getattr(self.phase, second), second_val)

            self.phase.TDX = 500, 2.5, 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'Y', (None, None, values['Y']))
            self.assertNear(getattr(self.phase, first), first_val)
            self.assertNear(getattr(self.phase, second), second_val)

    def test_setter_errors(self):
        with self.assertRaises(Exception):
            self.phase.TD = 400

        with self.assertRaises(AssertionError):
            self.phase.TP = 300, 101325, 'CH4:1.0'

        with self.assertRaises(AssertionError):
            self.phase.HPY = 1.2e6, 101325

        with self.assertRaises(AssertionError):
            self.phase.UVX = -4e5, 4.4, 'H2:1.0', -1

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

        S,V,X = self.phase.SVX
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(X, self.phase.X)

    def test_getState_mass(self):
        self.phase.TDY = 350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.check_getters()

    def test_getState_mole(self):
        self.phase.basis = 'molar'
        self.phase.TDX = 350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'
        self.check_getters()

    def test_getState(self):
        self.assertNear(self.phase.P, ct.one_atm)
        self.assertNear(self.phase.T, 300)

    def test_partial_molar(self):
        self.phase.TDY = 350.0, 0.6, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.assertNear(sum(self.phase.partial_molar_enthalpies * self.phase.X),
                        self.phase.enthalpy_mole)

        self.assertNear(sum(self.phase.partial_molar_entropies * self.phase.X),
                        self.phase.entropy_mole)

        self.assertNear(sum(self.phase.partial_molar_int_energies * self.phase.X),
                        self.phase.int_energy_mole)

        self.assertNear(sum(self.phase.chemical_potentials * self.phase.X),
                        self.phase.gibbs_mole)

        self.assertNear(sum(self.phase.partial_molar_cp * self.phase.X),
                        self.phase.cp_mole)

    def test_nondimensional(self):
        self.phase.TDY = 850.0, 0.2, 'H2:0.1, H2O:0.6, AR:0.3'
        H = (sum(self.phase.standard_enthalpies_RT * self.phase.X) *
             ct.gas_constant * self.phase.T)
        self.assertNear(H, self.phase.enthalpy_mole)

        U = (sum(self.phase.standard_int_energies_RT * self.phase.X) *
             ct.gas_constant * self.phase.T)
        self.assertNear(U, self.phase.int_energy_mole)

        cp = sum(self.phase.standard_cp_R * self.phase.X) * ct.gas_constant
        self.assertNear(cp, self.phase.cp_mole)

    def test_isothermal_compressibility(self):
        self.assertNear(self.phase.isothermal_compressibility, 1.0/self.phase.P)

    def test_thermal_expansion_coeff(self):
        self.assertNear(self.phase.thermal_expansion_coeff, 1.0/self.phase.T)

    def test_ref_info(self):
        self.assertNear(self.phase.reference_pressure, ct.one_atm)
        self.assertNear(self.phase.min_temp, 300.0)
        self.assertNear(self.phase.max_temp, 3500.0)


class TestThermo(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.ThermoPhase('h2o2.xml')

    def test_setSV_lowT(self):
        """
        Set state in terms of (s,v) when the end temperature is below the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        s1, v1 = self.gas.SV
        self.gas.SV = s1, 3 * v1

        self.assertNear(self.gas.s, s1)
        self.assertNear(self.gas.v, 3 * v1)
        self.assertTrue(self.gas.T < self.gas.min_temp)

    def test_setSV_highT(self):
        """
        Set state in terms of (s,v) when the end temperature is above the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 2900, 1e5, 'H2:1.0, O2:0.4, AR:3'
        s1, v1 = self.gas.SV
        self.gas.SV = s1, 0.3 * v1

        self.assertNear(self.gas.s, s1)
        self.assertNear(self.gas.v, 0.3 * v1)
        self.assertTrue(self.gas.T > self.gas.max_temp)

    def test_setHP_lowT(self):
        """
        Set state in terms of (s,v) when the end temperature is below the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 450, 1e5, 'H2:1.0, O2:0.4, AR:3'
        deltaH = 1.25e5
        h1, p1 = self.gas.HP
        self.gas.HP = h1 - deltaH, None

        self.assertNear(self.gas.h, h1 - deltaH)
        self.assertNear(self.gas.P, p1)
        self.assertTrue(self.gas.T < self.gas.min_temp)

    def test_setHP_highT(self):
        """
        Set state in terms of (s,v) when the end temperature is above the
        phase's nominal temperature limit.
        """

        self.gas.TPX = 2800, 1e5, 'H2:1.0, O2:0.4, AR:3'
        deltaH = 8.25e5
        h1, p1 = self.gas.HP
        self.gas.HP = h1 + deltaH, None

        self.assertNear(self.gas.h, h1 + deltaH)
        self.assertNear(self.gas.P, p1)
        self.assertTrue(self.gas.T > self.gas.max_temp)


class TestInterfacePhase(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('diamond.xml', 'gas')
        self.solid = ct.Solution('diamond.xml', 'diamond')
        self.interface = ct.Solution('diamond.xml', 'diamond_100',
                                      (self.gas, self.solid))

    def test_properties(self):
        self.interface.site_density = 100
        self.assertNear(self.interface.site_density, 100)
