from .utilities import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestThermoPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')

    def test_phases(self):
        self.assertEqual(self.phase.n_phases, 1)

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
        data = [(1, 'O', 'O'), (2, 'O', 'O2'), (1, b'H', b'OH'),
                (2, 'H', 'H2O'), (2, u'O', u'H2O2'), (1, 'Ar', 'AR'),
                (0, 'O', 'H'), (0, 'H', 'AR'), (0, 'Ar', 'HO2')]
        for (n, elem, species) in data:
            self.assertEqual(self.phase.n_atoms(species, elem), n)
            mElem = self.phase.element_index(elem)
            kSpec = self.phase.species_index(species)
            self.assertEqual(self.phase.n_atoms(kSpec, mElem), n)

        with self.assertRaises(ValueError):
            self.phase.n_atoms('C', 'H2')
        with self.assertRaises(ValueError):
            self.phase.n_atoms('H', 'CH4')

    def test_elemental_mass_fraction(self):
        self.phase.Y = 'H2O:0.5, O2:0.5'
        Zo = self.phase.elemental_mass_fraction('O')
        Zh = self.phase.elemental_mass_fraction('H')
        Zar = self.phase.elemental_mass_fraction('Ar')

        mO = self.phase.element_index('O')
        self.assertEqual(Zo, self.phase.elemental_mass_fraction(mO))
        self.assertNear(Zo, 0.5 + 0.5 * (15.9994 / 18.01528))
        self.assertNear(Zh, 0.5 * (2.01588 / 18.01528))
        self.assertEqual(Zar, 0.0)

        with self.assertRaises(ValueError):
            self.phase.elemental_mass_fraction('C')
        with self.assertRaises(ValueError):
            self.phase.elemental_mass_fraction(5)

    def test_elemental_mole_fraction(self):
        self.phase.X = 'H2O:0.5, O2:0.5'
        Zo = self.phase.elemental_mole_fraction('O')
        Zh = self.phase.elemental_mole_fraction('H')
        Zar = self.phase.elemental_mole_fraction('Ar')

        mO = self.phase.element_index('O')
        self.assertEqual(Zo, self.phase.elemental_mole_fraction(mO))
        self.assertNear(Zo, 0.5/3 + 0.5)
        self.assertNear(Zh, 0.5*2/3)
        self.assertEqual(Zar, 0.0)

        with self.assertRaises(ValueError):
            self.phase.elemental_mole_fraction('C')
        with self.assertRaises(ValueError):
            self.phase.elemental_mole_fraction(5)

    def test_weights(self):
        atomic_weights = self.phase.atomic_weights
        molecular_weights = self.phase.molecular_weights
        self.assertEqual(self.phase.n_elements, len(atomic_weights))
        self.assertEqual(self.phase.n_species, len(molecular_weights))

        for i,mw in enumerate(molecular_weights):
            test_weight = 0.0
            for j,aw in enumerate(atomic_weights):
                test_weight += aw * self.phase.n_atoms(i,j)
            self.assertAlmostEqual(test_weight, mw)

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

        with self.assertRaises(Exception):
            self.phase.X = 'H2:1.0, CO2:1.5'

    def test_setCompositionStringBad(self):
        X0 = self.phase.X
        with self.assertRaises(Exception):
            self.phase.X = 'H2:1.0, O2:asdf'
        self.assertArrayNear(X0, self.phase.X)

        with self.assertRaises(Exception):
            self.phase.X = 'H2:1e-x4'
        self.assertArrayNear(X0, self.phase.X)

        with self.assertRaises(Exception):
            self.phase.X = 'H2:1e-1.4'
        self.assertArrayNear(X0, self.phase.X)

    def test_setCompositionDict(self):
        self.phase.X = {b'H2':1.0, b'O2':3.0}
        X = self.phase.X
        self.assertNear(X[0], 0.25)
        self.assertNear(X[3], 0.75)

        self.phase.Y = {u'H2':1.0, u'O2':3.0}
        Y = self.phase.Y
        self.assertNear(Y[0], 0.25)
        self.assertNear(Y[3], 0.75)

    def test_getCompositionDict(self):
        self.phase.X = 'OH:1e-9, O2:0.4, AR:0.6'
        self.assertEqual(len(self.phase.mole_fraction_dict(1e-7)), 2)
        self.assertEqual(len(self.phase.mole_fraction_dict()), 3)

        self.phase.Y = 'O2:0.4, AR:0.6'
        Y1 = self.phase.mass_fraction_dict()
        self.assertNear(Y1['O2'], 0.4)
        self.assertNear(Y1['AR'], 0.6)

    def test_setCompositionNoNorm(self):
        X = np.zeros(self.phase.n_species)
        X[2] = 1.0
        X[0] = 0.01
        self.phase.set_unnormalized_mole_fractions(X)
        self.assertArrayNear(self.phase.X, X)
        self.assertNear(sum(X), 1.01)

        Y = np.zeros(self.phase.n_species)
        Y[2] = 1.0
        Y[0] = 0.01
        self.phase.set_unnormalized_mass_fractions(Y)
        self.assertArrayNear(self.phase.Y, Y)
        self.assertNear(sum(Y), 1.01)

    def test_setCompositionNoNormBad(self):
        X = np.zeros(self.phase.n_species - 1)
        with self.assertRaises(ValueError):
            self.phase.set_unnormalized_mole_fractions(X)

        with self.assertRaises(ValueError):
            self.phase.set_unnormalized_mass_fractions([1,2,3])

    @unittest.expectedFailure
    def test_setCompositionDict_bad1(self):
        # Non-existent species should raise an exception
        with self.assertRaises(Exception):
            self.phase.X = {'H2':1.0, 'HCl':3.0}

    def test_setCompositionDict_bad2(self):
        with self.assertRaises(Exception):
            self.phase.Y = {'H2':1.0, 'O2':'xx'}

    def test_setCompositionSlice(self):
        self.phase['H2', 'O2'].X = 0.1, 0.9
        X = self.phase.X
        self.assertNear(X[0], 0.1)
        self.assertNear(X[3], 0.9)

    def test_setCompositionSlice_bad(self):
        with self.assertRaises(ValueError):
            self.phase['H2','O2'].Y = [0.1, 0.2, 0.3]

    def test_full_report(self):
        report = self.phase.report(threshold=0.0)
        self.assertIn(self.phase.name, report)
        self.assertIn('temperature', report)
        self.assertNotIn('minor', report)
        for name in self.phase.species_names:
            self.assertIn(name, report)

    def test_default_report(self):
        self.phase.X = 'H2:0.1, O2:0.9, HO2:1e-10, H2O2:1e-20'
        report = self.phase.report()
        self.assertIn('minor', report)
        for name in (' H2 ', ' O2 ', ' HO2 '):
            self.assertIn(name, report)
        for name in (' H2O2 ', ' OH ', ' AR '):
            self.assertNotIn(name, report)

    def test_name(self):
        self.assertEqual(self.phase.name, 'ohmech')

        self.phase.name = 'something'
        self.assertEqual(self.phase.name, 'something')
        self.assertIn('something', self.phase.report())

    def test_ID(self):
        self.assertEqual(self.phase.ID, 'ohmech')

        self.phase.ID = 'something'
        self.assertEqual(self.phase.ID, 'something')
        self.assertEqual(self.phase.name, 'ohmech')

    def test_badLength(self):
        X = np.zeros(5)
        with self.assertRaises(ValueError):
            self.phase.X = X
        with self.assertRaises(ValueError):
            self.phase.Y = X

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

    def test_invalid_property(self):
        x = self.phase
        with self.assertRaises(AttributeError):
            x.foobar = 300
        with self.assertRaises(AttributeError):
            x.foobar

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

    def test_unpicklable(self):
        import pickle
        with self.assertRaises(NotImplementedError):
            pickle.dumps(self.phase)

    def test_uncopyable(self):
        import copy
        with self.assertRaises(NotImplementedError):
            copy.copy(self.phase)


class TestThermo(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.ThermoPhase('h2o2.xml')
        self.gas.TPX = 450, 2e5, 'H2:1.0, O2:0.4, AR:3, H2O:0.1'

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

    def test_volume(self):
        """ This phase should follow the ideal gas law """
        g = self.gas
        self.assertAlmostEqual(g.P, g.density_mole * ct.gas_constant * g.T)

        self.assertAlmostEqual(
            g.P / g.density,
            ct.gas_constant / g.mean_molecular_weight * g.T)

        self.assertAlmostEqual(g.density, 1.0 / g.volume_mass)

    def test_energy(self):
        g = self.gas
        mmw = g.mean_molecular_weight
        self.assertAlmostEqual(g.enthalpy_mass, g.enthalpy_mole / mmw)
        self.assertAlmostEqual(g.int_energy_mass, g.int_energy_mole / mmw)
        self.assertAlmostEqual(g.gibbs_mass, g.gibbs_mole / mmw)
        self.assertAlmostEqual(g.entropy_mass, g.entropy_mole / mmw)

        self.assertAlmostEqual(g.cv_mass, g.cv_mole / mmw)
        self.assertAlmostEqual(g.cp_mass, g.cp_mole / mmw)
        self.assertAlmostEqual(g.cv_mole + ct.gas_constant, g.cp_mole)

    def test_nondimensional(self):
        g = self.gas
        R = ct.gas_constant

        self.assertAlmostEqual(np.dot(g.standard_cp_R, g.X),
                               g.cp_mole / R)
        self.assertAlmostEqual(np.dot(g.standard_enthalpies_RT, g.X),
                               g.enthalpy_mole / (R*g.T))

        Smix_R = - np.dot(g.X, np.log(g.X+1e-20))
        self.assertAlmostEqual(np.dot(g.standard_entropies_R, g.X) + Smix_R,
                               g.entropy_mole / R)
        self.assertAlmostEqual(np.dot(g.standard_gibbs_RT, g.X) - Smix_R,
                               g.gibbs_mole / (R*g.T))


class TestInterfacePhase(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('diamond.xml', 'gas')
        self.solid = ct.Solution('diamond.xml', 'diamond')
        self.interface = ct.Interface('diamond.xml', 'diamond_100',
                                      (self.gas, self.solid))

    def test_properties(self):
        self.interface.site_density = 100
        self.assertNear(self.interface.site_density, 100)

    def test_coverages_array(self):
        C = np.zeros(self.interface.n_species)
        C[1] = 0.25
        C[3] = 0.125
        C[4] = 0.125
        self.interface.coverages = C
        C = self.interface.coverages
        # should now be normalized
        self.assertNear(C[1], 0.5)
        self.assertNear(C[3], 0.25)
        self.assertNear(C[4], 0.25)
        self.assertNear(sum(C), 1.0)

    def test_coverages_string(self):
        self.interface.coverages = 'c6HM:0.2, c6H*:0.8'
        C = self.interface.coverages
        self.assertNear(C[self.interface.species_index('c6HM')], 0.2)
        self.assertNear(C[self.interface.species_index('c6H*')], 0.8)

    def test_coverages_dict(self):
        self.interface.coverages = {'c6**':1.0, 'c6*M':3.0}
        C = self.interface.coverages
        self.assertNear(C[self.interface.species_index('c6**')], 0.25)
        self.assertNear(C[self.interface.species_index('c6*M')], 0.75)


class ImportTest(utilities.CanteraTest):
    """
    Test the various ways of creating a Solution object
    """
    def check(self, gas, name, T, P, nSpec, nElem):
        self.assertEqual(gas.name, name)
        self.assertAlmostEqual(gas.T, T)
        self.assertAlmostEqual(gas.P, P)
        self.assertEqual(gas.n_species, nSpec)
        self.assertEqual(gas.n_elements, nElem)

    def test_import_phase_cti(self):
        gas1 = ct.Solution('../data/air-no-reactions.cti', 'air')
        self.check(gas1, 'air', 300, 101325, 8, 3)

        gas2 = ct.Solution('../data/air-no-reactions.cti', 'notair')
        self.check(gas2, 'notair', 900, 5*101325, 7, 2)

    def test_import_phase_cti2(self):
        # This should import the first phase, i.e. 'air'
        gas = ct.Solution('../data/air-no-reactions.cti')
        self.check(gas, 'air', 300, 101325, 8, 3)

    def test_import_phase_xml(self):
        gas1 = ct.Solution('../data/air-no-reactions.xml', 'air')
        self.check(gas1, 'air', 300, 101325, 8, 3)

        gas2 = ct.Solution('../data/air-no-reactions.xml', 'notair')
        self.check(gas2, 'notair', 900, 5*101325, 7, 2)

    def test_import_phase_cti_text(self):
        cti_def = """
ideal_gas(name='spam', elements='O H',
          species='gri30: all',
          options='skip_undeclared_elements',
          initial_state=state(temperature=350, pressure=2e6))
"""
        gas = ct.Solution(source=cti_def)
        self.check(gas, 'spam', 350, 2e6, 8, 2)

    def test_import_phase_xml_text(self):
        xml_def = """
<?xml version="1.0"?>
<ctml>
  <validate reactions="yes" species="yes"/>
  <phase dim="3" id="spam">
    <elementArray datasrc="elements.xml">O</elementArray>
    <speciesArray datasrc="gri30.xml#species_data">all
      <skip element="undeclared"/>
    </speciesArray>
    <state>
      <temperature units="K">350.0</temperature>
      <pressure units="Pa">2000000.0</pressure>
    </state>
    <thermo model="IdealGas"/>
    <kinetics model="GasKinetics"/>
    <transport model="None"/>
  </phase>
</ctml>"""
        gas = ct.Solution(source=xml_def)
        self.check(gas, 'spam', 350, 2e6, 2, 1)

    def test_import_from_species(self):
        gas1 = ct.Solution('h2o2.xml')
        gas1.TPX = 350, 101325, 'H2:0.3, O2:0.7'
        gas1.equilibrate('HP')

        species = ct.Species.listFromFile('h2o2.xml')
        gas2 = ct.ThermoPhase(thermo='IdealGas', species=species)
        gas2.TPX = 350, 101325, 'H2:0.3, O2:0.7'
        gas2.equilibrate('HP')
        self.assertEqual(gas1.n_elements, gas2.n_elements)
        self.assertEqual(gas1.species_names, gas2.species_names)
        self.assertNear(gas1.T, gas2.T)
        self.assertArrayNear(gas1.X, gas2.X)

    def test_checkReactionBalance(self):
        with self.assertRaises(Exception):
            ct.Solution('../data/h2o2_unbalancedReaction.xml')


class TestSpecies(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('h2o2.xml')

    def test_standalone(self):
        s = ct.Species('CH4', {'C':1, 'H':4})

        self.assertEqual(s.name, 'CH4')
        c = s.composition
        self.assertEqual(len(c), 2)
        self.assertEqual(c['C'], 1)
        self.assertEqual(c['H'], 4)

    def test_defaults(self):
        s = ct.Species('H2')
        self.assertEqual(s.size, 1.0)
        self.assertEqual(s.charge, 0.0)

        self.assertIsNone(s.thermo)
        self.assertIsNone(s.transport)

    def test_index_accessor(self):
        for k in range(self.gas.n_species):
            s = self.gas.species(k)
            self.assertEqual(s.name, self.gas.species_name(k))

            for m,n in s.composition.items():
                self.assertEqual(n, self.gas.n_atoms(k,m))

    def test_species_noargs(self):
        for k,s in enumerate(self.gas.species()):
            self.assertEqual(s.name, self.gas.species_name(k))

    def test_name_accessor(self):
        for name in self.gas.species_names:
            s = self.gas.species(name)
            self.assertEqual(s.name, name)

    def test_fromCti(self):
        h2_cti = """
            species(
                name="H2",
                atoms="H:2",
                thermo=(
                    NASA([200.00, 1000.00],
                         [2.344331120E+00, 7.980520750E-03, -1.947815100E-05,
                          2.015720940E-08, -7.376117610E-12, -9.179351730E+02,
                          6.830102380E-01]),
                    NASA([1000.00, 3500.00],
                         [3.337279200E+00, -4.940247310E-05, 4.994567780E-07,
                         -1.795663940E-10, 2.002553760E-14, -9.501589220E+02,
                         -3.205023310E+00])
                ),
                transport=gas_transport(geom="linear",
                                        diam=2.92,
                                        well_depth=38.00,
                                        polar=0.79,
                                        rot_relax=280.00),
                note = "TPIS78"
            )"""
        s1 = self.gas.species('H2')
        s2 = ct.Species.fromCti(h2_cti)
        self.assertEqual(s2.name, 'H2')
        self.assertEqual(s1.composition, s2.composition)
        self.assertEqual(s1.thermo.cp(350), s2.thermo.cp(350))

    def test_fromXml(self):
        import xml.etree.ElementTree as ET
        root = ET.parse('../../build/data/h2o2.xml').getroot()
        h2_node = root.find('.//species[@name="H2"]')
        h2_string = ET.tostring(h2_node)

        s1 = self.gas.species('H2')
        s2 = ct.Species.fromXml(h2_string)

        self.assertEqual(s2.name, 'H2')
        self.assertEqual(s1.composition, s2.composition)
        self.assertEqual(s1.thermo.cp(350), s2.thermo.cp(350))

    def test_listFromFile_cti(self):
        S = ct.Species.listFromFile('h2o2.cti')
        self.assertEqual({sp.name for sp in S},
                         set(self.gas.species_names))

    def test_listFromFile_xml(self):
        S = ct.Species.listFromFile('h2o2.xml')
        self.assertEqual({sp.name for sp in S},
                         set(self.gas.species_names))

    def test_listFromCti(self):
        S = ct.Species.listFromCti(open('../../build/data/h2o2.cti').read())

        self.assertEqual({sp.name for sp in S},
                         set(self.gas.species_names))

    def test_listFromXml(self):
        S = ct.Species.listFromXml(open('../../build/data/h2o2.xml').read())

        self.assertEqual({sp.name for sp in S},
                         set(self.gas.species_names))



class TestSpeciesThermo(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('h2o2.xml')
        self.gas.X = 'H2O:1.0'

    def test_create(self):
        st = ct.NasaPoly2(300, 3500, 101325,
                [1000.0, 3.03399249E+00, 2.17691804E-03, -1.64072518E-07,
                 -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00,
                 4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09,
                 1.77197817E-12, -3.02937267E+04, -8.49032208E-01])

        for T in [300, 500, 900, 1200, 2000]:
            self.gas.TP = T, 101325
            self.assertAlmostEqual(st.cp(T), self.gas.cp_mole)
            self.assertAlmostEqual(st.h(T), self.gas.enthalpy_mole)
            self.assertAlmostEqual(st.s(T), self.gas.entropy_mole)

    def test_invalid(self):
        with self.assertRaises(ValueError):
            # not enough coefficients
            st = ct.NasaPoly2(300, 3500, 101325,
                              [1000.0, 3.03399249E+00, 2.17691804E-03])

    def test_wrap(self):
        st = self.gas.species('H2O').thermo

        self.assertTrue(isinstance(st, ct.NasaPoly2))

        for T in [300, 500, 900, 1200, 2000]:
            self.gas.TP = T, 101325
            self.assertAlmostEqual(st.cp(T), self.gas.cp_mole)
            self.assertAlmostEqual(st.h(T), self.gas.enthalpy_mole)
            self.assertAlmostEqual(st.s(T), self.gas.entropy_mole)
