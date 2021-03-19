import itertools
import numpy as np
import warnings
from pint import UnitRegistry
units = UnitRegistry()
Q_ = units.Quantity

import cantera.units as ct
from . import utilities

class TestSolutionUnits(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution("h2o2.yaml")

    def test_mass_basis(self):
        self.assertEqual(self.phase.basis_units, 'kg')
        self.assertQuantityNear(self.phase.density_mass, self.phase.density)
        self.assertQuantityNear(self.phase.enthalpy_mass, self.phase.h)
        self.assertQuantityNear(self.phase.entropy_mass, self.phase.s)
        self.assertQuantityNear(self.phase.int_energy_mass, self.phase.u)
        self.assertQuantityNear(self.phase.volume_mass, self.phase.v)
        self.assertQuantityNear(self.phase.gibbs_mass, self.phase.g)
        self.assertQuantityNear(self.phase.cp_mass, self.phase.cp)
        self.assertQuantityNear(self.phase.cv_mass, self.phase.cv)

    def test_molar_basis(self):
        self.assertEqual(self.phase.basis_units, 'kmol')
        self.assertQuantityNear(self.phase.density_mole, self.phase.density)
        self.assertQuantityNear(self.phase.enthalpy_mole, self.phase.h)
        self.assertQuantityNear(self.phase.entropy_mole, self.phase.s)
        self.assertQuantityNear(self.phase.int_energy_mole, self.phase.u)
        self.assertQuantityNear(self.phase.volume_mole, self.phase.v)
        self.assertQuantityNear(self.phase.gibbs_mole, self.phase.g)
        self.assertQuantityNear(self.phase.cp_mole, self.phase.cp)
        self.assertQuantityNear(self.phase.cv_mole, self.phase.cv)

    def test_dimensions(self):
        #basis-independent
        dims_T = self.phase.T.dimensionality
        self.assertIn("[temperature]", dims_T)
        self.assertEqual(dims_T["[temperature]"], 1.0)
        dims_P = self.phase.P.dimensionality
        self.assertIn("[pressure]", dims_P)
        self.assertEqual(dims_P["[pressure]"], 1.0)
        dims_X = self.phase.X.dimensionality
        self.assertIn("[dimensionless]", dims_X)
        self.assertEqual(dims_X["[dimensionless]"], 1.0)
        dims_Y = self.phase.Y.dimensionality
        self.assertIn("[dimensionless]", dims_Y)
        self.assertEqual(dims_Y["[dimensionless]"], 1.0)
        dims_Q = self.phase.Q.dimensionality
        self.assertIn("[dimensionless]", dims_Q)
        self.assertEqual(dims_Q["[dimensionless]"], 1.0)
        dims_T_sat = self.phase.T_sat.dimensionality
        self.assertIn("[temperature]", dims_T_sat)
        self.assertEqual(dims_T_sat["[temperature]"], 1.0)
        dims_P_sat = self.phase.P_sat.dimensionality
        self.assertIn("[pressure]", dims_P_sat)
        self.assertEqual(dims_P_sat["[pressure]"], 1.0)
        dims_atomic_weight = self.phase.atomic_weight.dimensionality
        self.assertIn("[mass/molar]", dims_atomic_weight)
        self.assertEqual(dims_atomic_weight["[mass]"], 1.0)
        self.assertEqual(dims_atomic_weight["[molar]"], -1.0)
        dims_chemical_potentials = self.phase.chemical_potentials.dimensionality
        self.assertIn("[energy/molar]", dims_chemical_potentials)
        self.assertEqual(dims_chemical_potentials["[energy]"], 1.0)
        self.assertEqual(dims_chemical_potentials["[molar]"], -1.0)
        dims_concentration = self.phase.concentration.dimensionality
        self.assertIn("[molar/length]", dims_concentration)
        self.assertEqual(dims_concentration["[molar]"], 1.0)
        self.assertEqual(dims_concentration["[length]"], -3.0)
        dims_critical_temperature = self.phase.critical_temperature.dimensionality
        self.assertIn("[temperature]", dims_critical_temperature)
        self.assertEqual(dims_critical_temperature["[temperature]"], 1.0)
        dims_critical_pressure = self.phase.critical_pressure.dimensionality
        self.assertIn("[pressure]", dims_critical_pressure)
        self.assertEqual(dims_critical_pressure["[pressure]"], 1.0)
        dims_electric_potential = self.phase.electric_potential.dimensionality
        self.assertIn("[volts]", dims_electric_potential)
        self.assertEqual(dims_electric_potential["[volts]"], 1.0)
        dims_electrochemical_potentials = self.phase.electrochemical_potentials.dimensionality
        self.assertIn("[energy/molar]", dims_electrochemical_potentials)
        self.assertEqual(dims_electrochemical_potentials["[energy]"], 1.0)
        self.assertEqual(dims_electrochemical_potentials["[molar]"], -1.0)
        dims_isothermal_compressibility = self.phase.isothermal_compressibility.dimensionality
        self.assertIn("[1/pressure]", dims_isothermal_compressibility)
        self.assertEqual(dims_isothermal_compressibility["[pressure]"], -1.0)
        dims_max_temp = self.phase.max_temp.dimensionality
        self.assertIn("[temperature]", dims_max_temp)
        self.assertEqual(dims_max_temp["[temperature]"], 1.0)
        dims_mean_molecular_weight = self.phase.mean_molecular_weight.dimensionality
        self.assertIn("[mass/molar]", dims_mean_molecular_weight)
        self.assertEqual(dims_mean_molecular_weight["[mass]"], 1.0)
        self.assertEqual(dims_mean_molecular_weight["[molar]"], -1.0)
        dims_min_temp = self.phase.min_temp.dimensionality
        self.assertIn("[temperature]", dims_min_temp)
        self.assertEqual(dims_min_temp["[temperature]"], 1.0)
        dims_molecular_weights = self.phase.molecular_weights.dimensionality
        self.assertIn("[mass/molar]", dims_molecular_weights)
        self.assertEqual(dims_molecular_weights["[mass]"], 1.0)
        self.assertEqual(dims_molecular_weights["[molar]"], -1.0)
        dims_partial_molar_cp = self.phase.partial_molar_cp.dimensionality
        self.assertIn("[energy/molar/temperature]", dims_partial_molar_cp)
        self.assertEqual(dims_partial_molar_cp["[energy]"], 1.0)
        self.assertEqual(dims_partial_molar_cp["[molar]"], -1.0)
        self.assertEqual(dims_partial_molar_cp["[temperature]"], -1.0)
        dims_partial_molar_enthalpies = self.phase.partial_molar_enthalpies.dimensionality
        self.assertIn("[energy/molar/temperature]", dims_partial_molar_enthalpies)
        self.assertEqual(dims_partial_molar_enthalpies["[energy]"], 1.0)
        self.assertEqual(dims_partial_molar_enthalpies["[molar]"], -1.0)
        self.assertEqual(dims_partial_molar_enthalpies["[temperature]"], -1.0)
        dims_partial_molar_entropies = self.phase.partial_molar_entropies.dimensionality
        self.assertIn("[energy/molar/temperature]", dims_partial_molar_entropies)
        self.assertEqual(dims_partial_molar_entropies["[energy]"], 1.0)
        self.assertEqual(dims_partial_molar_entropies["[molar]"], -1.0)
        self.assertEqual(dims_partial_molar_entropies["[temperature]"], -1.0)
        dims_partial_molar_int_energies = self.phase.partial_molar_int_energies.dimensionality
        self.assertIn("[energy/molar]", dims_partial_molar_int_energies)
        self.assertEqual(dims_partial_molar_int_energies["[energy]"], 1.0)
        self.assertEqual(dims_partial_molar_int_energies["[molar]"], -1.0)
        dims_partial_molar_volumes = self.phase.partial_molar_volumes.dimensionality
        self.assertIn("[length/molar]", dims_partial_molar_volumes)
        self.assertEqual(dims_partial_molar_volumes["[length]"], 3.0)
        self.assertEqual(dims_partial_molar_volumes["[molar]"], -1.0)
        dims_reference_pressure = self.phase.reference_pressure.dimensionality
        self.assertIn("[pressure]", dims_reference_pressure)
        self.assertEqual(dims_reference_pressure["[pressure]"], 1.0)
        dims_thermal_expansion_coeff = self.phase.thermal_expansion_coeff.dimensionality
        self.assertIn("[1/temperature]", dims_thermal_expansion_coeff)
        self.assertEqual(dims_thermal_expansion_coeff["[temperature]"], -1.0)
        #basis-dependent (mass)
        dims_density_mass = self.phase.density_mass.dimensionality
        self.assertIn("[mass/length]", dims_density_mass)
        self.assertEqual(dims_density_mass["[mass]"], 1.0)
        self.assertEqual(dims_density_mass["[length]"], -3.0)
        dims_enthalpy_mass = self.phase.enthalpy_mass.dimensionality
        self.assertIn("[energy/mass]", dims_enthalpy_mass)
        self.assertEqual(dims_enthalpy_mass["[energy]"], 1.0)
        self.assertEqual(dims_enthalpy_mass["[mass]"], -1.0)
        dims_entropy_mass = self.phase.entropy_mass.dimensionality
        self.assertIn("[energy/mass/temperature]", dims_entropy_mass)
        self.assertEqual(dims_entropy_mass["[energy]"], 1.0)
        self.assertEqual(dims_entropy_mass["[mass]"], -1.0)
        self.assertEqual(dims_entropy_mass["[temperature]"], -1.0)
        dims_int_energy_mass = self.phase.int_energy_mass.dimensionality
        self.assertIn("[energy/mass]", dims_int_energy_mass)
        self.assertEqual(dims_int_energy_mass["[energy]"], 1.0)
        self.assertEqual(dims_int_energy_mass["[mass]"], -1.0)
        dims_volume_mass = self.phase.volume_mass.dimensionality
        self.assertIn("[length/mass]", dims_volume_mass)
        self.assertEqual(dims_volume_mass["[length]"], 3.0)
        self.assertEqual(dims_volume_mass["[mass]"], 1.0)
        dims_gibbs_mass = self.phase.gibbs_mass.dimensionality
        self.assertIn("[energy/mass]", dims_gibbs_mass)
        self.assertEqual(dims_gibbs_mass["[energy]"], 1.0)
        self.assertEqual(dims_gibbs_mass["[mass]"], -1.0)
        dims_cp_mass = self.phase.cp_mass.dimensionality
        self.assertIn("[energy/mass/temperature]", dims_cp_mass)
        self.assertEqual(dims_cp_mass["[energy]"], 1.0)
        self.assertEqual(dims_cp_mass["[mass]"], -1.0)
        self.assertEqual(dims_cp_mass["[temperature]"], -1.0)
        dims_cv_mass = self.phase.cv.dimensionality
        self.assertIn("[energy/mass/temperature]", dims_cv_mass)
        self.assertEqual(dims_cv_mass["[energy]"], 1.0)
        self.assertEqual(dims_cv_mass["[mass]"], -1.0)
        self.assertEqual(dims_cv_mass["[temperature]"], -1.0)
        #basis-dependent (molar)
        dims_density_mole = self.phase.density_mole.dimensionality
        self.assertIn("[molar/length]", dims_density_mole)
        self.assertEqual(dims_density_mole["[molar]"], 1.0)
        self.assertEqual(dims_density_mole["[length]"], -3.0)
        dims_enthalpy_mole = self.phase.enthalpy_mole.dimensionality
        self.assertIn("[energy/molar]", dims_enthalpy_mole)
        self.assertEqual(dims_enthalpy_mole["[energy]"], 1.0)
        self.assertEqual(dims_enthalpy_mole["[molar]"], -1.0)
        dims_entropy_mole = self.phase.entropy_mole.dimensionality
        self.assertIn("[energy/molar/temperature]", dims_entropy_mole)
        self.assertEqual(dims_entropy_mole["[energy]"], 1.0)
        self.assertEqual(dims_entropy_mole["[molar]"], -1.0)
        self.assertEqual(dims_entropy_mole["[temperature]"], -1.0)
        dims_int_energy_mole = self.phase.int_energy_mole.dimensionality
        self.assertIn("[energy/molar]", dims_int_energy_mole)
        self.assertEqual(dims_int_energy_mole["[energy]"], 1.0)
        self.assertEqual(dims_int_energy_mole["[molar]"], -1.0)
        dims_volume_mole = self.phase.volume_mole.dimensionality
        self.assertIn("[length/molar]", dims_volume_mole)
        self.assertEqual(dims_volume_mole["[length]"], 3.0)
        self.assertEqual(dims_volume_mole["[molar]"], 1.0)
        dims_gibbs_mole = self.phase.gibbs_mole.dimensionality
        self.assertIn("[energy/molar]", dims_gibbs_mole)
        self.assertEqual(dims_gibbs_mole["[energy]"], 1.0)
        self.assertEqual(dims_gibbs_mole["[molar]"], -1.0)
        dims_cp_mole = self.phase.cp_mole.dimensionality
        self.assertIn("[energy/molar/temperature]", dims_cp_mole)
        self.assertEqual(dims_cp_mole["[energy]"], 1.0)
        self.assertEqual(dims_cp_mole["[molar]"], -1.0)
        self.assertEqual(dims_cp_mole["[temperature]"], -1.0)
        dims_cv_mole = self.phase.cv.dimensionality
        self.assertIn("[energy/molar/temperature]", dims_cv_mole)
        self.assertEqual(dims_cv_mole["[energy]"], 1.0)
        self.assertEqual(dims_cv_mole["[molar]"], -1.0)
        self.assertEqual(dims_cv_mole["[temperature]"], -1.0)

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

        self.phase.TDY = T0, rho0, Y0
        self.phase.DPX = rho1, P1, X1
        check_state(T1, rho1, Y1)

        self.phase.TDY = T0, rho0, Y0
        self.phase.DPY = rho1, P1, Y1
        check_state(T1, rho1, Y1)

    def test_setState_mass(self):
        self.check_setters(T1 = Q_(500.0, "K"), rho1 = Q_(1.5, 'self.basis+"/m**3"'),
                           Y1 = Q_([0.1, 0.0, 0.0, 0.1, 0.4, 0.2, 0.0, 0.0, 0.2], "dimensionless")

    def test_setState_mole(self):
        self.phase.basis = 'molar'
        self.check_setters(T1 = Q_(750.0, "K"), rho1 = Q_(0.02, 'self.basis+"/m**3"'),
                           Y1 = Q_([0.2, 0.1, 0.0, 0.3, 0.1, 0.0, 0.0, 0.2, 0.1], "dimensionless"))

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

    def check_getters(self):
        T,D,X = self.phase.TDX
        self.assertNear(T, self.phase.T)
        self.assertNear(D, self.phase.density)
        self.assertArrayNear(X, self.phase.X)

        T,D,Y = self.phase.TDY
        self.assertNear(T, self.phase.T)
        self.assertNear(D, self.phase.density)
        self.assertArrayNear(Y, self.phase.Y)

        T,D = self.phase.TD
        self.assertNear(T, self.phase.T)
        self.assertNear(D, self.phase.density)

        T,P,X = self.phase.TPX
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        T,P,Y = self.phase.TPY
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        T,P = self.phase.TP
        self.assertNear(T, self.phase.T)
        self.assertNear(P, self.phase.P)

        H,P,X = self.phase.HPX
        self.assertNear(H, self.phase.h)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        H,P,Y = self.phase.HPY
        self.assertNear(H, self.phase.h)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        H,P = self.phase.HP
        self.assertNear(H, self.phase.h)
        self.assertNear(P, self.phase.P)

        U,V,X = self.phase.UVX
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(X, self.phase.X)

        U,V,Y = self.phase.UVY
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(Y, self.phase.Y)

        U,V = self.phase.UV
        self.assertNear(U, self.phase.u)
        self.assertNear(V, self.phase.v)

        S,P,X = self.phase.SPX
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        S,P,Y = self.phase.SPY
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        S,P = self.phase.SP
        self.assertNear(S, self.phase.s)
        self.assertNear(P, self.phase.P)

        S,V,X = self.phase.SVX
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(X, self.phase.X)

        S,V,Y = self.phase.SVY
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)
        self.assertArrayNear(Y, self.phase.Y)

        S,V = self.phase.SV
        self.assertNear(S, self.phase.s)
        self.assertNear(V, self.phase.v)

        D,P,X = self.phase.DPX
        self.assertNear(D, self.phase.density)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(X, self.phase.X)

        D,P,Y = self.phase.DPY
        self.assertNear(D, self.phase.density)
        self.assertNear(P, self.phase.P)
        self.assertArrayNear(Y, self.phase.Y)

        D,P = self.phase.DP
        self.assertNear(D, self.phase.density)
        self.assertNear(P, self.phase.P)

    def test_getState_mass(self):
        self.phase.TDY = 350.0, 0.7, 'H2:0.1, H2O2:0.1, AR:0.8'
        self.check_getters()

    def test_getState_mole(self):
        self.phase.basis = 'molar'
        self.phase.TDX = 350.0, 0.01, 'H2:0.1, O2:0.3, AR:0.6'
        self.check_getters()

    def test_isothermal_compressibility(self):
        self.assertNear(self.phase.isothermal_compressibility, 1.0/self.phase.P)

class TestPureFluidUnits(utilities.CanteraTest):
    def setUp(self):
        self.water = ct.Water()

    def test_critical_properties(self):
        self.assertNear(self.water.critical_pressure, 22.089e6)
        self.assertNear(self.water.critical_temperature, 647.286)
        self.assertNear(self.water.critical_density, 317.0)

    def test_temperature_limits(self):
        co2 = ct.CarbonDioxide()
        self.assertNear(co2.min_temp, 216.54)
        self.assertNear(co2.max_temp, 1500.0)

    def test_set_state(self):
        self.water.PQ = 101325, 0.5
        self.assertNear(self.water.P, 101325)
        self.assertNear(self.water.Q, 0.5)

        self.water.TQ = 500, 0.8
        self.assertNear(self.water.T, 500)
        self.assertNear(self.water.Q, 0.8)

    def test_substance_set(self):
        self.water.TV = 400, 1.45
        self.assertNear(self.water.T, 400)
        self.assertNear(self.water.v, 1.45)
        with self.assertRaisesRegex(ct.CanteraError, 'Negative specific volume'):
            self.water.TV = 300, -1.

        self.water.PV = 101325, 1.45
        self.assertNear(self.water.P, 101325)
        self.assertNear(self.water.v, 1.45)

        self.water.UP = -1.45e7, 101325
        self.assertNear(self.water.u, -1.45e7)
        self.assertNear(self.water.P, 101325)

        self.water.VH = 1.45, -1.45e7
        self.assertNear(self.water.v, 1.45)
        self.assertNear(self.water.h, -1.45e7)

        self.water.TH = 400, -1.45e7
        self.assertNear(self.water.T, 400)
        self.assertNear(self.water.h, -1.45e7)

        self.water.SH = 5000, -1.45e7
        self.assertNear(self.water.s, 5000)
        self.assertNear(self.water.h, -1.45e7)

        self.water.ST = 5000, 400
        self.assertNear(self.water.s, 5000)
        self.assertNear(self.water.T, 400)

    def test_states(self):
        self.assertEqual(self.water._native_state, ('T', 'D'))
        self.assertNotIn('TPY', self.water._full_states.values())
        self.assertIn('TQ', self.water._partial_states.values())

    def test_set_Q(self):
        self.water.TQ = 500, 0.0
        p = self.water.P
        self.water.Q = 0.8
        self.assertNear(self.water.P, p)
        self.assertNear(self.water.T, 500)
        self.assertNear(self.water.Q, 0.8)

        self.water.TP = 650, 101325
        with self.assertRaises(ct.CanteraError):
            self.water.Q = 0.1

        self.water.TP = 300, 101325
        with self.assertRaises(ValueError):
            self.water.Q = 0.3

    def test_set_minmax(self):
        self.water.TP = self.water.min_temp, 101325
        self.assertNear(self.water.T, self.water.min_temp)

        self.water.TP = self.water.max_temp, 101325
        self.assertNear(self.water.T, self.water.max_temp)

    def check_fd_properties(self, T1, P1, T2, P2, tol):
        # Properties which are computed as finite differences
        self.water.TP = T1, P1
        h1a = self.water.enthalpy_mass
        cp1 = self.water.cp_mass
        cv1 = self.water.cv_mass
        k1 = self.water.isothermal_compressibility
        alpha1 = self.water.thermal_expansion_coeff
        h1b = self.water.enthalpy_mass

        self.water.TP = T2, P2
        h2a = self.water.enthalpy_mass
        cp2 = self.water.cp_mass
        cv2 = self.water.cv_mass
        k2 = self.water.isothermal_compressibility
        alpha2 = self.water.thermal_expansion_coeff
        h2b = self.water.enthalpy_mass

        self.assertNear(cp1, cp2, tol)
        self.assertNear(cv1, cv2, tol)
        self.assertNear(k1, k2, tol)
        self.assertNear(alpha1, alpha2, tol)

        # calculating these finite difference properties should not perturb the
        # state of the object (except for checks on edge cases)
        self.assertNear(h1a, h1b, 1e-9)
        self.assertNear(h2a, h2b, 1e-9)

    def test_properties_near_min(self):
        self.check_fd_properties(self.water.min_temp*(1+1e-5), 101325,
                                 self.water.min_temp*(1+1e-4), 101325, 1e-2)

    def test_properties_near_max(self):
        self.check_fd_properties(self.water.max_temp*(1-1e-5), 101325,
                                 self.water.max_temp*(1-1e-4), 101325, 1e-2)

    def test_properties_near_sat1(self):
        for T in [340,390,420]:
            self.water.TQ = T, 0.0
            P = self.water.P
            self.check_fd_properties(T, P+0.01, T, P+0.5, 1e-4)

    def test_properties_near_sat2(self):
        for T in [340,390,420]:
            self.water.TQ = T, 0.0
            P = self.water.P
            self.check_fd_properties(T, P-0.01, T, P-0.5, 1e-4)

    def test_isothermal_compressibility_lowP(self):
        # Low-pressure limit corresponds to ideal gas
        ref = ct.Solution('gri30.xml')
        ref.TPX = 450, 12, 'H2O:1.0'
        self.water.TP = 450, 12
        self.assertNear(ref.isothermal_compressibility,
                        self.water.isothermal_compressibility, 1e-5)

    def test_thermal_expansion_coeff_lowP(self):
        # Low-pressure limit corresponds to ideal gas
        ref = ct.Solution('gri30.xml')
        ref.TPX = 450, 12, 'H2O:1.0'
        self.water.TP = 450, 12
        self.assertNear(ref.thermal_expansion_coeff,
                        self.water.thermal_expansion_coeff, 1e-5)

    def test_thermal_expansion_coeff_TD(self):
        for T in [440, 550, 660]:
            self.water.TD = T, 0.1
            self.assertNear(T * self.water.thermal_expansion_coeff, 1.0, 1e-2)

    def test_pq_setter_triple_check(self):
        self.water.PQ = 101325, .2
        T = self.water.T
        # change T such that it would result in a Psat larger than P
        self.water.TP = 400, 101325
        # ensure that correct triple point pressure is recalculated
        # (necessary as this value is not stored by the C++ base class)
        self.water.PQ = 101325, .2
        self.assertNear(T, self.water.T, 1e-9)
        with self.assertRaisesRegex(ct.CanteraError, 'below triple point'):
            # min_temp is triple point temperature
            self.water.TP = self.water.min_temp, 101325
            P = self.water.P_sat # triple-point pressure
            self.water.PQ = .999*P, .2

    def test_quality_exceptions(self):
        # Critical point
        self.water.TP = 300, ct.one_atm
        self.water.TQ = self.water.critical_temperature, .5
        self.assertNear(self.water.P, self.water.critical_pressure)
        self.water.TP = 300, ct.one_atm
        self.water.PQ = self.water.critical_pressure, .5
        self.assertNear(self.water.T, self.water.critical_temperature)

        # Supercritical
        with self.assertRaisesRegex(ct.CanteraError, 'supercritical'):
            self.water.TQ = 1.001 * self.water.critical_temperature, 0.
        with self.assertRaisesRegex(ct.CanteraError, 'supercritical'):
            self.water.PQ = 1.001 * self.water.critical_pressure, 0.

        # Q negative
        with self.assertRaisesRegex(ct.CanteraError, 'Invalid vapor fraction'):
            self.water.TQ = 373.15, -.001
        with self.assertRaisesRegex(ct.CanteraError, 'Invalid vapor fraction'):
            self.water.PQ = ct.one_atm, -.001

        # Q larger than one
        with self.assertRaisesRegex(ct.CanteraError, 'Invalid vapor fraction'):
            self.water.TQ = 373.15, 1.001
        with self.assertRaisesRegex(ct.CanteraError, 'Invalid vapor fraction'):
            self.water.PQ = ct.one_atm, 1.001

    def test_saturated_mixture(self):
        self.water.TP = 300, ct.one_atm
        with self.assertRaisesRegex(ct.CanteraError, 'Saturated mixture detected'):
            self.water.TP = 300, self.water.P_sat

        w = ct.Water()

        # Saturated vapor
        self.water.TQ = 373.15, 1.
        self.assertEqual(self.water.phase_of_matter, 'liquid-gas-mix')
        w.TP = self.water.T, .999 * self.water.P_sat
        self.assertNear(self.water.cp, w.cp, 1.e-3)
        self.assertNear(self.water.cv, w.cv, 1.e-3)
        self.assertNear(self.water.thermal_expansion_coeff, w.thermal_expansion_coeff, 1.e-3)
        self.assertNear(self.water.isothermal_compressibility, w.isothermal_compressibility, 1.e-3)

        # Saturated mixture
        self.water.TQ = 373.15, .5
        self.assertEqual(self.water.phase_of_matter, 'liquid-gas-mix')
        self.assertEqual(self.water.cp, np.inf)
        self.assertTrue(np.isnan(self.water.cv))
        self.assertEqual(self.water.isothermal_compressibility, np.inf)
        self.assertEqual(self.water.thermal_expansion_coeff, np.inf)

        # Saturated liquid
        self.water.TQ = 373.15, 0.
        self.assertEqual(self.water.phase_of_matter, 'liquid-gas-mix')
        w.TP = self.water.T, 1.001 * self.water.P_sat
        self.assertNear(self.water.cp, w.cp, 1.e-3)
        self.assertNear(self.water.cv, w.cv, 1.e-3)
        self.assertNear(self.water.thermal_expansion_coeff, w.thermal_expansion_coeff, 1.e-3)
        self.assertNear(self.water.isothermal_compressibility, w.isothermal_compressibility, 1.e-3)

    def test_saturation_near_limits(self):
        # Low temperature limit (triple point)
        self.water.TP = 300, ct.one_atm
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = self.water.min_temp, ct.one_atm
        psat = self.water.P_sat
        self.water.TP = 300, ct.one_atm
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = 300, psat
        self.assertNear(self.water.T_sat, self.water.min_temp)

        # High temperature limit (critical point) - saturation temperature
        self.water.TP = 300, ct.one_atm
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = self.water.critical_temperature, self.water.critical_pressure
        self.assertNear(self.water.T_sat, self.water.critical_temperature)

        # High temperature limit (critical point) - saturation pressure
        self.water.TP = 300, ct.one_atm
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = self.water.critical_temperature, self.water.critical_pressure
        self.assertNear(self.water.P_sat, self.water.critical_pressure)

        # Supercricital
        with self.assertRaisesRegex(ct.CanteraError, 'Illegal temperature value'):
            self.water.TP = 1.001 * self.water.critical_temperature, self.water.critical_pressure
            self.water.P_sat
        with self.assertRaisesRegex(ct.CanteraError, 'Illegal pressure value'):
            self.water.TP = self.water.critical_temperature, 1.001 * self.water.critical_pressure
            self.water.T_sat

        # Below triple point
        with self.assertRaisesRegex(ct.CanteraError, 'Illegal temperature'):
            self.water.TP = .999 * self.water.min_temp, ct.one_atm
            self.water.P_sat
        # @TODO: test disabled pending fix of GitHub issue #605
        # with self.assertRaisesRegex(ct.CanteraError, 'Illegal pressure value'):
        #     self.water.TP = 300, .999 * psat
        #     self.water.T_sat

    def test_TPQ(self):
        self.water.TQ = 400, 0.8
        T, P, Q = self.water.TPQ
        self.assertNear(T, 400)
        self.assertNear(Q, 0.8)

        # a supercritical state
        self.water.TPQ = 800, 3e7, 1
        self.assertNear(self.water.T, 800)
        self.assertNear(self.water.P, 3e7)

        self.water.TPQ = T, P, Q
        self.assertNear(self.water.Q, 0.8)
        with self.assertRaisesRegex(ct.CanteraError, 'inconsistent'):
            self.water.TPQ = T, .999*P, Q
        with self.assertRaisesRegex(ct.CanteraError, 'inconsistent'):
            self.water.TPQ = T, 1.001*P, Q
        with self.assertRaises(TypeError):
            self.water.TPQ = T, P, 'spam'

        self.water.TPQ = 500, 1e5, 1  # superheated steam
        self.assertNear(self.water.P, 1e5)
        with self.assertRaisesRegex(ct.CanteraError, 'inconsistent'):
            self.water.TPQ = 500, 1e5, 0  # vapor fraction should be 1 (T < Tc)
        with self.assertRaisesRegex(ct.CanteraError, 'inconsistent'):
            self.water.TPQ = 700, 1e5, 0  # vapor fraction should be 1 (T > Tc)

    def test_phase_of_matter(self):
        self.water.TP = 300, 101325
        self.assertEqual(self.water.phase_of_matter, "liquid")
        self.water.TP = 500, 101325
        self.assertEqual(self.water.phase_of_matter, "gas")
        self.water.TP = self.water.critical_temperature*2, 101325
        self.assertEqual(self.water.phase_of_matter, "supercritical")
        self.water.TP = 300, self.water.critical_pressure*2
        self.assertEqual(self.water.phase_of_matter, "supercritical")
        self.water.TQ = 300, 0.4
        self.assertEqual(self.water.phase_of_matter, "liquid-gas-mix")

        # These cases work after fixing GH-786
        n2 = ct.Nitrogen()
        n2.TP = 100, 1000
        self.assertEqual(n2.phase_of_matter, "gas")

        co2 = ct.CarbonDioxide()
        self.assertEqual(co2.phase_of_matter, "gas")

    def test_water_backends(self):
        w = ct.Water(backend='Reynolds')
        self.assertEqual(w.thermo_model, 'PureFluid')
        w = ct.Water(backend='IAPWS95')
        self.assertEqual(w.thermo_model, 'liquid-water-IAPWS95')
        with self.assertRaisesRegex(KeyError, 'Unknown backend'):
            ct.Water('foobar')

    def test_water_iapws(self):
        w = ct.Water(backend='IAPWS95')
        self.assertNear(w.critical_density, 322.)
        self.assertNear(w.critical_temperature, 647.096)
        self.assertNear(w.critical_pressure, 22064000.0)

        # test internal TP setters (setters update temperature at constant
        # density before updating pressure)
        w.TP = 300, ct.one_atm
        dens = w.density
        w.TP = 2000, ct.one_atm # supercritical
        self.assertEqual(w.phase_of_matter, "supercritical")
        w.TP = 300, ct.one_atm # state goes from supercritical -> gas -> liquid
        self.assertNear(w.density, dens)
        self.assertEqual(w.phase_of_matter, "liquid")

        # test setters for critical conditions
        w.TP = w.critical_temperature, w.critical_pressure
        self.assertNear(w.density, 322.)
        w.TP = 2000, ct.one_atm # uses current density as initial guess
        w.TP = 273.16, ct.one_atm # uses fixed density as initial guess
        self.assertNear(w.density, 999.84376)
        self.assertEqual(w.phase_of_matter, "liquid")
        w.TP = w.T, w.P_sat
        self.assertEqual(w.phase_of_matter, "liquid")
        with self.assertRaisesRegex(ct.CanteraError, "assumes liquid phase"):
            w.TP = 273.1599999, ct.one_atm
        with self.assertRaisesRegex(ct.CanteraError, "assumes liquid phase"):
            w.TP = 500, ct.one_atm