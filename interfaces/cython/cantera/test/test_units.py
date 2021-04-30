import itertools
import numpy as np
import warnings
from pint import UnitRegistry
units = UnitRegistry()
Q_ = units.Quantity

import cantera.units as ct
import cantera as cn
from . import utilities

class CanteraUnitsTest(utilities.CanteraTest):
    def assertQuantityNear(self, a, b, atol=1.0E-12, rtol=1.0E-8):
        if not np.isclose(a, b, atol=atol, rtol=rtol):
            message = (f'AssertNear: {a:.14g} - {b:.14g} = {a-b:.14g}\n' +
                f'rtol = {rtol:10e}; atol = {atol:10e}')
            self.fail(message)

    def assertArrayQuantityNear(self, A, B, rtol=1e-8, atol=1e-12, msg=None):
        if A.shape != B.shape:
            self.fail(f"Arrays are of different lengths ({A.shape}, {B.shape})")
        isclose = np.isclose(A, B, atol=atol, rtol=rtol)
        if not isclose.all():
            bad_A = A[~isclose]
            bad_B = B[~isclose]
            message = (
                f"AssertNear: {bad_A:.14g} - {bad_B:.14g} = {bad_A - bad_B:.14g}\n" +
                f"Error for {(~isclose).sum()} element(s) exceeds rtol = {rtol:10e}," +
                f"atol = {atol:10e}."
            )
            if msg is not None:
                message = msg + "\n" + message
            self.fail(message)


class TestSolutionUnits(CanteraUnitsTest):
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
        self.phase._phase.basis = "molar"
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
        self.assertEqual(dims_T["[temperature]"], 1.0)
        dims_P = self.phase.P.dimensionality
        self.assertEqual(dims_P["[mass]"], 1.0)
        self.assertEqual(dims_P["[length]"], -1.0)
        self.assertEqual(dims_P["[time]"], -2.0)
        dims_X = self.phase.X.dimensionality
        #units container for dimensionless is empty
        self.assertEqual(len(dims_X), 0)
        dims_Y = self.phase.Y.dimensionality
        self.assertEqual(len(dims_Y), 0)
        # dims_T_sat = self.phase.T_sat.dimensionality
        # self.assertEqual(dims_T_sat["[temperature]"], 1.0)
        # dims_P_sat = self.phase.P_sat.dimensionality
        # self.assertEqual(dims_P_sat["[mass]"], 1.0)
        # self.assertEqual(dims_P_sat["[length]"], -1.0)
        # self.assertEqual(dims_P_sat["[time]"], -2.0)
        dims_atomic_weight = self.phase.atomic_weight.dimensionality
        self.assertEqual(dims_atomic_weight["[mass]"], 1.0)
        self.assertEqual(dims_atomic_weight["[substance]"], -1.0)
        dims_chemical_potentials = self.phase.chemical_potentials.dimensionality
        self.assertEqual(dims_chemical_potentials["[mass]"], 1.0)
        self.assertEqual(dims_chemical_potentials["[length]"], 2.0)
        self.assertEqual(dims_chemical_potentials["[time]"], -2.0)
        self.assertEqual(dims_chemical_potentials["[substance]"], -1.0)
        dims_concentrations = self.phase.concentrations.dimensionality
        self.assertEqual(dims_concentrations["[substance]"], 1.0)
        self.assertEqual(dims_concentrations["[length]"], -3.0)
        # dims_critical_temperature = self.phase.critical_temperature.dimensionality
        # self.assertEqual(dims_critical_temperature["[temperature]"], 1.0)
        # dims_critical_pressure = self.phase.critical_pressure.dimensionality
        # self.assertEqual(dims_critical_pressure["[mass]"], 1.0)
        # self.assertEqual(dims_critical_pressure["[length]"], -1.0)
        # self.assertEqual(dims_critical_pressure["[time]"], -2.0)
        dims_electric_potential = self.phase.electric_potential.dimensionality
        self.assertEqual(dims_electric_potential["[mass]"], 1.0)
        self.assertEqual(dims_electric_potential["[length]"], 2.0)
        self.assertEqual(dims_electric_potential["[time]"], -3.0)
        self.assertEqual(dims_electric_potential["[current]"], -1.0)
        dims_electrochemical_potentials = self.phase.electrochemical_potentials.dimensionality
        self.assertEqual(dims_electrochemical_potentials["[mass]"], 1.0)
        self.assertEqual(dims_electrochemical_potentials["[length]"], 2.0)
        self.assertEqual(dims_electrochemical_potentials["[time]"], -2.0)
        self.assertEqual(dims_electrochemical_potentials["[substance]"], -1.0)
        dims_isothermal_compressibility = self.phase.isothermal_compressibility.dimensionality
        self.assertEqual(dims_isothermal_compressibility["[mass]"], -1.0)
        self.assertEqual(dims_isothermal_compressibility["[length]"], 1.0)
        self.assertEqual(dims_isothermal_compressibility["[time]"], 2.0)
        dims_max_temp = self.phase.max_temp.dimensionality
        self.assertEqual(dims_max_temp["[temperature]"], 1.0)
        dims_mean_molecular_weight = self.phase.mean_molecular_weight.dimensionality
        self.assertEqual(dims_mean_molecular_weight["[mass]"], 1.0)
        self.assertEqual(dims_mean_molecular_weight["[substance]"], -1.0)
        dims_min_temp = self.phase.min_temp.dimensionality
        self.assertEqual(dims_min_temp["[temperature]"], 1.0)
        dims_molecular_weights = self.phase.molecular_weights.dimensionality
        self.assertEqual(dims_molecular_weights["[mass]"], 1.0)
        self.assertEqual(dims_molecular_weights["[substance]"], -1.0)
        dims_partial_molar_cp = self.phase.partial_molar_cp.dimensionality
        self.assertEqual(dims_partial_molar_cp["[mass]"], 1.0)
        self.assertEqual(dims_partial_molar_cp["[length]"], 2.0)
        self.assertEqual(dims_partial_molar_cp["[time]"], -2.0)
        self.assertEqual(dims_partial_molar_cp["[substance]"], -1.0)
        self.assertEqual(dims_partial_molar_cp["[temperature]"], -1.0)
        dims_partial_molar_enthalpies = self.phase.partial_molar_enthalpies.dimensionality
        self.assertEqual(dims_partial_molar_enthalpies["[mass]"], 1.0)
        self.assertEqual(dims_partial_molar_enthalpies["[length]"], 2.0)
        self.assertEqual(dims_partial_molar_enthalpies["[time]"], -2.0)
        self.assertEqual(dims_partial_molar_enthalpies["[substance]"], -1.0)
        dims_partial_molar_entropies = self.phase.partial_molar_entropies.dimensionality
        self.assertEqual(dims_partial_molar_entropies["[mass]"], 1.0)
        self.assertEqual(dims_partial_molar_entropies["[length]"], 2.0)
        self.assertEqual(dims_partial_molar_entropies["[time]"], -2.0)
        self.assertEqual(dims_partial_molar_entropies["[substance]"], -1.0)
        self.assertEqual(dims_partial_molar_entropies["[temperature]"], -1.0)
        dims_partial_molar_int_energies = self.phase.partial_molar_int_energies.dimensionality
        self.assertEqual(dims_partial_molar_int_energies["[mass]"], 1.0)
        self.assertEqual(dims_partial_molar_int_energies["[length]"], 2.0)
        self.assertEqual(dims_partial_molar_int_energies["[time]"], -2.0)
        self.assertEqual(dims_partial_molar_int_energies["[substance]"], -1.0)
        dims_partial_molar_volumes = self.phase.partial_molar_volumes.dimensionality
        self.assertEqual(dims_partial_molar_volumes["[length]"], 3.0)
        self.assertEqual(dims_partial_molar_volumes["[substance]"], -1.0)
        dims_reference_pressure = self.phase.reference_pressure.dimensionality
        self.assertEqual(dims_reference_pressure["[mass]"], 1.0)
        self.assertEqual(dims_reference_pressure["[length]"], -1.0)
        self.assertEqual(dims_reference_pressure["[time]"], -2.0)
        dims_thermal_expansion_coeff = self.phase.thermal_expansion_coeff.dimensionality
        self.assertEqual(dims_thermal_expansion_coeff["[temperature]"], -1.0)
        #basis-dependent (mass)
        dims_density_mass = self.phase.density_mass.dimensionality
        self.assertEqual(dims_density_mass["[mass]"], 1.0)
        self.assertEqual(dims_density_mass["[length]"], -3.0)
        dims_enthalpy_mass = self.phase.enthalpy_mass.dimensionality
        self.assertEqual(dims_enthalpy_mass["[length]"], 2.0)
        self.assertEqual(dims_enthalpy_mass["[time]"], -2.0)
        dims_entropy_mass = self.phase.entropy_mass.dimensionality
        self.assertEqual(dims_entropy_mass["[length]"], 2.0)
        self.assertEqual(dims_entropy_mass["[time]"], -2.0)
        self.assertEqual(dims_entropy_mass["[temperature]"], -1.0)
        dims_int_energy_mass = self.phase.int_energy_mass.dimensionality
        self.assertEqual(dims_int_energy_mass["[length]"], 2.0)
        self.assertEqual(dims_int_energy_mass["[time]"], -2.0)
        dims_volume_mass = self.phase.volume_mass.dimensionality
        self.assertEqual(dims_volume_mass["[length]"], 3.0)
        self.assertEqual(dims_volume_mass["[mass]"], -1.0)
        dims_gibbs_mass = self.phase.gibbs_mass.dimensionality
        self.assertEqual(dims_gibbs_mass["[length]"], 2.0)
        self.assertEqual(dims_gibbs_mass["[time]"], -2.0)
        dims_cp_mass = self.phase.cp_mass.dimensionality
        self.assertEqual(dims_cp_mass["[length]"], 2.0)
        self.assertEqual(dims_cp_mass["[time]"], -2.0)
        self.assertEqual(dims_cp_mass["[temperature]"], -1.0)
        dims_cv_mass = self.phase.cv.dimensionality
        self.assertEqual(dims_cv_mass["[length]"], 2.0)
        self.assertEqual(dims_cv_mass["[time]"], -2.0)
        self.assertEqual(dims_cv_mass["[temperature]"], -1.0)
        #basis-dependent (molar)
        dims_density_mole = self.phase.density_mole.dimensionality
        self.assertEqual(dims_density_mole["[substance]"], 1.0)
        self.assertEqual(dims_density_mole["[length]"], -3.0)
        dims_enthalpy_mole = self.phase.enthalpy_mole.dimensionality
        self.assertEqual(dims_enthalpy_mole["[mass]"], 1.0)
        self.assertEqual(dims_enthalpy_mole["[length]"], 2.0)
        self.assertEqual(dims_enthalpy_mole["[time]"], -2.0)
        self.assertEqual(dims_enthalpy_mole["[substance]"], -1.0)
        dims_entropy_mole = self.phase.entropy_mole.dimensionality
        self.assertEqual(dims_entropy_mole["[mass]"], 1.0)
        self.assertEqual(dims_entropy_mole["[length]"], 2.0)
        self.assertEqual(dims_entropy_mole["[time]"], -2.0)
        self.assertEqual(dims_entropy_mole["[substance]"], -1.0)
        self.assertEqual(dims_entropy_mole["[temperature]"], -1.0)
        dims_int_energy_mole = self.phase.int_energy_mole.dimensionality
        self.assertEqual(dims_int_energy_mole["[mass]"], 1.0)
        self.assertEqual(dims_int_energy_mole["[length]"], 2.0)
        self.assertEqual(dims_int_energy_mole["[time]"], -2.0)
        self.assertEqual(dims_int_energy_mole["[substance]"], -1.0)
        dims_volume_mole = self.phase.volume_mole.dimensionality
        self.assertEqual(dims_volume_mole["[length]"], 3.0)
        self.assertEqual(dims_volume_mole["[substance]"], -1.0)
        dims_gibbs_mole = self.phase.gibbs_mole.dimensionality
        self.assertEqual(dims_gibbs_mole["[mass]"], 1.0)
        self.assertEqual(dims_gibbs_mole["[length]"], 2.0)
        self.assertEqual(dims_gibbs_mole["[time]"], -2.0)
        self.assertEqual(dims_gibbs_mole["[substance]"], -1.0)
        dims_cp_mole = self.phase.cp_mole.dimensionality
        self.assertEqual(dims_cp_mole["[mass]"], 1.0)
        self.assertEqual(dims_cp_mole["[length]"], 2.0)
        self.assertEqual(dims_cp_mole["[time]"], -2.0)
        self.assertEqual(dims_cp_mole["[substance]"], -1.0)
        self.assertEqual(dims_cp_mole["[temperature]"], -1.0)
        dims_cv_mole = self.phase.cv_mole.dimensionality
        self.assertEqual(dims_cv_mole["[mass]"], 1.0)
        self.assertEqual(dims_cv_mole["[length]"], 2.0)
        self.assertEqual(dims_cv_mole["[time]"], -2.0)
        self.assertEqual(dims_cv_mole["[substance]"], -1.0)
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
            self.assertQuantityNear(self.phase.T, T)
            self.assertQuantityNear(self.phase.density, rho)
            self.assertArrayQuantityNear(self.phase.Y, Y)

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
        self.check_setters(T1 = Q_(500.0, "K"), rho1 = Q_(1.5, "kg/m**3"),
                           Y1 = Q_([0.1, 0.0, 0.0, 0.1, 0.4, 0.2, 0.0, 0.0, 0.2, 0.0], "dimensionless"))

    # def test_setState_mole(self):
    #     self.phase.basis = 'molar'
    #     self.check_setters(T1 = Q_(750.0, "K"), rho1 = Q_(0.02, "kmol/m**3"),
    #                        Y1 = Q_([0.2, 0.1, 0.0, 0.3, 0.1, 0.0, 0.0, 0.2, 0.1, 0.0], "dimensionless"))

    def test_setters_hold_constant(self):
        props = ('T','P','s','h','u','v','X','Y')
        pairs = [('TP', 'T', 'P'), ('SP', 's', 'P'),
                 ('UV', 'u', 'v')]

        self.phase.TDX = Q_(1000, "K"), Q_(1.5, "kg/m**3"), 'H2O:0.1, O2:0.95, AR:3.0'
        values = {}
        for p in props:
            values[p] = getattr(self.phase, p)

        for pair, first, second in pairs:
            self.phase.TDX = Q_(500, "K"), Q_(2.5, "kg/m**3"), 'H2:0.1, O2:1.0, AR:3.0'
            first_val = getattr(self.phase, first)
            second_val = getattr(self.phase, second)

            setattr(self.phase, pair, (values[first], None))
            self.assertQuantityNear(getattr(self.phase, first), values[first])
            self.assertQuantityNear(getattr(self.phase, second), second_val)

            self.phase.TDX = Q_(500, "K"), Q_(2.5, "kg/m**3"), 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair, (None, values[second]))
            self.assertQuantityNear(getattr(self.phase, first), first_val)
            self.assertQuantityNear(getattr(self.phase, second), values[second])

            self.phase.TDX = Q_(500, "K"), Q_(2.5, "kg/m**3"), 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'X', (None, None, values['X']))
            self.assertQuantityNear(getattr(self.phase, first), first_val)
            self.assertQuantityNear(getattr(self.phase, second), second_val)

            self.phase.TDX = Q_(500, "K"), Q_(2.5, "kg/m**3"), 'H2:0.1, O2:1.0, AR:3.0'
            setattr(self.phase, pair + 'Y', (None, None, values['Y']))
            self.assertQuantityNear(getattr(self.phase, first), first_val)
            self.assertQuantityNear(getattr(self.phase, second), second_val)

    def check_getters(self):
        T, D, X = self.phase.TDX
        self.assertQuantityNear(T, self.phase.T)
        self.assertQuantityNear(D, self.phase.density)
        self.assertArrayQuantityNear(X, self.phase.X)

        T, D, Y = self.phase.TDY
        self.assertQuantityNear(T, self.phase.T)
        self.assertQuantityNear(D, self.phase.density)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        T, D = self.phase.TD
        self.assertQuantityNear(T, self.phase.T)
        self.assertQuantityNear(D, self.phase.density)

        T, P, X = self.phase.TPX
        self.assertQuantityNear(T, self.phase.T)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(X, self.phase.X)

        T, P, Y = self.phase.TPY
        self.assertQuantityNear(T, self.phase.T)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        T, P = self.phase.TP
        self.assertQuantityNear(T, self.phase.T)
        self.assertQuantityNear(P, self.phase.P)

        H, P, X = self.phase.HPX
        self.assertQuantityNear(H, self.phase.h)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(X, self.phase.X)

        H, P, Y = self.phase.HPY
        self.assertQuantityNear(H, self.phase.h)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        H, P = self.phase.HP
        self.assertQuantityNear(H, self.phase.h)
        self.assertQuantityNear(P, self.phase.P)

        U, V, X = self.phase.UVX
        self.assertQuantityNear(U, self.phase.u)
        self.assertQuantityNear(V, self.phase.v)
        self.assertArrayQuantityNear(X, self.phase.X)

        U, V, Y = self.phase.UVY
        self.assertQuantityNear(U, self.phase.u)
        self.assertQuantityNear(V, self.phase.v)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        U, V = self.phase.UV
        self.assertQuantityNear(U, self.phase.u)
        self.assertQuantityNear(V, self.phase.v)

        S, P, X = self.phase.SPX
        self.assertQuantityNear(S, self.phase.s)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(X, self.phase.X)

        S, P, Y = self.phase.SPY
        self.assertQuantityNear(S, self.phase.s)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        S, P = self.phase.SP
        self.assertQuantityNear(S, self.phase.s)
        self.assertQuantityNear(P, self.phase.P)

        S, V, X = self.phase.SVX
        self.assertQuantityNear(S, self.phase.s)
        self.assertQuantityNear(V, self.phase.v)
        self.assertArrayQuantityNear(X, self.phase.X)

        S, V, Y = self.phase.SVY
        self.assertQuantityNear(S, self.phase.s)
        self.assertQuantityNear(V, self.phase.v)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        S, V = self.phase.SV
        self.assertQuantityNear(S, self.phase.s)
        self.assertQuantityNear(V, self.phase.v)

        D, P, X = self.phase.DPX
        self.assertQuantityNear(D, self.phase.density)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(X, self.phase.X)

        D, P, Y = self.phase.DPY
        self.assertQuantityNear(D, self.phase.density)
        self.assertQuantityNear(P, self.phase.P)
        self.assertArrayQuantityNear(Y, self.phase.Y)

        D, P = self.phase.DP
        self.assertQuantityNear(D, self.phase.density)
        self.assertQuantityNear(P, self.phase.P)

    def test_getState_mass(self):
        self.phase.TDY = Q_(350.0, "K"), Q_(0.7, "kg/m**3"), 'H2:0.1, H2O2:0.1, AR:0.8'
        self.check_getters()

    def test_getState_mole(self):
        self.phase.basis = 'molar'
        self.phase.TDX = Q_(350.0, "K"), Q_(0.01, "kg/m**3"), 'H2:0.1, O2:0.3, AR:0.6'
        self.check_getters()

    def test_isothermal_compressibility(self):
        self.assertQuantityNear(self.phase.isothermal_compressibility, 1.0/self.phase.P)

class TestPureFluidUnits(CanteraUnitsTest):
    def setUp(self):
        self.water = ct.Water()

    def test_critical_properties(self):
        self.assertQuantityNear(self.water.critical_pressure, 22.089e6 * ct.units.Pa)
        self.assertQuantityNear(self.water.critical_temperature, 647.286 * ct.units.K)
        self.assertQuantityNear(self.water.critical_density, 317.0 * ct.units.dimensionless)

    def test_temperature_limits(self):
        co2 = ct.CarbonDioxide()
        self.assertQuantityNear(co2.min_temp, 216.54 * ct.units.K)
        self.assertQuantityNear(co2.max_temp, 1500.0 * ct.units.K)

    def test_set_state(self):
        self.water.PQ = Q_(101325, "Pa"), Q_(0.5, "dimensionless")
        self.assertQuantityNear(self.water.P, 101325 * ct.units.Pa)
        self.assertQuantityNear(self.water.Q, 0.5 * ct.units.dimensionless)

        self.water.TQ = Q_(500, "K"), Q_(0.8, "dimensionless")
        self.assertQuantityNear(self.water.T, 500 * ct.units.K)
        self.assertQuantityNear(self.water.Q, 0.8 * ct.units.dimensionless)

    def test_substance_set(self):
        self.water.TV = Q_(400, "K"), Q_(1.45, "m**3/kg")
        self.assertQuantityNear(self.water.T, 400 * ct.units.K)
        self.assertQuantityNear(self.water.v, 1.45 * ct.units.m**3 / ct.units.kg)
        with self.assertRaisesRegex(cn.CanteraError, 'Negative specific volume'):
            self.water.TV = Q_(300, "K"), Q_(-1, "m**3/kg")

        self.water.PV = Q_(101325, "Pa"), Q_(1.45, "m**3/kg")
        self.assertQuantityNear(self.water.P, 101325 * ct.units.Pa)
        self.assertQuantityNear(self.water.v, 1.45 * ct.units.m**3 / ct.units.kg)

        self.water.UP = Q_(-1.45e7, "J/kg"), Q_(101325, "Pa")
        self.assertQuantityNear(self.water.u, -1.45e7 * ct.units.J / ct.units.kg)
        self.assertQuantityNear(self.water.P, 101325 * ct.units.Pa)

        self.water.VH = Q_(1.45, "m**3/kg"), Q_(-1.45e7, "J/kg")
        self.assertQuantityNear(self.water.v, 1.45 * ct.units.m**3 / ct.units.kg)
        self.assertQuantityNear(self.water.h, -1.45e7 * ct.units.J / ct.units.kg)

        self.water.TH = Q_(400, "K"), Q_(-1.45e7, "J/kg")
        self.assertQuantityNear(self.water.T, 400 * ct.units.K)
        self.assertQuantityNear(self.water.h, -1.45e7 * ct.units.J / ct.units.kg)

        self.water.SH = Q_(5000, "J/kg/K"), Q_(-1.45e7, "J/kg")
        self.assertQuantityNear(self.water.s, 5000 * ct.units.J / ct.units.kg / ct.units.K)
        self.assertQuantityNear(self.water.h, -1.45e7 * ct.units.J / ct.units.kg)

        self.water.ST = Q_(5000, "J/kg/K"), Q_(400, "K")
        self.assertQuantityNear(self.water.s, 5000 * ct.units.J / ct.units.kg / ct.units.K)
        self.assertQuantityNear(self.water.T, 400 * ct.units.K)

    def test_states(self):
        self.assertEqual(self.water._native_state, ('T', 'D'))
        self.assertNotIn('TPY', self.water._full_states.values())
        self.assertIn('TQ', self.water._partial_states.values())

    def test_set_Q(self):
        self.water.TQ = Q_(500, "K"), Q_(0.0, "dimensionless")
        p = self.water.P
        self.water.Q = Q_(0.8, "dimensionless")
        self.assertQuantityNear(self.water.P, p)
        self.assertQuantityNear(self.water.T, 500 * ct.units.K)
        self.assertQuantityNear(self.water.Q, 0.8 * ct.units.dimensionless)

        self.water.TP = Q_(650, "K"), Q_(101325, "Pa")
        with self.assertRaises(cn.CanteraError):
            self.water.Q = Q_(0.1, "dimensionless")

        self.water.TP = Q_(300, "K"), Q_(101325, "Pa")
        with self.assertRaises(ValueError):
            self.water.Q = Q_(0.3, "dimensionless")

    def test_set_minmax(self):
        self.water.TP = self.water.min_temp, Q_(101325, "Pa")
        self.assertQuantityNear(self.water.T, self.water.min_temp)

        self.water.TP = self.water.max_temp, Q_(101325, "Pa")
        self.assertQuantityNear(self.water.T, self.water.max_temp)

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

        self.assertQuantityNear(cp1, cp2, tol, tol)
        self.assertQuantityNear(cv1, cv2, tol, tol)
        self.assertQuantityNear(k1, k2, tol, tol)
        self.assertQuantityNear(alpha1, alpha2, tol, tol)

        # calculating these finite difference properties should not perturb the
        # state of the object (except for checks on edge cases)
        self.assertQuantityNear(h1a, h1b, 1e-9)
        self.assertQuantityNear(h2a, h2b, 1e-9)

    def test_properties_near_min(self):
        self.check_fd_properties(self.water.min_temp*(1+1e-5), 101325 * ct.units.Pa,
                                 self.water.min_temp*(1+1e-4), 101325 * ct.units.Pa, 1e-2)

    def test_properties_near_max(self):
        self.check_fd_properties(self.water.max_temp*(1-1e-5), 101325 * ct.units.Pa,
                                 self.water.max_temp*(1-1e-4), 101325 * ct.units.Pa, 1e-2)

    def test_properties_near_sat1(self):
        for T in [340 * ct.units.K, 390 * ct.units.K, 420 * ct.units.K]:
            self.water.TQ = T, Q_(0.0, "dimensionless")
            P = self.water.P
            self.check_fd_properties(T, P+0.01 * ct.units.Pa, T, P+0.5 * ct.units.Pa, 1e-4)

    def test_properties_near_sat2(self):
        for T in [340 * ct.units.K, 390 * ct.units.K, 420 * ct.units.K]:
            self.water.TQ = T, Q_(0.0, "dimensionless")
            P = self.water.P
            self.check_fd_properties(T, P-0.01 * ct.units.Pa, T, P-0.5 * ct.units.Pa, 1e-4)

    def test_isothermal_compressibility_lowP(self):
        # Low-pressure limit corresponds to ideal gas
        ref = ct.Solution('gri30.xml')
        ref.TPX = Q_(450, "K"), Q_(12, "Pa"), 'H2O:1.0'
        self.water.TP = Q_(450, "K"), Q_(12, "Pa")
        self.assertQuantityNear(ref.isothermal_compressibility,
                        self.water.isothermal_compressibility, 1e-5 * ct.units.dimensionless)

    def test_thermal_expansion_coeff_lowP(self):
        # Low-pressure limit corresponds to ideal gas
        ref = ct.Solution('gri30.xml')
        ref.TPX = Q_(450, "K"), Q_(12, "Pa"), 'H2O:1.0'
        self.water.TP = Q_(450, "K"), Q_(12, "Pa")
        self.assertQuantityNear(ref.thermal_expansion_coeff,
                        self.water.thermal_expansion_coeff, 1e-5 * ct.units.dimensionless)

    def test_thermal_expansion_coeff_TD(self):
        for T in [440 * ct.units.K, 550 * ct.units.K, 660 * ct.units.K]:
            self.water.TD = T, Q_(0.1, "kg/m**3")
            self.assertQuantityNear(T * self.water.thermal_expansion_coeff, 1.0, 1e-2)

    def test_pq_setter_triple_check(self):
        self.water.PQ = Q_(101325, "Pa"), Q_(.2, "dimensionless")
        T = self.water.T
        # change T such that it would result in a Psat larger than P
        self.water.TP = Q_(400, "K"), Q_(101325, "Pa")
        # ensure that correct triple point pressure is recalculated
        # (necessary as this value is not stored by the C++ base class)
        self.water.PQ = Q_(101325, "Pa"), Q_(.2, "dimensionless")
        self.assertQuantityNear(T, self.water.T, 1e-9)
        with self.assertRaisesRegex(cn.CanteraError, 'below triple point'):
            # min_temp is triple point temperature
            self.water.TP = self.water.min_temp, Q_(101325, "Pa")
            P = self.water.P_sat # triple-point pressure
            self.water.PQ = .999*P, Q_(.2, "dimensionless")

    def test_quality_exceptions(self):
        # Critical point
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        self.water.TQ = self.water.critical_temperature, Q_(.5, "dimensionless")
        self.assertQuantityNear(self.water.P, self.water.critical_pressure)
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        self.water.PQ = self.water.critical_pressure, Q_(.5, "dimensionless")
        self.assertQuantityNear(self.water.T, self.water.critical_temperature)

        # Supercritical
        with self.assertRaisesRegex(cn.CanteraError, 'supercritical'):
            self.water.TQ = 1.001 * self.water.critical_temperature, Q_(0, "dimensionless")
        with self.assertRaisesRegex(cn.CanteraError, 'supercritical'):
            self.water.PQ = 1.001 * self.water.critical_pressure, Q_(0, "dimensionless")

        # Q negative
        with self.assertRaisesRegex(cn.CanteraError, 'Invalid vapor fraction'):
            self.water.TQ = Q_(373.15, "K"), Q_(-.001, "dimensionless")
        with self.assertRaisesRegex(cn.CanteraError, 'Invalid vapor fraction'):
            self.water.PQ = Q_(cn.one_atm, "Pa"), Q_(-.001, "dimensionless")

        # Q larger than one
        with self.assertRaisesRegex(cn.CanteraError, 'Invalid vapor fraction'):
            self.water.TQ = Q_(373.15, "K"), Q_(1.001, "dimensionless")
        with self.assertRaisesRegex(cn.CanteraError, 'Invalid vapor fraction'):
            self.water.PQ = Q_(cn.one_atm, "Pa"), Q_(1.001, "dimensionless")

    def test_saturated_mixture(self):
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        with self.assertRaisesRegex(cn.CanteraError, 'Saturated mixture detected'):
            self.water.TP = Q_(300, "K"), self.water.P_sat

        w = ct.Water()

        # Saturated vapor
        self.water.TQ = Q_(373.15, "K"), Q_(1, "dimensionless")
        self.assertEqual(self.water.phase_of_matter, 'liquid-gas-mix')
        w.TP = self.water.T, .999 * self.water.P_sat
        self.assertQuantityNear(self.water.cp, w.cp, 1.e-1)
        self.assertQuantityNear(self.water.cv, w.cv, 1.e-1)
        self.assertQuantityNear(self.water.thermal_expansion_coeff, w.thermal_expansion_coeff, 1.e-1)
        self.assertQuantityNear(self.water.isothermal_compressibility, w.isothermal_compressibility, 1.e-1)

        # Saturated mixture
        self.water.TQ = Q_(373.15, "K"), Q_(.5, "dimensionless")
        self.assertEqual(self.water.phase_of_matter, 'liquid-gas-mix')
        self.assertEqual(self.water.cp, np.inf * ct.units.J / ct.units.kg / ct.units.K)
        self.assertTrue(np.isnan(self.water.cv))
        self.assertEqual(self.water.isothermal_compressibility, np.inf * 1 / ct.units.Pa)
        self.assertEqual(self.water.thermal_expansion_coeff, np.inf * 1 / ct.units.K)

        # Saturated liquid
        self.water.TQ = Q_(373.15, "K"), Q_(0, "dimensionless")
        self.assertEqual(self.water.phase_of_matter, 'liquid-gas-mix')
        w.TP = self.water.T, 1.001 * self.water.P_sat
        self.assertQuantityNear(self.water.cp, w.cp, 9.e-1)
        self.assertQuantityNear(self.water.cv, w.cv, 9.e-1)
        self.assertQuantityNear(self.water.thermal_expansion_coeff, w.thermal_expansion_coeff, 1.e-1)
        self.assertQuantityNear(self.water.isothermal_compressibility, w.isothermal_compressibility, 1.e-1)

    def test_saturation_near_limits(self):
        # Low temperature limit (triple point)
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = self.water.min_temp, Q_(cn.one_atm, "Pa")
        psat = self.water.P_sat
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = Q_(300, "K"), psat
        self.assertQuantityNear(self.water.T_sat, self.water.min_temp)

        # High temperature limit (critical point) - saturation temperature
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = self.water.critical_temperature, self.water.critical_pressure
        self.assertQuantityNear(self.water.T_sat, self.water.critical_temperature)

        # High temperature limit (critical point) - saturation pressure
        self.water.TP = Q_(300, "K"), Q_(cn.one_atm, "Pa")
        self.water.P_sat # ensure that solver buffers sufficiently different values
        self.water.TP = self.water.critical_temperature, self.water.critical_pressure
        self.assertQuantityNear(self.water.P_sat, self.water.critical_pressure)

        # Supercricital
        with self.assertRaisesRegex(cn.CanteraError, 'Illegal temperature value'):
            self.water.TP = 1.001 * self.water.critical_temperature, self.water.critical_pressure
            self.water.P_sat
        with self.assertRaisesRegex(cn.CanteraError, 'Illegal pressure value'):
            self.water.TP = self.water.critical_temperature, 1.001 * self.water.critical_pressure
            self.water.T_sat

        # Below triple point
        with self.assertRaisesRegex(cn.CanteraError, 'Illegal temperature'):
            self.water.TP = .999 * self.water.min_temp, Q_(cn.one_atm, "Pa")
            self.water.P_sat
        # @TODO: test disabled pending fix of GitHub issue #605
        # with self.assertRaisesRegex(ct.CanteraError, 'Illegal pressure value'):
        #     self.water.TP = 300, .999 * psat
        #     self.water.T_sat

    def test_TPQ(self):
        self.water.TQ = Q_(400, "K"), Q_(0.8, "dimensionless")
        T, P, Q = self.water.TPQ
        self.assertQuantityNear(T, 400 * ct.units.K)
        self.assertQuantityNear(Q, 0.8 * ct.units.dimensionless)

        # a supercritical state
        self.water.TPQ = Q_(800, "K"), Q_(3e7, "Pa"), Q_(1, "dimensionless")
        self.assertQuantityNear(self.water.T, 800 * ct.units.K)
        self.assertQuantityNear(self.water.P, 3e7 * ct.units.Pa)

        self.water.TPQ = T, P, Q
        self.assertQuantityNear(self.water.Q, 0.8 * ct.units.dimensionless)
        with self.assertRaisesRegex(cn.CanteraError, 'inconsistent'):
            self.water.TPQ = T, .999*P, Q
        with self.assertRaisesRegex(cn.CanteraError, 'inconsistent'):
            self.water.TPQ = T, 1.001*P, Q
        with self.assertRaises(TypeError):
            self.water.TPQ = T, P, 'spam'

        self.water.TPQ = Q_(500, "K"), Q_(1e5, "Pa"), Q_(1, "dimensionless")  # superheated steam
        self.assertQuantityNear(self.water.P, 1e5 * ct.units.Pa)
        with self.assertRaisesRegex(cn.CanteraError, 'inconsistent'):
            self.water.TPQ = Q_(500, "K"), Q_(1e5, "Pa"), Q_(0, "dimensionless")  # vapor fraction should be 1 (T < Tc)
        with self.assertRaisesRegex(cn.CanteraError, 'inconsistent'):
            self.water.TPQ = Q_(700, "K"), Q_(1e5, "Pa"), Q_(0, "dimensionless")  # vapor fraction should be 1 (T > Tc)

    def test_phase_of_matter(self):
        self.water.TP = Q_(300, "K"), Q_(101325, "Pa")
        self.assertEqual(self.water.phase_of_matter, "liquid")
        self.water.TP = Q_(500, "K"), Q_(101325, "Pa")
        self.assertEqual(self.water.phase_of_matter, "gas")
        self.water.TP = self.water.critical_temperature*2, Q_(101325, "Pa")
        self.assertEqual(self.water.phase_of_matter, "supercritical")
        self.water.TP = Q_(300, "K"), self.water.critical_pressure*2
        self.assertEqual(self.water.phase_of_matter, "supercritical")
        self.water.TQ = Q_(300, "K"), Q_(0.4, "dimensionless")
        self.assertEqual(self.water.phase_of_matter, "liquid-gas-mix")

        # These cases work after fixing GH-786
        n2 = ct.Nitrogen()
        n2.TP = Q_(100, "K"), Q_(1000, "Pa")
        self.assertEqual(n2.phase_of_matter, "gas")

        co2 = ct.CarbonDioxide()
        self.assertEqual(co2.phase_of_matter, "gas")
