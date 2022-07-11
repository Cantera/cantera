import math
import re

import numpy as np
import pytest
from .utilities import unittest

import cantera as ct
from . import utilities

class TestReactor(utilities.CanteraTest):
    reactorClass = ct.Reactor

    def make_reactors(self, independent=True, n_reactors=2,
                      T1=300, P1=101325, X1='O2:1.0',
                      T2=300, P2=101325, X2='O2:1.0'):

        self.net = ct.ReactorNet()

        self.gas1 = ct.Solution('h2o2.yaml', transport_model=None)
        self.gas1.TPX = T1, P1, X1
        self.r1 = self.reactorClass(self.gas1)
        self.net.add_reactor(self.r1)

        if independent:
            self.gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        else:
            self.gas2 = self.gas1

        if n_reactors >= 2:
            self.gas2.TPX = T2, P2, X2
            self.r2 = self.reactorClass(self.gas2)
            self.net.add_reactor(self.r2)

    def add_wall(self, **kwargs):
        self.w = ct.Wall(self.r1, self.r2, **kwargs)
        return self.w

    def test_verbose(self):
        self.make_reactors(independent=False, n_reactors=1)
        self.assertFalse(self.net.verbose)
        self.net.verbose = True
        self.assertTrue(self.net.verbose)

    def test_insert(self):
        R = self.reactorClass()
        with self.assertRaisesRegex(ct.CanteraError, 'No phase'):
            R.T
        with self.assertRaisesRegex(ct.CanteraError, 'No phase'):
            R.kinetics.net_production_rates

        g = ct.Solution('h2o2.yaml', transport_model=None)
        g.TP = 300, 101325
        R.insert(g)

        self.assertNear(R.T, 300)
        self.assertEqual(len(R.kinetics.net_production_rates), g.n_species)

    def test_volume(self):
        R = self.reactorClass(volume=11)
        self.assertEqual(R.volume, 11)

        R.volume = 9
        self.assertEqual(R.volume, 9)

    def test_names(self):
        self.make_reactors()

        pattern = re.compile(r'(\d+)')
        digits1 = pattern.search(self.r1.name).group(0)
        digits2 = pattern.search(self.r2.name).group(0)

        self.assertEqual(int(digits2), int(digits1) + 1)

        self.r1.name = 'hello'
        self.assertEqual(self.r1.name, 'hello')

    def test_types(self):
        self.make_reactors()
        self.assertEqual(self.r1.type, self.reactorClass.__name__)

    def test_component_index(self):
        self.make_reactors(n_reactors=1)
        self.net.step()

        N0 = self.net.n_vars - self.gas1.n_species
        for i, name in enumerate(self.gas1.species_names):
            self.assertEqual(i + N0, self.r1.component_index(name))

    def test_component_names(self):
        self.make_reactors(n_reactors=2)
        N = self.net.n_vars // 2
        for i in range(N):
            self.assertEqual(self.r1.component_index(self.r1.component_name(i)), i)
            self.assertEqual(self.net.component_name(i),
                '{}: {}'.format(self.r1.name, self.r1.component_name(i)))
            self.assertEqual(self.net.component_name(N+i),
                '{}: {}'.format(self.r2.name, self.r2.component_name(i)))

    def test_disjoint(self):
        T1, P1 = 300, 101325
        T2, P2 = 500, 300000

        self.make_reactors(T1=T1, T2=T2, P1=P1, P2=P2)
        self.net.advance(1.0)

        # Nothing should change from the initial condition
        self.assertNear(T1, self.gas1.T)
        self.assertNear(T2, self.gas2.T)
        self.assertNear(P1, self.gas1.P)
        self.assertNear(P2, self.gas2.P)

    def test_disjoint2(self):
        T1, P1 = 300, 101325
        T2, P2 = 500, 300000

        self.make_reactors(T1=T1, T2=T2, P1=P1, P2=P2, independent=False)
        self.net.advance(1.0)

        # Nothing should change from the initial condition
        self.assertNear(T1, self.r1.T)
        self.assertNear(T2, self.r2.T)
        self.assertNear(P1, self.r1.thermo.P)
        self.assertNear(P2, self.r2.thermo.P)

    def test_derivative(self):
        T1, P1 = 300, 101325

        self.make_reactors(n_reactors=1, T1=T1, P1=P1)
        self.net.advance(1.0)

        # compare cvode derivative to numerical derivative
        dydt = self.net.get_derivative(1)
        dt = -self.net.time
        dy = -self.net.get_state()
        self.net.step()
        dt += self.net.time
        dy += self.net.get_state()
        for i in range(self.net.n_vars):
            self.assertNear(dydt[i], dy[i]/dt)

    def test_timestepping(self):
        self.make_reactors()

        tStart = 0.3
        tEnd = 10.0
        dt_max = 0.07
        t = tStart

        self.net.max_time_step = dt_max
        self.assertEqual(self.net.max_time_step, dt_max)
        self.net.set_initial_time(tStart)
        self.assertNear(self.net.time, tStart)

        while t < tEnd:
            tPrev = t
            t = self.net.step()
            self.assertTrue(t - tPrev <= 1.0001 * dt_max)
            self.assertNear(t, self.net.time)

        #self.assertNear(self.net.time, tEnd)

    def test_maxsteps(self):
        self.make_reactors()

        # set the up a case where we can't take
        # enough time-steps to reach the endtime
        max_steps = 10
        max_step_size = 1e-07
        self.net.set_initial_time(0)
        self.net.max_time_step = max_step_size
        self.net.max_steps = max_steps
        with self.assertRaisesRegex(
                ct.CanteraError, 'mxstep steps taken before reaching tout'):
            self.net.advance(1e-04)
        self.assertLessEqual(self.net.time, max_steps * max_step_size)
        self.assertEqual(self.net.max_steps, max_steps)

    def test_wall_type(self):
        self.make_reactors(P1=101325, P2=300000)
        self.add_wall(K=0.1, A=1.0)
        self.assertEqual(self.w.type, "Wall")

    def test_equalize_pressure(self):
        self.make_reactors(P1=101325, P2=300000)
        self.add_wall(K=0.1, A=1.0)

        self.assertEqual(len(self.r1.walls), 1)
        self.assertEqual(len(self.r2.walls), 1)
        self.assertEqual(self.r1.walls[0], self.w)
        self.assertEqual(self.r2.walls[0], self.w)

        self.net.advance(1.0)

        self.assertNear(self.net.time, 1.0)
        self.assertNear(self.gas1.P, self.gas2.P)
        self.assertNotAlmostEqual(self.r1.T, self.r2.T)

    def test_tolerances(self):
        def integrate(atol, rtol):
            P0 = 10 * ct.one_atm
            T0 = 1100
            X0 = 'H2:1.0, O2:0.5, AR:8.0'
            self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)
            self.net.rtol = rtol
            self.net.atol = atol

            self.assertEqual(self.net.rtol, rtol)
            self.assertEqual(self.net.atol, atol)

            tEnd = 1.0
            nSteps = 0
            t = 0

            while t < tEnd:
                t = self.net.step()
                nSteps += 1

            return nSteps

        n_baseline = integrate(1e-10, 1e-20)
        n_rtol = integrate(5e-7, 1e-20)
        n_atol = integrate(1e-10, 1e-6)

        self.assertTrue(n_baseline > n_rtol)
        self.assertTrue(n_baseline > n_atol)

    def test_advance_limits(self):
        P0 = 10 * ct.one_atm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'
        self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)

        limit_H2 = .01
        ix = self.net.global_component_index('H2', 0)
        self.r1.set_advance_limit('H2', limit_H2)
        self.assertEqual(self.net.advance_limits[ix], limit_H2)

        self.r1.set_advance_limit('H2', None)
        self.assertEqual(self.net.advance_limits[ix], -1.)

        self.r1.set_advance_limit('H2', limit_H2)
        self.net.advance_limits = None
        self.assertEqual(self.net.advance_limits[ix], -1.)

        self.r1.set_advance_limit('H2', limit_H2)
        self.net.advance_limits = 0 * self.net.advance_limits - 1.
        self.assertEqual(self.net.advance_limits[ix], -1.)

    def test_advance_with_limits(self):
        def integrate(limit_H2 = None, apply=True):
            P0 = 10 * ct.one_atm
            T0 = 1100
            X0 = 'H2:1.0, O2:0.5, AR:8.0'
            self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)
            if limit_H2 is not None:
                self.r1.set_advance_limit('H2', limit_H2)
                ix = self.net.global_component_index('H2', 0)
                self.assertEqual(self.net.advance_limits[ix], limit_H2)

            tEnd = 1.0
            tStep = 1.e-3
            nSteps = 0

            t = tStep
            while t < tEnd:
                t_curr = self.net.advance(t, apply_limit=apply)
                nSteps += 1
                if t_curr == t:
                    t += tStep

            return nSteps

        n_baseline = integrate()
        n_advance_coarse = integrate(.01)
        n_advance_fine = integrate(.001)
        n_advance_negative = integrate(-1.0)
        n_advance_override = integrate(.001, False)

        self.assertGreater(n_advance_coarse, n_baseline)
        self.assertGreater(n_advance_fine, n_advance_coarse)
        self.assertEqual(n_advance_negative, n_baseline)
        self.assertEqual(n_advance_override, n_baseline)

    def test_heat_transfer1(self):
        # Connected reactors reach thermal equilibrium after some time
        self.make_reactors(T1=300, T2=1000)
        self.add_wall(U=500, A=1.0)

        self.net.advance(10.0)
        self.assertNear(self.net.time, 10.0)
        self.assertNear(self.r1.T, self.r2.T, 5e-7)
        self.assertNotAlmostEqual(self.r1.thermo.P, self.r2.thermo.P)

    def test_advance_limits_invalid(self):
        self.make_reactors(n_reactors=1)

        with pytest.raises(ct.CanteraError, match="No component named 'spam'"):
            self.r1.set_advance_limit("spam", 0.1)

    def test_heat_transfer2(self):
        # Result should be the same if (m * cp) / (U * A) is held constant
        self.make_reactors(T1=300, T2=1000)
        self.add_wall(U=200, A=1.0)

        self.net.advance(1.0)
        T1a = self.r1.T
        T2a = self.r2.T

        self.make_reactors(T1=300, T2=1000)
        self.r1.volume = 0.25
        self.r2.volume = 0.25
        w = self.add_wall(U=100, A=0.5)

        self.assertNear(w.heat_transfer_coeff * w.area * (self.r1.T - self.r2.T),
                        w.qdot(0))
        self.net.advance(1.0)
        self.assertNear(w.heat_transfer_coeff * w.area * (self.r1.T - self.r2.T),
                        w.qdot(1.0))
        T1b = self.r1.T
        T2b = self.r2.T

        self.assertNear(T1a, T1b)
        self.assertNear(T2a, T2b)

    def test_equilibrium_UV(self):
        # Adiabatic, constant volume combustion should proceed to equilibrium
        # at constant internal energy and volume.

        P0 = 10 * ct.one_atm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'
        self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)

        self.net.advance(1.0)

        gas = ct.Solution('h2o2.yaml', transport_model=None)
        gas.TPX = T0, P0, X0
        gas.equilibrate('UV')

        self.assertNear(self.r1.T, gas.T)
        self.assertNear(self.r1.thermo.density, gas.density)
        self.assertNear(self.r1.thermo.P, gas.P)
        self.assertArrayNear(self.r1.thermo.X, gas.X)

    def test_equilibrium_HP(self):
        # Adiabatic, constant pressure combustion should proceed to equilibrium
        # at constant enthalpy and pressure.

        P0 = 10 * ct.one_atm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'

        gas1 = ct.Solution('h2o2.yaml', transport_model=None)
        gas1.TPX = T0, P0, X0
        r1 = ct.IdealGasConstPressureReactor(gas1)

        net = ct.ReactorNet()
        net.add_reactor(r1)
        net.advance(1.0)

        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TPX = T0, P0, X0
        gas2.equilibrate('HP')

        self.assertNear(r1.T, gas2.T)
        self.assertNear(r1.thermo.P, P0)
        self.assertNear(r1.thermo.density, gas2.density)
        self.assertArrayNear(r1.thermo.X, gas2.X)

    def test_wall_velocity(self):
        self.make_reactors()
        A = 0.2

        V1 = 2.0
        V2 = 5.0
        self.r1.volume = V1
        self.r2.volume = V2

        self.add_wall(A=A)

        def v(t):
            if 0 < t <= 1:
                return t
            elif 1 <= t <= 2:
                return 2 - t
            else:
                return 0.0

        self.w.set_velocity(v)
        self.net.advance(1.0)
        self.assertNear(self.w.vdot(1.0), 1.0 * A, 1e-7)
        self.net.advance(2.0)
        self.assertNear(self.w.vdot(2.0), 0.0, 1e-7)

        self.assertNear(self.r1.volume, V1 + 1.0 * A, 1e-7)
        self.assertNear(self.r2.volume, V2 - 1.0 * A, 1e-7)

    def test_disable_energy(self):
        self.make_reactors(T1=500)
        self.r1.energy_enabled = False
        self.add_wall(A=1.0, U=2500)
        self.net.advance(11.0)

        self.assertNear(self.r1.T, 500)
        self.assertNear(self.r2.T, 500)

    def test_disable_chemistry(self):
        self.make_reactors(T1=1000, n_reactors=1, X1='H2:2.0,O2:1.0')
        self.r1.chemistry_enabled = False

        self.net.advance(11.0)

        self.assertNear(self.r1.T, 1000)
        self.assertNear(self.r1.thermo.X[self.r1.thermo.species_index('H2')], 2.0/3.0)
        self.assertNear(self.r1.thermo.X[self.r1.thermo.species_index('O2')], 1.0/3.0)

    def test_heat_flux_func(self):
        self.make_reactors(T1=500, T2=300)
        self.r1.volume = 0.5

        U1a = self.r1.volume * self.r1.density * self.r1.thermo.u
        U2a = self.r2.volume * self.r2.density * self.r2.thermo.u

        V1a = self.r1.volume
        V2a = self.r2.volume

        self.add_wall(A=0.3)
        self.w.set_heat_flux(lambda t: 90000 * (1 - t**2) if t <= 1.0 else 0.0)
        Q = 0.3 * 60000

        self.net.advance(1.1)
        U1b = self.r1.volume * self.r1.density * self.r1.thermo.u
        U2b = self.r2.volume * self.r2.density * self.r2.thermo.u

        self.assertNear(V1a, self.r1.volume)
        self.assertNear(V2a, self.r2.volume)
        self.assertNear(U1a - Q, U1b, 1e-6)
        self.assertNear(U2a + Q, U2b, 1e-6)

    def test_mass_flow_controller(self):
        self.make_reactors(n_reactors=1)
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TPX = 300, 10*101325, 'H2:1.0'
        reservoir = ct.Reservoir(gas2)

        # Apply a mass flow rate that is not a smooth function of time. This
        # demonstrates the accuracy that can be achieved by having the mass flow rate
        # evaluated simultaneously with the rest of the governing equations, as opposed
        # to doing only between integrator time steps. In the latter case, larger
        # errors would be introduced when the integrator steps past the times where the
        # function's behavior changes and can't dynamically reduce the steps size for
        # the already-completed steps.
        mfc = ct.MassFlowController(reservoir, self.r1)
        # Triangular pulse with area = 0.1
        def mdot(t):
            if 0.2 <= t < 1.2:
                return 0.2 - 0.4 * abs(t - 0.7)
            else:
                return 0.0
        mfc.mass_flow_rate = mdot
        self.assertEqual(mfc.mass_flow_coeff, 1.)

        self.assertEqual(mfc.type, type(mfc).__name__)
        self.assertEqual(len(reservoir.inlets), 0)
        self.assertEqual(len(reservoir.outlets), 1)
        self.assertEqual(reservoir.outlets[0], mfc)
        self.assertEqual(len(self.r1.outlets), 0)
        self.assertEqual(len(self.r1.inlets), 1)
        self.assertEqual(self.r1.inlets[0], mfc)

        ma = self.r1.volume * self.r1.density
        Ya = self.r1.Y

        self.net.rtol = 1e-11
        self.net.max_time_step = 0.05

        self.net.advance(0.1)
        self.assertNear(mfc.mass_flow_rate, 0.)
        self.net.advance(0.3)
        self.assertNear(mfc.mass_flow_rate, 0.04)
        self.net.advance(1.0)
        self.assertNear(mfc.mass_flow_rate, 0.08)
        self.net.advance(1.2)
        self.assertNear(mfc.mass_flow_rate, 0.)

        self.net.advance(2.5)

        mb = self.r1.volume * self.r1.density
        Yb = self.r1.Y

        self.assertNear(ma + 0.1, mb)
        self.assertArrayNear(ma * Ya + 0.1 * gas2.Y, mb * Yb)

    def test_mass_flow_controller_errors(self):
        # Make sure Python error message actually gets displayed
        self.make_reactors(n_reactors=2)
        mfc = ct.MassFlowController(self.r1, self.r2)
        mfc.mass_flow_rate = lambda t: eggs

        with self.assertRaisesRegex(Exception, 'eggs'):
            self.net.step()

        with self.assertRaisesRegex(ct.CanteraError, 'NotImplementedError'):
            mfc.set_pressure_function(lambda p: p**2)

    def test_valve1(self):
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.valve_coeff = k

        self.assertEqual(self.r1.outlets, self.r2.inlets)
        self.assertEqual(valve.valve_coeff, k)
        self.assertTrue(self.r1.energy_enabled)
        self.assertTrue(self.r2.energy_enabled)
        self.net.initialize()
        self.assertNear((self.r1.thermo.P - self.r2.thermo.P) * k,
                        valve.mass_flow_rate)

        m1a = self.r1.thermo.density * self.r1.volume
        m2a = self.r2.thermo.density * self.r2.volume
        Y1a = self.r1.thermo.Y
        Y2a = self.r2.thermo.Y

        self.net.advance(0.1)

        m1b = self.r1.thermo.density * self.r1.volume
        m2b = self.r2.thermo.density * self.r2.volume

        self.assertNear((self.r1.thermo.P - self.r2.thermo.P) * k,
                        valve.mass_flow_rate)
        self.assertNear(m1a+m2a, m1b+m2b)
        Y1b = self.r1.thermo.Y
        Y2b = self.r2.thermo.Y
        self.assertArrayNear(m1a*Y1a + m2a*Y2a, m1b*Y1b + m2b*Y2b, atol=1e-10)
        self.assertArrayNear(Y1a, Y1b)

    def test_valve2(self):
        # Similar to test_valve1, but by disabling the energy equation
        # (constant T) we can compare with an analytical solution for
        # the mass of each reactor as a function of time
        self.make_reactors(P1=10*ct.one_atm)
        self.net.rtol = 1e-11
        self.r1.energy_enabled = False
        self.r2.energy_enabled = False
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.valve_coeff = k
        self.assertEqual(valve.valve_coeff, k)

        self.assertFalse(self.r1.energy_enabled)
        self.assertFalse(self.r2.energy_enabled)

        m1a = self.r1.thermo.density * self.r1.volume
        m2a = self.r2.thermo.density * self.r2.volume
        P1a = self.r1.thermo.P
        P2a = self.r2.thermo.P
        Y1 = self.r1.Y

        A = k * P1a * (1 + m2a/m1a)
        B = k * (P1a/m1a + P2a/m2a)

        for t in np.linspace(1e-5, 0.5):
            self.net.advance(t)
            m1 = self.r1.thermo.density * self.r1.volume
            m2 = self.r2.thermo.density * self.r2.volume
            self.assertNear(m2, (m2a - A/B) * np.exp(-B * t) + A/B)
            self.assertNear(m1a+m2a, m1+m2)
            self.assertArrayNear(self.r1.Y, Y1)

    def test_valve3(self):
        # This case specifies a non-linear relationship between pressure drop
        # and flow rate.
        self.make_reactors(P1=10*ct.one_atm, X1='AR:0.5, O2:0.5',
                           X2='O2:1.0')
        self.net.rtol = 1e-12
        self.net.atol = 1e-20
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.set_pressure_function(mdot)
        self.assertEqual(valve.valve_coeff, 1.)

        Y1 = self.r1.Y
        kO2 = self.gas1.species_index('O2')
        kAr = self.gas1.species_index('AR')

        def speciesMass(k):
            return self.r1.Y[k] * self.r1.mass + self.r2.Y[k] * self.r2.mass
        mO2 = speciesMass(kO2)
        mAr = speciesMass(kAr)

        t = 0
        while t < 1.0:
            t = self.net.step()
            p1 = self.r1.thermo.P
            p2 = self.r2.thermo.P
            self.assertNear(mdot(p1-p2), valve.mass_flow_rate)
            self.assertArrayNear(Y1, self.r1.Y)
            self.assertNear(speciesMass(kAr), mAr)
            self.assertNear(speciesMass(kO2), mO2)

    def test_valve_timing(self):
        # test timed valve
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.valve_coeff = k
        valve.set_time_function(lambda t: t>.01)

        mdot = lambda: valve.valve_coeff * (self.r1.thermo.P - self.r2.thermo.P)
        self.net.initialize()
        self.assertEqual(valve.mass_flow_rate, 0.0)
        self.net.advance(0.01)
        self.assertEqual(valve.mass_flow_rate, 0.0)
        self.net.advance(0.01 + 1e-9)
        self.assertNear(valve.mass_flow_rate, mdot())
        self.net.advance(0.02)
        self.assertNear(valve.mass_flow_rate, mdot())

    def test_valve_errors(self):
        self.make_reactors()
        res = ct.Reservoir()

        with self.assertRaisesRegex(ct.CanteraError, 'contents not defined'):
            # Must assign contents of both reactors before creating Valve
            v = ct.Valve(self.r1, res)

        v = ct.Valve(self.r1, self.r2)
        with self.assertRaisesRegex(ct.CanteraError, 'Already installed'):
            # inlet and outlet cannot be reassigned
            v._install(self.r2, self.r1)

    def test_pressure_controller1(self):
        self.make_reactors(n_reactors=1)
        g = ct.Solution('h2o2.yaml', transport_model=None)
        g.TPX = 500, 2*101325, 'H2:1.0'
        inlet_reservoir = ct.Reservoir(g)
        g.TP = 300, 101325
        outlet_reservoir = ct.Reservoir(g)

        mfc = ct.MassFlowController(inlet_reservoir, self.r1)
        mdot = lambda t: np.exp(-100*(t-0.5)**2)
        mfc.mass_flow_coeff = 1.
        mfc.set_time_function(mdot)

        pc = ct.PressureController(self.r1, outlet_reservoir)
        pc.set_master(mfc)
        pc.pressure_coeff = 1e-5
        self.assertEqual(pc.pressure_coeff, 1e-5)

        t = 0
        while t < 1.0:
            t = self.net.step()
            self.assertNear(mdot(t), mfc.mass_flow_rate)
            dP = self.r1.thermo.P - outlet_reservoir.thermo.P
            self.assertNear(mdot(t) + 1e-5 * dP, pc.mass_flow_rate)

    def test_pressure_controller2(self):
        self.make_reactors(n_reactors=1)
        g = ct.Solution('h2o2.yaml', transport_model=None)
        g.TPX = 500, 2*101325, 'H2:1.0'
        inlet_reservoir = ct.Reservoir(g)
        g.TP = 300, 101325
        outlet_reservoir = ct.Reservoir(g)

        mfc = ct.MassFlowController(inlet_reservoir, self.r1)
        mdot = lambda t: np.exp(-100*(t-0.5)**2)
        mfc.mass_flow_coeff = 1.
        mfc.set_time_function(mdot)

        pc = ct.PressureController(self.r1, outlet_reservoir)
        pc.set_master(mfc)
        pfunc = lambda dp: 1.e-5 * abs(dp)**.5
        pc.set_pressure_function(pfunc)
        self.assertEqual(pc.pressure_coeff, 1.)

        t = 0
        while t < 1.0:
            t = self.net.step()
            self.assertNear(mdot(t), mfc.mass_flow_rate)
            dP = self.r1.thermo.P - outlet_reservoir.thermo.P
            self.assertNear(mdot(t) + pfunc(dP), pc.mass_flow_rate)

    def test_pressure_controller_errors(self):
        self.make_reactors()
        res = ct.Reservoir(self.gas1)
        mfc = ct.MassFlowController(res, self.r1, mdot=0.6)

        p = ct.PressureController(self.r1, self.r2, master=mfc, K=0.5)

        with self.assertRaisesRegex(ct.CanteraError, 'is not ready'):
            p = ct.PressureController(self.r1, self.r2, K=0.5)
            p.mass_flow_rate

        with self.assertRaisesRegex(ct.CanteraError, 'is not ready'):
            p = ct.PressureController(self.r1, self.r2)
            p.mass_flow_rate

        with self.assertRaisesRegex(ct.CanteraError, 'NotImplementedError'):
            p = ct.PressureController(self.r1, self.r2)
            p.set_time_function(lambda t: t>1.)

    def test_set_initial_time(self):
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.set_pressure_function(mdot)

        t0 = 0.0
        tf = t0 + 0.5
        self.net.advance(tf)
        self.assertNear(self.net.time, tf)
        p1a = self.r1.thermo.P
        p2a = self.r2.thermo.P

        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.set_pressure_function(mdot)

        t0 = 0.2
        self.net.set_initial_time(t0)
        tf = t0 + 0.5
        self.net.advance(tf)
        self.assertNear(self.net.time, tf)
        p1b = self.r1.thermo.P
        p2b = self.r2.thermo.P

        self.assertNear(p1a, p1b)
        self.assertNear(p2a, p2b)

    def test_reinitialize(self):
        self.make_reactors(T1=300, T2=1000, independent=False)
        self.add_wall(U=200, A=1.0)
        self.net.advance(1.0)
        T1a = self.r1.T
        T2a = self.r2.T

        self.r1.thermo.TD = 300, None
        self.r1.syncState()

        self.r2.thermo.TD = 1000, None
        self.r2.syncState()

        self.assertNear(self.r1.T, 300)
        self.assertNear(self.r2.T, 1000)
        self.net.advance(2.0)
        T1b = self.r1.T
        T2b = self.r2.T

        self.assertNear(T1a, T1b)
        self.assertNear(T2a, T2b)

    def test_unpicklable(self):
        self.make_reactors()
        import pickle
        with self.assertRaises(NotImplementedError):
            pickle.dumps(self.r1)
        with self.assertRaises(NotImplementedError):
            pickle.dumps(self.net)

    def test_uncopyable(self):
        self.make_reactors()
        import copy
        with self.assertRaises(NotImplementedError):
            copy.copy(self.r1)
        with self.assertRaises(NotImplementedError):
            copy.copy(self.net)

    def test_invalid_property(self):
        self.make_reactors()
        for x in (self.r1, self.net):
            with self.assertRaises(AttributeError):
                x.foobar = 300
            with self.assertRaises(AttributeError):
                x.foobar

    def test_bad_kwarg(self):
        self.reactorClass(name='ok')
        with self.assertRaises(TypeError):
            r1 = self.reactorClass(foobar=3.14)


class TestIdealGasReactor(TestReactor):
    reactorClass = ct.IdealGasReactor


class TestWellStirredReactorIgnition(utilities.CanteraTest):
    """ Ignition (or not) of a well-stirred reactor """
    def setup(self, T0, P0, mdot_fuel, mdot_ox):
        gas_def = """
        phases:
        - name: gas
          species:
          - gri30.yaml/species: [H2, H, O, O2, OH, H2O, HO2, H2O2, CH2, CH2(S), CH3,
              CH4, CO, CO2, HCO, CH2O, CH2OH, CH3O, CH3OH, C2H4, C2H5, C2H6, N2, AR]
          thermo: ideal-gas
          kinetics: gas
          reactions:
          - gri30.yaml/reactions: declared-species
          skip-undeclared-third-bodies: true
        """

        self.gas = ct.Solution(yaml=gas_def)

        # fuel inlet
        self.gas.TPX = T0, P0, "CH4:1.0"
        self.fuel_in = ct.Reservoir(self.gas)

        # oxidizer inlet
        self.gas.TPX = T0, P0, "N2:3.76, O2:1.0"
        self.oxidizer_in = ct.Reservoir(self.gas)

        # reactor, initially filled with N2
        self.gas.TPX = T0, P0, "N2:1.0"
        self.combustor = ct.IdealGasReactor(self.gas)
        self.combustor.volume = 1.0

        # outlet
        self.exhaust = ct.Reservoir(self.gas)

        # connect the reactor to the reservoirs
        self.fuel_mfc = ct.MassFlowController(self.fuel_in, self.combustor)
        self.fuel_mfc.mass_flow_rate = mdot_fuel
        self.oxidizer_mfc = ct.MassFlowController(self.oxidizer_in, self.combustor)
        self.oxidizer_mfc.mass_flow_rate = mdot_ox
        self.valve = ct.Valve(self.combustor, self.exhaust)
        self.valve.valve_coeff = 1.0

        self.net = ct.ReactorNet()
        self.net.add_reactor(self.combustor)
        self.net.max_err_test_fails = 10

    def integrate(self, tf):
        t = 0.0
        times = []
        T = []
        i = 0
        while t < tf:
            i += 1
            t = self.net.step()
            times.append(t)
            T.append(self.combustor.T)
        return times, T

    def test_nonreacting(self):
        mdot_f = 1.0
        mdot_o = 5.0
        T0 = 900.0
        self.setup(T0, 10*ct.one_atm, mdot_f, mdot_o)
        self.gas.set_multiplier(0.0)
        t,T = self.integrate(100.0)

        for i in range(len(t)):
            self.assertNear(T[i], T0, rtol=1e-5)

        self.assertNear(self.combustor.thermo['CH4'].Y,
                        mdot_f / (mdot_o + mdot_f))

    def test_ignition1(self):
        self.setup(900.0, 10*ct.one_atm, 1.0, 5.0)
        t,T = self.integrate(10.0)

        self.assertTrue(T[-1] > 1200) # mixture ignited
        for i in range(len(t)):
            if T[i] > 0.5 * (T[0] + T[-1]):
                tIg = t[i]
                break

        # regression test; no external basis for this result
        self.assertNear(tIg, 2.2249, 1e-3)

    def test_ignition2(self):
        self.setup(900.0, 10*ct.one_atm, 1.0, 20.0)
        t,T = self.integrate(10.0)

        self.assertTrue(T[-1] > 1200) # mixture ignited
        for i in range(len(t)):
            if T[i] > 0.5 * (T[0] + T[-1]):
                tIg = t[i]
                break

        # regression test; no external basis for this result
        self.assertNear(tIg, 1.4856, 1e-3)

    def test_ignition3(self):
        self.setup(900.0, 10*ct.one_atm, 1.0, 80.0)
        self.net.max_time_step = 0.5
        t,T = self.integrate(100.0)
        self.assertTrue(T[-1] < 910) # mixture did not ignite

    def test_steady_state(self):
        self.setup(900.0, 10*ct.one_atm, 1.0, 20.0)
        residuals = self.net.advance_to_steady_state(return_residuals=True)
        # test if steady state is reached
        self.assertTrue(residuals[-1] < 10. * self.net.rtol)
        # regression test; no external basis for these results
        self.assertNear(self.combustor.T, 2498.94, 1e-5)
        self.assertNear(self.combustor.thermo['H2O'].Y[0], 0.103658, 1e-5)
        self.assertNear(self.combustor.thermo['HO2'].Y[0], 8.734515e-06, 1e-5)


class TestConstPressureReactor(utilities.CanteraTest):
    """
    The constant pressure reactor should give essentially the same results as
    as a regular "Reactor" with a wall with a very high expansion rate
    coefficient.
    """

    reactorClass = ct.ConstPressureReactor

    def create_reactors(self, add_Q=False, add_mdot=False, add_surf=False):
        gas_def = """
        phases:
        - name: gas
          species:
          - gri30.yaml/species: [H2, H, O, O2, OH, H2O, HO2, H2O2, CH3, CH4, CO, CO2,
              HCO, CH2O, CH3O, CH3OH, N2, AR]
          thermo: ideal-gas
          kinetics: gas
          reactions:
          - gri30.yaml/reactions: declared-species
          skip-undeclared-third-bodies: true
        """
        self.gas = ct.Solution(yaml=gas_def)
        self.gas.TPX = 900, 25*ct.one_atm, 'CO:0.5, H2O:0.2'

        self.gas1 = ct.Solution(yaml=gas_def)
        self.gas2 = ct.Solution(yaml=gas_def)
        resGas = ct.Solution(yaml=gas_def)
        solid = ct.Solution('diamond.yaml', 'diamond')

        T0 = 1200
        P0 = 25*ct.one_atm
        X0 = 'CH4:0.5, H2O:0.2, CO:0.3'

        self.gas1.TPX = T0, P0, X0
        self.gas2.TPX = T0, P0, X0

        self.r1 = ct.IdealGasReactor(self.gas1)
        self.r2 = self.reactorClass(self.gas2)

        self.r1.volume = 0.2
        self.r2.volume = 0.2

        resGas.TP = T0 - 300, P0
        env = ct.Reservoir(resGas)

        U = 300 if add_Q else 0

        self.w1 = ct.Wall(self.r1, env, K=1e3, A=0.1, U=U)
        self.w2 = ct.Wall(self.r2, env, A=0.1, U=U)

        if add_mdot:
            mfc1 = ct.MassFlowController(env, self.r1, mdot=0.05)
            mfc2 = ct.MassFlowController(env, self.r2, mdot=0.05)

        if add_surf:
            self.interface1 = ct.Interface('diamond.yaml', 'diamond_100',
                                      (self.gas1, solid))
            self.interface2 = ct.Interface('diamond.yaml', 'diamond_100',
                                      (self.gas2, solid))

            C = np.zeros(self.interface1.n_species)
            C[0] = 0.3
            C[4] = 0.7
            self.surf1 = ct.ReactorSurface(self.interface1, A=0.2)
            self.surf2 = ct.ReactorSurface(self.interface2, A=0.2)
            self.surf1.coverages = C
            self.surf2.coverages = C
            self.surf1.install(self.r1)
            self.surf2.install(self.r2)

        self.net1 = ct.ReactorNet([self.r1])
        self.net2 = ct.ReactorNet([self.r2])
        self.net1.max_time_step = 0.05
        self.net2.max_time_step = 0.05
        self.net2.max_err_test_fails = 10

    def test_component_index(self):
        self.create_reactors(add_surf=True)
        for (gas,net,iface,r) in ((self.gas1, self.net1, self.interface1, self.r1),
                                  (self.gas2, self.net2, self.interface2, self.r2)):
            net.step()

            N0 = net.n_vars - gas.n_species - iface.n_species
            N1 = net.n_vars - iface.n_species
            for i, name in enumerate(gas.species_names):
                self.assertEqual(i + N0, r.component_index(name))
            for i, name in enumerate(iface.species_names):
                self.assertEqual(i + N1, r.component_index(name))

    def test_component_names(self):
        self.create_reactors(add_surf=True)
        for i in range(self.net1.n_vars):
            self.assertEqual(self.r1.component_index(self.r1.component_name(i)), i)
            self.assertEqual(self.net1.component_name(i),
                '{}: {}'.format(self.r1.name, self.r1.component_name(i)))

    def integrate(self, surf=False):
        for t in np.arange(0.5, 50, 1.0):
            self.net1.advance(t)
            self.net2.advance(t)
            self.assertArrayNear(self.r1.thermo.Y, self.r2.thermo.Y,
                                 rtol=5e-4, atol=1e-6)
            self.assertNear(self.r1.T, self.r2.T, rtol=5e-5)
            self.assertNear(self.r1.thermo.P, self.r2.thermo.P, rtol=1e-6)
            if surf:
                self.assertArrayNear(self.surf1.coverages, self.surf2.coverages,
                                     rtol=1e-4, atol=1e-8)

    def test_closed(self):
        self.create_reactors()
        self.integrate()

    def test_with_heat_transfer(self):
        self.create_reactors(add_Q=True)
        self.integrate()

    def test_with_mdot(self):
        self.create_reactors(add_mdot=True)
        self.integrate()

    def test_with_surface_reactions(self):
        self.create_reactors(add_surf=True)
        self.net1.atol = self.net2.atol = 1e-18
        self.net1.rtol = self.net2.rtol = 1e-9
        self.integrate(surf=True)


class TestIdealGasConstPressureReactor(TestConstPressureReactor):
    reactorClass = ct.IdealGasConstPressureReactor


class TestIdealGasConstPressureMoleReactor(TestIdealGasConstPressureReactor):
    reactorClass = ct.IdealGasConstPressureMoleReactor

    def create_reactors(self, **kwargs):
        super().create_reactors(**kwargs)
        if "add_surf" not in kwargs.keys():
            self.net2.preconditioner = ct.AdaptivePreconditioner()
            self.net2.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}

    def test_get_solver_type(self):
        self.create_reactors()
        self.assertEqual(self.net2.linear_solver_type, "GMRES")

class TestIdealGasMoleReactor(TestReactor):
    reactorClass = ct.IdealGasMoleReactor

    def test_adaptive_precon_integration(self):
        # Network one with non-mole reactor
        net1 = ct.ReactorNet()
        gas1 = ct.Solution('h2o2.yaml', transport_model=None)
        gas1.TP = 300, ct.one_atm
        gas1.set_equivalence_ratio(1, "H2", "O2:1, N2:3.76")
        r1 = ct.IdealGasReactor(gas1)
        net1.add_reactor(r1)
        # Network two with mole reactor and preconditioner
        net2 = ct.ReactorNet()
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TP = 300, ct.one_atm
        gas2.set_equivalence_ratio(1, "H2", "O2:1, N2:3.76")
        r2 = ct.IdealGasMoleReactor(gas2)
        net2.add_reactor(r2)
        # add preconditioner
        net2.preconditioner = ct.AdaptivePreconditioner()
        net2.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
        # integrate
        for i in range(1, 11, 1):
            adv_time = i * 0.1
            net1.advance(adv_time)
            net2.advance(adv_time)
            self.assertNear(r1.T, r2.T)
            self.assertNear(r1.thermo.density, r2.thermo.density)
            self.assertNear(r1.thermo.P, r1.thermo.P)
            self.assertArrayNear(r1.thermo.X, r1.thermo.X)


class TestFlowReactor(utilities.CanteraTest):
    gas_def = """
    phases:
    - name: gas
      species:
      - gri30.yaml/species: [H2, H, O, O2, OH, H2O, HO2, H2O2, CH3, CH4, CO, CO2,
          HCO, CH2O, CH3O, CH3OH, AR, N2]
      thermo: ideal-gas
      kinetics: gas
      reactions:
      - gri30.yaml/reactions: declared-species
      skip-undeclared-third-bodies: true
    """

    def test_nonreacting(self):
        g = ct.Solution(yaml=self.gas_def)
        g.TPX = 300, 101325, 'O2:1.0'
        r = ct.FlowReactor(g)
        r.mass_flow_rate = 10

        net = ct.ReactorNet()
        net.add_reactor(r)

        t = 0
        v0 = r.speed
        self.assertNear(v0, 10 / r.density)
        while t < 10.0:
            t = net.step()

            self.assertNear(v0, r.speed)
            self.assertNear(r.distance, v0 * t)

    def test_reacting(self):
        g = ct.Solution(yaml=self.gas_def)
        g.TPX = 1400, 20*101325, 'CO:1.0, H2O:1.0'

        r = ct.FlowReactor(g)
        r.mass_flow_rate = 10

        net = ct.ReactorNet()
        net.add_reactor(r)
        net.atol = 1e-18
        net.rtol = 1e-9
        net.max_err_test_fails = 10

        t = 0
        self.assertNear(r.speed, 10 / r.density)
        while t < 1.0:
            t1 = net.time
            x1 = r.distance

            t = net.step()

            v = (r.distance - x1) / (net.time - t1)
            self.assertNear(r.speed, v, 1e-3)


class TestSurfaceKinetics(utilities.CanteraTest):
    def make_reactors(self):
        self.net = ct.ReactorNet()

        self.gas = ct.Solution('diamond.yaml', 'gas')
        self.solid = ct.Solution('diamond.yaml', 'diamond')
        self.interface = ct.Interface('diamond.yaml', 'diamond_100',
                                      (self.gas, self.solid))
        self.gas.TPX = None, 1.0e3, 'H:0.002, H2:1, CH4:0.01, CH3:0.0002'
        self.r1 = ct.IdealGasReactor(self.gas)
        self.r1.volume = 0.01
        self.net.add_reactor(self.r1)

        self.r2 = ct.IdealGasReactor(self.gas)
        self.r2.volume = 0.01
        self.net.add_reactor(self.r2)

    def test_coverages(self):
        self.make_reactors()
        surf1 = ct.ReactorSurface(self.interface, self.r1)

        surf1.coverages = {'c6HH':0.3, 'c6HM':0.7}
        self.assertNear(surf1.coverages[0], 0.3)
        self.assertNear(surf1.coverages[1], 0.0)
        self.assertNear(surf1.coverages[4], 0.7)
        self.net.advance(1e-5)
        C_left = surf1.coverages

        self.make_reactors()
        surf2 = ct.ReactorSurface(self.interface, self.r2)
        surf2.coverages = 'c6HH:0.3, c6HM:0.7'
        self.assertNear(surf2.coverages[0], 0.3)
        self.assertNear(surf2.coverages[4], 0.7)
        self.net.advance(1e-5)
        C_right = surf2.coverages

        self.assertNear(sum(C_left), 1.0)
        self.assertArrayNear(C_left, C_right)

    def test_coverages_regression1(self):
        # Test with energy equation disabled
        self.make_reactors()
        self.r1.energy_enabled = False
        self.r2.energy_enabled = False
        surf1 = ct.ReactorSurface(self.interface, self.r1)

        C = np.zeros(self.interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        surf1.coverages = C
        self.assertArrayNear(surf1.coverages, C)
        data = []
        test_file = self.test_work_path / "test_coverages_regression1.csv"
        reference_file = self.test_data_path / "WallKinetics-coverages-regression1.csv"
        data = []
        for t in np.linspace(1e-6, 1e-3):
            self.net.advance(t)
            data.append([t, self.r1.T, self.r1.thermo.P, self.r1.mass] +
                        list(self.r1.thermo.X) + list(surf1.coverages))
        np.savetxt(test_file, data, delimiter=',')

        bad = utilities.compareProfiles(reference_file, test_file,
                                        rtol=1e-5, atol=1e-9, xtol=1e-12)
        self.assertFalse(bool(bad), bad)

    def test_coverages_regression2(self):
        # Test with energy equation enabled
        self.make_reactors()
        surf = ct.ReactorSurface(self.interface, self.r1)

        C = np.zeros(self.interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        surf.coverages = C
        self.assertArrayNear(surf.coverages, C)
        data = []
        test_file = self.test_work_path / "test_coverages_regression2.csv"
        reference_file = self.test_data_path / "WallKinetics-coverages-regression2.csv"
        data = []
        for t in np.linspace(1e-6, 1e-3):
            self.net.advance(t)
            data.append([t, self.r1.T, self.r1.thermo.P, self.r1.mass] +
                        list(self.r1.thermo.X) + list(surf.coverages))
        np.savetxt(test_file, data, delimiter=',')

        bad = utilities.compareProfiles(reference_file, test_file,
                                        rtol=1e-5, atol=1e-9, xtol=1e-12)
        self.assertFalse(bool(bad), bad)


class TestReactorSensitivities(utilities.CanteraTest):
    def test_sensitivities1(self):
        net = ct.ReactorNet()
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 1300, 20*101325, 'CO:1.0, H2:0.1, CH4:0.1, H2O:0.5'
        r1 = ct.IdealGasReactor(gas)
        net.add_reactor(r1)

        self.assertEqual(net.n_sensitivity_params, 0)
        r1.add_sensitivity_reaction(40)
        r1.add_sensitivity_reaction(41)

        net.advance(0.1)

        self.assertEqual(net.n_sensitivity_params, 2)
        self.assertEqual(net.n_vars,
                         gas.n_species + r1.component_index(gas.species_name(0)))
        S = net.sensitivities()
        self.assertEqual(S.shape, (net.n_vars, net.n_sensitivity_params))

    def test_sensitivities2(self):
        net = ct.ReactorNet()

        gas1 = ct.Solution("diamond.yaml", "gas")
        solid = ct.Solution("diamond.yaml", "diamond")
        interface = ct.Interface("diamond.yaml", "diamond_100", (gas1, solid))
        r1 = ct.IdealGasReactor(gas1)
        net.add_reactor(r1)
        net.atol_sensitivity = 1e-10
        net.rtol_sensitivity = 1e-8

        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'
        r2 = ct.IdealGasReactor(gas2)
        net.add_reactor(r2)

        surf = ct.ReactorSurface(interface, r1, A=1.5)

        C = np.zeros(interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        surf.coverages = C
        surf.add_sensitivity_reaction(2)
        r2.add_sensitivity_reaction(18)

        for T in (901, 905, 910, 950, 1500):
            while r2.T < T:
                net.step()

            S = net.sensitivities()

            # number of non-species variables in each reactor
            Ns = r1.component_index(gas1.species_name(0))

            # Index of first variable corresponding to r2
            K2 = Ns + gas1.n_species + interface.n_species

            # Constant volume should generate zero sensitivity coefficient
            self.assertArrayNear(S[1,:], np.zeros(2))
            self.assertArrayNear(S[K2+1,:], np.zeros(2))

            # Sensitivity coefficients for the disjoint reactors should be zero
            self.assertNear(np.linalg.norm(S[Ns:K2,1]), 0.0, atol=1e-5)
            self.assertNear(np.linalg.norm(S[K2+Ns:,0]), 0.0, atol=1e-5)

    def _test_parameter_order1(self, reactorClass):
        # Single reactor, changing the order in which parameters are added
        gas = ct.Solution('h2o2.yaml', transport_model=None)

        def setup(params):
            net = ct.ReactorNet()
            gas.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'

            r = reactorClass(gas)
            net.add_reactor(r)

            for kind, p in params:
                if kind == 'r':
                    r.add_sensitivity_reaction(p)
                elif kind == 's':
                    r.add_sensitivity_species_enthalpy(p)
            return r, net

        def integrate(r, net):
            while r.T < 910:
                net.step()
            return net.sensitivities()

        def check_names(reactor, net, params):
            for i,(kind,p) in enumerate(params):
                rname, comp = net.sensitivity_parameter_name(i).split(': ')
                self.assertEqual(reactor.name, rname)
                if kind == 'r':
                    self.assertEqual(gas.reaction(p).equation, comp)
                elif kind == 's':
                    self.assertEqual(p + ' enthalpy', comp)

        params1 = [('r', 2), ('r', 10), ('r', 18), ('r', 19), ('s', 'O2'),
                   ('s', 'OH'), ('s', 'H2O2')]
        r1,net1 = setup(params1)
        S1 = integrate(r1, net1)
        check_names(r1, net1, params1)

        params2 = [('r', 19), ('s', 'H2O2'), ('s', 'OH'), ('r', 10),
                   ('s', 'O2'), ('r', 2), ('r', 18)]
        r2,net2 = setup(params2)
        S2 = integrate(r2, net2)
        check_names(r2, net2, params2)

        for i,j in enumerate((5,3,6,0,4,2,1)):
            self.assertArrayNear(S1[:,i], S2[:,j])

    def test_parameter_order1a(self):
        self._test_parameter_order1(ct.IdealGasReactor)

    @utilities.slow_test
    def test_parameter_order1b(self):
        self._test_parameter_order1(ct.IdealGasConstPressureReactor)

    @utilities.slow_test
    def test_parameter_order2(self):
        # Multiple reactors, changing the order in which parameters are added
        gas = ct.Solution('h2o2.yaml', transport_model=None)

        def setup(reverse=False):
            net = ct.ReactorNet()
            gas1 = ct.Solution('h2o2.yaml', transport_model=None)
            gas1.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'
            rA = ct.IdealGasReactor(gas1)

            gas2 = ct.Solution('h2o2.yaml', transport_model=None)
            gas2.TPX = 920, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:0.5'
            rB = ct.IdealGasReactor(gas2)
            if reverse:
                net.add_reactor(rB)
                net.add_reactor(rA)
            else:
                net.add_reactor(rA)
                net.add_reactor(rB)

            return rA, rB, net

        def integrate(r, net):
            net.advance(1e-4)
            return net.sensitivities()

        S = []

        for reverse in (True,False):
            rA1,rB1,net1 = setup(reverse)
            params1 = [(rA1,2),(rA1,19),(rB1,10),(rB1,18)]
            for r,p in params1:
                r.add_sensitivity_reaction(p)
            S.append(integrate(rA1, net1))

            pname = lambda r,i: '%s: %s' % (r.name, gas.reaction(i).equation)
            for i,(r,p) in enumerate(params1):
                self.assertEqual(pname(r,p), net1.sensitivity_parameter_name(i))

            rA2,rB2,net2 = setup(reverse)
            params2 = [(rB2,10),(rA2,19),(rB2,18),(rA2,2)]
            for r,p in params2:
                r.add_sensitivity_reaction(p)
            S.append(integrate(rA2, net2))

            for i,(r,p) in enumerate(params2):
                self.assertEqual(pname(r,p), net2.sensitivity_parameter_name(i))

        # Check that the results reflect the changed parameter ordering
        for a,b in ((0,1), (2,3)):
            for i,j in enumerate((3,1,0,2)):
                self.assertArrayNear(S[a][:,i], S[b][:,j])

        # Check that results are consistent after changing the order that
        # reactors are added to the network
        N = gas.n_species + r.component_index(gas.species_name(0))
        self.assertArrayNear(S[0][:N], S[2][N:], 1e-5, 1e-5)
        self.assertArrayNear(S[0][N:], S[2][:N], 1e-5, 1e-5)
        self.assertArrayNear(S[1][:N], S[3][N:], 1e-5, 1e-5)
        self.assertArrayNear(S[1][N:], S[3][:N], 1e-5, 1e-5)

    @utilities.slow_test
    def test_parameter_order3(self):
        # Test including reacting surfaces
        gas1 = ct.Solution("diamond.yaml", "gas")
        solid = ct.Solution("diamond.yaml", "diamond")
        interface = ct.Interface("diamond.yaml", "diamond_100", (gas1, solid))

        gas2 = ct.Solution('h2o2.yaml', transport_model=None)

        def setup(order):
            gas1.TPX = 1200, 1e3, 'H:0.002, H2:1, CH4:0.01, CH3:0.0002'
            gas2.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'
            net = ct.ReactorNet()
            rA = ct.IdealGasReactor(gas1)
            rB = ct.IdealGasReactor(gas2)

            if order % 2 == 0:
                surfX = ct.ReactorSurface(interface, rA, A=0.1)
                surfY = ct.ReactorSurface(interface, rA, A=10)
            else:
                surfY = ct.ReactorSurface(interface, rA, A=10)
                surfX = ct.ReactorSurface(interface, rA, A=0.1)

            C1 = np.zeros(interface.n_species)
            C2 = np.zeros(interface.n_species)
            C1[0] = 0.3
            C1[4] = 0.7

            C2[0] = 0.9
            C2[4] = 0.1
            surfX.coverages = C1
            surfY.coverages = C2

            if order // 2 == 0:
                net.add_reactor(rA)
                net.add_reactor(rB)
            else:
                net.add_reactor(rB)
                net.add_reactor(rA)

            return rA, rB, surfX, surfY, net

        def integrate(r, net):
            net.advance(1e-4)
            return net.sensitivities()

        S = []

        for order in range(4):
            rA, rB, surfX, surfY, net = setup(order)
            for (obj,k) in [(rB,2), (rB,18), (surfX,2),
                            (surfY,0), (surfY,2)]:
                obj.add_sensitivity_reaction(k)
            integrate(rB, net)
            S.append(net.sensitivities())

            rA, rB, surfX, surfY, net = setup(order)
            for (obj,k) in [(surfY,2), (surfX,2), (rB,18),
                            (surfX,0), (rB,2)]:
                obj.add_sensitivity_reaction(k)

            integrate(rB, net)
            S.append(net.sensitivities())

        for a,b in [(0,1),(2,3),(4,5),(6,7)]:
            for i,j in enumerate((4,2,1,3,0)):
                self.assertArrayNear(S[a][:,i], S[b][:,j], 1e-2, 1e-3)

    def setup_ignition_delay(self):
        gas = ct.Solution('h2o2.yaml', transport_model=None)
        gas.TP = 900, 5*ct.one_atm
        gas.set_equivalence_ratio(0.4, 'H2', 'O2:1.0, AR:4.0')
        r = ct.IdealGasReactor(gas)
        net = ct.ReactorNet([r])
        net.rtol_sensitivity = 2e-5
        return gas, r, net

    def calc_tig(self, species, dH):
        gas, r, net = self.setup_ignition_delay()

        S = gas.species(species)
        st = S.thermo
        coeffs = st.coeffs
        coeffs[[6, 13]] += dH / ct.gas_constant
        snew = ct.NasaPoly2(st.min_temp, st.max_temp, st.reference_pressure, coeffs)
        S.thermo = snew
        gas.modify_species(gas.species_index(species), S)
        t = []
        T = []
        while net.time < 0.6:
            t.append(net.time)
            T.append(r.thermo.T)
            net.step()
        T = np.array(T)
        t = np.array(t)
        To = T[0]
        Tf = T[-1]

        return (t[-1]*T[-1] - np.trapz(T,t)) / (T[-1] - T[0])

    def calc_dtdh(self, species):
        gas, r, net = self.setup_ignition_delay()
        for s in species:
            r.add_sensitivity_species_enthalpy(s)

        t = [0.0]
        T = [r.T]
        S = [[0.0]*len(species)]
        iTemp = r.component_index('temperature')
        while net.time < 0.6:
            net.step()
            t.append(net.time)
            T.append(r.thermo.T)
            S.append(net.sensitivities()[iTemp])

        T = np.array(T)
        t = np.array(t)
        S = np.array(S)

        To = T[0]
        Tf = T[-1]
        tig = (t[-1]*Tf - np.trapz(T,t))/(Tf-To)
        dtdp = ((t[-1] - tig)*S[-1,:]*Tf - np.trapz(S*T[:,None], t, axis=0))/(Tf-To)
        return dtdp

    # See https://github.com/Cantera/enhancements/issues/55
    @unittest.skip("Integration of sensitivity ODEs is unreliable")
    def test_ignition_delay_sensitivity(self):
        species = ('H2', 'H', 'O2', 'H2O2', 'H2O', 'OH', 'HO2')
        dtigdh_cvodes = self.calc_dtdh(species)
        tig0 = self.calc_tig('H2', 0)
        dH = 1e4
        for i,s in enumerate(species):
            dtigdh = (self.calc_tig(s, dH) - tig0) / dH
            self.assertNear(dtigdh_cvodes[i], dtigdh, atol=1e-14, rtol=5e-2)


class CombustorTestImplementation:
    """
    These tests are based on the sample:

        interfaces/cython/cantera/examples/reactors/combustor.py

    with some simplifications so that they run faster and produce more
    consistent output.
    """

    def setUp(self):
        self.referenceFile = utilities.TEST_DATA_PATH / "CombustorTest-integrateWithAdvance.csv"
        self.gas = ct.Solution('h2o2.yaml', transport_model=None)

        # create a reservoir for the fuel inlet, and set to pure methane.
        self.gas.TPX = 300.0, ct.one_atm, 'H2:1.0'
        fuel_in = ct.Reservoir(self.gas)
        fuel_mw = self.gas.mean_molecular_weight

        # Oxidizer inlet
        self.gas.TPX = 300.0, ct.one_atm, 'O2:1.0, AR:3.0'
        oxidizer_in = ct.Reservoir(self.gas)
        oxidizer_mw = self.gas.mean_molecular_weight

        # to ignite the fuel/air mixture, we'll introduce a pulse of radicals.
        # The steady-state behavior is independent of how we do this, so we'll
        # just use a stream of pure atomic hydrogen.
        self.gas.TPX = 300.0, ct.one_atm, 'H:1.0'
        self.igniter = ct.Reservoir(self.gas)

        # create the combustor, and fill it in initially with a diluent
        self.gas.TPX = 300.0, ct.one_atm, 'AR:1.0'
        self.combustor = ct.IdealGasReactor(self.gas)

        # create a reservoir for the exhaust
        self.exhaust = ct.Reservoir(self.gas)

        # compute fuel and air mass flow rates
        factor = 0.1
        oxidizer_mdot = 4 * factor*oxidizer_mw
        fuel_mdot = factor*fuel_mw

        # The igniter will use a time-dependent igniter mass flow rate.
        def igniter_mdot(t, t0=0.1, fwhm=0.05, amplitude=0.1):
            return amplitude * math.exp(-(t-t0)**2 * 4 * math.log(2) / fwhm**2)

        # create and install the mass flow controllers. Controllers
        # m1 and m2 provide constant mass flow rates, and m3 provides
        # a short Gaussian pulse only to ignite the mixture
        self.m1 = ct.MassFlowController(fuel_in, self.combustor, mdot=fuel_mdot)
        self.m2 = ct.MassFlowController(oxidizer_in, self.combustor, mdot=oxidizer_mdot)
        self.m3 = ct.MassFlowController(self.igniter, self.combustor, mdot=igniter_mdot)

        # put a valve on the exhaust line to regulate the pressure
        self.v = ct.Valve(self.combustor, self.exhaust, K=1.0)

        # the simulation only contains one reactor
        self.sim = ct.ReactorNet([self.combustor])

    def test_integrateWithStep(self):
        tnow = 0.0
        tfinal = 0.25
        self.data = []
        while tnow < tfinal:
            tnow = self.sim.step()
            self.data.append([tnow, self.combustor.T] +
                             list(self.combustor.thermo.X))

        self.assertTrue(tnow >= tfinal)
        bad = utilities.compareProfiles(self.referenceFile, self.data,
                                        rtol=1e-3, atol=1e-9)
        self.assertFalse(bad, bad)

    def test_integrateWithAdvance(self, saveReference=False):
        self.data = []
        for t in np.linspace(0, 0.25, 101)[1:]:
            self.sim.advance(t)
            self.data.append([t, self.combustor.T] +
                             list(self.combustor.thermo.X))

        if saveReference:
            np.savetxt(self.referenceFile, np.array(self.data), '%11.6e', ', ')
        else:
            bad = utilities.compareProfiles(self.referenceFile, self.data,
                                            rtol=1e-6, atol=1e-12)
            self.assertFalse(bad, bad)

    def test_invasive_mdot_function(self):
        def igniter_mdot(t, t0=0.1, fwhm=0.05, amplitude=0.1):
            # Querying properties of the igniter changes the state of the
            # underlying ThermoPhase object, but shouldn't affect the
            # integration
            self.igniter.density
            return amplitude * math.exp(-(t-t0)**2 * 4 * math.log(2) / fwhm**2)
        self.m3.mass_flow_rate = igniter_mdot

        self.data = []
        for t in np.linspace(0, 0.25, 101)[1:]:
            self.sim.advance(t)
            self.data.append([t, self.combustor.T] +
                             list(self.combustor.thermo.X))

        bad = utilities.compareProfiles(self.referenceFile, self.data,
                                        rtol=1e-6, atol=1e-12)
        self.assertFalse(bad, bad)

class WallTestImplementation:
    """
    These tests are based on the sample:

        interfaces/cython/cantera/examples/reactors/reactor2.py

    with some simplifications so that they run faster and produce more
    consistent output.
    """

    def setUp(self):
        self.referenceFile = utilities.TEST_DATA_PATH / "WallTest-integrateWithAdvance.csv"
        # reservoir to represent the environment
        self.gas0 = ct.Solution("air.yaml")
        self.gas0.TP = 300, ct.one_atm
        self.env = ct.Reservoir(self.gas0)

        # reactor to represent the side filled with Argon
        self.gas1 = ct.Solution("air.yaml")
        self.gas1.TPX = 1000.0, 30*ct.one_atm, 'AR:1.0'
        self.r1 = ct.Reactor(self.gas1)

        # reactor to represent the combustible mixture
        self.gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        self.gas2.TPX = 500.0, 1.5*ct.one_atm, 'H2:0.5, O2:1.0, AR:10.0'
        self.r2 = ct.Reactor(self.gas2)

        # Wall between the two reactors
        self.w1 = ct.Wall(self.r2, self.r1, A=1.0, K=2e-4, U=400.0)

        # Wall to represent heat loss to the environment
        self.w2 = ct.Wall(self.r2, self.env, A=1.0, U=2000.0)

        # Create the reactor network
        self.sim = ct.ReactorNet([self.r1, self.r2])

    def test_integrateWithStep(self):
        tnow = 0.0
        tfinal = 0.01
        self.data = []
        while tnow < tfinal:
            tnow = self.sim.step()
            self.data.append([tnow,
                              self.r1.T, self.r2.T,
                              self.r1.thermo.P, self.r2.thermo.P,
                              self.r1.volume, self.r2.volume])

        self.assertTrue(tnow >= tfinal)
        bad = utilities.compareProfiles(self.referenceFile, self.data,
                                        rtol=1e-3, atol=1e-8)
        self.assertFalse(bad, bad)

    def test_integrateWithAdvance(self, saveReference=False):
        self.data = []
        for t in np.linspace(0, 0.01, 200)[1:]:
            self.sim.advance(t)
            self.data.append([t,
                              self.r1.T, self.r2.T,
                              self.r1.thermo.P, self.r2.thermo.P,
                              self.r1.volume, self.r2.volume])

        if saveReference:
            np.savetxt(self.referenceFile, np.array(self.data), '%11.6e', ', ')
        else:
            bad = utilities.compareProfiles(self.referenceFile, self.data,
                                            rtol=2e-5, atol=1e-9)
            self.assertFalse(bad, bad)


# Keep the implementations separate from the unittest-derived class
# so that they can be run independently to generate the reference data files.
class CombustorTest(CombustorTestImplementation, unittest.TestCase): pass
class WallTest(WallTestImplementation, unittest.TestCase): pass


class PureFluidReactorTest(utilities.CanteraTest):
    def test_Reactor(self):
        phase = ct.PureFluid("liquidvapor.yaml", "oxygen")
        air = ct.Solution("air.yaml")

        phase.TP = 93, 4e5
        r1 = ct.Reactor(phase)
        r1.volume = 0.1

        air.TP = 300, 4e5
        r2 = ct.Reactor(air)
        r2.volume = 10.0

        air.TP = 500, 4e5
        env = ct.Reservoir(air)

        w1 = ct.Wall(r1,r2)
        w1.expansion_rate_coeff = 1e-3
        w2 = ct.Wall(env,r1, Q=500000, A=1)
        net = ct.ReactorNet([r1,r2])
        net.atol = 1e-10
        net.rtol = 1e-6

        states = ct.SolutionArray(phase, extra='t')
        for t in np.arange(0.0, 60.0, 1):
            net.advance(t)
            states.append(TD=r1.thermo.TD, t=net.time)

        self.assertEqual(states.Q[0], 0)
        self.assertEqual(states.Q[-1], 1)
        self.assertNear(states.Q[30], 0.54806, 1e-4)

    def test_Reactor_2(self):
        phase = ct.PureFluid("liquidvapor.yaml", "carbon-dioxide")
        air = ct.Solution("air.yaml")

        phase.TP = 218, 5e6
        r1 = ct.Reactor(phase)
        r1.volume = 0.1

        air.TP = 500, 5e6
        r2 = ct.Reactor(air)
        r2.volume = 10.0

        w1 = ct.Wall(r1, r2, U=10000, A=1)
        w1.expansion_rate_coeff = 1e-3
        net = ct.ReactorNet([r1,r2])

        states = ct.SolutionArray(phase, extra='t')
        for t in np.arange(0.0, 60.0, 1):
            net.advance(t)
            states.append(TD=r1.thermo.TD, t=net.time)

        self.assertEqual(states.Q[0], 0)
        self.assertEqual(states.Q[-1], 1)
        self.assertNear(states.Q[20], 0.644865, 1e-4)


    def test_ConstPressureReactor(self):
        phase = ct.Nitrogen()
        air = ct.Solution("air.yaml")

        phase.TP = 75, 4e5
        r1 = ct.ConstPressureReactor(phase)
        r1.volume = 0.1

        air.TP = 500, 4e5
        env = ct.Reservoir(air)

        w2 = ct.Wall(env,r1, Q=250000, A=1)
        net = ct.ReactorNet([r1])

        states = ct.SolutionArray(phase, extra='t')
        for t in np.arange(0.0, 100.0, 10):
            net.advance(t)
            states.append(TD=r1.thermo.TD, t=t)

        self.assertEqual(states.Q[1], 0)
        self.assertEqual(states.Q[-2], 1)
        for i in range(3,7):
            self.assertNear(states.T[i], states.T[2])


class AdvanceCoveragesTest(utilities.CanteraTest):
    def setup(self, model="ptcombust.yaml", gas_phase="gas", interface_phase="Pt_surf"):
        # create gas and interface
        self.gas = ct.Solution(model, gas_phase)
        self.surf = ct.Interface(model, interface_phase, [self.gas])

    def test_advance_coverages_parameters(self):
        # create gas and interface
        self.setup()

        # first, test max step size & max steps
        dt = 1.0
        max_steps = 10
        max_step_size = dt / (max_steps + 1)
        # this should throw an error, as we can't reach dt
        with self.assertRaises(ct.CanteraError):
            self.surf.advance_coverages(
                dt=dt, max_step_size=max_step_size, max_steps=max_steps)

        # next, run with different tolerances
        self.setup()
        self.surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        self.gas.TP = self.surf.TP

        self.surf.advance_coverages(dt=dt, rtol=1e-5, atol=1e-12)
        cov = self.surf.coverages[:]

        self.surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        self.gas.TP = self.surf.TP
        self.surf.advance_coverages(dt=dt, rtol=1e-7, atol=1e-14)

        # check that the solutions are similar, but not identical
        self.assertArrayNear(cov, self.surf.coverages)
        self.assertTrue(any(cov != self.surf.coverages))


class ExtensibleReactorTest(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution("h2o2.yaml")

    def test_extra_variable(self):
        class InertialWallReactor(ct.ExtensibleIdealGasReactor):
            def __init__(self, *args, neighbor, **kwargs):
                super().__init__(*args, **kwargs)
                self.v_wall = 0
                self.k_wall = 1e-5
                self.neighbor = neighbor

            def after_initialize(self, t0):
                self.n_vars += 1
                self.i_wall = self.n_vars - 1

            def after_get_state(self, y):
                y[self.i_wall] = self.v_wall

            def after_update_state(self, y):
                self.v_wall = y[self.i_wall]
                self.walls[0].set_velocity(self.v_wall)

            def after_eval(self, t, LHS, RHS):
                # Extra equation is d(v_wall)/dt = k * delta P
                a = self.k_wall * (self.thermo.P - self.neighbor.thermo.P)
                RHS[self.i_wall] = a

            def before_component_index(self, name):
                if name == 'v_wall':
                    return self.i_wall

            def before_component_name(self, i):
                if i == self.i_wall:
                    return 'v_wall'

        self.gas.TP = 300, ct.one_atm
        res = ct.Reservoir(self.gas)
        self.gas.TP = 300, 2 * ct.one_atm
        r = InertialWallReactor(self.gas, neighbor=res)
        w = ct.Wall(r, res)
        net = ct.ReactorNet([r])

        V = []
        for i in range(20):
            net.advance(0.05 * i)
            V.append(r.volume)

        # Wall is accelerating
        self.assertTrue((np.diff(V, 2) > 0).all())

        self.assertIn('v_wall', net.component_name(self.gas.n_species + 3))
        self.assertEqual(r.component_index('volume'), 1)
        self.assertEqual(r.component_name(self.gas.n_species + 3), 'v_wall')
        self.assertEqual(r.component_name(2), 'temperature')

    def test_replace_equations(self):
        nsp = self.gas.n_species
        tau = np.linspace(0.5, 2, nsp + 3)
        class DummyReactor(ct.ExtensibleReactor):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.y = np.ones(nsp + 3)

            def replace_get_state(self, y):
                y[:] = self.y

            def replace_update_state(self, y):
                self.y[:] = y

            def replace_eval(self, t, LHS, RHS):
                RHS[:] = - self.y / tau

        r = DummyReactor(self.gas)
        net = ct.ReactorNet([r])
        net.rtol *= 0.1

        while(net.time < 1):
            net.step()
            self.assertArrayNear(r.get_state(), np.exp(- net.time / tau))

    def test_error_handling(self):
        class DummyReactor1(ct.ExtensibleReactor):
            def replace_eval(self, t): # wrong number of arguments
                pass

        with self.assertRaisesRegex(ValueError, "right number of arguments"):
            DummyReactor1(self.gas)

        class DummyReactor2(ct.ExtensibleReactor):
            def replace_component_index(self, name):
                if name == "succeed":
                    return 0
                elif name == "wrong-type":
                    return "spam"
                # Otherwise, does not return a value

        r2 = DummyReactor2(self.gas)
        self.assertEqual(r2.component_index("succeed"), 0)
        with self.assertRaises(TypeError):
            r2.component_index("wrong-type")
        # Error information should have been reset
        self.assertEqual(r2.component_index("succeed"), 0)
        with self.assertRaisesRegex(ct.CanteraError, "did not return a value"):
            r2.component_index("H2")
        self.assertEqual(r2.component_index("succeed"), 0)

    def test_delegate_throws(self):
        class TestException(Exception):
            pass

        class DummyReactor(ct.ExtensibleConstPressureReactor):
            def before_eval(self, t, LHS, RHS):
                if t > 0.1:
                    raise TestException()

            def before_component_index(self, name):
                if name == "fail":
                    raise TestException()

        r = DummyReactor(self.gas)
        net = ct.ReactorNet([r])

        # Because the TestException is raised inside code called by CVODES, the actual
        # error raised will be a CanteraError
        with self.assertRaises(ct.CanteraError):
            net.advance(0.2)

        self.assertEqual(r.component_index("enthalpy"), 1)
        with self.assertRaises(TestException):
            r.component_index("fail")

    def test_misc(self):
        class DummyReactor(ct.ExtensibleReactor):
            def __init__(self, gas):
                super().__init__(gas)
                self.sync_calls = 0

            def after_species_index(self, name):
                # This will cause returned species indices to be higher by 5 than they
                # would be otherwise
                return 5

            def before_sync_state(self):
                self.sync_calls += 1

        r = DummyReactor(self.gas)
        net = ct.ReactorNet([r])
        self.assertEqual(r.component_index("H2"), 5 + 3 + self.gas.species_index("H2"))
        r.syncState()
        net.advance(1)
        r.syncState()
        self.assertEqual(r.sync_calls, 2)

    def test_RHS_LHS(self):
        # set initial state
        self.gas.TPX = 500, ct.one_atm, 'H2:2,O2:1,N2:4'
        gas_initial_enthalpy = self.gas.enthalpy_mass

        # define properties of gas and solid
        mass_gas = 20  # [kg]
        Q = 100  # [J/s]
        mass_lump = 10  # [kg]
        cp_lump = 1.0  # [J/kg/K]

        # initialize time at zero
        time = 0  # [s]
        n_steps = 300

        # define a class representing reactor with a solid mass and gas inside of it
        class DummyReactor(ct.ExtensibleIdealGasConstPressureReactor):
            # modify energy equation to include solid mass in reactor
            def after_eval(self,t,LHS,RHS):
                self.m_mass = mass_gas
                LHS[1] = mass_lump * cp_lump + self.m_mass * self.thermo.cp_mass
                RHS[1] = Q

        r1 = DummyReactor(self.gas)
        r1_net = ct.ReactorNet([r1])

        for n in range(n_steps):
            time += 4.e-4
            r1_net.advance(time)

        # compare heat added (add_heat) to the equivalent energy contained by the solid
        # and gaseous mass in the reactor
        r1_heat = (mass_lump * cp_lump * (r1.thermo.T - 500) +
                   mass_gas * (self.gas.enthalpy_mass - gas_initial_enthalpy))
        add_heat = Q * time
        self.assertNear(add_heat, r1_heat, atol=1e-5)

    def test_heat_addition(self):
        # Applying heat via 'qdot' property should be equivalent to adding it via a wall
        Qext = 100
        Qwall = -66
        class HeatedReactor(ct.ExtensibleReactor):
            def after_eval_walls(self, y):
                self.qdot += Qext

        self.gas.TPX = 300, ct.one_atm, "N2:1.0"
        r1 = HeatedReactor(self.gas)
        res = ct.Reservoir(self.gas)
        wall = ct.Wall(res, r1, Q=Qwall, A=1)
        net = ct.ReactorNet([r1])
        U0 = r1.thermo.int_energy_mass * r1.mass
        for t in np.linspace(0.1, 5, 10):
            net.advance(t)
            U = r1.thermo.int_energy_mass * r1.mass
            self.assertNear(U - U0, (Qext + Qwall) * t)
            self.assertNear(r1.qdot, Qext + Qwall)

    def test_with_surface(self):
        phase_defs = """
            units: {length: cm, quantity: mol, activation-energy: J/mol}
            phases:
            - name: gas
              thermo: ideal-gas
              species:
              - gri30.yaml/species: [H2, H2O, H2O2, O2]
              kinetics: gas
              reactions:
              - gri30.yaml/reactions: declared-species
              skip-undeclared-third-bodies: true
            - name: Pt_surf
              thermo: ideal-surface
              species:
              - ptcombust.yaml/species: [PT(S), H(S), H2O(S), OH(S), O(S)]
              kinetics: surface
              reactions: [ptcombust.yaml/reactions: declared-species]
              site-density: 3e-09
        """

        gas = ct.Solution(yaml=phase_defs, name="gas")
        surf = ct.Interface(yaml=phase_defs, name="Pt_surf", adjacent=[gas])
        gas.TPX = 800, 0.01*ct.one_atm, "H2:2.0, O2:1.0"
        surf.TP = 800, 0.01*ct.one_atm
        surf.coverages = {"H(S)": 1.0}  # dummy values to be replaced by delegate

        kHs = surf.species_index("H(S)")
        kPts = surf.species_index("PT(S)")
        kH2 = gas.species_index("H2")
        kO2 = gas.species_index("O2")
        class SurfReactor(ct.ExtensibleIdealGasReactor):
            def replace_eval_surfaces(self, LHS, RHS, sdot):
                site_density = self.surfaces[0].kinetics.site_density
                sdot[:] = 0.0
                LHS[:] = 1.0
                RHS[:] = 0.0
                C = self.thermo.concentrations
                theta = self.surfaces[0].coverages

                # Replace actual reactions with simple absorption of H2 -> H(S)
                rop = 1e-4 * C[kH2] * theta[kPts]
                RHS[kHs] = 2 * rop / site_density
                RHS[kPts] = - 2 * rop / site_density
                sdot[kH2] = - rop * self.surfaces[0].area

            def replace_get_surface_initial_conditions(self, y):
                y[:] = 0
                y[kPts] = 1

            def replace_update_surface_state(self, y):
                # this is the same thing the original method does
                self.surfaces[0].coverages = y

        r1 = SurfReactor(gas)
        r1.volume = 1e-6 # 1 cm^3
        r1.energy_enabled = False
        rsurf = ct.ReactorSurface(surf, A=0.01)
        rsurf.install(r1)
        net = ct.ReactorNet([r1])

        Hweight = ct.Element("H").weight
        total_sites = rsurf.area * surf.site_density
        def masses():
            mass_H = (gas.elemental_mass_fraction("H") * r1.mass +
                      total_sites * r1.surfaces[0].kinetics["H(s)"].X * Hweight)
            mass_O = gas.elemental_mass_fraction("O") * r1.mass
            return mass_H, mass_O

        net.step()
        mH0, mO0 = masses()

        for t in np.linspace(0.1, 1, 12):
            net.advance(t)
            mH, mO = masses()
            self.assertNear(mH, mH0)
            self.assertNear(mO, mO0)

        # Regression test values
        self.assertNear(r1.thermo.P, 647.56016304)
        self.assertNear(r1.thermo.X[kH2], 0.4784268406)
        self.assertNear(r1.thermo.X[kO2], 0.5215731594)
        self.assertNear(r1.surfaces[0].kinetics.X[kHs], 0.3665198138)
        self.assertNear(r1.surfaces[0].kinetics.X[kPts], 0.6334801862)

    def test_interactions(self):
        # Reactors connected by a movable, H2-permeable surface
        kH2 = self.gas.species_index("H2")
        class TestReactor(ct.ExtensibleIdealGasReactor):
            def __init__(self, gas):
                super().__init__(gas)
                self.neighbor = None
                self.h2coeff = 12  # mass transfer coeff
                self.p_coeff = 5 # expansion coeff

            def after_update_connected(self, do_pressure):
                self.conc_H2 = self.thermo.concentrations[kH2]
                self.P = self.thermo.P

            def replace_eval_walls(self, t):
                if self.neighbor:
                    self.vdot = np.clip(self.p_coeff * (self.P - self.neighbor.P),
                                        -1.7, 1.7)

            def after_eval(self, t, LHS, RHS):
                if self.neighbor:
                    mdot_H2 = self.h2coeff * (self.neighbor.conc_H2 - self.conc_H2)
                    RHS[kH2+3] += mdot_H2
                    RHS[3:] -= self.thermo.Y * mdot_H2
                    RHS[0] += mdot_H2
                    # enthalpy flux is neglected, so energy isn't properly conserved

        self.gas.TPX = 300, 2*101325, "H2:0.8, N2:0.2"
        r1 = TestReactor(self.gas)
        self.gas.TPX = 300, 101325, "H2:0.1, O2:0.9"
        r2 = TestReactor(self.gas)
        r1.neighbor = r2
        r2.neighbor = r1
        net = ct.ReactorNet([r1, r2])
        states1 = ct.SolutionArray(self.gas, extra=["t", "mass", "vdot"])
        states2 = ct.SolutionArray(self.gas, extra=["t", "mass", "vdot"])
        net.step()
        M0 = r1.mass * r1.thermo.Y + r2.mass * r2.thermo.Y

        def deltaC():
            return r1.thermo["H2"].concentrations[0] - r2.thermo["H2"].concentrations[0]
        deltaCprev = deltaC()

        V0 = r1.volume + r2.volume
        for t in np.linspace(0.01, 0.2, 50):
            net.advance(t)
            self.assertNear(r1.vdot, -r2.vdot)
            self.assertNear(V0, r1.volume + r2.volume)
            deltaCnow = deltaC()
            self.assertLess(deltaCnow, deltaCprev) # difference is always decreasing
            deltaCprev = deltaCnow
            self.assertArrayNear(M0, r1.mass * r1.thermo.Y + r2.mass * r2.thermo.Y, rtol=2e-8)
            states1.append(r1.thermo.state, t=net.time, mass=r1.mass, vdot=r1.vdot)
            states2.append(r2.thermo.state, t=net.time, mass=r2.mass, vdot=r2.vdot)

        # Regression test values
        self.assertNear(r1.thermo.P, 151561.15, rtol=1e-6)
        self.assertNear(r1.thermo["H2"].Y[0], 0.13765976, rtol=1e-6)
        self.assertNear(r2.thermo["O2"].Y[0], 0.94617029, rtol=1e-6)
