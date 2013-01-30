import numpy as np
import unittest
import re

import cantera as ct
from . import utilities


class TestReactor(utilities.CanteraTest):
    def make_reactors(self, independent=True, n_reactors=2,
                      T1=300, P1=101325, X1='O2:1.0',
                      T2=300, P2=101325, X2='O2:1.0'):

        self.net = ct.ReactorNet()

        self.gas1 = ct.Solution('h2o2.xml')
        self.gas1.TPX = T1, P1, X1
        self.r1 = ct.Reactor(self.gas1)
        self.net.add_reactor(self.r1)

        if independent:
            self.gas2 = ct.Solution('h2o2.xml')
        else:
            self.gas2 = self.gas1

        if n_reactors >= 2:
            self.gas2.TPX = T2, P2, X2
            self.r2 = ct.Reactor(self.gas2)
            self.net.add_reactor(self.r2)

    def add_wall(self, **kwargs):
        self.w = ct.Wall(self.r1, self.r2, **kwargs)
        return self.w

    def test_insert(self):
        R = ct.Reactor()
        f1 = lambda r: r.T
        f2 = lambda r: r.kinetics.net_production_rates
        self.assertRaises(Exception, f1, R)
        self.assertRaises(Exception, f2, R)

        g = ct.Solution('h2o2.xml')
        g.TP = 300, 101325
        R.insert(g)

        self.assertNear(R.T, 300)
        self.assertEqual(len(R.kinetics.net_production_rates), g.n_species)

    def test_names(self):
        self.make_reactors()

        pattern = re.compile(r'(\d+)')
        digits1 = pattern.search(self.r1.name).group(0)
        digits2 = pattern.search(self.r2.name).group(0)

        self.assertEqual(int(digits2), int(digits1) + 1)

        self.r1.name = 'hello'
        self.assertEqual(self.r1.name, 'hello')

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

    def test_timestepping(self):
        self.make_reactors()

        tStart = 0.3
        tEnd = 10.0
        dt_max = 0.07
        t = tStart
        self.net.set_max_time_step(dt_max)
        self.net.set_initial_time(tStart)
        self.assertNear(self.net.time, tStart)

        while t < tEnd:
            tPrev = t
            t = self.net.step(tEnd)
            self.assertTrue(t - tPrev <= 1.0001 * dt_max)
            self.assertNear(t, self.net.time)

        #self.assertNear(self.net.time, tEnd)

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
                t = self.net.step(tEnd)
                nSteps += 1

            return nSteps

        n_baseline = integrate(1e-10, 1e-20)
        n_rtol = integrate(1e-6, 1e-20)
        n_atol = integrate(1e-10, 1e-6)

        self.assertTrue(n_baseline > n_rtol)
        self.assertTrue(n_baseline > n_atol)

    def test_heat_transfer1(self):
        # Connected reactors reach thermal equilibrium after some time
        self.make_reactors(T1=300, T2=1000)
        self.add_wall(U=500, A=1.0)

        self.net.advance(10.0)
        self.assertNear(self.net.time, 10.0)
        self.assertNear(self.r1.T, self.r2.T, 1e-7)
        self.assertNotAlmostEqual(self.r1.thermo.P, self.r2.thermo.P)

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
        # Adiabatic, constant volume combustion should proceed to equilibrum
        # at constant internal energy and volume.

        P0 = 10 * ct.one_atm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'
        self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)

        self.net.advance(1.0)

        gas = ct.Solution('h2o2.xml')
        gas.TPX = T0, P0, X0
        gas.equilibrate('UV')

        self.assertNear(self.r1.T, gas.T)
        self.assertNear(self.r1.thermo.density, gas.density)
        self.assertNear(self.r1.thermo.P, gas.P)
        self.assertArrayNear(self.r1.thermo.X, gas.X)

    def test_equilibrium_HP(self):
        # Adiabatic, constant pressure combustion should proceed to equilibrum
        # at constant enthalpy and pressure.

        P0 = 10 * ct.one_atm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'

        gas1 = ct.Solution('h2o2.xml')
        gas1.TPX = T0, P0, X0
        r1 = ct.ConstPressureReactor(gas1)

        net = ct.ReactorNet()
        net.add_reactor(r1)
        net.advance(1.0)

        gas2 = ct.Solution('h2o2.xml')
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

        self.assertNear(self.r1.volume, V1 + 1.0 * A)
        self.assertNear(self.r2.volume, V2 - 1.0 * A)

    def test_disable_energy(self):
        self.make_reactors(T1=500)
        self.r1.energy_enabled = False
        self.add_wall(A=1.0, U=2500)

        self.net.advance(11.0)

        self.assertNear(self.r1.T, 500)
        self.assertNear(self.r2.T, 500)

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
        gas2 = ct.Solution('h2o2.xml')
        gas2.TPX = 300, 10*101325, 'H2:1.0'
        reservoir = ct.Reservoir(gas2)

        mfc = ct.MassFlowController(reservoir, self.r1)
        mfc.set_mass_flow_rate(lambda t: 0.1 if 0.2 <= t < 1.2 else 0.0)

        self.assertEqual(len(reservoir.inlets), 0)
        self.assertEqual(len(reservoir.outlets), 1)
        self.assertEqual(reservoir.outlets[0], mfc)
        self.assertEqual(len(self.r1.outlets), 0)
        self.assertEqual(len(self.r1.inlets), 1)
        self.assertEqual(self.r1.inlets[0], mfc)

        ma = self.r1.volume * self.r1.density
        Ya = self.r1.Y

        self.net.advance(2.5)

        mb = self.r1.volume * self.r1.density
        Yb = self.r1.Y

        self.assertNear(ma + 0.1, mb)
        self.assertArrayNear(ma * Ya + 0.1 * gas2.Y, mb * Yb)

    def test_valve1(self):
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.set_valve_coeff(k)

        self.assertEqual(self.r1.outlets, self.r2.inlets)
        self.assertTrue(self.r1.energy_enabled)
        self.assertTrue(self.r2.energy_enabled)
        self.assertTrue((self.r1.thermo.P - self.r2.thermo.P) * k,
                        valve.mdot(0))

        m1a = self.r1.thermo.density * self.r1.volume
        m2a = self.r2.thermo.density * self.r2.volume
        Y1a = self.r1.thermo.Y
        Y2a = self.r2.thermo.Y

        self.net.advance(0.1)

        m1b = self.r1.thermo.density * self.r1.volume
        m2b = self.r2.thermo.density * self.r2.volume

        self.assertTrue((self.r1.thermo.P - self.r2.thermo.P) * k,
                        valve.mdot(0.1))
        self.assertNear(m1a+m2a, m1b+m2b)
        Y1b = self.r1.thermo.Y
        Y2b = self.r2.thermo.Y
        self.assertArrayNear(m1a*Y1a + m2a*Y2a, m1b*Y1b + m2b*Y2b)

    def test_valve2(self):
        # Similar to test_valve1, but by disabling the energy equation
        # (constant T) we can compare with an analytical solution for
        # the mass of each reactor as a function of time
        self.make_reactors(P1=10*ct.one_atm)
        self.r1.energy_enabled = False
        self.r2.energy_enabled = False
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.set_valve_coeff(k)

        self.assertFalse(self.r1.energy_enabled)
        self.assertFalse(self.r2.energy_enabled)

        m1a = self.r1.thermo.density * self.r1.volume
        m2a = self.r2.thermo.density * self.r2.volume
        P1a = self.r1.thermo.P
        P2a = self.r2.thermo.P

        A = k * P1a * (1 + m2a/m1a)
        B = k * (P1a/m1a + P2a/m2a)

        for t in np.linspace(1e-5, 0.5):
            self.net.advance(t)
            m1 = self.r1.thermo.density * self.r1.volume
            m2 = self.r2.thermo.density * self.r2.volume
            self.assertNear(m2, (m2a - A/B) * np.exp(-B * t) + A/B)
            self.assertNear(m1a+m2a, m1+m2)

    def test_valve3(self):
        # This case specifies a non-linear relationship between pressure drop
        # and flow rate.
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.set_valve_coeff(mdot)

        t = 0
        while t < 1.0:
            t = self.net.step(1.0)
            p1 = self.r1.thermo.P
            p2 = self.r2.thermo.P
            self.assertNear(mdot(p1-p2), valve.mdot(t))

    def test_pressure_controller(self):
        self.make_reactors(n_reactors=1)
        g = ct.Solution('h2o2.xml')
        g.TPX = 500, 2*101325, 'H2:1.0'
        inlet_reservoir = ct.Reservoir(g)
        g.TP = 300, 101325
        outlet_reservoir = ct.Reservoir(g)

        mfc = ct.MassFlowController(inlet_reservoir, self.r1)
        mdot = lambda t: np.exp(-100*(t-0.5)**2)
        mfc.set_mass_flow_rate(mdot)

        pc = ct.PressureController(self.r1, outlet_reservoir)
        pc.set_master(mfc)
        pc.set_pressure_coeff(1e-5)

        t = 0
        while t < 1.0:
            t = self.net.step(1.0)
            self.assertNear(mdot(t), mfc.mdot(t))
            dP = self.r1.thermo.P - outlet_reservoir.thermo.P
            self.assertNear(mdot(t) + 1e-5 * dP, pc.mdot(t))

    def test_set_initial_time(self):
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.set_valve_coeff(mdot)

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
        valve.set_valve_coeff(mdot)

        t0 = 0.2
        self.net.set_initial_time(t0)
        tf = t0 + 0.5
        self.net.advance(tf)
        self.assertNear(self.net.time, tf)
        p1b = self.r1.thermo.P
        p2b = self.r2.thermo.P

        self.assertNear(p1a, p1b)
        self.assertNear(p2a, p2b)


class TestFlowReactor(utilities.CanteraTest):
    def test_nonreacting(self):
        g = ct.Solution('h2o2.xml')
        g.TPX = 300, 101325, 'O2:1.0'
        r = ct.FlowReactor(g)
        r.mass_flow_rate = 10

        net = ct.ReactorNet()
        net.add_reactor(r)

        t = 0
        v0 = r.speed
        self.assertNear(v0, 10 / r.density)
        while t < 10.0:
            t = net.step(10.0)

            self.assertNear(v0, r.speed)
            self.assertNear(r.distance, v0 * t)

    def test_reacting(self):
        g = ct.Solution('gri30.xml')
        g.TPX = 1400, 20*101325, 'CO:1.0, H2O:1.0'

        r = ct.FlowReactor(g)
        r.mass_flow_rate = 10

        net = ct.ReactorNet()
        net.add_reactor(r)
        net.atol = 1e-18
        net.rtol = 1e-9

        t = 0
        self.assertNear(r.speed, 10 / r.density)
        while t < 1.0:
            t1 = net.time
            x1 = r.distance

            t = net.step(1.0)

            v = (r.distance - x1) / (net.time - t1)
            self.assertNear(r.speed, v, 1e-3)


class TestWallKinetics(utilities.CanteraTest):
    def make_reactors(self):
        self.net = ct.ReactorNet()

        self.gas = ct.Solution('diamond.xml', 'gas')
        self.solid = ct.Solution('diamond.xml', 'diamond')
        self.interface = ct.Interface('diamond.xml', 'diamond_100',
                                      (self.gas, self.solid))
        self.r1 = ct.Reactor(self.gas)
        self.net.add_reactor(self.r1)

        self.r2 = ct.Reactor(self.gas)
        self.net.add_reactor(self.r2)

        self.w = ct.Wall(self.r1, self.r2)

    def test_coverages(self):
        self.make_reactors()
        self.w.left.kinetics = self.interface

        C = np.zeros(self.interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        self.w.left.coverages = C
        self.assertArrayNear(self.w.left.coverages, C)
        self.net.advance(1e-5)
        C_left = self.w.left.coverages

        self.assertEqual(self.w.right.kinetics, None)
        self.assertRaises(Exception, lambda: self.w.right.coverages)

        self.make_reactors()
        self.w.right.kinetics = self.interface
        self.w.right.coverages = C
        self.assertArrayNear(self.w.right.coverages, C)
        self.assertEqual(self.w.left.kinetics, None)
        self.assertRaises(Exception, lambda: self.w.left.coverages)
        self.net.advance(1e-5)
        C_right = self.w.right.coverages

        self.assertNear(sum(C_left), 1.0)
        self.assertArrayNear(C_left, C_right)


@unittest.skipUnless(ct._have_sundials(),
                     "Sensitivity calculations require Sundials")
class TestReactorSensitivities(utilities.CanteraTest):
    def test_sensitivities1(self):
        net = ct.ReactorNet()
        gas = ct.Solution('gri30.xml')
        gas.TPX = 1300, 20*101325, 'CO:1.0, H2:0.1, CH4:0.1, H2O:0.5'
        r1 = ct.Reactor(gas)
        net.add_reactor(r1)

        self.assertEqual(net.n_sensitivity_params, 0)
        r1.add_sensitivity_reaction(40)
        r1.add_sensitivity_reaction(41)

        net.advance(0.1)

        self.assertEqual(net.n_sensitivity_params, 2)
        self.assertEqual(net.n_vars, gas.n_species + 2)
        S = net.sensitivities()
        self.assertEqual(S.shape, (net.n_vars, net.n_sensitivity_params))

    def test_sensitivities2(self):
        net = ct.ReactorNet()

        gas1 = ct.Solution('diamond.xml', 'gas')
        solid = ct.Solution('diamond.xml', 'diamond')
        interface = ct.Interface('diamond.xml', 'diamond_100',
                                 (gas1, solid))
        r1 = ct.Reactor(gas1)
        net.add_reactor(r1)

        gas2 = ct.Solution('h2o2.xml')
        gas2.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'
        r2 = ct.Reactor(gas2)
        net.add_reactor(r2)

        w = ct.Wall(r1, r2)
        w.left.kinetics = interface

        C = np.zeros(interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        w.left.coverages = C
        w.left.add_sensitivity_reaction(2)
        r2.add_sensitivity_reaction(18)

        for T in (901, 905, 910, 950, 1500):
            while r2.T < T:
                net.step(1.0)

            S = net.sensitivities()
            K1 = gas1.n_species + interface.n_species

            # Constant internal energy and volume should generate zero
            # sensitivity coefficients
            self.assertArrayNear(S[0:2,:], np.zeros((2,2)))
            self.assertArrayNear(S[K1+2:K1+4,:], np.zeros((2,2)))

            S11 = np.linalg.norm(S[2:K1+2,0])
            S21 = np.linalg.norm(S[2:K1+2,1])
            S12 = np.linalg.norm(S[K1+4:,0])
            S22 = np.linalg.norm(S[K1+4:,1])

            self.assertTrue(S11 > 1e5 * S12)
            self.assertTrue(S22 > 1e5 * S21)

    def test_parameter_order1(self):
        # Single reactor, changing the order in which parameters are added
        gas = ct.Solution('h2o2.xml')

        def setup():
            net = ct.ReactorNet()
            gas.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'

            r = ct.Reactor(gas)
            net.add_reactor(r)
            return r, net

        def integrate(r, net):
            while r.T < 910:
                net.step(1.0)
            return net.sensitivities()

        r1,net1 = setup()
        params1 = [2,10,18,19]
        for p in params1:
            r1.add_sensitivity_reaction(p)
        S1 = integrate(r1, net1)

        pname = lambda r,i: '%s: %s' % (r.name, gas.reaction_equation(i))
        for i,p in enumerate(params1):
            self.assertEqual(pname(r1,p), net1.sensitivity_parameter_name(i))

        r2,net2 = setup()
        params2 = [19,10,2,18]
        for p in params2:
            r2.add_sensitivity_reaction(p)
        S2 = integrate(r2, net2)

        for i,p in enumerate(params2):
            self.assertEqual(pname(r2,p), net2.sensitivity_parameter_name(i))

        for i,j in enumerate((2,1,3,0)):
            self.assertArrayNear(S1[:,i], S2[:,j])

    def test_parameter_order2(self):
        # Multiple reactors, changing the order in which parameters are added
        gas = ct.Solution('h2o2.xml')

        def setup(reverse=False):
            net = ct.ReactorNet()
            gas1 = ct.Solution('h2o2.xml')
            gas1.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'
            rA = ct.Reactor(gas1)

            gas2 = ct.Solution('h2o2.xml')
            gas2.TPX = 920, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:0.5'
            rB = ct.Reactor(gas2)
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

            pname = lambda r,i: '%s: %s' % (r.name, gas.reaction_equation(i))
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
        N = gas.n_species + 2
        self.assertArrayNear(S[0][:N], S[2][N:], 1e-5, 1e-5)
        self.assertArrayNear(S[0][N:], S[2][:N], 1e-5, 1e-5)
        self.assertArrayNear(S[1][:N], S[3][N:], 1e-5, 1e-5)
        self.assertArrayNear(S[1][N:], S[3][:N], 1e-5, 1e-5)

    def test_parameter_order3(self):
        # Test including reacting surfaces
        gas1 = ct.Solution('diamond.xml', 'gas')
        solid = ct.Solution('diamond.xml', 'diamond')
        interface = ct.Interface('diamond.xml', 'diamond_100',
                                 (gas1, solid))

        gas2 = ct.Solution('h2o2.xml')

        def setup(order):
            gas1.TPX = 1200, 1e3, 'H:0.002, H2:1, CH4:0.01, CH3:0.0002'
            gas2.TPX = 900, 101325, 'H2:0.1, OH:1e-7, O2:0.1, AR:1e-5'
            net = ct.ReactorNet()
            rA = ct.Reactor(gas1)
            rB = ct.Reactor(gas2)

            if order % 2 == 0:
                wA = ct.Wall(rA, rB)
                wB = ct.Wall(rB, rA)
            else:
                wB = ct.Wall(rB, rA)
                wA = ct.Wall(rA, rB)

            wA.left.kinetics = interface
            wB.right.kinetics = interface

            wA.area = 0.1
            wB.area = 10

            C1 = np.zeros(interface.n_species)
            C2 = np.zeros(interface.n_species)
            C1[0] = 0.3
            C1[4] = 0.7

            C2[0] = 0.9
            C2[4] = 0.1
            wA.left.coverages = C1
            wB.right.coverages = C2

            if order // 2 == 0:
                net.add_reactor(rA)
                net.add_reactor(rB)
            else:
                net.add_reactor(rB)
                net.add_reactor(rA)

            return rA,rB,wA,wB,net

        def integrate(r, net):
            net.advance(1e-4)
            return net.sensitivities()

        S = []

        for order in range(4):
            rA,rB,wA,wB,net = setup(order)
            for (obj,k) in [(rB,2), (rB,18), (wA.left,2),
                            (wA.left,0), (wB.right,2)]:
                obj.add_sensitivity_reaction(k)
            integrate(rB, net)
            S.append(net.sensitivities())

            rA,rB,wA,wB,net = setup(order)
            for (obj,k) in [(wB.right,2), (wA.left,2), (rB,18),
                            (wA.left,0), (rB,2)]:
                obj.add_sensitivity_reaction(k)

            integrate(rB, net)
            S.append(net.sensitivities())

        for a,b in [(0,1),(2,3),(4,5),(6,7)]:
            for i,j in enumerate((4,2,1,3,0)):
                self.assertArrayNear(S[a][:,i], S[b][:,j], 1e-2, 1e-3)
