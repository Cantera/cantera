import numpy as np
import re

import cantera as ct
from . import utilities


class TestReactor(utilities.CanteraTest):
    def makeReactors(self, independent=True, nReactors=2,
                      T1=300, P1=101325, X1='O2:1.0',
                      T2=300, P2=101325, X2='O2:1.0'):

        self.net = ct.ReactorNet()

        self.gas1 = ct.Solution('h2o2.xml')
        self.gas1.TPX = T1, P1, X1
        self.r1 = ct.Reactor(self.gas1)
        self.net.addReactor(self.r1)

        if independent:
            self.gas2 = ct.Solution('h2o2.xml')
        else:
            self.gas2 = self.gas1

        if nReactors >= 2:
            self.gas2.TPX = T2, P2, X2
            self.r2 = ct.Reactor(self.gas2)
            self.net.addReactor(self.r2)

    def addWall(self, **kwargs):
        self.w = ct.Wall(self.r1, self.r2, **kwargs)
        return self.w

    def test_insert(self):
        R = ct.Reactor()
        f1 = lambda r: r.T
        f2 = lambda r: r.kinetics.netProductionRates
        self.assertRaises(Exception, f1, R)
        self.assertRaises(Exception, f2, R)

        g = ct.Solution('h2o2.xml')
        g.TP = 300, 101325
        R.insert(g)

        self.assertNear(R.T, 300)
        self.assertEqual(len(R.kinetics.netProductionRates), g.nSpecies)

    def test_names(self):
        self.makeReactors()

        pattern = re.compile(r'(\d+)')
        digits1 = pattern.search(self.r1.name).group(0)
        digits2 = pattern.search(self.r2.name).group(0)

        self.assertEqual(int(digits2), int(digits1) + 1)

        self.r1.name = 'hello'
        self.assertEqual(self.r1.name, 'hello')

    def test_disjoint(self):
        T1, P1 = 300, 101325
        T2, P2 = 500, 300000

        self.makeReactors(T1=T1, T2=T2, P1=P1, P2=P2)
        self.net.advance(1.0)

        # Nothing should change from the initial condition
        self.assertNear(T1, self.gas1.T)
        self.assertNear(T2, self.gas2.T)
        self.assertNear(P1, self.gas1.P)
        self.assertNear(P2, self.gas2.P)

    def test_disjoint2(self):
        T1, P1 = 300, 101325
        T2, P2 = 500, 300000

        self.makeReactors(T1=T1, T2=T2, P1=P1, P2=P2, independent=False)
        self.net.advance(1.0)

        # Nothing should change from the initial condition
        self.assertNear(T1, self.r1.T)
        self.assertNear(T2, self.r2.T)
        self.assertNear(P1, self.r1.thermo.P)
        self.assertNear(P2, self.r2.thermo.P)

    def test_timestepping(self):
        self.makeReactors()

        tStart = 0.3
        tEnd = 10.0
        dt_max = 0.07
        t = tStart
        self.net.setMaxTimeStep(dt_max)
        self.net.setInitialTime(tStart)
        self.assertNear(self.net.time, tStart)

        while t < tEnd:
            tPrev = t
            t = self.net.step(tEnd)
            self.assertTrue(t - tPrev <= 1.0001 * dt_max)
            self.assertNear(t, self.net.time)

        #self.assertNear(self.net.time, tEnd)

    def test_equalizePressure(self):
        self.makeReactors(P1=101325, P2=300000)
        self.addWall(K=0.1, A=1.0)

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
            P0 = 10 * ct.OneAtm
            T0 = 1100
            X0 = 'H2:1.0, O2:0.5, AR:8.0'
            self.makeReactors(nReactors=1, T1=T0, P1=P0, X1=X0)
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

        n_baseline = integrate(1e-6, 1e-10)
        n_rtol = integrate(1e-10, 1e-10)
        n_atol = integrate(1e-6, 1e-18)

        self.assertTrue(n_baseline < n_rtol)
        self.assertTrue(n_baseline < n_atol)

    def test_heatTransfer1(self):
        # Connected reactors reach thermal equilibrium after some time
        self.makeReactors(T1=300, T2=1000)
        self.addWall(U=500, A=1.0)

        self.net.advance(10.0)
        self.assertNear(self.net.time, 10.0)
        self.assertNear(self.r1.T, self.r2.T)
        self.assertNotAlmostEqual(self.r1.thermo.P, self.r2.thermo.P)

    def test_heatTransfer2(self):
        # Result should be the same if (m * cp) / (U * A) is held constant
        self.makeReactors(T1=300, T2=1000)
        self.addWall(U=200, A=1.0)

        self.net.advance(1.0)
        T1a = self.r1.T
        T2a = self.r2.T

        self.makeReactors(T1=300, T2=1000)
        self.r1.volume = 0.25
        self.r2.volume = 0.25
        w = self.addWall(U=100, A=0.5)

        self.assertNear(w.heatTransferCoeff * w.area * (self.r1.T - self.r2.T),
                        w.qdot(0))
        self.net.advance(1.0)
        self.assertNear(w.heatTransferCoeff * w.area * (self.r1.T - self.r2.T),
                        w.qdot(1.0))
        T1b = self.r1.T
        T2b = self.r2.T

        self.assertNear(T1a, T1b)
        self.assertNear(T2a, T2b)

    def test_equilibrium_UV(self):
        # Adiabatic, constant volume combustion should proceed to equilibrum
        # at constant internal energy and volume.

        P0 = 10 * ct.OneAtm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'
        self.makeReactors(nReactors=1, T1=T0, P1=P0, X1=X0)

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

        P0 = 10 * ct.OneAtm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'

        gas1 = ct.Solution('h2o2.xml')
        gas1.TPX = T0, P0, X0
        r1 = ct.ConstPressureReactor(gas1)

        net = ct.ReactorNet()
        net.addReactor(r1)
        net.advance(1.0)

        gas2 = ct.Solution('h2o2.xml')
        gas2.TPX = T0, P0, X0
        gas2.equilibrate('HP')

        self.assertNear(r1.T, gas2.T)
        self.assertNear(r1.thermo.P, P0)
        self.assertNear(r1.thermo.density, gas2.density)
        self.assertArrayNear(r1.thermo.X, gas2.X)

    def test_wall_velocity(self):
        self.makeReactors()
        A = 0.2

        V1 = 2.0
        V2 = 5.0
        self.r1.volume = V1
        self.r2.volume = V2

        self.addWall(A=A)

        def v(t):
            if 0 < t <= 1:
                return t
            elif 1 <= t <= 2:
                return 2 - t
            else:
                return 0.0

        self.w.setVelocity(v)
        self.net.advance(1.0)
        self.assertNear(self.w.vdot(1.0), 1.0 * A, 1e-7)
        self.net.advance(2.0)
        self.assertNear(self.w.vdot(2.0), 0.0, 1e-7)

        self.assertNear(self.r1.volume, V1 + 1.0 * A)
        self.assertNear(self.r2.volume, V2 - 1.0 * A)

    def test_disable_energy(self):
        self.makeReactors(T1=500)
        self.r1.energyEnabled = False
        self.addWall(A=1.0, U=2500)

        self.net.advance(11.0)

        self.assertNear(self.r1.T, 500)
        self.assertNear(self.r2.T, 500)

    def test_heat_flux_func(self):
        self.makeReactors(T1=500, T2=300)
        self.r1.volume = 0.5

        U1a = self.r1.volume * self.r1.density * self.r1.thermo.u
        U2a = self.r2.volume * self.r2.density * self.r2.thermo.u

        V1a = self.r1.volume
        V2a = self.r2.volume

        self.addWall(A=0.3)
        self.w.setHeatFlux(lambda t: 90000 * (1 - t**2) if t <= 1.0 else 0.0)
        Q = 0.3 * 60000

        self.net.advance(1.1)
        U1b = self.r1.volume * self.r1.density * self.r1.thermo.u
        U2b = self.r2.volume * self.r2.density * self.r2.thermo.u

        self.assertNear(V1a, self.r1.volume)
        self.assertNear(V2a, self.r2.volume)
        self.assertNear(U1a - Q, U1b, 1e-6)
        self.assertNear(U2a + Q, U2b, 1e-6)

    def test_mass_flow_controller(self):
        self.makeReactors(nReactors=1)
        gas2 = ct.Solution('h2o2.xml')
        gas2.TPX = 300, 10*101325, 'H2:1.0'
        reservoir = ct.Reservoir(gas2)

        mfc = ct.MassFlowController(reservoir, self.r1)
        mfc.setMassFlowRate(lambda t: 0.1 if 0.2 <= t < 1.2 else 0.0)

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
        self.makeReactors(P1=10*ct.OneAtm, X1='AR:1.0', X2='O2:1.0')
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.setValveCoeff(k)

        self.assertEqual(self.r1.outlets, self.r2.inlets)
        self.assertTrue(self.r1.energyEnabled)
        self.assertTrue(self.r2.energyEnabled)
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
        self.makeReactors(P1=10*ct.OneAtm)
        self.r1.energyEnabled = False
        self.r2.energyEnabled = False
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.setValveCoeff(k)

        self.assertFalse(self.r1.energyEnabled)
        self.assertFalse(self.r2.energyEnabled)

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
        self.makeReactors(P1=10*ct.OneAtm, X1='AR:1.0', X2='O2:1.0')
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.setValveCoeff(mdot)

        t = 0
        while t < 1.0:
            t = self.net.step(1.0)
            p1 = self.r1.thermo.P
            p2 = self.r2.thermo.P
            self.assertNear(mdot(p1-p2), valve.mdot(t))

    def test_pressure_controller(self):
        self.makeReactors(nReactors=1)
        g = ct.Solution('h2o2.xml')
        g.TPX = 500, 2*101325, 'H2:1.0'
        inletReservoir = ct.Reservoir(g)
        g.TP = 300, 101325
        outletReservoir = ct.Reservoir(g)

        mfc = ct.MassFlowController(inletReservoir, self.r1)
        mdot = lambda t: np.exp(-100*(t-0.5)**2)
        mfc.setMassFlowRate(mdot)

        pc = ct.PressureController(self.r1, outletReservoir)
        pc.setMaster(mfc)
        pc.setPressureCoeff(1e-5)

        t = 0
        while t < 1.0:
            t = self.net.step(1.0)
            self.assertNear(mdot(t), mfc.mdot(t))
            dP = self.r1.thermo.P - outletReservoir.thermo.P
            self.assertNear(mdot(t) + 1e-5 * dP, pc.mdot(t))

    def test_setInitialTime(self):
        self.makeReactors(P1=10*ct.OneAtm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.setValveCoeff(mdot)

        t0 = 0.0
        tf = t0 + 0.5
        self.net.advance(tf)
        self.assertNear(self.net.time, tf)
        p1a = self.r1.thermo.P
        p2a = self.r2.thermo.P

        self.makeReactors(P1=10*ct.OneAtm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.setValveCoeff(mdot)

        t0 = 0.2
        self.net.setInitialTime(t0)
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
        r.massFlowRate = 10

        net = ct.ReactorNet()
        net.addReactor(r)

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
        r.massFlowRate = 10

        net = ct.ReactorNet()
        net.addReactor(r)
        net.atol = 1e-22
        net.rtol = 1e-8

        t = 0
        self.assertNear(r.speed, 10 / r.density)
        while t < 1.0:
            t1 = net.time
            x1 = r.distance

            t = net.step(1.0)

            v = (r.distance - x1) / (net.time - t1)
            self.assertNear(r.speed, v, 1e-3)


class TestWallKinetics(utilities.CanteraTest):
    def makeReactors(self):

        self.net = ct.ReactorNet()

        self.gas = ct.Solution('diamond.xml', 'gas')
        self.solid = ct.Solution('diamond.xml', 'diamond')
        self.interface = ct.Interface('diamond.xml', 'diamond_100',
                                      (self.gas, self.solid))
        self.r1 = ct.Reactor(self.gas)
        self.net.addReactor(self.r1)

        self.r2 = ct.Reactor(self.gas)
        self.net.addReactor(self.r2)

        self.w = ct.Wall(self.r1, self.r2)

    def test_coverages(self):
        self.makeReactors()
        self.w.leftKinetics = self.interface

        C = np.zeros(self.interface.nSpecies)
        C[0] = 0.3
        C[4] = 0.7

        self.w.leftCoverages = C
        self.assertArrayNear(self.w.leftCoverages, C)
        self.net.advance(1e-5)
        C_left = self.w.leftCoverages

        self.assertEqual(self.w.rightKinetics, None)
        self.assertRaises(Exception, lambda: self.w.rightCoverages)

        self.makeReactors()
        self.w.rightKinetics = self.interface
        self.w.rightCoverages = C
        self.assertArrayNear(self.w.rightCoverages, C)
        self.assertEqual(self.w.leftKinetics, None)
        self.assertRaises(Exception, lambda: self.w.leftCoverages)
        self.net.advance(1e-5)
        C_right = self.w.rightCoverages

        self.assertNear(sum(C_left), 1.0)
        self.assertArrayNear(C_left, C_right)
