import math
import numpy as np
import pytest
from pytest import approx
import re

import cantera as ct

try:
    ct.drawnetwork._import_graphviz()
except ImportError:
    pass

from cantera.drawnetwork import _graphviz

from .utilities import (
    compareProfiles
)

class TestReactor:
    reactorClass = ct.Reactor

    def make_reactors(self, independent=True, n_reactors=2,
                      T1=300, P1=101325, X1='O2:1.0',
                      T2=300, P2=101325, X2='O2:1.0'):

        self.net = ct.ReactorNet()
        assert self.net.initial_time == 0.

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
        assert not self.net.verbose
        self.net.verbose = True
        assert self.net.verbose

    def test_volume(self):
        g = ct.Solution('h2o2.yaml', transport_model=None)
        R = self.reactorClass(g, volume=11)
        assert R.volume == 11

        R.volume = 9
        assert R.volume ==  9

    def test_names(self):
        self.make_reactors()

        pattern = re.compile(r'(\d+)')
        digits1 = pattern.search(self.r1.name).group(0)
        digits2 = pattern.search(self.r2.name).group(0)

        assert int(digits2) == int(digits1) + 1

        self.r1.name = 'hello'
        assert self.r1.name == 'hello'

    def test_types(self):
        self.make_reactors()
        assert self.r1.type == self.reactorClass.__name__

    def test_component_index(self):
        self.make_reactors(n_reactors=1)
        self.net.step()

        N0 = self.net.n_vars - self.gas1.n_species
        for i, name in enumerate(self.gas1.species_names):
            assert i + N0 == self.r1.component_index(name)

    def test_component_names(self):
        self.make_reactors(n_reactors=2)
        self.net.initialize()
        N = self.net.n_vars // 2
        for i in range(N):
            assert self.r1.component_index(self.r1.component_name(i)) ==  i
            assert (self.net.component_name(i)
                    == '{}: {}'.format(self.r1.name, self.r1.component_name(i)))
            assert (self.net.component_name(N+i)
                    == '{}: {}'.format(self.r2.name, self.r2.component_name(i)))

    def test_independent_variable(self):
        self.make_reactors(independent=False, n_reactors=1)

        with pytest.raises(ct.CanteraError, match="independent variable"):
            self.net.distance

        assert self.net.time == 0.0

    def test_disjoint(self):
        T1, P1 = 300, 101325
        T2, P2 = 500, 300000

        self.make_reactors(T1=T1, T2=T2, P1=P1, P2=P2)
        self.net.advance(1.0)

        # Nothing should change from the initial condition
        assert T1 == approx(self.gas1.T)
        assert T2 == approx(self.gas2.T)
        assert P1 == approx(self.gas1.P)
        assert P2 == approx(self.gas2.P)

    def test_disjoint2(self):
        T1, P1 = 300, 101325
        T2, P2 = 500, 300000

        self.make_reactors(T1=T1, T2=T2, P1=P1, P2=P2, independent=False)
        self.net.advance(1.0)

        # Nothing should change from the initial condition
        assert T1 == approx(self.r1.T)
        assert T2 == approx(self.r2.T)
        assert P1 == approx(self.r1.thermo.P)
        assert P2 == approx(self.r2.thermo.P)

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
            assert dydt[i] == approx(dy[i]/dt)

    def test_finite_difference_jacobian(self):
        self.make_reactors(n_reactors=1, T1=900, P1=101325, X1="H2:0.4, O2:0.4, N2:0.2")
        kH2 = self.gas1.species_index("H2")
        while self.r1.thermo.X[kH2] > 0.3:
            self.net.step()

        J = self.r1.finite_difference_jacobian
        assert J.shape == (self.r1.n_vars, self.r1.n_vars)

        # state variables that should be constant, depending on reactor type
        constant = {"mass", "volume", "int_energy", "enthalpy", "pressure"}
        variable = {"temperature"}
        for i in range(3):
            name = self.r1.component_name(i)
            if name in constant:
                assert all(J[i,:] == 0), (i, name)
            elif name in variable:
                assert any(J[i,:] != 0)

        # Disabling energy equation should zero these terms
        self.r1.energy_enabled = False
        J = self.r1.finite_difference_jacobian
        for i in range(3):
            name = self.r1.component_name(i)
            if name == "temperature":
                assert all(J[i,:] == 0)

        # Disabling species equations should zero these terms
        self.r1.energy_enabled = True
        self.r1.chemistry_enabled = False
        J = self.r1.finite_difference_jacobian
        constant = set(self.gas1.species_names)
        species_start = self.r1.component_index(self.gas1.species_name(0))
        for i in range(self.r1.n_vars):
            name = self.r1.component_name(i)
            if name in constant:
                assert all(J[i, species_start:] == 0), (i, name)

    def test_timestepping(self):
        self.make_reactors()

        tStart = 0.3
        tEnd = 10.0
        dt_max = 0.07
        t = tStart

        self.net.max_time_step = dt_max
        assert self.net.max_time_step == dt_max
        self.net.initial_time = tStart
        assert self.net.initial_time == tStart
        assert self.net.time == approx(tStart)

        while t < tEnd:
            tPrev = t
            t = self.net.step()
            assert t - tPrev <= 1.0001 * dt_max
            assert t == approx(self.net.time)

    def test_maxsteps(self):
        self.make_reactors()

        # set the up a case where we can't take
        # enough time-steps to reach the endtime
        max_steps = 10
        max_step_size = 1e-07
        self.net.initial_time = 0.
        self.net.max_time_step = max_step_size
        self.net.max_steps = max_steps
        with pytest.raises(
                ct.CanteraError, match='Maximum number of timesteps'):
            self.net.advance(1e-04)
        assert self.net.time <= max_steps * max_step_size
        assert self.net.max_steps == max_steps

    def test_wall_type1(self):
        self.make_reactors(P1=101325, P2=300000)
        self.add_wall(K=0.1, A=1.0)
        net = ct.ReactorNet([self.r1, self.r2])  # assigns default names
        assert self.r1.name.startswith(f"{self.r1.type}_")  # default name
        assert self.r2.name.startswith(f"{self.r2.type}_")  # default name
        assert self.w.type == "Wall"
        assert self.w.name.startswith("Wall_")  # default name
        self.w.name = "name-of-wall"
        assert self.w.name == "name-of-wall"

    def test_wall_type2(self):
        self.make_reactors(n_reactors=1)
        res = ct.Reservoir(self.gas1)
        w = ct.Wall(self.r1, res)
        net = ct.ReactorNet([self.r1])  # assigns default names
        assert self.r1.name.startswith(f"{self.r1.type}_")  # default name
        assert w.type == "Wall"
        assert w.name.startswith("Wall_")  # default name

    def test_wall_type3(self):
        self.make_reactors(n_reactors=1)
        res = ct.Reservoir(self.gas1)
        w = ct.Wall(res, self.r1)
        net = ct.ReactorNet([self.r1])  # assigns default names
        assert self.r1.name.startswith(f"{self.r1.type}_")  # default name
        assert w.type == "Wall"
        assert w.name.startswith("Wall_")  # default name

    def test_equalize_pressure(self):
        self.make_reactors(P1=101325, P2=300000)
        self.add_wall(K=0.1, A=1.0)

        assert len(self.r1.walls) == 1
        assert len(self.r2.walls) == 1
        assert self.r1.walls[0] == self.w
        assert self.r2.walls[0] == self.w

        self.net.advance(1.0)

        assert self.net.time == approx(1.0)
        assert self.gas1.P == approx(self.gas2.P)
        assert self.r1.T != approx(self.r2.T)

    def test_tolerances(self, rtol_lim=1e-10, atol_lim=1e-20):
        def integrate(atol, rtol):
            P0 = 10 * ct.one_atm
            T0 = 1100
            X0 = 'H2:1.0, O2:0.5, AR:8.0'
            self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)
            self.net.rtol = rtol
            self.net.atol = atol

            assert self.net.rtol == rtol
            assert self.net.atol == atol

            tEnd = 1.0
            nSteps = 0
            t = 0

            while t < tEnd:
                t = self.net.step()
                nSteps += 1

            return nSteps

        n_baseline = integrate(rtol_lim, atol_lim)
        n_rtol = integrate(rtol_lim * 1e2, atol_lim)
        n_atol = integrate(rtol_lim, atol_lim * 1e15)
        assert n_baseline > n_rtol
        assert n_baseline > n_atol

    def test_advance_limits(self):
        P0 = 10 * ct.one_atm
        T0 = 1100
        X0 = 'H2:1.0, O2:0.5, AR:8.0'
        self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)

        limit_H2 = .01
        ix = self.net.global_component_index('H2', 0)
        self.r1.set_advance_limit('H2', limit_H2)
        assert self.net.advance_limits[ix] == limit_H2

        self.r1.set_advance_limit('H2', None)
        assert self.net.advance_limits[ix] == -1.

        self.r1.set_advance_limit('H2', limit_H2)
        self.net.advance_limits = None
        assert self.net.advance_limits[ix] == -1.

        self.r1.set_advance_limit('H2', limit_H2)
        self.net.advance_limits = 0 * self.net.advance_limits - 1.
        assert self.net.advance_limits[ix] == -1.

    @pytest.mark.xfail(reason="See GitHub Issue #1453")
    def test_advance_with_limits(self):
        def integrate(limit_H2 = None, apply=True):
            P0 = 10 * ct.one_atm
            T0 = 1100
            X0 = 'H2:1.0, O2:0.5, AR:8.0'
            self.make_reactors(n_reactors=1, T1=T0, P1=P0, X1=X0)
            if limit_H2 is not None:
                self.r1.set_advance_limit('H2', limit_H2)
                ix = self.net.global_component_index('H2', 0)
                assert self.net.advance_limits[ix] == limit_H2

            tEnd = 0.1
            tStep = 7e-4
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

        assert n_advance_coarse > n_baseline
        assert n_advance_fine > n_advance_coarse
        assert n_advance_negative == n_baseline
        assert n_advance_override == n_baseline

    def test_heat_transfer1(self):
        # Connected reactors reach thermal equilibrium after some time
        self.make_reactors(T1=300, T2=1000)
        self.add_wall(U=500, A=1.0)

        self.net.advance(10.0)
        assert self.net.time == approx(10.0)
        assert self.r1.T == approx(self.r2.T, rel=5e-7)
        assert self.r1.thermo.P != approx(self.r2.thermo.P)

    def test_advance_limits_invalid(self):
        self.make_reactors(n_reactors=1)

        with pytest.raises(ct.CanteraError, match="No component named 'spam'"):
            self.r1.set_advance_limit("spam", 0.1)

    def test_advance_reverse(self):
        self.make_reactors(n_reactors=1)
        self.net.advance(0.1)
        with pytest.raises(ct.CanteraError, match="backwards in time"):
            self.net.advance(0.09)

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

        assert (w.heat_transfer_coeff * w.area * (self.r1.T - self.r2.T)
                == approx(w.heat_rate))
        self.net.advance(1.0)
        assert (w.heat_transfer_coeff * w.area * (self.r1.T - self.r2.T)
                == approx(w.heat_rate))
        T1b = self.r1.T
        T2b = self.r2.T

        assert T1a == approx(T1b)
        assert T2a == approx(T2b)

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

        assert self.r1.T == approx(gas.T)
        assert self.r1.thermo.density == approx(gas.density)
        assert self.r1.thermo.P == approx(gas.P)
        assert self.r1.thermo.X == approx(gas.X)

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

        assert r1.T == approx(gas2.T)
        assert r1.thermo.P == approx(P0)
        assert r1.thermo.density == approx(gas2.density)
        assert r1.thermo.X == approx(gas2.X)

    def test_wall_velocity(self):
        self.make_reactors()
        A = 0.2

        V1 = 2.0
        V2 = 5.0
        self.r1.volume = V1
        self.r2.volume = V2

        self.add_wall(A=A)

        v = ct.Tabulated1([0.0, 1.0, 2.0], [0.0, 1.0, 0.0])

        self.w.velocity = v
        self.net.advance(1.0)
        assert self.w.velocity == approx(v(1.0))
        assert self.w.expansion_rate == approx(1.0 * A, rel=1e-7)
        self.net.advance(2.0)
        assert self.w.expansion_rate == approx(0.0, rel=1e-7)

        assert self.r1.volume == approx(V1 + 1.0 * A, rel=1e-7)
        assert self.r2.volume == approx(V2 - 1.0 * A, rel=1e-7)

    def test_disable_energy(self):
        self.make_reactors(T1=500)
        self.r1.energy_enabled = False
        self.add_wall(A=1.0, U=2500)
        self.net.advance(11.0)

        assert self.r1.T == approx(500)
        assert self.r2.T == approx(500)

    def test_disable_chemistry(self):
        self.make_reactors(T1=1000, n_reactors=1, X1='H2:2.0,O2:1.0')
        self.r1.chemistry_enabled = False

        self.net.advance(11.0)

        assert self.r1.T == approx(1000)
        assert self.r1.thermo.X[self.r1.thermo.species_index('H2')] == approx(2.0/3.0)
        assert self.r1.thermo.X[self.r1.thermo.species_index('O2')] == approx(1.0/3.0)

    def test_heat_flux_func(self):
        self.make_reactors(T1=500, T2=300)
        self.r1.volume = 0.5

        U1a = self.r1.volume * self.r1.density * self.r1.thermo.u
        U2a = self.r2.volume * self.r2.density * self.r2.thermo.u

        V1a = self.r1.volume
        V2a = self.r2.volume

        self.add_wall(A=0.3)
        hfunc = lambda t: 90000 * (1 - t**2) if t <= 1.0 else 0.0
        self.w.heat_flux = hfunc
        Q = 0.3 * 60000

        self.net.advance(1.1)
        assert self.w.heat_flux == hfunc(1.1)
        U1b = self.r1.volume * self.r1.density * self.r1.thermo.u
        U2b = self.r2.volume * self.r2.density * self.r2.thermo.u

        assert V1a == approx(self.r1.volume)
        assert V2a == approx(self.r2.volume)
        assert U1a - Q == approx(U1b, rel=1e-6)
        assert U2a + Q == approx(U2b, rel=1e-6)

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
        assert mfc.mass_flow_coeff == 1.

        assert mfc.type == type(mfc).__name__
        assert len(reservoir.inlets) == 0
        assert len(reservoir.outlets) == 1
        assert reservoir.outlets[0] == mfc
        assert len(self.r1.outlets) == 0
        assert len(self.r1.inlets) == 1
        assert self.r1.inlets[0] == mfc

        ma = self.r1.volume * self.r1.density
        Ya = self.r1.Y

        self.net.rtol = 1e-11
        self.net.max_time_step = 0.05

        self.net.advance(0.1)
        assert mfc.mass_flow_rate == approx(0.)
        self.net.advance(0.3)
        assert mfc.mass_flow_rate == approx(0.04)
        self.net.advance(1.0)
        assert mfc.mass_flow_rate == approx(0.08)
        self.net.advance(1.2)
        assert mfc.mass_flow_rate == approx(0.)

        self.net.advance(2.5)

        mb = self.r1.volume * self.r1.density
        Yb = self.r1.Y

        assert ma + 0.1 == approx(mb)
        assert ma * Ya + 0.1 * gas2.Y == approx(mb * Yb)

    def test_mass_flow_controller_type(self):
        self.make_reactors(n_reactors=2)
        mfc = ct.MassFlowController(self.r1, self.r2)
        net = ct.ReactorNet([self.r1, self.r2])  # assigns default names
        assert mfc.type == "MassFlowController"
        assert mfc.name.startswith("MassFlowController_")  # default name
        mfc.name = "name-of-mfc"
        assert mfc.name == "name-of-mfc"

    def test_mass_flow_controller_errors(self):
        # Make sure Python error message actually gets displayed
        self.make_reactors(n_reactors=2)
        mfc = ct.MassFlowController(self.r1, self.r2)
        mfc.mass_flow_rate = lambda t: eggs

        with pytest.raises(Exception, match='eggs'):
            self.net.step()

        with pytest.raises(NotImplementedError):
            mfc.pressure_function = lambda p: p**2

    def test_valve1(self):
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.valve_coeff = k

        assert self.r1.outlets == self.r2.inlets
        assert valve.valve_coeff == k
        assert self.r1.energy_enabled
        assert self.r2.energy_enabled
        self.net.initialize()
        assert (self.r1.thermo.P - self.r2.thermo.P) * k == approx(
                valve.mass_flow_rate)

        m1a = self.r1.thermo.density * self.r1.volume
        m2a = self.r2.thermo.density * self.r2.volume
        Y1a = self.r1.thermo.Y
        Y2a = self.r2.thermo.Y

        self.net.advance(0.1)

        m1b = self.r1.thermo.density * self.r1.volume
        m2b = self.r2.thermo.density * self.r2.volume

        assert (self.r1.thermo.P - self.r2.thermo.P) * k == approx(
                valve.mass_flow_rate)
        assert m1a + m2a == approx(m1b + m2b)
        Y1b = self.r1.thermo.Y
        Y2b = self.r2.thermo.Y
        assert m1a*Y1a + m2a*Y2a == approx(m1b*Y1b + m2b*Y2b, abs=1e-10)
        assert Y1a == approx(Y1b)

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
        assert valve.valve_coeff == k

        assert not self.r1.energy_enabled
        assert not self.r2.energy_enabled

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
            assert m2 == approx((m2a - A/B) * np.exp(-B * t) + A/B)
            assert m1a + m2a == approx(m1 + m2)
            assert self.r1.Y == approx(Y1)

    def test_valve3(self):
        # This case specifies a non-linear relationship between pressure drop
        # and flow rate.
        self.make_reactors(P1=10*ct.one_atm, X1='AR:0.5, O2:0.5',
                           X2='O2:1.0')
        self.net.rtol = 1e-12
        self.net.atol = 1e-20
        valve = ct.Valve(self.r1, self.r2)
        mdot = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.pressure_function = mdot
        assert valve.valve_coeff == 1.

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
            assert mdot(p1-p2) == approx(valve.mass_flow_rate)
            assert Y1 == approx(self.r1.Y)
            assert speciesMass(kAr) == approx(mAr)
            assert speciesMass(kO2) == approx(mO2)

    def test_valve_timing(self):
        # test timed valve
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        k = 2e-5
        valve.valve_coeff = k
        valve.time_function = lambda t: t > .01

        delta_p = lambda: self.r1.thermo.P - self.r2.thermo.P
        mdot = lambda: valve.valve_coeff * (self.r1.thermo.P - self.r2.thermo.P)
        self.net.initialize()
        assert valve.time_function == 0.0
        assert valve.pressure_function == approx(delta_p())
        assert valve.mass_flow_rate == 0.0
        self.net.advance(0.01)
        assert valve.time_function == 0.0
        assert valve.pressure_function == approx(delta_p())
        assert valve.mass_flow_rate == 0.0
        self.net.advance(0.01 + 1e-9)
        assert valve.time_function == 1.0
        assert valve.pressure_function == approx(delta_p())
        assert valve.mass_flow_rate == approx(mdot())
        self.net.advance(0.02)
        assert valve.time_function == 1.0
        assert valve.pressure_function == approx(delta_p())
        assert valve.mass_flow_rate == approx(mdot())

    def test_valve_type1(self):
        self.make_reactors()
        res = ct.Reservoir(self.gas1)
        v = ct.Valve(self.r1, res)
        ct.ReactorNet([self.r1])  # assigns default names
        assert self.r1.name.startswith(f"{self.r1.type}_")  # default name
        assert res.name.startswith(f"{res.type}_")  # default name
        assert v.type == "Valve"
        assert v.name.startswith("Valve_")  # default name
        v.name = "name-of-valve"
        assert v.name == "name-of-valve"

    def test_valve_type2(self):
        self.make_reactors()
        res = ct.Reservoir(self.gas1)
        ct.Valve(res, self.r1)
        ct.ReactorNet([self.r1])  # assigns default names
        assert self.r1.name.startswith(f"{self.r1.type}_")  # default name
        assert res.name.startswith(f"{res.type}_")  # default name

    def test_valve_errors(self):
        self.make_reactors()
        v = ct.Valve(self.r1, self.r2)
        with pytest.raises(ct.CanteraError, match='Already installed'):
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
        mfc.time_function = mdot

        pc = ct.PressureController(self.r1, outlet_reservoir)
        pc.primary = mfc
        pc.pressure_coeff = 1e-5
        assert pc.pressure_coeff == 1e-5

        t = 0
        while t < 1.0:
            t = self.net.step()
            assert mdot(t) == approx(mfc.mass_flow_rate)
            dP = self.r1.thermo.P - outlet_reservoir.thermo.P
            assert mdot(t) + 1e-5 * dP == approx(pc.mass_flow_rate)

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
        mfc.time_function = mdot

        pc = ct.PressureController(self.r1, outlet_reservoir)
        pc.primary = mfc
        pfunc = lambda dp: 1.e-5 * abs(dp)**.5
        pc.pressure_function = pfunc
        assert pc.pressure_coeff == 1.

        t = 0
        while t < 1.0:
            t = self.net.step()
            assert mdot(t) == approx(mfc.mass_flow_rate)
            dP = self.r1.thermo.P - outlet_reservoir.thermo.P
            assert mdot(t) + pfunc(dP) == approx(pc.mass_flow_rate)

    def test_pressure_controller_type(self):
        self.make_reactors()
        res = ct.Reservoir(self.gas1)
        mfc = ct.MassFlowController(res, self.r1, mdot=0.6)
        p = ct.PressureController(self.r1, self.r2, primary=mfc, K=0.5)
        net = ct.ReactorNet([self.r1, self.r2])  # assigns default names
        assert p.type == "PressureController"
        assert p.name.startswith("PressureController_")  # default name
        p.name = "name-of-pressure-controller"
        assert p.name == "name-of-pressure-controller"

    def test_pressure_controller_errors(self):
        self.make_reactors()
        res = ct.Reservoir(self.gas1)
        mfc = ct.MassFlowController(res, self.r1, mdot=0.6)

        p = ct.PressureController(self.r1, self.r2, primary=mfc, K=0.5)

        with pytest.raises(ct.CanteraError, match='is not ready'):
            p = ct.PressureController(self.r1, self.r2, K=0.5)
            p.mass_flow_rate

        with pytest.raises(ct.CanteraError, match='is not ready'):
            p = ct.PressureController(self.r1, self.r2)
            p.mass_flow_rate

        with pytest.raises(NotImplementedError):
            p = ct.PressureController(self.r1, self.r2)
            p.time_function = lambda t: t>1.

    def test_set_initial_time(self):
        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        pfunc_a = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.pressure_function = pfunc_a

        t0 = 0.0
        tf = t0 + 0.5
        self.net.advance(tf)
        assert self.net.time == approx(tf)
        p1a = self.r1.thermo.P
        p2a = self.r2.thermo.P
        assert valve.pressure_function == approx(pfunc_a(p1a - p2a))

        self.make_reactors(P1=10*ct.one_atm, X1='AR:1.0', X2='O2:1.0')
        self.net.rtol = 1e-12
        valve = ct.Valve(self.r1, self.r2)
        pfunc_b = lambda dP: 5e-3 * np.sqrt(dP) if dP > 0 else 0.0
        valve.pressure_function = pfunc_b

        t0 = 0.2
        self.net.initial_time = t0
        tf = t0 + 0.5
        self.net.advance(tf)
        assert self.net.time == approx(tf)
        p1b = self.r1.thermo.P
        p2b = self.r2.thermo.P
        assert valve.pressure_function == approx(pfunc_b(p1b - p2b))

        assert p1a == approx(p1b)
        assert p2a == approx(p2b)

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

        assert self.r1.T == approx(300)
        assert self.r2.T == approx(1000)
        self.net.advance(2.0)
        T1b = self.r1.T
        T2b = self.r2.T

        assert T1a == approx(T1b)
        assert T2a == approx(T2b)

    def test_unpicklable(self):
        self.make_reactors()
        import pickle
        with pytest.raises(NotImplementedError):
            pickle.dumps(self.r1)
        with pytest.raises(NotImplementedError):
            pickle.dumps(self.net)

    def test_uncopyable(self):
        self.make_reactors()
        import copy
        with pytest.raises(NotImplementedError):
            copy.copy(self.r1)
        with pytest.raises(NotImplementedError):
            copy.copy(self.net)

    def test_invalid_property(self):
        self.make_reactors()
        for x in (self.r1, self.net):
            with pytest.raises(AttributeError):
                x.foobar = 300
            with pytest.raises(AttributeError):
                x.foobar

    def test_bad_kwarg(self):
        g = ct.Solution('h2o2.yaml', transport_model=None)
        self.reactorClass(g, name='ok')
        with pytest.raises(TypeError):
            self.reactorClass(g, foobar=3.14)

    def test_preconditioner_unsupported(self):
        self.make_reactors()
        self.net.preconditioner = ct.AdaptivePreconditioner()
        # initialize should throw an error because the mass fraction
        # reactors do not support preconditioning
        with pytest.raises(ct.CanteraError):
            self.net.initialize()

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_reactor(self):
        self.make_reactors()
        T1, P1, X1 = 300, 101325, 'O2:1.0'
        self.gas1.TPX = T1, P1, X1
        # set attributes during creation
        r1 = self.reactorClass(self.gas1, node_attr={'fillcolor': 'red'})
        r1.name = "Name"
        # overwrite fillcolor in object attributes
        r1.node_attr = {'style': 'filled', 'fillcolor': 'green'}
        graph = r1.draw()
        expected = ['\tName [fillcolor=green style=filled]\n']
        assert graph.body == expected

        # overwrite style during call to draw
        expected = [('\tName [label="{Name|{T (K)\\n300.00|P (bar)\\n1.013}}" '
                     'color=blue fillcolor=green shape=Mrecord style=""]\n')]
        graph = r1.draw(print_state=True, node_attr={"style": "", "color": "blue"})
        assert graph.body == expected

        # print state with mole fractions
        r1.node_attr = {}
        expected = [('\tName [label="{Name|{{T (K)\\n300.00|P (bar)\\n1.013}|X (%)'
                     '\\nO2: 100.00}}" shape=Mrecord]\n')]
        graph = r1.draw(print_state=True, species="X")
        assert graph.body == expected

        # print state with mass fractions
        expected = [('\tName [label="{Name|{{T (K)\\n300.00|P (bar)\\n1.013}|Y (%)'
                     '\\nO2: 100.00}}" shape=Mrecord]\n')]
        graph = r1.draw(print_state=True, species="Y")
        assert graph.body == expected

        # print state with specified species
        expected = [('\tName [label="{Name|{{T (K)\\n300.00|P (bar)\\n1.013}|X (%)'
                     '\\nH2: 0.00}}" shape=Mrecord]\n')]
        graph = r1.draw(print_state=True, species=["H2"])
        assert graph.body == expected

        # print state with specified species and specified unit
        expected = [('\tName [label="{Name|{{T (K)\\n300.00|P (bar)\\n1.013}|X (ppm)'
                     '\\nO2: 1000000.0}}" shape=Mrecord]\n')]
        graph = r1.draw(print_state=True, species=["O2"], species_units="ppm")
        assert graph.body == expected

        # add reactor to existing graph
        graph = _graphviz.Digraph()
        r1.draw(graph)
        expected = ['\tName\n']
        assert graph.body == expected

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_reactors_same_name(self):
        self.make_reactors()
        self.r1.name = 'Reactor'
        self.r2.name = 'Reactor'
        with pytest.raises(AssertionError, match="unique names"):
            self.net.draw()

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_grouped_reactors(self):
        self.make_reactors()
        self.r1.name = "Reactor 1"
        self.r2.name = "Reactor 2"
        self.r1.group_name = "Group 1"
        self.r2.group_name = "Group 2"
        graph = self.net.draw()
        expected = ['\tsubgraph "cluster_Group 1" {\n',
                    '\t\t"Reactor 1"\n',
                    '\t\tlabel="Group 1"\n',
                    '\t}\n',
                    '\tsubgraph "cluster_Group 2" {\n',
                    '\t\t"Reactor 2"\n',
                    '\t\tlabel="Group 2"\n',
                    '\t}\n']
        assert set(graph.body) == set(expected)

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_wall(self):
        T1, P1, X1 = 300, 101325, 'O2:1.0'
        T2, P2, X2 = 600, 101325, 'O2:1.0'
        self.make_reactors(T1=T1, P1=P1, X1=X1, T2=T2, P2=P2, X2=X2)
        self.r1.name = "Name 1"
        self.r2.name = "Name 2"
        w = ct.Wall(self.r1, self.r2, U=0.1, edge_attr={'style': 'dotted'}, name="wall")
        w.edge_attr = {'color': 'green'}
        graph = w.draw(node_attr={'style': 'filled'},
                       edge_attr={'style': 'dashed', 'color': 'blue'})
        expected = [('\t"Name 2" -> "Name 1" [label="wall\\nQ̇ = 30 W" '
                     'color=blue style=dashed]\n')]
        assert graph.body == expected

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_moving_wall(self):
        T1, P1, X1 = 300, 101325, 'O2:1.0'
        T2, P2, X2 = 600, 101325, 'O2:1.0'
        self.make_reactors(T1=T1, P1=P1, X1=X1, T2=T2, P2=P2, X2=X2)
        self.r1.name = "Name 1"
        self.r2.name = "Name 2"
        w = ct.Wall(self.r1, self.r2, U=0.1, velocity=1, name="wall")
        graph = w.draw()
        expected = [('\t"Name 1" -> "Name 2" [label="wall\\nv = 1 m/s" '
                     'arrowhead=icurveteecurve arrowtail=icurveteecurve '
                     'dir=both style=dotted]\n'),
                    ('\t"Name 2" -> "Name 1" [label="wall\\nQ̇ = 30 W" '
                     'color=red style=dashed]\n')]
        assert graph.body == expected

        # omit heat flow if zero
        w.heat_transfer_coeff = 0
        graph = w.draw()
        expected = [('\t"Name 1" -> "Name 2" [label="wall\\nv = 1 m/s" '
                     'arrowhead=icurveteecurve arrowtail=icurveteecurve '
                     'dir=both style=dotted]\n')]
        assert graph.body == expected

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_flow_controller(self):
        self.make_reactors(n_reactors=1)
        self.r1.name = "Reactor"
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TPX = 300, 10*101325, 'H2:1.0'
        inlet_reservoir = ct.Reservoir(gas2, name='Inlet')
        gas2.TPX = 300, 101325, 'H2:1.0'
        outlet_reservoir = ct.Reservoir(gas2, name='Outlet')
        mfc = ct.MassFlowController(inlet_reservoir, self.r1, mdot=2,
                                    edge_attr={'xlabel': 'MFC'}, name="MFC")
        ct.PressureController(self.r1, outlet_reservoir, primary=mfc, name="PC")
        mfc.edge_attr.update({'color': 'purple'})
        self.net.advance_to_steady_state()
        graph = mfc.draw(node_attr={'style': 'filled'},
                         edge_attr={'style': 'dotted', 'color': 'blue'})
        expected = [('\tInlet -> Reactor [label="MFC\\nṁ = 2 kg/s" color=blue '
                     'style=dotted xlabel=MFC]\n')]
        assert graph.body == expected

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_network(self):
        self.make_reactors()
        self.r1.name = "RH"
        self.r2.name = "RC"
        self.add_wall(U=10, name="wall")
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TPX = 600, 101325, 'O2:1.0'
        hot_inlet = ct.Reservoir(gas2, name='InH')
        gas2.TPX = 200, 101325, 'O2:1.0'
        cold_inlet = ct.Reservoir(gas2, name='InC')
        outlet = ct.Reservoir(gas2, name='Out')
        mfc_hot1 = ct.MassFlowController(hot_inlet, self.r1, mdot=1.5, name='mfc_h1')
        mfc_hot2 = ct.MassFlowController(hot_inlet, self.r1, mdot=1, name='mfc_h2')
        mfc_cold = ct.MassFlowController(cold_inlet, self.r2, mdot=2, name='mfc_c')
        ct.PressureController(self.r1, outlet, primary=mfc_hot1, name='pc_h1')
        ct.PressureController(self.r1, outlet, primary=mfc_hot2, name='pc_h2')
        ct.PressureController(self.r2, outlet, primary=mfc_cold, name='pc_c')
        self.net.advance_to_steady_state()
        graph = self.net.draw(mass_flow_attr={'color': 'green'},
                              heat_flow_attr={'color': 'orange'},
                              print_state=True)
        expected = {
            '\tRC [label="{RC|{T (K)\\n202.18|P (bar)\\n1.013}}" shape=Mrecord]\n',
            '\tOut [label="{Out|{T (K)\\n200.00|P (bar)\\n1.013}}" shape=Mrecord]\n',
            '\tInH [label="{InH|{T (K)\\n600.00|P (bar)\\n1.013}}" shape=Mrecord]\n',
            '\tInC [label="{InC|{T (K)\\n200.00|P (bar)\\n1.013}}" shape=Mrecord]\n',
            '\tRH [label="{RH|{T (K)\\n598.42|P (bar)\\n1.013}}" shape=Mrecord]\n',
            '\tInH -> RH [label="mfc_h1\\nṁ = 2.5 kg/s" color=green]\n',
            '\tInH -> RH [label="mfc_h2\\nṁ = 2.5 kg/s" color=green]\n',
            '\tRH -> Out [label="pc_h1\\nṁ = 2.5 kg/s" color=green]\n',
            '\tRH -> Out [label="pc_h2\\nṁ = 2.5 kg/s" color=green]\n',
            '\tRC -> Out [label="pc_c\\nṁ = 2 kg/s" color=green]\n',
            '\tInC -> RC [label="mfc_c\\nṁ = 2 kg/s" color=green]\n',
            '\tRH -> RC [label="wall\\nQ̇ = 4e+03 W" color=orange style=dashed]\n'
        }
        # use sets because order can be random
        # expected defines two alternatives for inH -> RH and RH -> Out as test defines
        # pairs of redundant mass flow controllers and pressure controllers; dependent
        # on type of reactor, either first or second is chosen for the visualization.
        assert not set(graph.body) - expected


class TestMoleReactor(TestReactor):
    reactorClass = ct.MoleReactor

    def test_mole_reactor_surface_chem(self):
        model = "ptcombust.yaml"
        # initial conditions
        T0 = 1500
        P0 = ct.one_atm
        equiv_ratio = 1
        surf_area = 0.527
        fuel = "CH4"
        air = "O2:1.0, N2:3.773"
        # reactor 1
        gas1 = ct.Solution(model, transport_model=None)
        gas1.TP = T0, P0
        gas1.set_equivalence_ratio(equiv_ratio, fuel, air)
        r1 = self.reactorClass(gas1)
        # comparison reactor
        gas2 = ct.Solution(model, transport_model=None)
        gas2.TP = T0, P0
        gas2.set_equivalence_ratio(equiv_ratio, fuel, air)
        if "ConstPressure" in r1.type:
            r2 = ct.ConstPressureReactor(gas2)
        else:
            r2 = ct.Reactor(gas2)
        # surf 1
        surf1 = ct.Interface(model, "Pt_surf", [gas1])
        surf1.TP = T0, P0
        surf1.coverages = {"PT(S)":1}
        rsurf1 = ct.ReactorSurface(surf1, r1, A=surf_area)
        # surf 2
        surf2 = ct.Interface(model, "Pt_surf", [gas2])
        surf2.TP = T0, P0
        surf2.coverages = {"PT(S)":1}
        rsurf2 = ct.ReactorSurface(surf2, r2, A=surf_area)
        # reactor network setup
        net1 = ct.ReactorNet([r1,])
        net2 = ct.ReactorNet([r2,])
        net1.rtol = net2.rtol = 1e-9
        net1.atol = net2.atol = 1e-18
        # steady state occurs at ~0.002 seconds
        for i in np.linspace(0, 0.0025, 50)[1:]:
            net1.advance(i)
            net2.advance(i)
            assert r1.thermo.Y == approx(r2.thermo.Y, rel=5e-4, abs=1e-6)
            assert r1.T == approx(r2.T, rel=5e-5)
            assert r1.thermo.P == approx(r2.thermo.P, rel=1e-6)
            assert rsurf1.coverages == approx(rsurf2.coverages, rel=1e-4, abs=1e-8)

    def test_tolerances(self, rtol_lim=1e-8, atol_lim=1e-18):
        super().test_tolerances(rtol_lim, atol_lim)

class TestIdealGasReactor(TestReactor):
    reactorClass = ct.IdealGasReactor


class TestWellStirredReactorIgnition:
    """ Ignition (or not) of a well-stirred reactor """

    def setup_reactor(self, T0, P0, mdot_fuel, mdot_ox):
        """ Runs before tests """
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
        fuel_in = ct.Reservoir(self.gas)

        # oxidizer inlet
        self.gas.TPX = T0, P0, "N2:3.76, O2:1.0"
        oxidizer_in = ct.Reservoir(self.gas)

        # reactor, initially filled with N2
        self.gas.TPX = T0, P0, "N2:1.0"
        self.combustor = ct.IdealGasReactor(self.gas)
        self.combustor.volume = 1.0

        # outlet
        exhaust = ct.Reservoir(self.gas)

        # connect the reactor to the reservoirs
        fuel_mfc = ct.MassFlowController(fuel_in, self.combustor)
        fuel_mfc.mass_flow_rate = mdot_fuel
        oxidizer_mfc = ct.MassFlowController(oxidizer_in, self.combustor)
        oxidizer_mfc.mass_flow_rate = mdot_ox
        valve = ct.Valve(self.combustor, exhaust)
        valve.valve_coeff = 1.0

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
        self.setup_reactor(T0, 10*ct.one_atm, mdot_f, mdot_o)
        self.gas.set_multiplier(0.0)
        t,T = self.integrate(100.0)

        for i in range(len(t)):
            assert T[i] == approx(T0, rel=1e-5)

        assert self.combustor.thermo['CH4'].Y == approx(mdot_f / (mdot_o + mdot_f))

    def test_ignition1(self):
        self.setup_reactor(900.0, 10*ct.one_atm, 1.0, 5.0)
        t,T = self.integrate(10.0)

        assert T[-1] > 1200 # mixture ignited
        for i in range(len(t)):
            if T[i] > 0.5 * (T[0] + T[-1]):
                tIg = t[i]
                break

        # regression test; no external basis for this result
        assert tIg == approx(2.2249, rel=1e-3)

    def test_ignition2(self):
        self.setup_reactor(900.0, 10*ct.one_atm, 1.0, 20.0)
        t,T = self.integrate(10.0)

        assert T[-1] > 1200 # mixture ignited
        for i in range(len(t)):
            if T[i] > 0.5 * (T[0] + T[-1]):
                tIg = t[i]
                break

        # regression test; no external basis for this result
        assert tIg == approx(1.4856, rel=1e-3)

    def test_ignition3(self):
        self.setup_reactor(900.0, 10*ct.one_atm, 1.0, 80.0)
        self.net.max_time_step = 0.5
        t,T = self.integrate(100.0)
        assert T[-1] < 910 # mixture did not ignite

    def test_steady_state(self):
        self.setup_reactor(900.0, 10*ct.one_atm, 1.0, 20.0)
        residuals = self.net.advance_to_steady_state(return_residuals=True)
        # test if steady state is reached
        assert residuals[-1] < 10. * self.net.rtol
        # regression test; no external basis for these results
        assert self.combustor.T == approx(2498.94, rel=1e-5)
        assert self.combustor.thermo['H2O'].Y[0] == approx(0.103658, rel=1e-5)
        assert self.combustor.thermo['HO2'].Y[0] == approx(8.734515e-06, rel=1e-5)


class TestConstPressureReactor:
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

    def test_reactor_surface_type(self):
        self.create_reactors(add_surf=True)
        assert self.surf1.type == "ReactorSurface"
        assert self.surf1.name.startswith("ReactorSurface_")  # default name
        self.surf1.name = "name-of-reactor-surface"
        assert self.surf1.name == "name-of-reactor-surface"

    def test_component_index(self):
        self.create_reactors(add_surf=True)
        for (gas,net,iface,r) in ((self.gas1, self.net1, self.interface1, self.r1),
                                  (self.gas2, self.net2, self.interface2, self.r2)):
            net.step()

            N0 = net.n_vars - gas.n_species - iface.n_species
            N1 = net.n_vars - iface.n_species
            for i, name in enumerate(gas.species_names):
                assert i + N0 == r.component_index(name)
            for i, name in enumerate(iface.species_names):
                assert i + N1 == r.component_index(name)

    def test_component_names(self):
        self.create_reactors(add_surf=True)
        for i in range(self.net1.n_vars):
            assert self.r1.component_index(self.r1.component_name(i)) == i
            assert self.net1.component_name(i) == '{}: {}'.format(self.r1.name,
                                                                  self.r1.component_name(i))

    def integrate(self, surf=False):
        for t in np.arange(0.5, 50, 1.0):
            self.net1.advance(t)
            self.net2.advance(t)
            assert self.r1.thermo.Y == approx(self.r2.thermo.Y, rel=5e-4, abs=1e-6)
            assert self.r1.T == approx(self.r2.T, rel=5e-5)
            assert self.r1.thermo.P == approx(self.r2.thermo.P, rel=1e-6)
            if surf:
                assert self.surf1.coverages == approx(self.surf2.coverages,
                                                      rel=1e-4, abs=1e-8)

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

    def test_preconditioner_unsupported(self):
        self.create_reactors()
        self.net2.preconditioner = ct.AdaptivePreconditioner()
        # initialize should throw an error because the mass fraction
        # reactors do not support preconditioning
        with pytest.raises(ct.CanteraError):
            self.net2.initialize()

class TestConstPressureMoleReactor(TestConstPressureReactor):
    """
    The constant pressure reactor should give the same results as
    as a regular "Reactor" with a wall with a very high expansion rate
    coefficient.
    """
    reactorClass = ct.ConstPressureMoleReactor
    test_mole_reactor_surface_chem = TestMoleReactor.test_mole_reactor_surface_chem


class TestIdealGasConstPressureReactor(TestConstPressureReactor):
    reactorClass = ct.IdealGasConstPressureReactor


class TestIdealGasConstPressureMoleReactor(TestConstPressureMoleReactor):
    reactorClass = ct.IdealGasConstPressureMoleReactor
    test_preconditioner_unsupported = None

    def create_reactors(self, **kwargs):
        super().create_reactors(**kwargs)
        self.precon = ct.AdaptivePreconditioner()
        self.net2.preconditioner = self.precon
        self.net2.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True,
            "skip-coverage-dependence":True}

    def test_get_solver_type(self):
        self.create_reactors()
        assert self.precon.side == "right"
        self.net2.initialize()
        assert self.net2.linear_solver_type == "GMRES"


class TestIdealGasMoleReactor(TestMoleReactor):
    reactorClass = ct.IdealGasMoleReactor
    test_preconditioner_unsupported = None

    def test_adaptive_precon_integration(self):
        # Network one with non-mole reactor
        net1 = ct.ReactorNet()
        T0 = 900
        P0 = ct.one_atm
        gas1 = ct.Solution("gri30.yaml")
        gas1.TP = T0, P0
        gas1.set_equivalence_ratio(1, "CH4", "O2:1, N2:3.76")
        r1 = ct.IdealGasMoleReactor(gas1)
        net1.add_reactor(r1)
        # Network two with mole reactor and preconditioner
        net2 = ct.ReactorNet()
        gas2 = ct.Solution("gri30.yaml")
        gas2.TP = T0, P0
        gas2.set_equivalence_ratio(1, "CH4", "O2:1, N2:3.76")
        r2 = ct.IdealGasMoleReactor(gas2)
        net2.add_reactor(r2)
        # add preconditioner
        net2.preconditioner = ct.AdaptivePreconditioner()
        net2.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
        # tolerances
        net1.atol = net2.atol = 1e-16
        net1.rtol = net1.rtol = 1e-8
        # integrate
        for t in np.arange(0.5, 5, 0.5):
            net1.advance(t)
            net2.advance(t)
            assert r1.thermo.Y == approx(r2.thermo.Y, rel=5e-4, abs=1e-6)
            assert r1.T == approx(r2.T, rel=1e-5)
            assert r1.thermo.P == approx(r2.thermo.P, rel=1e-5)


class TestReactorJacobians:

    def test_multi_surface_simple(self):
        # conditions for simulation
        yml = "simple_surface.yaml"
        fuel = "A:1.0, B:1.0"
        # gas kinetics
        gas = ct.Solution(yml, "gas")
        gas.TPX = 1000, 2e5, fuel
        gas.set_multiplier(0)
        # surface kinetics for the simulation
        surf = ct.Interface(yml, 'surf', [gas])
        surf2 = ct.Interface(yml, 'surf', [gas])
        surf.coverages = 'A(S):0.1, B(S):0.2, C(S):0.3, D(S):0.2, (S):0.2'
        surf2.coverages = 'A(S):0.1, D(S):0.2, (S):0.2'
        # create reactor
        r = ct.IdealGasMoleReactor(gas)
        r.volume = 3
        # create surfaces
        rsurf1 = ct.ReactorSurface(surf, r, A=9e-4)
        rsurf2 = ct.ReactorSurface(surf2, r, A=5e-4)
        # create network
        net = ct.ReactorNet([r])
        net.step()
        # get jacobians
        jacobian = r.jacobian
        fd_jacobian = r.finite_difference_jacobian
        # the volume row is not considered in comparisons because it is presently
        # not calculated.
        # check first row is near, terms which are generally on the order of 1e5 to 1e7
        assert jacobian[0, 2:] == approx(fd_jacobian[0, 2:], rel=1e-1, abs=1e-2)
        # check first col is near, these are finite difference terms and should be close
        assert jacobian[2:, 0] == approx(fd_jacobian[2:, 0], rel=1e-3, abs=1e-4)
        # check all species are near, these terms are usually ~ 1e2
        assert jacobian[2:, 2:] == approx(fd_jacobian[2:, 2:], rel=1e-3, abs=1e-4)

    def test_gas_simple(self):
        # conditions for simulation
        yml = "simple_surface.yaml"
        fuel = "A:1.0, B:1.0, C:1.0, D:1.0"
        # gas kinetics
        gas = ct.Solution(yml, "gas")
        gas.TPX = 1000, 2e5, fuel
        # create reactor
        r = ct.IdealGasMoleReactor(gas)
        r.volume = 1
        # create network
        net = ct.ReactorNet([r])
        net.initialize()
        # compare jacobians
        assert r.jacobian == approx(r.finite_difference_jacobian, rel=1e-6)

    def test_const_volume_hydrogen_single(self):
        # conditions for simulation
        yml = "h2o2.yaml"
        # gas kinetics
        gas = ct.Solution(yml, "ohmech")
        gas.TPX = 1000, ct.one_atm, "H2:1.0, H:2.0, H2O:2.0"
        gas.set_multiplier(0)
        gas.set_multiplier(1, 14)
        # create reactor
        r = ct.IdealGasMoleReactor(gas)
        r.volume = 2
        # create network
        net = ct.ReactorNet([r])
        net.step()
        # get jacobians
        jacobian = r.jacobian
        fd_jacobian = r.finite_difference_jacobian
        # the volume row is not considered in comparisons because it is presently
        # not calculated.
        # check first row is near, terms which are generally on the order of 1e5 to 1e7
        assert jacobian[0, 2:] == approx(fd_jacobian[0, 2:], rel=1e-2, abs=1e-3)
        # check first col is near, these are finite difference terms and should be close
        assert jacobian[2:, 0] == approx(fd_jacobian[2:, 0], rel=1e-3, abs=1e-4)
        # check all species are near, these terms are usually ~ 1e2
        assert jacobian[2:, 2:] == approx(fd_jacobian[2:, 2:], rel=1e-3, abs=1e-4)

    def test_const_pressure_hydrogen_single(self):
        # conditions for simulation
        yml = "h2o2.yaml"
        # gas kinetics
        gas = ct.Solution(yml, "ohmech")
        gas.TPX = 1000, ct.one_atm, "H2:1.0, H:2.0, H2O:2.0"
        gas.set_multiplier(0)
        gas.set_multiplier(1, 14)
        # create reactor
        r = ct.IdealGasConstPressureMoleReactor(gas)
        r.volume = 2
        # create network
        net = ct.ReactorNet([r])
        net.step()
        # compare analytical jacobian with finite difference
        assert r.jacobian == approx(r.finite_difference_jacobian, rel=1e-3, abs=1e-4)

    def test_phase_order_surf_jacobian(self):
        # create gas phase
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
        gas = ct.Solution(yaml=gas_def)
        # set gas phase conditions
        T0 = 1200
        P0 = 25*ct.one_atm
        X0 = 'CH4:0.5, H2O:0.2, CO:0.3'
        gas.TPX = T0, P0, X0
        # create reactors
        r1 = ct.IdealGasMoleReactor(gas)
        r2 = ct.IdealGasMoleReactor(gas)
        r1.volume = 0.1
        r2.volume = 0.1
        # create solid and interfaces
        solid = ct.Solution('diamond.yaml', 'diamond')
        # first interface
        interface1 = ct.Interface('diamond.yaml', 'diamond_100', (gas, solid))
        # interface with reversed phase order
        interface2 = ct.Interface('diamond.yaml', 'diamond_100', (solid, gas))
        # creating initial coverages
        C = np.zeros(interface1.n_species)
        C[0] = 0.3
        C[4] = 0.7
        # creating reactor surfaces
        surf1 = ct.ReactorSurface(interface1, A=1e-3)
        surf2 = ct.ReactorSurface(interface2, A=1e-3)
        surf1.coverages = C
        surf2.coverages = C
        surf1.install(r1)
        surf2.install(r2)
        # create reactor network
        net = ct.ReactorNet([r1, r2])
        # set derivative settings
        net.derivative_settings = {"skip-coverage-dependence":True}
        net.initialize()
        # check that they two arrays are the same
        assert r1.jacobian == approx(r2.jacobian, rel=2e-6, abs=1e-8)

# A rate type used for testing integrator error handling
class FailRateData(ct.ExtensibleRateData):
    def __init__(self):
        self.fail = False

    def update(self, gas):
        self.T = gas.T
        if self.T < 1400:
            self.fail = True
        return True

@ct.extension(name="fail-rate", data=FailRateData)
class FailRate(ct.ExtensibleRate):
    def __init__(self, *args, recoverable, **kwargs):
        super().__init__(*args, **kwargs)
        self.recoverable = recoverable
        self.count = 0

    def eval(self, data):
        if data.fail:
            self.count += 1
            if self.count < 3 or not self.recoverable:
                raise ValueError("spam")
        return 0.0


class TestFlowReactor:
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

        x = 0
        v0 = r.speed
        assert v0 == approx(10 / r.density)
        while x < 10.0:
            x = net.step()
            assert v0 == approx(r.speed)

    def test_reacting(self):
        g = ct.Solution(yaml=self.gas_def)
        g.TPX = 1400, 20*101325, 'CO:1.0, H2O:1.0'

        r = ct.FlowReactor(g)
        r.mass_flow_rate = 10

        net = ct.ReactorNet()
        net.add_reactor(r)

        i = 0
        while net.distance < 1.0:
            net.step()
            i += 1
            assert r.speed * r.density * r.area == approx(10)

        stats = net.solver_stats
        assert stats['steps'] == i
        assert 'err_tests_fails' in stats

        # advancing to the current distance should be a no-op
        x_now = net.distance
        net.advance(x_now)
        assert net.solver_stats['steps'] == i

    def test_catalytic_surface(self):
        # Regression test based roughly on surf_pfr.py
        T0 = 1073.15
        P0 = ct.one_atm
        X0 = 'CH4:1, O2:1.5, AR:0.1'

        surf = ct.Interface('methane_pox_on_pt.yaml', 'Pt_surf')
        gas = surf.adjacent['gas']
        gas.TPX = T0, P0, X0
        surf.TP = T0, P0

        r = ct.FlowReactor(gas)
        r.area = 1e-4
        porosity = 0.3
        velocity = 0.4 / 60
        mdot = velocity * gas.density * r.area * porosity
        r.surface_area_to_volume_ratio = porosity * 1e5
        r.mass_flow_rate = mdot
        r.energy_enabled = False

        rsurf = ct.ReactorSurface(surf, r)

        sim = ct.ReactorNet([r])
        kCH4 = gas.species_index('CH4')
        kH2 = gas.species_index('H2')
        kCO = gas.species_index('CO')

        sim.advance(1e-7)
        X = r.thermo['CH4', 'H2', 'CO'].X
        assert X == approx([0.10578801, 0.001654415, 0.012103974])
        assert r.thermo.density * r.area * r.speed == approx(mdot)

        sim.advance(1e-5)
        X = r.thermo['CH4', 'H2', 'CO'].X
        assert X == approx([0.07748481, 0.048165072, 0.01446654])
        assert r.thermo.density * r.area * r.speed == approx(mdot)

        sim.advance(1e-3)
        X = r.thermo['CH4', 'H2', 'CO'].X
        assert X == approx([0.01815402, 0.21603645, 0.045640395])
        assert r.thermo.density * r.area * r.speed == approx(mdot)

    def test_component_names(self):
        surf = ct.Interface('methane_pox_on_pt.yaml', 'Pt_surf')
        gas = surf.adjacent['gas']
        r = ct.FlowReactor(gas)
        r.mass_flow_rate = 0.1
        rsurf = ct.ReactorSurface(surf, r)
        sim = ct.ReactorNet([r])
        sim.initialize()

        assert r.n_vars == 4 + gas.n_species + surf.n_species
        assert sim.n_vars == r.n_vars

        for i in range(r.n_vars):
            name = r.component_name(i)
            assert r.component_index(name) == i

        with pytest.raises(IndexError, match="No such component: 'spam'"):
            r.component_index('spam')

        with pytest.raises(ct.CanteraError, match='out of bounds'):
            r.component_name(200)

class TestFlowReactor2:
    def import_phases(self):
        surf = ct.Interface('SiF4_NH3_mec.yaml', 'SI3N4')
        return surf, surf.adjacent['gas']

    def make_reactors(self, gas, surf):
        r = ct.FlowReactor(gas)
        r.area = 1e-4
        r.surface_area_to_volume_ratio = 5000
        r.mass_flow_rate = 0.02
        rsurf = ct.ReactorSurface(surf, r)
        sim = ct.ReactorNet([r])
        return r, rsurf, sim

    def test_advance_reverse(self):
        surf, gas = self.import_phases()
        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        r, rsurf, sim = self.make_reactors(gas, surf)

        sim.advance(0.1)
        with pytest.raises(ct.CanteraError, match="backwards in time"):
            sim.advance(0.09)

    def test_no_mass_flow_rate(self):
        surf, gas = self.import_phases()
        r = ct.FlowReactor(gas)
        rsurf = ct.ReactorSurface(surf, r)
        sim = ct.ReactorNet([r])
        with pytest.raises(ct.CanteraError, match="mass flow rate"):
            sim.initialize()

    def test_mixed_reactor_types(self):
        surf, gas = self.import_phases()
        r1 = ct.FlowReactor(gas)
        r2 = ct.IdealGasReactor(gas)
        with pytest.raises(ct.CanteraError, match="Cannot mix Reactor types"):
            ct.ReactorNet([r1, r2])

    def test_unrecoverable_integrator_errors(self):
        surf, gas = self.import_phases()

        # To cause integrator failures, add a reaction that will fail under
        # certain conditions (T < 1400)
        fail = ct.Reaction(equation='NH3 <=> NH3', rate=FailRate(recoverable=False))
        gas.add_reaction(fail)

        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP
        r, rsurf, sim = self.make_reactors(gas, surf)

        with pytest.raises(ct.CanteraError,
                           match="repeated recoverable residual errors"):
            while r.thermo.T > 1300:
                sim.step()

    def test_integrator_errors_advance(self):
        surf, gas = self.import_phases()

        # To cause integrator failures, add a reaction that will fail under
        # certain conditions (T < 1400)
        fail = ct.Reaction(equation='NH3 <=> NH3', rate=FailRate(recoverable=False))
        gas.add_reaction(fail)

        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP
        r, rsurf, sim = self.make_reactors(gas, surf)

        with pytest.raises(ct.CanteraError,
                           match="repeated recoverable residual errors"):
            while r.thermo.T > 1300:
                sim.advance(sim.distance + 0.01)

    def test_recoverable_integrator_errors(self):
        surf, gas = self.import_phases()

        # Test integrator behavior on "recoverable" errors that are resolved by
        # calling taking a different step size
        fail = ct.Reaction(equation='NH3 <=> NH3', rate=FailRate(recoverable=True))
        gas.add_reaction(fail)

        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP
        r, rsurf, sim = self.make_reactors(gas, surf)

        while r.thermo.T > 1300:
            sim.step()

        # At least some "recoverable" errors occurred
        assert fail.rate.count > 0

    def test_max_steps(self):
        surf, gas = self.import_phases()
        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP
        r, rsurf, sim = self.make_reactors(gas, surf)

        sim.max_steps = 13
        assert sim.max_steps == 13

        with pytest.raises(ct.CanteraError, match="Maximum number of timesteps"):
            sim.advance(0.1)

        assert sim.solver_stats['steps'] == 13

    def test_independent_variable(self):
        surf, gas = self.import_phases()
        r, rsurf, sim = self.make_reactors(gas, surf)
        with pytest.raises(ct.CanteraError, match="independent variable"):
            sim.time

        assert sim.distance == 0.0

    def test_max_time_step(self):
        surf, gas = self.import_phases()
        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP
        r, rsurf, sim = self.make_reactors(gas, surf)

        sim.max_time_step = 0.002
        sim.advance(0.01)
        sim.step()
        x1 = sim.distance
        x2 = sim.step()
        dx_limit = 0.1 * (x2-x1)

        sim.max_time_step = dx_limit
        assert sim.max_time_step == dx_limit

        # Setting a step size limit seems to take one additional step before it's
        # fully enforced
        sim.step()

        for i in range(20):
            tPrev = sim.distance
            tNow = sim.step()
            assert tNow - tPrev <= 1.0001 * dx_limit

    def test_tolerances(self):
        surf, gas = self.import_phases()
        gas.TPX = 1500, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP

        r0, rsurf0, sim0 = self.make_reactors(gas, surf)
        r1, rsurf1, sim1 = self.make_reactors(gas, surf)
        r2, rsurf2, sim2 = self.make_reactors(gas, surf)

        sim0.advance(0.05)
        baseline = sim0.solver_stats

        # Expect that satisfying tighter tolerances will require more time steps
        sim1.atol = 0.001 * sim0.atol
        sim1.advance(0.05)
        tight_atol = sim1.solver_stats
        assert tight_atol['steps'] > baseline['steps']

        sim2.rtol = 0.001 * sim0.rtol
        sim2.advance(0.05)
        tight_rtol = sim2.solver_stats
        assert tight_rtol['steps'] > baseline['steps']

    def test_iteration_limits(self):
        surf, gas = self.import_phases()
        gas.TPX = 1700, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP

        r, rsurf, sim = self.make_reactors(gas, surf)

        # Too restrictive limits should cause integration errors:
        sim.advance(0.1)
        sim.max_nonlinear_iterations = 1
        sim.max_nonlinear_convergence_failures = 1
        sim.include_algebraic_in_error_test = True
        sim.max_err_test_fails = 2
        sim.rtol = 1e-12
        with pytest.raises(ct.CanteraError, match="corrector convergence"):
            sim.advance(0.2)

    def test_solver_order(self):
        surf, gas = self.import_phases()
        gas.TPX = 1700, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP

        r, rsurf, sim = self.make_reactors(gas, surf)

        sim.max_order = 7  # Invalid, will be caught later
        with pytest.raises(ct.CanteraError, match="IDA_ILL_INPUT"):
            sim.initialize()

        sim.max_order = 2
        sim.advance(0.1)
        assert sim.solver_stats['last_order'] == 2

        sim.max_order = 4
        sim.advance(0.4)
        assert sim.solver_stats['last_order'] == 4

        with pytest.raises(ct.CanteraError, match="IDA_ILL_INPUT"):
            sim.max_order = -1

    def test_reinitialization(self):
        surf, gas = self.import_phases()
        gas.TPX = 1700, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP

        r, rsurf, sim = self.make_reactors(gas, surf)
        r.mass_flow_rate = 0.01
        sim.advance(0.6)
        Ygas1 = r.thermo.Y
        cov1 = rsurf.kinetics.coverages

        # Reset the reactor to the same initial state
        gas.TPX = 1700, 4000, 'NH3:1.0, SiF4:0.4'
        surf.TP = gas.TP
        r.mass_flow_rate = 0.01
        r.syncState()

        sim.initial_time = 0.
        sim.advance(0.6)
        Ygas2 = r.thermo.Y
        cov2 = rsurf.kinetics.coverages

        assert Ygas1 == approx(Ygas2)
        assert cov1 == approx(cov2)

    def test_initial_condition_tolerances(self):
        surf, gas = self.import_phases()
        gas.TPX = 1700, 4000, 'NH3:1.0, SiF4:0.4'
        surf.coverages = np.ones(surf.n_species)
        surf.TP = gas.TP
        r, rsurf, sim = self.make_reactors(gas, surf)

        # With tight tolerances, some error test failures should be expected
        r.inlet_surface_atol = 1e-14
        r.inlet_surface_rtol = 1e-20
        r.inlet_surface_max_error_failures = 1
        with pytest.raises(ct.CanteraError, match="error test failed repeatedly"):
            sim.initialize()

        # With few steps allowed, won't be able to reach steady-state
        r.inlet_surface_max_error_failures = 10
        r.inlet_surface_max_steps = 200
        with pytest.raises(ct.CanteraError, match="Maximum number of timesteps"):
            sim.initialize()

        # Relaxing the tolerances will allow the integrator to reach steady-state
        # in fewer timesteps
        surf.coverages = np.ones(surf.n_species)
        r.inlet_surface_atol = 0.001
        r.inlet_surface_rtol = 0.001
        sim.initialize()


class TestSurfaceKinetics:
    def make_reactors(self):
        self.net = ct.ReactorNet()

        self.interface = ct.Interface('diamond.yaml', 'diamond_100')
        self.gas = self.interface.adjacent['gas']
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
        assert surf1.coverages[0] == approx(0.3)
        assert surf1.coverages[1] == approx(0.0)
        assert surf1.coverages[4] == approx(0.7)
        self.net.advance(1e-5)
        C_left = surf1.coverages

        self.make_reactors()
        surf2 = ct.ReactorSurface(self.interface, self.r2)
        surf2.coverages = 'c6HH:0.3, c6HM:0.7'
        assert surf2.coverages[0] == approx(0.3)
        assert surf2.coverages[4] == approx(0.7)
        self.net.advance(1e-5)
        C_right = surf2.coverages

        assert sum(C_left) == approx(1.0)
        assert C_left == approx(C_right)

        with pytest.raises(ValueError):
            surf2.coverages = np.ones(self.interface.n_species + 1)

    def test_coverages_regression1(self, test_data_path):
        # Test with energy equation disabled
        self.make_reactors()
        self.r1.energy_enabled = False
        self.r2.energy_enabled = False
        surf1 = ct.ReactorSurface(self.interface, self.r1)

        C = np.zeros(self.interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        surf1.coverages = C
        assert surf1.coverages == approx(C)
        data = []
        test_file = self.test_work_path / "test_coverages_regression1.csv"
        reference_file = test_data_path / "WallKinetics-coverages-regression1.csv"
        data = []
        for t in np.linspace(1e-6, 1e-3):
            self.net.advance(t)
            data.append([t, self.r1.T, self.r1.thermo.P, self.r1.mass] +
                        list(self.r1.thermo.X) + list(surf1.coverages))
        np.savetxt(test_file, data, delimiter=',')

        bad = compareProfiles(reference_file, test_file,
                                        rtol=1e-5, atol=1e-9, xtol=1e-12)
        assert not bool(bad), bad

    def test_coverages_regression2(self, test_data_path):
        # Test with energy equation enabled
        self.make_reactors()
        surf = ct.ReactorSurface(self.interface, self.r1)

        C = np.zeros(self.interface.n_species)
        C[0] = 0.3
        C[4] = 0.7

        surf.coverages = C
        assert surf.coverages == approx(C)
        data = []
        test_file = self.test_work_path / "test_coverages_regression2.csv"
        reference_file = test_data_path / "WallKinetics-coverages-regression2.csv"
        data = []
        for t in np.linspace(1e-6, 1e-3):
            self.net.advance(t)
            data.append([t, self.r1.T, self.r1.thermo.P, self.r1.mass] +
                        list(self.r1.thermo.X) + list(surf.coverages))
        np.savetxt(test_file, data, delimiter=',')

        bad = compareProfiles(reference_file, test_file,
                                        rtol=1e-5, atol=1e-9, xtol=1e-12)
        assert not bool(bad), bad

    @pytest.mark.skipif(_graphviz is None, reason="graphviz is not installed")
    def test_draw_ReactorSurface(self):
        self.make_reactors()
        surf = ct.ReactorSurface(self.interface, self.r1)
        self.r1.name = "Reactor"

        graph = surf.draw(node_attr={'style': 'filled'},
                          surface_edge_attr={'color': 'red'}, print_state=True)
        expected = [('\tReactor [label="{Reactor|{T (K)\\n1200.00|P (bar)\\n0.010}}" '
                     'shape=Mrecord style=filled]\n'),
                    '\t"Reactor surface" [shape=underline style=filled]\n',
                    ('\tReactor -> "Reactor surface" '
                     '[arrowhead=none color=red style=dotted]\n')]
        assert graph.body == expected


class TestReactorSensitivities:
    def test_sensitivities1(self):
        net = ct.ReactorNet()
        gas = ct.Solution('gri30.yaml', transport_model=None)
        gas.TPX = 1300, 20*101325, 'CO:1.0, H2:0.1, CH4:0.1, H2O:0.5'
        r1 = ct.IdealGasReactor(gas)
        net.add_reactor(r1)

        assert net.n_sensitivity_params == 0
        r1.add_sensitivity_reaction(40)
        r1.add_sensitivity_reaction(41)

        net.advance(0.1)

        assert net.n_sensitivity_params == 2
        assert net.n_vars == gas.n_species + r1.component_index(gas.species_name(0))
        S = net.sensitivities()
        assert S.shape == (net.n_vars, net.n_sensitivity_params)

    def test_sensitivities2(self):
        net = ct.ReactorNet()

        interface = ct.Interface("diamond.yaml", "diamond_100")
        gas1 = interface.adjacent["gas"]
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
            assert S[1,:] == approx(np.zeros(2))
            assert S[K2+1,:] == approx(np.zeros(2))

            # Sensitivity coefficients for the disjoint reactors should be zero
            assert np.linalg.norm(S[Ns:K2,1]) == approx(0.0, abs=1e-5)
            assert np.linalg.norm(S[K2+Ns:,0]) == approx(0.0, abs=1e-5)

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
                assert reactor.name == rname
                if kind == 'r':
                    assert gas.reaction(p).equation == comp
                elif kind == 's':
                    assert p + ' enthalpy' == comp

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
            assert S1[:,i] == approx(S2[:,j])

    def test_parameter_order1a(self):
        self._test_parameter_order1(ct.IdealGasReactor)

    @pytest.mark.slow_test
    def test_parameter_order1b(self):
        self._test_parameter_order1(ct.IdealGasConstPressureReactor)

    @pytest.mark.slow_test
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
                assert pname(r,p) == net1.sensitivity_parameter_name(i)

            rA2,rB2,net2 = setup(reverse)
            params2 = [(rB2,10),(rA2,19),(rB2,18),(rA2,2)]
            for r,p in params2:
                r.add_sensitivity_reaction(p)
            S.append(integrate(rA2, net2))

            for i,(r,p) in enumerate(params2):
                assert pname(r,p) == net2.sensitivity_parameter_name(i)

        # Check that the results reflect the changed parameter ordering
        for a,b in ((0,1), (2,3)):
            for i,j in enumerate((3,1,0,2)):
                assert S[a][:,i] == approx(S[b][:,j])

        # Check that results are consistent after changing the order that
        # reactors are added to the network
        N = gas.n_species + r.component_index(gas.species_name(0))
        assert S[0][:N] == approx(S[2][N:], rel=1e-5, abs=1e-5)
        assert S[0][N:] == approx(S[2][:N], rel=1e-5, abs=1e-5)
        assert S[1][:N] == approx(S[3][N:], rel=1e-5, abs=1e-5)
        assert S[1][N:] == approx(S[3][:N], rel=1e-5, abs=1e-5)

    @pytest.mark.slow_test
    def test_parameter_order3(self):
        # Test including reacting surfaces
        interface = ct.Interface("diamond.yaml", "diamond_100")
        gas1 = interface.adjacent["gas"]

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
                assert S[a][:,i] == approx(S[b][:,j], rel=1e-2, abs=1e-3)

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

    # @todo: replace np.trapz with np.trapezoid when dropping support for NumPy 1.x
    @pytest.mark.skip(reason="Integration of sensitivity ODEs is unreliable, "
                              "see: https://github.com/Cantera/enhancements/issues/55")
    @pytest.mark.filterwarnings("ignore:`trapz` is deprecated")
    def test_ignition_delay_sensitivity(self):
        species = ('H2', 'H', 'O2', 'H2O2', 'H2O', 'OH', 'HO2')
        dtigdh_cvodes = self.calc_dtdh(species)
        tig0 = self.calc_tig('H2', 0)
        dH = 1e4
        for i,s in enumerate(species):
            dtigdh = (self.calc_tig(s, dH) - tig0) / dH
            assert dtigdh_cvodes[i] == approx(dtigdh, rel=5e-2, abs=1e-14)


class TestCombustor:
    """
    These tests are based on the sample:

        samples/python/reactors/combustor.py

    with some simplifications so that they run faster and produce more
    consistent output.

    Note: to re-create the reference file:
    (1) set PYTHONPATH to build/python.
    (2) go into test/python directory and run:
        pytest --save-reference=combustor test_reactor.py::TestCombustor::test_integrateWithAdvance
    (3) Compare the reference files created in the current working directory with
        the ones in test/data and replace them if needed.
    """
    referenceFile = "CombustorTest-integrateWithAdvance.csv"

    @pytest.fixture
    def setup_combustor_tests(self):
        gas = ct.Solution('h2o2.yaml', transport_model=None)

        # create a reservoir for the fuel inlet, and set to pure methane.
        gas.TPX = 300.0, ct.one_atm, 'H2:1.0'
        fuel_in = ct.Reservoir(gas)
        fuel_mw = gas.mean_molecular_weight

        # Oxidizer inlet
        gas.TPX = 300.0, ct.one_atm, 'O2:1.0, AR:3.0'
        oxidizer_in = ct.Reservoir(gas)
        oxidizer_mw = gas.mean_molecular_weight

        # to ignite the fuel/air mixture, we'll introduce a pulse of radicals.
        # The steady-state behavior is independent of how we do this, so we'll
        # just use a stream of pure atomic hydrogen.
        gas.TPX = 300.0, ct.one_atm, 'H:1.0'
        igniter = ct.Reservoir(gas)

        # create the combustor, and fill it with a diluent
        gas.TPX = 300.0, ct.one_atm, 'AR:1.0'
        combustor = ct.IdealGasReactor(gas)

        # create a reservoir for the exhaust
        exhaust = ct.Reservoir(gas)

        # compute fuel and air mass flow rates
        factor = 0.1
        oxidizer_mdot = 4 * factor * oxidizer_mw
        fuel_mdot = factor * fuel_mw

        # The igniter will use a time-dependent igniter mass flow rate.
        def igniter_mdot(t, t0=0.1, fwhm=0.05, amplitude=0.1):
            return amplitude * math.exp(-(t - t0) ** 2 * 4 * math.log(2) / fwhm ** 2)

        # create and install the mass flow controllers. Controllers m1 and m2 provide
        # constant mass flow rates, and m3 provides a short Gaussian pulse only to ignite
        # the mixture
        m1 = ct.MassFlowController(fuel_in, combustor, mdot=fuel_mdot)
        m2 = ct.MassFlowController(oxidizer_in, combustor, mdot=oxidizer_mdot)
        m3 = ct.MassFlowController(igniter, combustor, mdot=igniter_mdot)

        # put a valve on the exhaust line to regulate the pressure
        valve = ct.Valve(combustor, exhaust, K=1.0)

        # the simulation only contains one reactor
        sim = ct.ReactorNet([combustor])

        return {
            'fuel_in': fuel_in,
            'oxidizer_in': oxidizer_in,
            'igniter': igniter,
            'combustor': combustor,
            'exhaust': exhaust,
            'm1': m1,
            'm2': m2,
            'm3': m3,
            'valve': valve,
            'sim': sim
        }

    def test_integrateWithStep(self, test_data_path, setup_combustor_tests):
        sim = setup_combustor_tests['sim']
        combustor = setup_combustor_tests['combustor']

        tnow = 0.0
        tfinal = 0.25
        data = []
        while tnow < tfinal:
            tnow = sim.step()
            data.append([tnow, combustor.T] + list(combustor.thermo.X))

        assert tnow >= tfinal
        bad = compareProfiles(test_data_path / self.referenceFile, data,
                              rtol=1e-3, atol=1e-9)
        assert not bad, bad

    def test_integrateWithAdvance(self, request, test_data_path, setup_combustor_tests):
        sim = setup_combustor_tests['sim']
        combustor = setup_combustor_tests['combustor']

        data = []
        for t in np.linspace(0, 0.25, 101)[1:]:
            sim.advance(t)
            data.append([t, combustor.T] + list(combustor.thermo.X))

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'combustor':
            np.savetxt(self.referenceFile, np.array(data), '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / self.referenceFile, data,
                                  rtol=1e-6, atol=1e-12)
            assert not bad, bad

    def test_invasive_mdot_function(self, test_data_path, setup_combustor_tests):
        igniter = setup_combustor_tests['igniter']
        m3 = setup_combustor_tests['m3']
        sim = setup_combustor_tests['sim']
        combustor = setup_combustor_tests['combustor']

        def igniter_mdot(t, t0=0.1, fwhm=0.05, amplitude=0.1):
            # Querying properties of the igniter changes the state of the
            # underlying ThermoPhase object, but shouldn't affect the
            # integration
            igniter.density
            return amplitude * math.exp(-(t-t0)**2 * 4 * math.log(2) / fwhm**2)
        m3.mass_flow_rate = igniter_mdot

        data = []
        for t in np.linspace(0, 0.25, 101)[1:]:
            sim.advance(t)
            data.append([t, combustor.T] + list(combustor.thermo.X))

        bad = compareProfiles(test_data_path / self.referenceFile, data,
                              rtol=1e-6, atol=1e-12)
        assert not bad, bad


class TestWall:
    """
    These tests are based on the sample:

        samples/python/reactors/reactor2.py

    with some simplifications so that they run faster and produce more
    consistent output.

    Note: to re-create the reference file:
    (1) set PYTHONPATH to build/python.
    (2) go into test/python directory and run:
        pytest --save-reference=wall test_reactor.py::TestWall::test_integrateWithAdvance
    (3) Compare the reference files created in the current working directory with
        the ones in test/data and replace them if needed.
    """
    referenceFile = "WallTest-integrateWithAdvance.csv"

    @pytest.fixture
    def setup_wall_tests(self):
        # reservoir to represent the environment
        gas0 = ct.Solution("air.yaml")
        gas0.TP = 300, ct.one_atm
        env = ct.Reservoir(gas0)

        # reactor to represent the side filled with Argon
        gas1 = ct.Solution("air.yaml")
        gas1.TPX = 1000.0, 30*ct.one_atm, 'AR:1.0'
        r1 = ct.Reactor(gas1)

        # reactor to represent the combustible mixture
        gas2 = ct.Solution('h2o2.yaml', transport_model=None)
        gas2.TPX = 500.0, 1.5*ct.one_atm, 'H2:0.5, O2:1.0, AR:10.0'
        r2 = ct.Reactor(gas2)

        # Wall between the two reactors
        w1 = ct.Wall(r2, r1, A=1.0, K=2e-4, U=400.0)

        # Wall to represent heat loss to the environment
        w2 = ct.Wall(r2, env, A=1.0, U=2000.0)

        # Create the reactor network
        sim = ct.ReactorNet([r1, r2])

        return sim, r1, r2

    def test_integrateWithStep(self, test_data_path, setup_wall_tests):
        sim, r1, r2 = setup_wall_tests
        tnow = 0.0
        tfinal = 0.01
        data = []
        while tnow < tfinal:
            tnow = sim.step()
            data.append([tnow, r1.T, r2.T, r1.thermo.P,
                        r2.thermo.P, r1.volume, r2.volume])

        assert tnow >= tfinal
        bad = compareProfiles(test_data_path / self.referenceFile, data,
                              rtol=1e-3, atol=1e-8)
        assert not bad, bad

    def test_integrateWithAdvance(self, request, test_data_path, setup_wall_tests):
        sim, r1, r2 = setup_wall_tests
        data = []
        for t in np.linspace(0, 0.01, 200)[1:]:
            sim.advance(t)
            data.append([t, r1.T, r2.T, r1.thermo.P,
                        r2.thermo.P, r1.volume, r2.volume])

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'wall':
            np.savetxt(self.referenceFile, np.array(self.data), '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / self.referenceFile, data,
                                  rtol=2e-5, atol=1e-9)
            assert not bad, bad


class TestPureFluidReactor:
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

        assert states.Q[0] == 0
        assert states.Q[-1] == 1
        assert states.Q[30] == approx(0.54806, rel=1e-4)

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

        assert states.Q[0] == 0
        assert states.Q[-1] == 1
        assert states.Q[20] == approx(0.644865, rel=1e-4)


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

        assert states.Q[1] == 0
        assert states.Q[-2] == 1
        for i in range(3,7):
            assert states.T[i] == approx(states.T[2])


@pytest.fixture(scope='function')
def setup_advance_converages_data(request):
    mechanism_file = 'ptcombust.yaml'
    interface_phase = 'Pt_surf'
    request.cls.surf = ct.Interface(mechanism_file, interface_phase)
    request.cls.gas = request.cls.surf.adjacent["gas"]

@pytest.mark.usefixtures('setup_advance_converages_data')
class TestAdvanceCoverages:

    def test_bad_timestep_specification(self):
        # first, test max step size & max steps
        dt = 1.0
        max_steps = 10
        max_step_size = dt / (max_steps + 1)
        # this should throw an error, as we can't reach dt
        with pytest.raises(ct.CanteraError):
            self.surf.advance_coverages(
                dt=dt, max_step_size=max_step_size, max_steps=max_steps)

    def test_different_tolerances(self):
        dt = 1.0

        #Run with different tolerances
        self.surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        self.gas.TP = self.surf.TP

        self.surf.advance_coverages(dt=dt, rtol=1e-5, atol=1e-12)
        cov = self.surf.coverages[:]

        self.surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        self.gas.TP = self.surf.TP
        self.surf.advance_coverages(dt=dt, rtol=1e-7, atol=1e-14)

        # check that the solutions are similar, but not identical
        assert cov == approx(self.surf.coverages)
        assert any(cov != self.surf.coverages)


@pytest.fixture(scope='function')
def setup_extensible_reactor_data(request):
    request.cls.gas = ct.Solution("h2o2.yaml")

@pytest.mark.usefixtures('setup_extensible_reactor_data')
class TestExtensibleReactor:

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
                self.walls[0].velocity = self.v_wall

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
        assert (np.diff(V, 2) > 0).all()

        assert 'v_wall' in net.component_name(self.gas.n_species + 3)
        assert r.component_index('volume') == 1
        assert r.component_name(self.gas.n_species + 3) == 'v_wall'
        assert r.component_name(2) == 'temperature'

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
            assert r.get_state() == approx(np.exp(- net.time / tau))

    def test_error_handling(self):
        class DummyReactor1(ct.ExtensibleReactor):
            def replace_eval(self, t): # wrong number of arguments
                pass

        with pytest.raises(ValueError, match="right number of arguments"):
            DummyReactor1(self.gas)

        class DummyReactor2(ct.ExtensibleReactor):
            def replace_component_index(self, name):
                if name == "succeed":
                    return 0
                elif name == "wrong-type":
                    return "spam"
                # Otherwise, does not return a value

        r2 = DummyReactor2(self.gas)
        assert r2.component_index("succeed") == 0
        with pytest.raises(TypeError):
            r2.component_index("wrong-type")
        # Error information should have been reset
        assert r2.component_index("succeed") == 0
        with pytest.raises(ct.CanteraError, match="did not return a value"):
            r2.component_index("H2")
        assert r2.component_index("succeed") == 0

    def test_delegate_throws(self):
        class TestException(Exception):
            pass

        class DummyReactor(ct.ExtensibleConstPressureReactor):
            def before_eval(self, t, LHS, RHS):
                if t > 0.1:
                    raise TestException("spam")

            def before_component_index(self, name):
                if name == "fail":
                    raise TestException()

        r = DummyReactor(self.gas)
        net = ct.ReactorNet([r])
        net.max_steps = 10

        # Because the TestException is raised inside code called by CVODES, the actual
        # error raised will be a CanteraError
        with pytest.raises(ct.CanteraError, match="TestException: spam"):
            net.advance(0.2)

        assert r.component_index("enthalpy") == 1
        with pytest.raises(TestException):
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
        assert r.component_index("H2") == 5 + 3 + self.gas.species_index("H2")
        r.syncState()
        net.advance(1)
        r.syncState()
        assert r.sync_calls == 2

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
        assert add_heat == approx(r1_heat, abs=1e-5)

    def test_heat_addition(self):
        # Applying heat via 'heat_rate' property should be equivalent to adding it via
        # a wall
        Qext = 100
        Qwall = -66
        class HeatedReactor(ct.ExtensibleReactor):
            def after_eval_walls(self, y):
                self.heat_rate += Qext

        self.gas.TPX = 300, ct.one_atm, "N2:1.0"
        r1 = HeatedReactor(self.gas)
        res = ct.Reservoir(self.gas)
        wall = ct.Wall(res, r1, Q=Qwall, A=1)
        net = ct.ReactorNet([r1])
        U0 = r1.thermo.int_energy_mass * r1.mass
        for t in np.linspace(0.1, 5, 10):
            net.advance(t)
            U = r1.thermo.int_energy_mass * r1.mass
            assert U - U0 == approx((Qext + Qwall) * t)
            assert r1.heat_rate == approx(Qext + Qwall)

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
            assert mH == approx(mH0)
            assert mO == approx(mO0)

        # Regression test values
        assert r1.thermo.P == approx(647.56016304)
        assert r1.thermo.X[kH2] == approx(0.4784268406)
        assert r1.thermo.X[kO2] == approx(0.5215731594)
        assert r1.surfaces[0].kinetics.X[kHs] == approx(0.3665198138)
        assert r1.surfaces[0].kinetics.X[kPts] == approx(0.6334801862)

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
                    self.expansion_rate = np.clip(
                        self.p_coeff * (self.P - self.neighbor.P), -1.7, 1.7)

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
            assert r1.expansion_rate == approx(-r2.expansion_rate)
            assert V0 == approx(r1.volume + r2.volume)
            deltaCnow = deltaC()
            assert deltaCnow < deltaCprev # difference is always decreasing
            deltaCprev = deltaCnow
            assert M0 == approx(r1.mass * r1.thermo.Y + r2.mass * r2.thermo.Y,rel=2e-8)
            states1.append(r1.thermo.state, t=net.time, mass=r1.mass,
                           vdot=r1.expansion_rate)
            states2.append(r2.thermo.state, t=net.time, mass=r2.mass,
                           vdot=r2.expansion_rate)

        # Regression test values
        assert r1.thermo.P == approx(151561.15, rel=1e-6)
        assert r1.thermo["H2"].Y[0] == approx(0.13765976, rel=1e-6)
        assert r2.thermo["O2"].Y[0] == approx(0.94617029, rel=1e-6)
