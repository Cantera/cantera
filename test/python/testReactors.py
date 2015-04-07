import unittest
import pprint

import numpy as np

import Cantera as ct
from Cantera import Reactor as reactors
from Cantera.Func import Gaussian

import utilities

class CombustorTestImplementation(object):
    """
    These tests are based on the sample:

        python/reactors/combustor_sim/combustor.py

    with some simplifications so that they run faster and produce more
    consistent output.
    """

    referenceFile = '../data/CombustorTest-integrateWithAdvance.csv'
    def setUp(self):
        self.gas = ct.importPhase('h2o2.cti')

        # create a reservoir for the fuel inlet, and set to pure methane.
        self.gas.set(T=300.0, P=ct.OneAtm, X='H2:1.0')
        fuel_in = reactors.Reservoir(self.gas)
        fuel_mw = self.gas.meanMolarMass()

        # Oxidizer inlet
        self.gas.set(T=300.0, P=ct.OneAtm, X='O2:1.0, AR:3.0')
        oxidizer_in = reactors.Reservoir(self.gas)
        oxidizer_mw = self.gas.meanMolarMass()

        # to ignite the fuel/air mixture, we'll introduce a pulse of radicals.
        # The steady-state behavior is independent of how we do this, so we'll
        # just use a stream of pure atomic hydrogen.
        self.gas.set(T=300.0, P=ct.OneAtm, X='H:1.0')
        self.igniter = reactors.Reservoir(self.gas)

        # create the combustor, and fill it in initially with a diluent
        self.gas.set(T=300.0, P=ct.OneAtm, X='AR:1.0')
        self.combustor = reactors.Reactor(contents=self.gas, volume=1.0)

        # create a reservoir for the exhaust
        self.exhaust = reactors.Reservoir(self.gas)

        # compute fuel and air mass flow rates
        factor = 0.1
        oxidizer_mdot = 4 * factor*oxidizer_mw
        fuel_mdot = factor*fuel_mw

        # create and install the mass flow controllers. Controllers
        # m1 and m2 provide constant mass flow rates, and m3 provides
        # a short Gaussian pulse only to ignite the mixture
        m1 = reactors.MassFlowController(upstream=fuel_in,
                                         downstream=self.combustor,
                                         mdot=fuel_mdot)

        m2 = reactors.MassFlowController(upstream=oxidizer_in,
                                         downstream=self.combustor,
                                         mdot=oxidizer_mdot)

        # The igniter will use a Gaussian 'functor' object to specify the
        # time-dependent igniter mass flow rate.
        igniter_mdot = Gaussian(t0=0.1, FWHM=0.05, A=0.1)
        m3 = reactors.MassFlowController(upstream=self.igniter,
                                         downstream=self.combustor,
                                         mdot=igniter_mdot)

        # put a valve on the exhaust line to regulate the pressure
        self.v = reactors.Valve(upstream=self.combustor,
                                downstream=self.exhaust, Kv=1.0)

        # the simulation only contains one reactor
        self.sim = reactors.ReactorNet([self.combustor])
        #self.sim.setTolerances(1e-8, 1e-12)

    def test_integrateWithStep(self):
        tnow = 0.0
        tfinal = 0.25
        self.data = []
        while tnow < tfinal:
            tnow = self.sim.step(tfinal)
            self.data.append([tnow, self.combustor.temperature()] +
                             list(self.combustor.moleFractions()))

        self.assertTrue(tnow >= tfinal)
        bad = utilities.compareTimeSeries(self.referenceFile, self.data,
                                          rtol=1e-3, atol=1e-9)
        self.assertFalse(bad, bad)

    def test_integrateWithAdvance(self, saveReference=False):
        times = np.linspace(0, 0.25, 101)
        self.data = []
        for t in times[1:]:
            self.sim.advance(t)
            self.data.append([t, self.combustor.temperature()] +
                             list(self.combustor.moleFractions()))

        if saveReference:
            np.savetxt(self.referenceFile, np.array(self.data), '%11.6e', ', ')
        else:
            bad = utilities.compareTimeSeries(self.referenceFile, self.data,
                                              rtol=1e-6, atol=1e-12)
            self.assertFalse(bad, bad)


class WallTestImplementation(object):
    """
    These tests are based on the sample:

        python/reactors/reactor2_sim/reactor2.py

    with some simplifications so that they run faster and produce more
    consistent output.
    """

    referenceFile = '../data/WallTest-integrateWithAdvance.csv'
    def setUp(self):
        # reservoir to represent the environment
        self.gas0 = ct.importPhase('air.cti')
        self.gas0.set(T=300, P=ct.OneAtm)
        self.env = reactors.Reservoir(self.gas0)

        # reactor to represent the side filled with Argon
        self.gas1 = ct.importPhase('air.cti')
        self.gas1.set(T=1000.0, P=30*ct.OneAtm, X='AR:1.0')
        self.r1 = reactors.Reactor(self.gas1)

        # reactor to represent the combustible mixture
        self.gas2 = ct.importPhase('h2o2.cti')
        self.gas2.set(T=500.0, P=1.5*ct.OneAtm, X='H2:0.5, O2:1.0, AR:10.0')
        self.r2 = reactors.Reactor(self.gas2)

        # Wall between the two reactors
        self.w1 = reactors.Wall(self.r2, self.r1)
        self.w1.set(area=1.0, K=2e-4, U=400.0)

        # Wall to represent heat loss to the environment
        self.w2 = reactors.Wall(self.r2, self.env)
        self.w2.set(area=1.0, U=2000.0)

        # Create the reactor network
        self.sim = reactors.ReactorNet([self.r1, self.r2])

    def test_integrateWithStep(self):
        tnow = 0.0
        tfinal = 0.01
        self.data = []
        while tnow < tfinal:
            tnow = self.sim.step(tfinal)
            self.data.append([tnow,
                              self.r1.temperature(),
                              self.r2.temperature(),
                              self.r1.pressure(),
                              self.r2.pressure(),
                              self.r1.volume(),
                              self.r2.volume()])

        self.assertTrue(tnow >= tfinal)
        bad = utilities.compareTimeSeries(self.referenceFile, self.data,
                                          rtol=1e-3, atol=1e-8)
        self.assertFalse(bad, bad)

    def test_integrateWithAdvance(self, saveReference=False):
        times = np.linspace(0, 0.01, 200)
        self.data = []
        for t in times[1:]:
            self.sim.advance(t)
            self.data.append([t,
                              self.r1.temperature(),
                              self.r2.temperature(),
                              self.r1.pressure(),
                              self.r2.pressure(),
                              self.r1.volume(),
                              self.r2.volume()])

        if saveReference:
            np.savetxt(self.referenceFile, np.array(self.data), '%11.6e', ', ')
        else:
            bad = utilities.compareTimeSeries(self.referenceFile, self.data,
                                              rtol=2e-5, atol=1e-9)
            self.assertFalse(bad, bad)


# Keep the implementations separate from the unittest-derived class
# so that they can be run independently to generate the reference data files.
class CombustorTest(CombustorTestImplementation, unittest.TestCase): pass
class WallTest(WallTestImplementation, unittest.TestCase): pass
