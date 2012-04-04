""" A combustor. Two separate stream - one pure methane and the other
 air, both at 300 K and 1 atm flow into an adiabatic combustor where
 they mix. We are interested in the steady-state burning
 solution. Since at 300 K no reaction will occur between methane and
 air, we need to use an 'igniter' to initiate the chemistry. A simple
 igniter is a pulsed flow of atomic hydrogen. After the igniter is
 turned off, the system approaches the steady burning solution."""

from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *

# use reaction mechanism GRI-Mech 3.0
gas = GRI30()

# create a reservoir for the fuel inlet, and set to pure methane.
gas.set(T = 300.0, P = OneAtm, X = 'CH4:1.0')
fuel_in = Reservoir(gas)
fuel_mw = gas.meanMolarMass()

# use predefined function Air() for the air inlet
air = Air()
air_in = Reservoir(air)
air_mw = air.meanMolarMass()

# to ignite the fuel/air mixture, we'll introduce a pulse of radicals.
# The steady-state behavior is independent of how we do this, so we'll
# just use a stream of pure atomic hydrogen.
gas.set(T = 300.0, P = OneAtm, X = 'H:1.0')
igniter = Reservoir(gas)


# create the combustor, and fill it in initially with N2
gas.set(T = 300.0, P = OneAtm, X = 'N2:1.0')
combustor = Reactor(contents = gas, volume = 1.0)

# create a reservoir for the exhaust
exhaust = Reservoir(gas)

# lean combustion, phi = 0.5
equiv_ratio = 0.5

# compute fuel and air mass flow rates
factor = 0.1
air_mdot = factor*9.52*air_mw
fuel_mdot = factor*equiv_ratio*fuel_mw

# create and install the mass flow controllers. Controllers
# m1 and m2 provide constant mass flow rates, and m3 provides
# a short Gaussian pulse only to ignite the mixture
m1 = MassFlowController(upstream = fuel_in,
                        downstream = combustor, mdot = fuel_mdot)

# note that this connects two reactors with different reaction
# mechanisms and different numbers of species. Downstream and upstream
# species are matched by name.
m2 = MassFlowController(upstream = air_in,
                        downstream = combustor, mdot = air_mdot)

# The igniter will use a Gaussian 'functor' object to specify the
# time-dependent igniter mass flow rate.
igniter_mdot = Gaussian(t0 = 1.0, FWHM = 0.2, A = 0.1)
m3 = MassFlowController(upstream = igniter,
                        downstream = combustor, mdot = igniter_mdot)

# put a valve on the exhaust line to regulate the pressure
v = Valve(upstream = combustor, downstream = exhaust, Kv = 1.0)

# the simulation only contains one reactor
sim = ReactorNet([combustor])

# take single steps to 6 s, writing the results to a CSV file
# for later plotting.
tfinal = 6.0
tnow = 0.0
f = open('combustor.csv','w')
while tnow < tfinal:
    tnow = sim.step(tfinal)
    tres = combustor.mass()/v.massFlowRate()
    writeCSV(f, [tnow, combustor.temperature(), tres]
             +list(combustor.moleFractions()))
f.close()
