# Mixing two streams.

# Since reactors can have multiple inlets and outlets, they can be
# used to implement mixers, splitters, etc. In this example, air and
# methane are mixed in stoichiometric proportions. Due to the low
# temperature, no reactions occur. Note that the air stream and the
# methane stream use *different* reaction mechanisms, with different
# numbers of species and reactions. When gas flows from one reactor or
# reservoir to another one with a different reaction mechanism,
# species are matched by name. If the upstream reactor contains a
# species that is not present in the downstream reaction mechanism, it
# will be ignored. In general, reaction mechanisms for downstream
# reactors should contain all species that might be present in any
# upstream reactor.
#
#-----------------------------------------------------------------------

from Cantera import *
from Cantera.Reactor import *


# Use air for stream a. Note that the Air() function does not set the
# composition correctly; thus, we need to explicitly set the
# composition to that of air.
gas_a = Air()
gas_a.set(T = 300.0, P = OneAtm, X = 'O2:0.21, N2:0.78, AR:0.01')
rho_a = gas_a.density()


# Use GRI-Mech 3.0 for stream b (methane) and for the mixer. If it is
# desired to have a pure mixer, with no chemistry, use instead a
# reaction mechanism for gas_b that has no reactions.
gas_b = GRI30()
gas_b.set(T = 300.0, P = OneAtm, X = 'CH4:1')
rho_b = gas_b.density()


# Create reservoirs for the two inlet streams and for the outlet
# stream.  The upsteam reservoirs could be replaced by reactors, which
# might themselves be connected to reactors further upstream. The
# outlet reservoir could be replaced with a reactor with no outlet, if
# it is desired to integrate the composition leaving the mixer in
# time, or by an arbitrary network of downstream reactors.
res_a = Reservoir(gas_a)
res_b = Reservoir(gas_b)
downstream   = Reservoir(gas_b)


# Create a reactor for the mixer. A reactor is required instead of a
# reservoir, since the state will change with time if the inlet mass
# flow rates change or if there is chemistry occurring.
mixer = Reactor(gas_b)


# create two mass flow controllers connecting the upstream reservoirs
# to the mixer, and set their mass flow rates to values corresponding
# to stoichiometric combustion.
mfc1 = MassFlowController(upstream = res_a, downstream = mixer,
                          mdot = rho_a*2.5/0.21)

mfc2 = MassFlowController(upstream = res_b, downstream = mixer,
                          mdot = rho_b*1.0)


# connect the mixer to the downstream reservoir with a valve.
outlet = Valve(upstream = mixer, downstream = downstream, Kv = 1.0)

sim = ReactorNet([mixer])

# Since the mixer is a reactor, we need to integrate in time to reach
# steady state. A few residence times should be enough.
t = 0.0
for n in range(30):
    tres = mixer.mass()/(mfc1.massFlowRate() + mfc2.massFlowRate())
    t += 0.5*tres
    sim.advance(t)
    print '%14.5g %14.5g %14.5g  %14.5g  %14.5g' % (t, mixer.temperature(),
                                                    mixer.enthalpy_mass(),
                                                    mixer.pressure(),
                                                    mixer.massFraction('CH4'))

# view the state of the gas in the mixer
print mixer.contents()
