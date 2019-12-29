"""
Mixing two streams.

Since reactors can have multiple inlets and outlets, they can be used to
implement mixers, splitters, etc. In this example, air and methane are mixed
in stoichiometric proportions. Due to the low temperature, no reactions occur.
Note that the air stream and the methane stream use *different* reaction
mechanisms, with different numbers of species and reactions. When gas flows
from one reactor or reservoir to another one with a different reaction
mechanism, species are matched by name. If the upstream reactor contains a
species that is not present in the downstream reaction mechanism, it will be
ignored. In general, reaction mechanisms for downstream reactors should
contain all species that might be present in any upstream reactor.

Compare this approach for the transient problem to the method used for the
steady-state problem in thermo/mixing.py.

Requires: cantera >= 2.5.0
"""

import cantera as ct

# Use air for stream a.
gas_a = ct.Solution('air.yaml')
gas_a.TPX = 300.0, ct.one_atm, 'O2:0.21, N2:0.78, AR:0.01'
rho_a = gas_a.density


# Use GRI-Mech 3.0 for stream b (methane) and for the mixer. If it is desired
# to have a pure mixer, with no chemistry, use instead a reaction mechanism
# for gas_b that has no reactions.
gas_b = ct.Solution('gri30.yaml')
gas_b.TPX = 300.0, ct.one_atm, 'CH4:1'
rho_b = gas_b.density

# Create reservoirs for the two inlet streams and for the outlet stream.  The
# upsteam reservoirs could be replaced by reactors, which might themselves be
# connected to reactors further upstream. The outlet reservoir could be
# replaced with a reactor with no outlet, if it is desired to integrate the
# composition leaving the mixer in time, or by an arbitrary network of
# downstream reactors.
res_a = ct.Reservoir(gas_a)
res_b = ct.Reservoir(gas_b)
downstream = ct.Reservoir(gas_b)

# Create a reactor for the mixer. A reactor is required instead of a
# reservoir, since the state will change with time if the inlet mass flow
# rates change or if there is chemistry occurring.
gas_b.TPX = 300.0, ct.one_atm, 'O2:0.21, N2:0.78, AR:0.01'
mixer = ct.IdealGasReactor(gas_b)

# create two mass flow controllers connecting the upstream reservoirs to the
# mixer, and set their mass flow rates to values corresponding to
# stoichiometric combustion.
mfc1 = ct.MassFlowController(res_a, mixer, mdot=rho_a*2.5/0.21)
mfc2 = ct.MassFlowController(res_b, mixer, mdot=rho_b*1.0)

# connect the mixer to the downstream reservoir with a valve.
outlet = ct.Valve(mixer, downstream, K=10.0)

sim = ct.ReactorNet([mixer])

# Since the mixer is a reactor, we need to integrate in time to reach steady
# state
sim.advance_to_steady_state()

# view the state of the gas in the mixer
print(mixer.thermo.report())
