"""
An example demonstrating how to use custom reaction objects.

For benchmark purposes, an ignition test is run to compare simulation times.

Requires: cantera >= 2.6.0
"""

from timeit import default_timer
import numpy as np
from math import exp

import cantera as ct

mech = 'gri30.yaml'
fuel = 'CH4'
gas0 = ct.Solution(mech)

species = gas0.species()
reactions = gas0.reactions()

# construct custom reactions: replace 2nd reaction with equivalent custom reaction
custom_reactions = [r for r in reactions]
custom_reactions[2] = ct.CustomReaction(
    equation='H2 + O <=> H + OH',
    rate=lambda T: 38.7 * T**2.7 * exp(-3150.15428/T),
    kinetics=gas0)

gas1 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                   species=species, reactions=custom_reactions)

# old framework - use xml input

gas2 = ct.Solution(mech.replace('.yaml', '.xml'))

# construct test case - simulate ignition

def ignition(gas):
    # set up reactor
    gas.TP = 1000., 5 * ct.one_atm
    gas.set_equivalence_ratio(0.8, fuel, 'O2:1.0, N2:3.773')
    r = ct.IdealGasReactor(gas)
    net = ct.ReactorNet([r])
    net.rtol_sensitivity = 2.e-5

    # time reactor integration
    t1 = default_timer()
    net.advance(.5)
    t2 = default_timer()

    return 1000 * (t2 - t1)

# output results

repeat = 100
print("Average time of {} simulation runs for '{}' "
      "({})".format(repeat, mech, fuel))

sim0 = 0
for i in range(repeat):
    sim0 += ignition(gas0)
sim0 /= repeat
print('- New framework (YAML): '
      '{0:.2f} ms (T_final={1:.2f})'.format(sim0, gas0.T))

sim1 = 0
for i in range(repeat):
    sim1 += ignition(gas1)
sim1 /= repeat
print('- One Python reaction: '
      '{0:.2f} ms (T_final={1:.2f}) ... '
      '{2:+.2f}%'.format(sim1, gas1.T, 100 * sim1 / sim0 - 100))

sim2 = 0
for i in range(repeat):
    sim2 += ignition(gas2)
sim2 /= repeat
print('- Old framework (XML): '
      '{0:.2f} ms (T_final={1:.2f}) ... '
      '{2:+.2f}%'.format(sim2, gas2.T, 100 * sim2 / sim0 - 100))
