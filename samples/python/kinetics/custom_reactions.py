"""
An example demonstrating how to use custom reaction objects.

For benchmark purposes, an ignition test is run to compare simulation times.

Requires: cantera >= 3.0.0
Keywords: kinetics, benchmarking, user-defined model
"""

from timeit import default_timer
import numpy as np
from math import exp
import warnings

import cantera as ct

mech = 'gri30.yaml'
fuel = 'CH4'
gas0 = ct.Solution(mech)

species = gas0.species()
reactions = gas0.reactions()

# construct custom reactions: replace 2nd reaction with equivalent custom reaction
custom_reactions = [r for r in reactions]
custom_reactions[2] = ct.Reaction(
    equation='H2 + O <=> H + OH',
    rate=lambda T: 38.7 * T**2.7 * exp(-3150.1542797022735/T))
custom_reactions[4] = ct.Reaction(
    equation='H2O2 + O <=> HO2 + OH',
    rate=lambda T: 9630.0 * T**2.0 * exp(-2012.8781339950629/T))

gas1 = ct.Solution(thermo='ideal-gas', kinetics='gas',
                   species=species, reactions=custom_reactions)

# construct reactions based on ExtensibleRate: replace 2nd reaction with equivalent
# ExtensibleRate
class ExtensibleArrheniusData(ct.ExtensibleRateData):
    __slots__ = ("T",)
    def __init__(self):
        self.T = None

    def update(self, gas):
        T = gas.T
        if self.T != T:
            self.T = T
            return True
        else:
            return False

@ct.extension(name="extensible-Arrhenius", data=ExtensibleArrheniusData)
class ExtensibleArrhenius(ct.ExtensibleRate):
    __slots__ = ("A", "b", "Ea_R")
    def set_parameters(self, params, units):
        self.A = params.convert_rate_coeff("A", units)
        self.b = params["b"]
        self.Ea_R = params.convert_activation_energy("Ea", "K")

    def get_parameters(self, params):
        params.set_quantity("A", self.A, self.conversion_units)
        params["b"] = self.b
        params.set_activation_energy("Ea", self.Ea_R, "K")

    def validate(self, equation, soln):
        if self.A < 0:
            raise ValueError(f"Found negative 'A' for reaction {equation}")

    def eval(self, data):
        return self.A * data.T**self.b * exp(-self.Ea_R/data.T)

extensible_yaml2 = """
    equation: H2 + O <=> H + OH
    type: extensible-Arrhenius
    units: {length: cm, quantity: mol, activation-energy: cal/mol}
    A: 3.87e+04
    b: 2.7
    Ea: 6260.0
    """

extensible_yaml4 = """
    equation: H2O2 + O <=> HO2 + OH
    type: extensible-Arrhenius
    units: {length: cm, quantity: mol, activation-energy: cal/mol}
    A: 9.63e+06
    b: 2
    Ea: 4000
    """

extensible_reactions = gas0.reactions()
extensible_reactions[2] = ct.Reaction.from_yaml(extensible_yaml2, gas0)
extensible_reactions[4] = ct.Reaction.from_yaml(extensible_yaml4, gas0)
gas2 = ct.Solution(thermo="ideal-gas", kinetics="gas",
                   species=species, reactions=extensible_reactions)

# construct test case - simulate ignition

def ignition(gas, dT=0):
    # set up reactor
    gas.TP = 1000 + dT, 5 * ct.one_atm
    gas.set_equivalence_ratio(0.8, fuel, 'O2:1.0, N2:3.773')
    r = ct.IdealGasReactor(gas)
    net = ct.ReactorNet([r])
    net.rtol_sensitivity = 2.e-5

    # time reactor integration
    t1 = default_timer()
    net.advance(.5)
    t2 = default_timer()

    return t2 - t1, net.solver_stats['steps']

# output results

repeat = 100
print("Average time of {} simulation runs for '{}' "
      "({})".format(repeat, mech, fuel))

sim0 = 0
sim0_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas0, dT=i)
    sim0 += elapsed
    sim0_steps += steps
sim0 *= 1e6 / sim0_steps
print('- New framework (YAML): '
      '{0:.2f} μs/step (T_final={1:.2f})'.format(sim0, gas0.T))

sim1 = 0
sim1_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas1, dT=i)
    sim1 += elapsed
    sim1_steps += steps
sim1 *= 1e6 / sim1_steps
print('- Two Custom reactions: '
      '{0:.2f} μs/step (T_final={1:.2f}) ...'
      '{2:+.2f}%'.format(sim1, gas1.T, 100 * sim1 / sim0 - 100))

sim2 = 0
sim2_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas2, dT=i)
    sim2 += elapsed
    sim2_steps += steps
sim2 *= 1e6 / sim2_steps
print('- Two Extensible reactions: '
      '{0:.2f} μs/step (T_final={1:.2f}) ...'
      '{2:+.2f}%'.format(sim2, gas1.T, 100 * sim2 / sim0 - 100))
