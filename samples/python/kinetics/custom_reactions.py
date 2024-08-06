"""
Custom reaction rates
=====================

An example demonstrating how to use custom reaction objects.

For benchmark purposes, an ignition test is run to compare simulation times.

Requires: cantera >= 3.1.0

.. tags:: Python, kinetics, benchmarking, user-defined model
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

# construct reactions based on CustomRate: replace two reactions with equivalent
# custom versions (uses Func1 functors defined by Python lambda functions)
reactions[2] = ct.Reaction(
    equation='H2 + O <=> H + OH',
    rate=lambda T: 38.7 * T**2.7 * exp(-3150.1542797022735/T))
reactions[4] = ct.Reaction(
    equation='H2O2 + O <=> HO2 + OH',
    rate=lambda T: 9630.0 * T**2.0 * exp(-2012.8781339950629/T))

gas1a = ct.Solution(thermo='ideal-gas', kinetics='gas',
                    species=species, reactions=reactions)

# construct reactions based on CustomRate: replace two reactions with equivalent
# custom versions (uses Func1 functors defined by C++ class Arrhenius1)
reactions[2] = ct.Reaction(
    equation='H2 + O <=> H + OH',
    rate=ct.Func1("Arrhenius", [38.7, 2.7, 3150.1542797022735]))
reactions[4] = ct.Reaction(
    equation='H2O2 + O <=> HO2 + OH',
    rate=ct.Func1("Arrhenius", [9630.0, 2.0, 2012.8781339950629]))

gas1b = ct.Solution(thermo='ideal-gas', kinetics='gas',
                    species=species, reactions=reactions)

# construct reactions based on ExtensibleRate: replace two reactions with equivalent
# extensible versions
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

reactions[2] = ct.Reaction.from_yaml(extensible_yaml2, gas0)
reactions[4] = ct.Reaction.from_yaml(extensible_yaml4, gas0)
gas2 = ct.Solution(thermo="ideal-gas", kinetics="gas",
                   species=species, reactions=reactions)

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
print(f"Average time of {repeat} simulation runs for '{mech}' ({fuel})")

sim0 = 0
sim0_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas0, dT=i)
    sim0 += elapsed
    sim0_steps += steps
sim0 *= 1e6 / sim0_steps
print(f'- Built-in rate parameterizations: {sim0:.2f} μs/step (T_final={gas0.T:.2f})')

sim1a = 0
sim1a_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas1a, dT=i)
    sim1a += elapsed
    sim1a_steps += steps
sim1a *= 1e6 / sim1a_steps
print('- Two Custom reactions (Python): '
      f'{sim1a:.2f} μs/step (T_final={gas1a.T:.2f}) ... '
      f'{100 * sim1a / sim0 - 100:+.2f}%')

sim1b = 0
sim1b_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas1b, dT=i)
    sim1b += elapsed
    sim1b_steps += steps
sim1b *= 1e6 / sim1b_steps
print('- Two Custom reactions (C++): '
      f'{sim1b:.2f} μs/step (T_final={gas1b.T:.2f}) ... '
      f'{100 * sim1b / sim0 - 100:+.2f}%')

sim2 = 0
sim2_steps = 0
for i in range(repeat):
    elapsed, steps = ignition(gas2, dT=i)
    sim2 += elapsed
    sim2_steps += steps
sim2 *= 1e6 / sim2_steps
print('- Two Extensible reactions: '
      f'{sim2:.2f} μs/step (T_final={gas2.T:.2f}) ... '
      f'{100 * sim2 / sim0 - 100:+.2f}%')
