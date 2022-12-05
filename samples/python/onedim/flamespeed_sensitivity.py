"""
Sensitivity analysis for a freely-propagating, premixed methane-air
flame. Computes the sensitivity of the laminar flame speed with respect
to each reaction rate constant.

Requires: cantera >= 2.5.0
Keywords: combustion, 1D flow, flame speed, premixed flame, sensitivity analysis
"""

import cantera as ct

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:0.45, O2:1.0, N2:3.76'

width = 0.03  # m

# Solution object used to compute mixture properties
gas = ct.Solution('gri30.yaml')
gas.TPX = Tin, p, reactants

# Flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)

f.solve(loglevel=1, auto=True)
print(f"\nmixture-averaged flamespeed = {f.velocity[0]:7f} m/s\n")

# Use the adjoint method to calculate sensitivities
sens = f.get_flame_speed_reaction_sensitivities()

print()
print('Rxn #   k/S*dS/dk    Reaction Equation')
print('-----   ----------   ----------------------------------')
for m in range(gas.n_reactions):
    print(f"{m: 5d}   {sens[m]: 10.3e}   {gas.reaction(m).equation}")
