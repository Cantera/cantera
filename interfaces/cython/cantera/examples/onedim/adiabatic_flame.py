"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.
"""

import cantera as ct
import numpy as np

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'H2:1.1, O2:1, AR:5'  # premixed gas composition

width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)
refine_grid = True  # 'True' to enable refinement, 'False' to disable

# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('h2o2.xml')
gas.TPX = Tin, p, reactants

# Flame object
f = ct.FreeFlame(gas, width=width)
f.show_solution()

# Solve with the energy equation disabled
f.energy_enabled = False
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, refine_grid=False)
f.save('h2_adiabatic.xml', 'no_energy',
       'solution with the energy equation disabled')

# Solve with the energy equation enabled
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.energy_enabled = True
f.solve(loglevel=loglevel, refine_grid=refine_grid)
f.save('h2_adiabatic.xml', 'energy',
       'solution with mixture-averaged transport')
f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))

# Solve with multi-component transport properties
f.transport_model = 'Multi'
f.solve(loglevel, refine_grid)
f.show_solution()
print('multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))
f.save('h2_adiabatic.xml','energy_multi',
       'solution with multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('h2_adiabatic.csv', quiet=False)
