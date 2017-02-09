"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.
"""

import cantera as ct
import numpy as np

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:1, O2:2, N2:7.52'  # premixed gas composition
width = 0.05  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30_ion.xml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.IonFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.show_solution()

# phase one
f.solve(loglevel=loglevel, auto=True)

# phase two
f.solve(loglevel=loglevel, stage=2, enable_energy=False)
f.solve(loglevel=loglevel, stage=2, enable_energy=True)

# phase three
f.solve(loglevel=loglevel, stage=3, enable_energy=True)

f.save('CH4_adiabatic.xml', 'mix', 'solution with mixture-averaged transport')
f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('CH4_adiabatic.csv', quiet=False)

