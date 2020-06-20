"""
A freely-propagating, premixed methane-air flat flame with charged species.

Requires: cantera >= 2.5.0
"""

import cantera as ct

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:1, O2:2, N2:7.52'  # premixed gas composition
width = 0.05  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30_ion.yaml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.IonFreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)
f.show_solution()

# stage one
f.solve(loglevel=loglevel, auto=True)

# stage two
f.solve(loglevel=loglevel, stage=2, enable_energy=True)

try:
    # save to HDF container file if h5py is installed
    f.write_hdf('ion_free_flame.h5', group='ion', mode='w',
                description='solution with ionized gas transport')
except ImportError:
    f.save('ion_free_flame.xml', 'ion', 'solution with ionized gas transport')

f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.velocity[0]))

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('ion_free_flame.csv', quiet=False)
