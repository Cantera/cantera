"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 2.5.0
"""

import cantera as ct

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'H2:1.1, O2:1, AR:5'  # premixed gas composition
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('h2o2.yaml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.show_solution()

# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)

# Solve with the energy equation enabled
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('adiabatic_flame.h5', group='mix', mode='w',
                description='solution with mixture-averaged transport')
except ImportError:
    f.save('adiabatic_flame.xml', 'mix',
           'solution with mixture-averaged transport')

f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.velocity[0]))

# Solve with multi-component transport properties
f.transport_model = 'Multi'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show_solution()
print('multicomponent flamespeed = {0:7f} m/s'.format(f.velocity[0]))
try:
    f.write_hdf('adiabatic_flame.h5', group='multi',
                description='solution with multicomponent transport')
except ImportError:
    f.save('adiabatic_flame.xml', 'multi',
           'solution with multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('adiabatic_flame.csv', quiet=False)
