"""
A burner-stabilized, premixed methane/air flat flame with multicomponent
transport properties and a specified temperature profile.

Requires: cantera >= 2.5.0
"""

import cantera as ct
import numpy as np

################################################################
# parameter values
p = ct.one_atm  # pressure
tburner = 373.7  # burner temperature
mdot = 0.04  # kg/m^2/s
comp = 'CH4:0.65, O2:1, N2:3.76'  # premixed gas composition

# The solution domain is chosen to be 1 cm
width = 0.01  # m

loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = True  # 'True' to enable refinement

################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic, and
# transport properties
gas = ct.Solution('gri30.yaml')

# set its state to that of the unburned gas at the burner
gas.TPX = tburner, p, comp

# create the BurnerFlame object.
f = ct.BurnerFlame(gas=gas, width=width)

# set the mass flow rate at the burner
f.burner.mdot = mdot

# read temperature vs. position data from a file.
# The file is assumed to have one z, T pair per line, separated by a comma.
zloc, tvalues = np.genfromtxt('tdata.dat', delimiter=',', comments='#').T
zloc /= max(zloc)

# set the temperature profile to the values read in
f.flame.set_fixed_temp_profile(zloc, tvalues)

# show the initial estimate for the solution
f.show_solution()

# don't solve the energy equation
f.energy_enabled = False

# first solve the flame with mixture-averaged transport properties
f.transport_model = 'Mix'
f.set_refine_criteria(ratio=3.0, slope=0.3, curve=1)

f.solve(loglevel, refine_grid)
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('flame_fixed_T.h5', group='mix', mode='w',
                description='solution with mixture-averaged transport')
except ImportError:
    f.save('flame_fixed_T.xml','mixav',
           'solution with mixture-averaged transport')

print('\n\n switching to multicomponent transport...\n\n')
f.transport_model = 'Multi'

f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2)
f.solve(loglevel, refine_grid)
try:
    f.write_hdf('flame_fixed_T.h5', group='multi',
                description='solution with multicomponent transport')
except ImportError:
    f.save('flame_fixed_T.xml','multi',
           'solution with  multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('flame_fixed_T.csv', quiet=False)
f.show_stats()
