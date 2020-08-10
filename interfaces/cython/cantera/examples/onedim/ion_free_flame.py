"""
A freely-propagating, premixed methane-air flat flame with charged species.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

import sys
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
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('ion_free_flame.h5', group='frozenIon', mode='w',
                description='solution with ionized gas transport - stage 1')
except ImportError:
    pass

# stage two
f.solve(loglevel=loglevel, stage=2, enable_energy=True)

try:
    f.write_hdf('ion_free_flame.h5', group='ion',
                description='solution with ionized gas transport - stage 2')
except ImportError:
    f.save('ion_free_flame.xml', 'ion', 'solution with ionized gas transport')

f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.velocity[0]))

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('ion_free_flame.csv', quiet=False)

if '--plot' in sys.argv:
    import matplotlib.pyplot as plt

    arr = f.to_solution_array()

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    ax1.plot(f.grid, f.heat_release_rate)
    ax1.set_ylabel(r'$\dot{q}$ [W/m${^3}$]')
    for s in ['H3O+', 'HCO+', 'E']:
        ax2.plot(f.grid, arr(s).X, label=s)
    ax2.legend()
    ax2.set_ylabel('X [-]')
    ax3.plot(f.grid, f.E)
    ax3.set_xlabel('Axial location [m]')
    ax3.set_ylabel('E [V/m]')
    ax3.set_xlim([.016, .022])

    fig.tight_layout()
    plt.show()
