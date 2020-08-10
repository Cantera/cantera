"""
A burner-stabilized premixed methane-air flame with charged species.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

import sys
import cantera as ct

p = ct.one_atm
tburner = 600.0
reactants = 'CH4:1.0, O2:2.0, N2:7.52'  # premixed gas composition
width = 0.5  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('gri30_ion.yaml')
gas.TPX = tburner, p, reactants
mdot = 0.15 * gas.density

f = ct.IonBurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.solve(loglevel, auto=True)
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('ion_burner_flame.h5', group='frozenIon', mode='w',
                description='solution with ionized gas transport - stage 1')
except ImportError:
    pass

f.solve(loglevel=loglevel, stage=2, enable_energy=True)
try:
    f.write_hdf('ion_burner_flame.h5', group='ion',
                description='solution with ionized gas transport - stage 2')
except ImportError:
    f.save('ion_burner_flame.xml', 'mix', 'solution with ionized gas transport')

f.write_csv('ion_burner_flame.csv', quiet=False)

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
    ax3.set_xlabel('Distance from burner [m]')
    ax3.set_ylabel('E [V/m]')
    ax3.set_xlim([0, .005])

    fig.tight_layout()
    plt.show()
