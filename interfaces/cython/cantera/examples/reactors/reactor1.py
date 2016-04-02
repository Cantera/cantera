"""
Constant-pressure, adiabatic kinetics simulation.
"""

import sys
import numpy as np

import cantera as ct

gri3 = ct.Solution('gri30.xml')
gri3.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r = ct.IdealGasConstPressureReactor(gri3)

sim = ct.ReactorNet([r])
time = 0.0
times = np.zeros(100)
states = ct.SolutionArray(gri3, 100)

print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
for n in range(100):
    time += 1.e-5
    sim.advance(time)
    times[n] = time * 1e3  # time in ms
    states[n].TDY = r.thermo.TDY
    print('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T,
                                           r.thermo.P, r.thermo.u))

# Plot the results if matplotlib is installed.
# See http://matplotlib.org/ to get it.
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(times, states.T)
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2, 2, 2)
    plt.plot(times, states.X[:,gri3.species_index('OH')])
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(times, states.X[:,gri3.species_index('H')])
    plt.xlabel('Time (ms)')
    plt.ylabel('H Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(times, states.X[:,gri3.species_index('H2')])
    plt.xlabel('Time (ms)')
    plt.ylabel('H2 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")
