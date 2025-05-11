"""

Requires: cantera >= XX.

.. tags:: Python, plasma
"""


import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import sys


# solution
gas = ct.Solution('two-step-air-plasma_Phelps.yaml')
gas.TPX = 300., 101325., 'N2:0.79, O2:0.21, O2+:1E-9, Electron:1E-9'
gas.EN = 200.0 * 1e-21 # Reduced electric field [V.m^2]
gas.update_EEDF()

# plasma reactor
r = ct.PlasmaReactor(gas)
r.dis_vol = 5e-3*np.pi*(1e-3)**2/4
print(r.dis_vol)
print(r.dis_power)

sim = ct.ReactorNet([r])
sim.verbose = True

# constant EN for 5 ns
t_end = 5e-9
dt_max = 1e-10
states = ct.SolutionArray(gas, extra=['t'])

print('{:10s} {:10s} {:10s} {:14s}'.format(
    't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
while sim.time < t_end:
    sim.advance(sim.time + dt_max)
    states.append(r.thermo.state, t=sim.time)
    print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(
            sim.time, r.T, r.thermo.P, r.thermo.u))

gas.EN = 1e-21
gas.update_EEDF()

# reinitialize cvode because of the steep change in reaction rates
sim.reinitialize()

# stop EN
t_end = 20.0e-9
dt_max = 1e-10
while sim.time < t_end:
    sim.advance(sim.time + dt_max)
    states.append(r.thermo.state, t=sim.time)
    print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(
            sim.time, r.T, r.thermo.P, r.thermo.u))

# figure
fig, ax = plt.subplots(2)

ax[0].plot(states.t, states.X[:, gas.species_index('Electron')], label='Electron')
ax[0].plot(states.t, states.X[:, gas.species_index('O2+')], label='O2+')
ax[0].plot(states.t, states.X[:, gas.species_index('O')], label='O')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(B3)')], label='N2(B3)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(C3)')], label='N2(C3)')
ax[0].set_yscale('log')
ax[0].set_ylim([1e-9, 1])

ax[1].plot(states.t, states.T, label='T')

for axx in ax:
    axx.legend(loc='lower right')
    axx.set_xlabel('Time [s]')

ax[0].set_ylabel('Mole fraction [-]')
ax[1].set_ylabel('Temperature [K]')

plt.tight_layout()
plt.show()
