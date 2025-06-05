"""
Nanosecond Pulse Plasma Simulation
==================================

This example simulates a nanosecond-scale pulse discharge in a reactor.
A Gaussian-shaped electric field pulse is applied over a short timescale.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, plasma
"""

import cantera as ct
ct.CanteraError.set_stack_trace_depth(10)

import numpy as np
import matplotlib.pyplot as plt

# Gaussian pulse parameters
EN_peak = 190 * 1e-21  # Td
pulse_center = 24e-9
pulse_width = 3e-9    # standard deviation in ns

def gaussian_EN(t):
    return EN_peak * np.exp(-((t - pulse_center)**2) / (2 * pulse_width**2))

# setup
gas = ct.Solution('example_data/gri30_plasma_cpavan.yaml')
gas.TPX = 300., 101325., 'CH4:0.095, O2:0.19, N2:0.715, e:1E-11'
gas.EN = gaussian_EN(0)
gas.update_EEDF()

r = ct.ConstPressureReactor(gas, energy="off")
#r.dis_vol = 5e-3 * np.pi * (1e-3)**2 / 4

sim = ct.ReactorNet([r])
sim.verbose = False

# simulation parameters
t_total = 90e-9
dt_max = 1e-10
dt_chunk = 1e-9  # 1 ns chunk
states = ct.SolutionArray(gas, extra=['t'])

print('{:>10} {:>10} {:>10} {:>14}'.format('t [s]', 'T [K]', 'P [Pa]', 'h [J/kg]'))

# simulate in 1 ns chunks
t = 0.0
while t < t_total:

    # integrate over the next chunk
    t_end = min(t + dt_chunk, t_total)
    while sim.time < t_end:
        sim.advance(sim.time + dt_max) #use sim.step
        states.append(r.thermo.state, t=sim.time)
        print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(
            sim.time, r.T, r.thermo.P, r.thermo.h))

    EN_t = gaussian_EN(t)
    gas.EN = EN_t
    gas.update_EEDF()

    # reinitialize integrator with new source terms
    sim.reinitialize()

    t = t_end

# Plotting
fig, ax = plt.subplots(2)

ax[0].plot(states.t, states.X[:, gas.species_index('e')], label='e')
ax[0].plot(states.t, states.X[:, gas.species_index('O2+')], label='O2+')
ax[0].plot(states.t, states.X[:, gas.species_index('N2+')], label='N2+')
ax[0].plot(states.t, states.X[:, gas.species_index('H2O+')], label='H2O+')
ax[0].plot(states.t, states.X[:, gas.species_index('CH4+')], label='CH4+')
ax[0].plot(states.t, states.X[:, gas.species_index('O')], label='O')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(A)')], label='N2(A)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(B)')], label='N2(B)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(C)')], label='N2(C)')
ax[0].plot(states.t, states.X[:, gas.species_index("N2(a')")], label="N2(a')")
ax[0].plot(states.t, states.X[:, gas.species_index('CH3')], label='CH3', linestyle='--')
ax[0].plot(states.t, states.X[:, gas.species_index('CO2')], label='CO2', linestyle='--')
ax[0].plot(states.t, states.X[:, gas.species_index('CO')], label='CO', linestyle='--')
ax[0].plot(states.t, states.X[:, gas.species_index('H2O')], label='H2O', linestyle='--')
ax[0].plot(states.t, states.X[:, gas.species_index('H')], label='H', linestyle='--')
ax[0].plot(states.t, states.X[:, gas.species_index('OH')], label='OH', linestyle='--')
# N2 vibrational states
""" ax[0].plot(states.t, states.X[:, gas.species_index('N2(v1)')], label='N2(v1)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v2)')], label='N2(v2)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v3)')], label='N2(v3)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v4)')], label='N2(v4)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v5)')], label='N2(v5)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v6)')], label='N2(v6)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v7)')], label='N2(v7)')
ax[0].plot(states.t, states.X[:, gas.species_index('N2(v8)')], label='N2(v8)') """

ax[0].set_yscale('log')
ax[0].set_ylim([1e-14, 1e-3])

ax[1].plot(states.t, states.T, label='T')
ax2 = ax[1].twinx()
EN_values = [gaussian_EN(t) for t in states.t]
ax2.plot(states.t, EN_values, label='E/N', color='tab:red', linestyle='--')
ax2.set_ylabel('E/N', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')

for axx in ax:
    axx.legend(loc='lower right')
    axx.set_xlabel('Time [s]')

ax[0].set_ylabel('Mole fraction [-]')
ax[1].set_ylabel('Temperature [K]')

plt.tight_layout()
plt.show()