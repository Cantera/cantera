"""
Nanosecond Pulse Plasma Simulation
==================================

This example simulates a nanosecond-scale pulse discharge in a reactor.
A Gaussian-shaped electric field pulse is applied over a short timescale.
The plasma reaction mechanism used is based on Colin Pavan's mechanism
for methane-air plasmas which is described in his Ph.D. dissertation and the
corresponding AIAA SciTech conference papers:

C. A. Pavan, "Nanosecond Pulsed Plasmas in Dynamic Combustion Environments,"
Ph.D. Thesis, Massachusetts Institute of Technology, 2023. Chapter 5.

C. A. Pavan and C. Guerra-Garcia, "Modelling the Impact of a Repetitively
Pulsed Nanosecond DBD Plasma on a Mesoscale Flame," in AIAA SCITECH 2022
Forum, Reston, Virginia: American Institute of Aeronautics and
Astronautics, Jan. 2022, pp. 1-15. DOI: 10.2514/6.2022-0975.

C. A. Pavan and C. Guerra-Garcia, "Modeling Flame Speed Modification by
Nanosecond Pulsed Discharges to Inform Experimental Design," in AIAA
SCITECH 2023 Forum, Reston, Virginia: American Institute of Aeronautics
and Astronautics, Jan. 2023, pp. 1-15. DOI: 10.2514/6.2023-2056.

Requires: cantera >= 3.2, matplotlib >= 2.0

.. tags:: Python, plasma, reactor network
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Gaussian pulse parameters
EN_peak = 190 * 1e-21  # 190 Td
pulse_center = 24e-9  # 24 ns
pulse_width = 3e-9  # standard deviation (3 ns)
pulse_fwhm = pulse_width * 2 * (2 * np.log(2))**.5
gaussian_EN = ct.Func1("Gaussian", [EN_peak, pulse_center, pulse_fwhm])

# setup
gas = ct.Solution('example_data/methane-plasma-pavan-2023.yaml')
gas.TPX = 300., 101325., 'CH4:0.095, O2:0.19, N2:0.715, e:1E-11'
gas.reduced_electric_field = gaussian_EN(0)
gas.update_electron_energy_distribution()

r = ct.ConstPressureReactor(gas, energy="off", clone=False)

sim = ct.ReactorNet([r])
sim.verbose = False

# simulation parameters
t_total = 90e-9
dt_max = 1e-10
dt_chunk = 1e-9  # 1 ns chunk
states = ct.SolutionArray(gas, extra=['t'])

print(f"{'t [s]':>10} {'T [K]':>10} {'P [Pa]':>10} {'h [J/kg]':>14}")

# simulate in 1 ns chunks
t = 0.0
while t < t_total:

    # integrate over the next chunk
    t_end = min(t + dt_chunk, t_total)
    while sim.time < t_end:
        sim.advance(sim.time + dt_max) #use sim.step
        states.append(r.contents.state, t=sim.time)
        print(f"{sim.time:10.3e} {r.T:10.3f} {r.contents.P:10.3f} {r.contents.h:14.6f}")

    EN_t = gaussian_EN(t)
    gas.reduced_electric_field = EN_t
    gas.update_electron_energy_distribution()

    # reinitialize integrator with new source terms
    sim.reinitialize()

    t = t_end

# Plotting
fig, ax = plt.subplots(2, layout="constrained")

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
    axx.legend(loc='lower right', ncol=2)
    axx.set_xlabel('Time [s]')

ax[0].set_ylabel('Mole fraction [-]')
ax[1].set_ylabel('Temperature [K]')

plt.show()
