"""
Catalytic combustion of methane on platinum
===========================================

This script solves a catalytic combustion problem. A stagnation flow is set
up, with a gas inlet 10 cm from a platinum surface at 900 K. The lean,
premixed methane/air mixture enters at ~6 cm/s (0.06 kg/m2/s), and burns
catalytically on the platinum surface. Gas-phase chemistry is included too,
and has some effect very near the surface.

The catalytic combustion mechanism is from Deutschmann et al., 26th
Symp. (Intl.) on Combustion,1996 pp. 1747-1754

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, catalysis, combustion, 1D flow, surface chemistry
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# %%
# Problem Definition
# ------------------
#
# Parameter values are collected here to make it easier to modify them
p = ct.one_atm  # pressure
tinlet = 300.0  # inlet temperature
tsurf = 900.0  # surface temperature
mdot = 0.06  # kg/m^2/s
transport = 'mixture-averaged'  # transport model

# %%
# We will solve first for a hydrogen/air case to use as the initial estimate
# for the methane/air case

# composition of the inlet premixed gas for the hydrogen/air case
comp1 = 'H2:0.05, O2:0.21, N2:0.78, AR:0.01'

# composition of the inlet premixed gas for the methane/air case
comp2 = 'CH4:0.095, O2:0.21, N2:0.78, AR:0.01'

# The inlet/surface separation is 10 cm.
width = 0.1  # m

loglevel = 1  # amount of diagnostic output (0 to 5)

# %%
# Create the phase objects
# ------------------------
#
# The ``surf_phase`` object will be used to evaluate all surface chemical production
# rates. It will be created from the interface definition ``Pt_surf`` in input file
# ``ptcombust.yaml``, which implements the reaction mechanism of Deutschmann et
# al., 1995 for catalytic combustion on platinum.
#
# This phase definition also references the phase ``gas`` in the same input file,
# which will be created and used to evaluate all thermodynamic, kinetic, and
# transport properties. It is a stripped-down version of GRI-Mech 3.0.
surf_phase = ct.Interface("ptcombust.yaml", "Pt_surf")
surf_phase.TP = tsurf, p
gas = surf_phase.adjacent["gas"]
gas.TPX = tinlet, p, comp1

# integrate the coverage equations in time for 1 s, holding the gas
# composition fixed to generate a good starting estimate for the coverages.
surf_phase.advance_coverages(1.0)

# create the object that simulates the stagnation flow, and specify an initial
# grid
sim = ct.ImpingingJet(gas=gas, width=width, surface=surf_phase)

# Objects of class ImpingingJet have members that represent the gas inlet
# ('inlet') and the surface ('surface'). Set some parameters of these objects.
sim.inlet.mdot = mdot
sim.inlet.T = tinlet
sim.inlet.X = comp1
sim.surface.T = tsurf

# %%
# Show the initial solution estimate
sim.show()

# %%
# Solving problems with stiff chemistry coupled to flow can require a
# sequential approach where solutions are first obtained for simpler problems
# and used as the initial guess for more difficult problems.

# disable the surface coverage equations, and turn off all gas and surface
# chemistry.
sim.surface.coverage_enabled = False
surf_phase.set_multiplier(0.0)
gas.set_multiplier(0.0)

# %%
# solve the problem, refining the grid if needed, to determine the non-
# reacting velocity and temperature distributions
sim.solve(loglevel, auto=True)

# %%
# now turn on the surface coverage equations, and turn the chemistry on slowly
sim.surface.coverage_enabled = True
for mult in np.logspace(-5, 0, 6):
    surf_phase.set_multiplier(mult)
    gas.set_multiplier(mult)
    print('Multiplier =', mult)
    sim.solve(loglevel)

# %%
# At this point, we should have the solution for the hydrogen/air problem.
sim.show()

# %%
# Now switch the inlet to the methane/air composition.
sim.inlet.X = comp2

# set more stringent grid refinement criteria
sim.set_refine_criteria(100.0, 0.15, 0.2, 0.0)

# solve the problem for the final time
sim.solve(loglevel)

# %%
# show the solution
sim.show()

# %%
# Save the full solution to HDF or YAML container files. The ``restore`` method can be
# used to restore or restart a simulation from a solution stored in this form.
if "native" in ct.hdf_support():
    filename = "catalytic_combustion.h5"
else:
    filename = "catalytic_combustion.yaml"
sim.save(filename, "soln1", description="catalytic combustion example",
         overwrite=True)

# save selected solution components in a CSV file for plotting in Excel or MATLAB.
sim.save('catalytic_combustion.csv', basis="mole", overwrite=True)

sim.show_stats(0)

# %%
# Temperature Profile
# -------------------
fig, ax = plt.subplots()
ax.plot(sim.grid, sim.T, color='C3')
ax.set_ylabel('heat release rate [MW/m³]')
ax.set(xlabel='distance from inlet [m]')
plt.show()

# %%
# Major Species Profiles
# ----------------------
fig, ax = plt.subplots()
major = ('O2', 'CH4', 'H2O', 'CO2')
states = sim.to_array()
ax.plot(states.grid, states(*major).X, label=major)
ax.set(xlabel='distance from inlet [m]', ylabel='mole fractions')
ax.legend()
plt.show()
