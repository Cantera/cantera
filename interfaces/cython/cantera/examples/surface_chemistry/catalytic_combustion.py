"""
CATCOMB -- Catalytic combustion of methane on platinum.

This script solves a catalytic combustion problem. A stagnation flow is set
up, with a gas inlet 10 cm from a platinum surface at 900 K. The lean,
premixed methane/air mixture enters at ~6 cm/s (0.06 kg/m2/s), and burns
catalytically on the platinum surface. Gas-phase chemistry is included too,
and has some effect very near the surface.

The catalytic combustion mechanism is from Deutschman et al., 26th
Symp. (Intl.) on Combustion,1996 pp. 1747-1754

Requires: cantera >= 2.5.0
"""

import numpy as np
import cantera as ct

#  Parameter values are collected here to make it easier to modify them
p = ct.one_atm  # pressure
tinlet = 300.0  # inlet temperature
tsurf = 900.0  # surface temperature
mdot = 0.06  # kg/m^2/s
transport = 'Mix'  # transport model

# We will solve first for a hydrogen/air case to use as the initial estimate
# for the methane/air case

# composition of the inlet premixed gas for the hydrogen/air case
comp1 = 'H2:0.05, O2:0.21, N2:0.78, AR:0.01'

# composition of the inlet premixed gas for the methane/air case
comp2 = 'CH4:0.095, O2:0.21, N2:0.78, AR:0.01'

# The inlet/surface separation is 10 cm.
width = 0.1  # m

loglevel = 1  # amount of diagnostic output (0 to 5)

################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic, and
# transport properties. The gas phase will be taken from the definition of
# phase 'gas' in input file 'ptcombust.yaml,' which is a stripped-down version
# of GRI-Mech 3.0.
gas = ct.Solution('ptcombust.yaml', 'gas', transport_model=transport)
gas.TPX = tinlet, p, comp1

################ create the interface object ##################
#
# This object will be used to evaluate all surface chemical production rates.
# It will be created from the interface definition 'Pt_surf' in input file
# 'ptcombust.yaml,' which implements the reaction mechanism of Deutschmann et
# al., 1995 for catalytic combustion on platinum.
#
surf_phase = ct.Interface('ptcombust.yaml', 'Pt_surf', [gas])
surf_phase.TP = tsurf, p

# integrate the coverage equations in time for 1 s, holding the gas
# composition fixed to generate a good starting estimate for the coverages.
surf_phase.advance_coverages(1.0)

# create the object that simulates the stagnation flow, and specify an initial
# grid
sim = ct.ImpingingJet(gas=gas, width=width, surface=surf_phase)

# Objects of class StagnationFlow have members that represent the gas inlet
# ('inlet') and the surface ('surface'). Set some parameters of these objects.
sim.inlet.mdot = mdot
sim.inlet.T = tinlet
sim.inlet.X = comp1
sim.surface.T = tsurf

# Show the initial solution estimate
sim.show_solution()

# Solving problems with stiff chemistry coulpled to flow can require a
# sequential approach where solutions are first obtained for simpler problems
# and used as the initial guess for more difficult problems.

# disable the surface coverage equations, and turn off all gas and surface
# chemistry.
sim.surface.coverage_enabled = False
surf_phase.set_multiplier(0.0)
gas.set_multiplier(0.0)

# solve the problem, refining the grid if needed, to determine the non-
# reacting velocity and temperature distributions
sim.solve(loglevel, auto=True)

# now turn on the surface coverage equations, and turn the chemistry on slowly
sim.surface.coverage_enabled = True
for mult in np.logspace(-5, 0, 6):
    surf_phase.set_multiplier(mult)
    gas.set_multiplier(mult)
    print('Multiplier =', mult)
    sim.solve(loglevel)

# At this point, we should have the solution for the hydrogen/air problem.
sim.show_solution()

# Now switch the inlet to the methane/air composition.
sim.inlet.X = comp2

# set more stringent grid refinement criteria
sim.set_refine_criteria(100.0, 0.15, 0.2, 0.0)

# solve the problem for the final time
sim.solve(loglevel)

# show the solution
sim.show_solution()

# save the solution in XML format. The 'restore' method can be used to restart
# a simulation from a solution stored in this form.
try:
    sim.write_hdf('catalytic_combustion.h5', group='soln1', mode='w',
                  description='catalytic combustion example')
except ImportError:
    sim.save("catalytic_combustion.xml", "soln1")

# save selected solution components in a CSV file for plotting in
# Excel or MATLAB.
sim.write_csv('catalytic_combustion.csv', quiet=False)

sim.show_stats(0)
