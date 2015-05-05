"""
An opposed-flow premixed strained flame

This script simulates a lean hydrogen-oxygen flame stabilized in a strained
flowfield, with an opposed flow consisting of equilibrium products.
"""

import cantera as ct
import numpy as np
import os

# parameter values
p = 0.05 * ct.one_atm  # pressure
T_in = 373.0  # inlet temperature
mdot_reactants = 0.12  # kg/m^2/s
mdot_products = 0.06  # kg/m^2/s
rxnmech = 'h2o2.cti'  # reaction mechanism file
comp = 'H2:1.6, O2:1, AR:7'  # premixed gas composition

initial_grid = np.linspace(0.0, 0.2, 12)  # m
tol_ss = [1.0e-7, 1.0e-13]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-7, 1.0e-11]  # [rtol atol] for time stepping
loglevel = 1  # amount of diagnostic output (0 to 5)

# Grid refinement parameters
ratio = 3
slope = 0.1
curve = 0.2
prune = 0.02

# Set up the problem
gas = ct.Solution(rxnmech)

# set state to that of the unburned gas at the burner
gas.TPX = T_in, p, comp

# Create the flame simulation object
sim = ct.CounterflowPremixedFlame(gas=gas, grid=initial_grid)

# set the boundary flow rates
sim.reactants.mdot = mdot_reactants
sim.products.mdot = mdot_products

sim.flame.set_steady_tolerances(default=tol_ss)
sim.flame.set_transient_tolerances(default=tol_ts)
sim.set_initial_guess()  # assume adiabatic equilibrium products
sim.show_solution()

sim.energy_enabled = False
sim.solve(loglevel, False)

sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
sim.energy_enabled = True
sim.solve(loglevel)

# write the velocity, temperature, and mole fractions to a CSV file
sim.write_csv('premixed_counterflow.csv', quiet=False)
sim.show_stats()
sim.show_solution()
