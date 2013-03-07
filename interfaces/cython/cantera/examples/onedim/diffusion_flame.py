"""
An opposed-flow ethane/air diffusion flame
"""

import cantera as ct
import numpy as np

# Input parameters
p = ct.one_atm  # pressure
tin_f = 300.0  # fuel inlet temperature
tin_o = 300.0  # oxidizer inlet temperature
mdot_o = 0.72  # kg/m^2/s
mdot_f = 0.24  # kg/m^2/s

comp_o = 'O2:0.21, N2:0.78, AR:0.01'  # air composition
comp_f = 'C2H6:1'  # fuel composition

# Distance between inlets is 2 cm.
# Start with an evenly-spaced 6-point grid.
initial_grid = np.linspace(0, 0.02, 6)

tol_ss = [1.0e-5, 1.0e-12]  # [rtol, atol] for steady-state problem
tol_ts = [5.0e-4, 1.0e-11]  # [rtol, atol] for time stepping

loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = 1   # 1 to enable refinement, 0 to disable

# Create the gas object used to evaluate all thermodynamic, kinetic, and
# transport properties.
gas = ct.Solution('gri30.xml', 'gri30_mix')
gas.TP = gas.T, p

# Create an object representing the counterflow flame configuration,
# which consists of a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right.
f = ct.CounterflowDiffusionFlame(gas, initial_grid)

# Set the state of the two inlets
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.X = comp_f
f.fuel_inlet.T = tin_f

f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.X = comp_o
f.oxidizer_inlet.T = tin_o

# Set error tolerances
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# construct the initial solution estimate. To do so, it is necessary
# to specify the fuel species. If a fuel mixture is being used,
# specify a representative species here for the purpose of
# constructing an initial guess.
f.set_initial_guess(fuel='C2H6')

# First disable the energy equation and solve the problem without
# refining the grid
f.energy_enabled = False
f.solve(loglevel, refine_grid=False)

# Now specify grid refinement criteria, turn on the energy equation,
# and solve the problem again.
f.energy_enabled = True
f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)
f.solve(loglevel, refine_grid=refine_grid)
f.show_solution()
f.save('c2h6_diffusion.xml')

# write the velocity, temperature, and mole fractions to a CSV file
f.write_csv('c2h6_diffusion.csv', quiet=False)

f.show_stats(0)
