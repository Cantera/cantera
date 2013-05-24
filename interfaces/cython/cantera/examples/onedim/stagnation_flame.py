"""
A detached flat flame stabilized at a stagnation point

This script simulates a lean hydrogen-oxygen flame stabilized in a strained
flowfield at an axisymmetric stagnation point on a non-reacting surface. The
solution begins with a flame attached to the inlet (burner), and the mass flow
rate is progressively increased, causing the flame to detach and move closer
to the surface.

This example illustrates use of the new 'prune' grid refinement parameter,
which allows grid points to be removed if they are no longer required to
resolve the solution. This is important here, since the flamefront moves as
the mass flowrate is increased. Without using 'prune', a large number of grid
points would be concentrated upsteam of the flame, where the flamefront had
been previously. (To see this, try setting prune to zero.)
"""

import cantera as ct
import numpy as np
import os

# parameter values
p = 0.05 * ct.one_atm  # pressure
tburner = 373.0  # burner temperature
tsurf = 500.0

# each mdot value will be solved to convergence, with grid refinement, and
# then that solution will be used for the next mdot
mdot = [0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12]  # kg/m^2/s

rxnmech = 'h2o2.cti'  # reaction mechanism file
comp = 'H2:1.8, O2:1, AR:7'  # premixed gas composition

# The solution domain is chosen to be 50 cm, and a point very near the
# downstream boundary is added to help with the zero-gradient boundary
# condition at this boundary.
initial_grid = np.linspace(0.0, 0.2, 12)  # m

tol_ss = [1.0e-5, 1.0e-13]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-4, 1.0e-9]  # [rtol atol] for time stepping
loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = True

# Grid refinement parameters
ratio = 3
slope = 0.1
curve = 0.2
prune = 0.06

# Set up the problem
gas = ct.Solution(rxnmech)

# set state to that of the unburned gas at the burner
gas.TPX = tburner, p, comp

# Create the stagnation flow object with a non-reactive surface.  (To make the
# surface reactive, supply a surface reaction mechanism. See example
# catalytic_combustion.py for how to do this.)
sim = ct.ImpingingJet(gas=gas, grid=initial_grid)

# set the properties at the inlet
sim.inlet.mdot = mdot[0]
sim.inlet.X = comp
sim.inlet.T = tburner

# set the surface state
sim.surface.T = tsurf

sim.flame.set_steady_tolerances(default=tol_ss)
sim.flame.set_transient_tolerances(default=tol_ts)
sim.set_grid_min(1e-4)
sim.energy_enabled = False

sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products
sim.show_solution()

sim.solve(loglevel, refine_grid)

sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
sim.energy_enabled = True

outfile = 'stflame1.xml'
if os.path.exists(outfile):
    os.remove(outfile)

for m,md in enumerate(mdot):
    sim.inlet.mdot = md
    sim.solve(loglevel,refine_grid)
    sim.save(outfile, 'mdot{0}'.format(m), 'mdot = {0} kg/m2/s'.format(md))

    # write the velocity, temperature, and mole fractions to a CSV file
    sim.write_csv('stflame1_{0}.csv'.format(m), quiet=False)

sim.show_stats()
