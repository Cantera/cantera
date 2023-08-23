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
points would be concentrated upstream of the flame, where the flamefront had
been previously. (To see this, try setting prune to zero.)

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, premixed flame, strained flame
"""

from pathlib import Path
import cantera as ct


# parameter values
p = 0.05 * ct.one_atm  # pressure
tburner = 373.0  # burner temperature
tsurf = 500.0

# each mdot value will be solved to convergence, with grid refinement, and
# then that solution will be used for the next mdot
mdot = [0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12]  # kg/m^2/s

rxnmech = 'h2o2.yaml'  # reaction mechanism file
comp = 'H2:1.8, O2:1, AR:7'  # premixed gas composition

# The solution domain is chosen to be 20 cm
width = 0.2  # m

loglevel = 1  # amount of diagnostic output (0 to 5)

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
sim = ct.ImpingingJet(gas=gas, width=width)

# set the mass flow rate at the inlet
sim.inlet.mdot = mdot[0]

# set the surface state
sim.surface.T = tsurf

sim.set_grid_min(1e-4)
sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products
sim.show()

sim.solve(loglevel, auto=True)

output_path = Path() / "stagnation_flame_data"
output_path.mkdir(parents=True, exist_ok=True)

if "native" in ct.hdf_support():
    output = output_path / "stagnation_flame.h5"
else:
    output = output_path / "stagnation_flame.yaml"
output.unlink(missing_ok=True)

for m, md in enumerate(mdot):
    sim.inlet.mdot = md
    sim.solve(loglevel)
    sim.save(output, name=f"mdot-{m}", description=f"mdot = {md} kg/m2/s")

    # write the velocity, temperature, and mole fractions to a CSV file
    sim.save(output_path / f"stagnation_flame_{m}.csv", basis="mole", overwrite=True)

sim.show_stats()
