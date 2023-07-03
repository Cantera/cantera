"""
An opposed-flow premixed strained flame

This script simulates a lean hydrogen-oxygen flame stabilized in a strained
flowfield, with an opposed flow consisting of equilibrium products.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, premixed flame, strained flame
"""

from pathlib import Path
import cantera as ct


# parameter values
p = 0.05 * ct.one_atm  # pressure
T_in = 373.0  # inlet temperature
mdot_reactants = 0.12  # kg/m^2/s
mdot_products = 0.06  # kg/m^2/s
rxnmech = 'h2o2.yaml'  # reaction mechanism file
comp = 'H2:1.6, O2:1, AR:7'  # premixed gas composition

width = 0.2  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

# Set up the problem
gas = ct.Solution(rxnmech)

# set state to that of the unburned gas at the burner
gas.TPX = T_in, p, comp

# Create the flame simulation object
sim = ct.CounterflowPremixedFlame(gas=gas, width=width)

# Set grid refinement parameters
sim.set_refine_criteria(ratio=3, slope=0.1, curve=0.2, prune=0.02)

# set the boundary flow rates
sim.reactants.mdot = mdot_reactants
sim.products.mdot = mdot_products

sim.set_initial_guess()  # assume adiabatic equilibrium products
sim.show()

sim.solve(loglevel, auto=True)

if "native" in ct.hdf_support():
    output = Path() / "premixed_counterflow_flame.h5"
else:
    output = Path() / "premixed_counterflow_flame.yaml"
output.unlink(missing_ok=True)

sim.save(output, name="mix", description="solution with mixture-averaged transport")

# write the velocity, temperature, and mole fractions to a CSV file
sim.save("premixed_counterflow_flame.csv", basis="mole", overwrite=True)
sim.show_stats()
sim.show()
